"""
Spatial reference-query mismatch correction for TissueResolve.

This module is relevant for **spatial (Visium) deconvolution only**.

Conceptual distinction from bulk protocol risk
----------------------------------------------
``protocol.risk`` (bulk)
    Identifies gene-level biases introduced by *known* protocol differences
    between bulk RNA-seq and the single-cell reference.  These biases are
    assessed *before* deconvolution and used to exclude or down-weight
    problem genes from the marker panel.

``protocol.mismatch`` (spatial — this module)
    Estimates per-gene multiplicative scale factors ``d_g`` that absorb the
    *empirical* expression difference between a Visium dataset and the
    pseudo-bulk reference.  These factors are unknown before fitting and are
    *estimated jointly with the cell-type proportions during deconvolution*.
    They correct for a wide range of unmeasured differences (library
    preparation, ambient RNA, batch effects) without requiring explicit
    protocol labels.

The class name :class:`SpatialMismatch` (and its alias
:data:`ProtocolMismatch`) is deliberately prefixed with ``Spatial`` to
prevent confusion with :class:`~protocol.risk.ProtocolRiskReport`.

Algorithm (from SpatCAR v1.1)
------------------------------
1. **Initialisation** — compute the mean library-size-normalised expression of
   each marker gene across all Visium spots (:func:`compute_spatial_discordance`).
   Divide by the mean CPM expression in the reference to get initial ``d_g``.

2. **Refinement** — every few deconvolution iterations, update ``d_g`` using
   the closed-form ratio estimator (:func:`update_mismatch_factors`):

   .. math::
       d_g^{\\text{new}} =
       \\frac{\\sum_s Y_{sg}}{\\sum_s \\hat{\\mu}_{sg}^{\\text{no-}d}}

   Highly discordant genes receive lower ``gene_weight`` and are updated more
   slowly towards their empirical ratio.

References
----------
SpatCAR v1.1 spatial deconvolution algorithm.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np

from tissueresolve.results import ReferenceSignature

__all__ = [
    "SpatialMismatch",
    "ProtocolMismatch",   # alias for backward-compat and user familiarity
    "compute_spatial_discordance",
    "compute_discordance",   # alias
    "update_mismatch_factors",
]

logger = logging.getLogger("tissueresolve.protocol.mismatch")


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------


@dataclass
class SpatialMismatch:
    """Per-gene spatial reference-query mismatch state.

    Carries the multiplicative scale factors ``d_g`` that correct for
    systematic expression differences between a Visium dataset and the
    pseudo-bulk reference used for deconvolution.

    **NOT** a protocol risk assessment.  Compare with
    :class:`~tissueresolve.protocol.risk.ProtocolRiskReport` which is a
    bulk-specific diagnostic computed before deconvolution.

    Attributes
    ----------
    d_g:
        Multiplicative scale factors per marker gene, shape ``(G_m,)``.
        Values close to 1 indicate no mismatch; > 1 means the gene is
        under-represented in the reference relative to Visium.
        Clipped to [0.1, 10.0].
    discord_score:
        Per-gene discordance ``|log2(mean_V + 1) − log2(mean_R + 1)|``,
        shape ``(G_m,)``.  Higher = more discordant.
    gene_weight:
        Down-weighting factor for discordant genes, shape ``(G_m,)``.
        Computed as ``exp(−discord_score / tau)``, clipped to [0.05, 1.0].
        Used to slow down ``d_g`` updates for highly discordant genes.
    gene_names:
        Ordered marker gene names, length ``G_m``.
    """

    d_g: np.ndarray            # (G_m,) float32
    discord_score: np.ndarray  # (G_m,) float32
    gene_weight: np.ndarray    # (G_m,) float32
    gene_names: list[str]

    def __post_init__(self) -> None:
        G_m = len(self.gene_names)
        if self.d_g.shape != (G_m,):
            raise ValueError(
                f"d_g shape {self.d_g.shape} does not match "
                f"len(gene_names) = {G_m}."
            )
        if self.discord_score.shape != (G_m,):
            raise ValueError(
                f"discord_score shape {self.discord_score.shape} does not "
                f"match len(gene_names) = {G_m}."
            )
        if self.gene_weight.shape != (G_m,):
            raise ValueError(
                f"gene_weight shape {self.gene_weight.shape} does not match "
                f"len(gene_names) = {G_m}."
            )

    @property
    def G_m(self) -> int:
        """Number of marker genes."""
        return len(self.gene_names)

    def summary(self) -> dict:
        """Return a summary dict for logging / saving in model metadata."""
        high_discord = self.discord_score > 2.0
        return {
            "n_marker_genes": self.G_m,
            "n_high_discord": int(high_discord.sum()),
            "mean_discord": float(np.mean(self.discord_score)),
            "d_g_median": float(np.median(self.d_g)),
            "d_g_min": float(self.d_g.min()),
            "d_g_max": float(self.d_g.max()),
        }


# Alias for backward-compat and to match the user task specification.
ProtocolMismatch = SpatialMismatch


# ---------------------------------------------------------------------------
# Discordance computation (initialisation)
# ---------------------------------------------------------------------------


def compute_spatial_discordance(
    visium_expr: "np.ndarray | anndata.AnnData",  # type: ignore[name-defined]
    ref: ReferenceSignature,
    marker_genes: list[str],
    *,
    lib_sizes: Optional[np.ndarray] = None,
    tau: float = 1.0,
    d_clip: tuple[float, float] = (0.1, 10.0),
    batch_size: int = 500,
) -> SpatialMismatch:
    """Compute per-gene spatial discordance and initialise mismatch factors.

    Accepts either an AnnData object (full Visium dataset) or a pre-subsetted
    ``(N, G_m)`` dense count matrix for testing and non-AnnData workflows.

    Parameters
    ----------
    visium_expr:
        Either:

        - ``anndata.AnnData`` — full Visium dataset.  Gene names are looked
          up in ``adata.var_names``.  Library sizes are read from
          ``adata.obs["total_counts"]``.
        - ``np.ndarray`` of shape ``(N, G_m)`` — pre-subsetted marker-gene
          count matrix (N spots × G_m marker genes, raw counts).  In this
          case *lib_sizes* must be provided.

    ref:
        :class:`~tissueresolve.results.ReferenceSignature` **already subsetted**
        to *marker_genes*.  ``ref.gene_names`` must equal *marker_genes*.
    marker_genes:
        Ordered list of marker gene names (length ``G_m``).
    lib_sizes:
        Per-spot library sizes, shape ``(N,)``.  Required when *visium_expr*
        is a numpy array.  Ignored when *visium_expr* is AnnData (read from
        ``adata.obs["total_counts"]``).
    tau:
        Down-weighting temperature.  Higher tau → less aggressive penalisation
        of discordant genes.
    d_clip:
        ``(lo, hi)`` clip range for the initial scale factors.
    batch_size:
        Batch size when computing Visium mean expression from AnnData.

    Returns
    -------
    SpatialMismatch
        Initialised mismatch object ready for the deconvolution model.

    Raises
    ------
    ValueError
        When *visium_expr* is a numpy array and *lib_sizes* is not provided,
        or when the marker gene lists are inconsistent.
    """
    G_m = len(marker_genes)

    # Verify ref gene alignment
    if list(ref.gene_names) != list(marker_genes):
        if set(ref.gene_names) != set(marker_genes):
            raise ValueError(
                "ref.gene_names and marker_genes do not match.  "
                "Call ref.subset_genes(marker_genes) before calling "
                "compute_spatial_discordance."
            )
        ref = ref.subset_genes(marker_genes)

    # --- Visium mean expression ----------------------------------------------
    try:
        import anndata as _ad
        _has_anndata = True
    except ImportError:
        _has_anndata = False

    if _has_anndata and hasattr(visium_expr, "obs"):
        # AnnData path
        mean_visium_cpm = _mean_from_adata(visium_expr, marker_genes, batch_size)
    elif isinstance(visium_expr, np.ndarray):
        if lib_sizes is None:
            raise ValueError(
                "lib_sizes must be provided when visium_expr is a numpy array."
            )
        mean_visium_cpm = _mean_from_array(visium_expr, lib_sizes)
    else:
        raise TypeError(
            f"visium_expr must be an np.ndarray or anndata.AnnData, "
            f"got {type(visium_expr).__name__}."
        )

    # --- Reference mean expression -------------------------------------------
    # Simple mean of CPM values across cell types (unweighted)
    R_cpm = ref.as_R_cpm().astype(np.float64)        # (K, G_m)
    mean_ref_cpm = R_cpm.mean(axis=0)                 # (G_m,)

    # --- Discordance score ---------------------------------------------------
    eps = 1.0   # 1 CPM floor — guards against log2(0)
    log2_vis = np.log2(mean_visium_cpm + eps)
    log2_ref = np.log2(mean_ref_cpm + eps)
    discord = np.abs(log2_vis - log2_ref).astype(np.float32)

    # --- Gene weight ---------------------------------------------------------
    gene_weight = np.exp(-discord / tau).clip(0.05, 1.0).astype(np.float32)

    # --- Initial scale factors -----------------------------------------------
    d_g = (mean_visium_cpm / np.maximum(mean_ref_cpm, eps)).astype(np.float32)
    d_g = np.clip(d_g, d_clip[0], d_clip[1])

    n_high = int((discord > 2.0).sum())
    logger.info(
        "Spatial discordance: %d / %d marker genes have |log2FC| > 2 "
        "(tau=%.1f, d_g median=%.3f).",
        n_high, G_m, tau, float(np.median(d_g)),
    )

    return SpatialMismatch(
        d_g=d_g,
        discord_score=discord,
        gene_weight=gene_weight,
        gene_names=list(marker_genes),
    )


# Alias matching the original SpatCAR function name.
compute_discordance = compute_spatial_discordance


# ---------------------------------------------------------------------------
# Scale-factor update (refinement during deconvolution)
# ---------------------------------------------------------------------------


def update_mismatch_factors(
    Y_marker: np.ndarray,     # (N, G_m) dense float32 observed counts
    Pi: np.ndarray,           # (N, K) current proportions
    R_lin: np.ndarray,        # (K, G_m) reference in proportion scale (CPM / 1e6)
    lib_sizes: np.ndarray,    # (N,) per-spot library sizes
    mismatch: SpatialMismatch,
    *,
    d_clip: tuple[float, float] = (0.1, 10.0),
) -> SpatialMismatch:
    """Update per-gene scale factors using the closed-form ratio estimator.

    Called every few deconvolution iterations.  Returns a *new*
    :class:`SpatialMismatch` object (the input is not mutated).

    Algorithm
    ---------
    For gene *g*:

    .. math::
        d_g^{\\text{new}} =
        \\frac{\\sum_s Y_{sg}}{\\sum_s \\hat{\\mu}_{sg}^{\\text{no-}d}}

    where :math:`\\hat{\\mu}_{sg}^{\\text{no-}d} = \\text{lib}_s \\cdot
    (\\pi_s \\cdot R_{\\cdot g})`.

    The update is then damped by ``gene_weight`` so that highly discordant
    genes change more slowly:

    .. math::
        d_g^{\\text{updated}} =
        w_g \\cdot d_g^{\\text{new}} + (1 - w_g) \\cdot d_g^{\\text{old}}

    Parameters
    ----------
    Y_marker:
        Dense observed counts, shape ``(N, G_m)``.
    Pi:
        Current proportion estimates, shape ``(N, K)``.
    R_lin:
        Reference in proportion scale ``(CPM / 1e6)``, shape ``(K, G_m)``.
    lib_sizes:
        Library sizes, shape ``(N,)``.
    mismatch:
        Current :class:`SpatialMismatch` (not mutated).
    d_clip:
        ``(lo, hi)`` clip range.

    Returns
    -------
    SpatialMismatch
        New object with updated ``d_g`` values.
    """
    # Expected expression without d_g factor
    Mu_base = lib_sizes[:, None] * (Pi @ R_lin)    # (N, G_m)

    numer = Y_marker.sum(axis=0).astype(np.float64)           # (G_m,)
    denom = np.maximum(Mu_base.sum(axis=0), 1e-6).astype(np.float64)  # (G_m,)

    d_g_new = (numer / denom).astype(np.float32)
    d_g_new = np.clip(d_g_new, d_clip[0], d_clip[1])

    # Damp with gene weight: trusted genes update quickly, discordant slowly
    w = mismatch.gene_weight
    d_g_updated = w * d_g_new + (1.0 - w) * mismatch.d_g

    logger.debug(
        "d_g update: median %.3f → %.3f (range [%.3f, %.3f]).",
        float(np.median(mismatch.d_g)),
        float(np.median(d_g_updated)),
        float(d_g_updated.min()),
        float(d_g_updated.max()),
    )

    return SpatialMismatch(
        d_g=d_g_updated,
        discord_score=mismatch.discord_score,
        gene_weight=mismatch.gene_weight,
        gene_names=mismatch.gene_names,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _mean_from_adata(
    adata: "anndata.AnnData",  # type: ignore[name-defined]
    marker_genes: list[str],
    batch_size: int,
) -> np.ndarray:
    """Compute mean library-size-normalised expression × 1e6 (CPM) for
    marker genes from an AnnData object."""
    import scipy.sparse as sp

    visium_genes = adata.var_names.tolist()
    gene_to_col = {g: i for i, g in enumerate(visium_genes)}
    marker_col_idx = np.array(
        [gene_to_col[g] for g in marker_genes if g in gene_to_col], dtype=np.intp
    )
    if len(marker_col_idx) != len(marker_genes):
        missing = [g for g in marker_genes if g not in gene_to_col]
        raise ValueError(
            f"{len(missing)} marker genes missing from Visium data: {missing[:5]}…"
        )

    lib_sizes = adata.obs["total_counts"].values.astype(np.float64)
    n_spots = adata.n_obs
    G_m = len(marker_genes)

    norm_sum = np.zeros(G_m, dtype=np.float64)
    for start in range(0, n_spots, batch_size):
        end = min(start + batch_size, n_spots)
        X_batch = adata.X[start:end]
        if sp.issparse(X_batch):
            X_batch = X_batch[:, marker_col_idx].toarray()
        else:
            X_batch = np.asarray(X_batch[:, marker_col_idx], dtype=np.float64)
        lib_batch = lib_sizes[start:end]
        norm_sum += (X_batch / np.maximum(lib_batch[:, None], 1.0)).sum(axis=0)

    mean_norm = norm_sum / max(n_spots, 1)
    return (mean_norm * 1e6).astype(np.float64)


def _mean_from_array(
    Y: np.ndarray,        # (N, G_m) raw counts
    lib_sizes: np.ndarray,  # (N,)
) -> np.ndarray:
    """Compute mean library-size-normalised expression × 1e6 (CPM) from
    a pre-subsetted dense count matrix."""
    Y = np.asarray(Y, dtype=np.float64)
    lib_safe = np.maximum(lib_sizes, 1.0).astype(np.float64)
    mean_norm = (Y / lib_safe[:, None]).mean(axis=0)
    return (mean_norm * 1e6).astype(np.float64)
