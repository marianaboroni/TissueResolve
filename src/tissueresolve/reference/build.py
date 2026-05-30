"""
Unified reference matrix construction for TissueResolve.

Combines:
- CHIMERA's donor-aware aggregation and cross-donor CV computation.
- SpatCAR's batched streaming, CPM normalisation, and NB overdispersion
  estimation.

The output :class:`~tissueresolve.results.ReferenceSignature` carries both
the ``phi`` (L1-normalised, G × K) representation required by the bulk solver
and the ``R_cpm`` / ``R_log`` (K × G) representation required by the spatial
NB model.

Algorithm
---------
For each cell type k:

1. Filter cells assigned to k.
2. If donors are provided and ≥ 2 donors are present:
   - Compute per-donor mean expression vectors.
   - Average donor means → mean expression profile (donor-aware).
   - Compute cross-donor coefficient of variation (CV) per gene.
3. Otherwise: simple cell mean.
4. CPM-normalise: ``cpm = mean / Σ_g(mean) × 10⁶``.
5. Store ``cpm + pseudocount`` as ``R_cpm[k]``.
6. L1-normalise ``R_cpm[k]`` to form ``phi[:, k]``.

Optional NB overdispersion estimation (for spatial models):
    Uses Welford's online algorithm (memory-efficient) when the count matrix
    is loaded batch-wise, or vectorised variance for in-memory matrices.

Input validation
----------------
- Duplicate gene names → ``ValueError``.
- Duplicate cell IDs → ``ValueError``.
- Missing cell-type column → ``KeyError``.
- Cell types with fewer than ``min_cells`` cells → ``UserWarning``; dropped.
- All cell types dropped → ``ValueError``.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

from tissueresolve.config import ReferenceConfig
from tissueresolve.io.validation import (
    validate_cell_type_column,
    validate_counts_matrix,
    validate_no_duplicate_genes,
)
from tissueresolve.results import ReferenceSignature

__all__ = ["ReferenceBuilder"]

logger = logging.getLogger("tissueresolve.reference.build")

_DEFAULT_PSEUDOCOUNT: float = 1.0   # CPM pseudocount added to every gene-type pair
_PHI_G_MIN: float = 0.5
_PHI_G_MAX: float = 100.0


class ReferenceBuilder:
    """Build a :class:`~tissueresolve.results.ReferenceSignature`.

    Parameters
    ----------
    config:
        :class:`~tissueresolve.config.ReferenceConfig`.
    pseudocount:
        CPM pseudocount added before log-transform (default 1.0).
    """

    def __init__(
        self,
        config: Optional[ReferenceConfig] = None,
        *,
        pseudocount: float = _DEFAULT_PSEUDOCOUNT,
    ) -> None:
        self.cfg = config or ReferenceConfig()
        self.pseudocount = pseudocount

    # ------------------------------------------------------------------
    # Entry points
    # ------------------------------------------------------------------

    def build_from_df(
        self,
        counts: pd.DataFrame,
        metadata: Union[pd.DataFrame, pd.Series],
        *,
        estimate_overdispersion: bool = False,
    ) -> ReferenceSignature:
        """Build from a genes × cells count matrix.

        Parameters
        ----------
        counts:
            Genes as rows, cells as columns.  Raw counts (integer or float).
            Index = gene identifiers, columns = cell IDs.
        metadata:
            Cell metadata indexed by cell ID.  Must contain at minimum
            the column named by ``config.celltype_col``.
            May be a Series (treated as single column = celltype_col).
        estimate_overdispersion:
            When True, estimate per-gene NB overdispersion (phi_g) from the
            count matrix using the method-of-moments estimator.  Required
            for the spatial deconvolution model; optional for bulk.

        Returns
        -------
        ReferenceSignature
        """
        if isinstance(metadata, pd.Series):
            metadata = metadata.rename(self.cfg.celltype_col).to_frame()

        # --- input validation ---
        validate_counts_matrix(counts)

        if counts.columns.duplicated().any():
            dupes = counts.columns[counts.columns.duplicated()].tolist()
            raise ValueError(
                f"Duplicate cell IDs in counts columns: {dupes[:5]}…  "
                "Each cell must have a unique identifier."
            )

        validate_cell_type_column(metadata, self.cfg.celltype_col)

        # --- cell alignment ---
        common = counts.columns.intersection(metadata.index)
        if len(common) == 0:
            raise ValueError(
                "No shared cell IDs between counts.columns and metadata.index.  "
                "Ensure that cell identifiers use the same format in both."
            )
        if len(common) < len(counts.columns):
            logger.warning(
                "Cell alignment: %d / %d cells in counts are in metadata.  "
                "Proceeding with the intersection.",
                len(common), len(counts.columns),
            )
        counts = counts[common]
        metadata = metadata.loc[common]

        # --- resolve donors ---
        donors: Optional[list[str]] = None
        if self.cfg.donor_col and self.cfg.donor_col in metadata.columns:
            donors = metadata[self.cfg.donor_col].astype(str).tolist()

        X = counts.to_numpy(dtype=np.float64)   # (G, C)
        cell_types = metadata[self.cfg.celltype_col].astype(str).tolist()

        return self._aggregate(
            X=X,
            gene_names=list(counts.index),
            cell_types=cell_types,
            donors=donors,
            estimate_overdispersion=estimate_overdispersion,
        )

    def build_from_csv(
        self,
        counts_path: Union[Path, str],
        meta_path: Union[Path, str],
        *,
        estimate_overdispersion: bool = False,
    ) -> ReferenceSignature:
        """Build from TSV/CSV files.

        Parameters
        ----------
        counts_path:
            Path to genes × cells count matrix (TSV or CSV).
        meta_path:
            Path to cell metadata file.
        estimate_overdispersion:
            See :meth:`build_from_df`.
        """
        from tissueresolve.io.reference import load_reference_csv
        counts, metadata = load_reference_csv(
            counts_path, meta_path, cell_type_col=self.cfg.celltype_col
        )
        return self.build_from_df(
            counts.astype(float),
            metadata,
            estimate_overdispersion=estimate_overdispersion,
        )

    def build_from_adata(
        self,
        adata: "anndata.AnnData",  # type: ignore[name-defined]
        *,
        estimate_overdispersion: bool = False,
    ) -> ReferenceSignature:
        """Build from an in-memory AnnData object.

        Parameters
        ----------
        adata:
            Single-cell AnnData.  ``adata.X`` must contain raw (integer) UMI
            counts.  ``adata.obs[celltype_col]`` must be present.
        estimate_overdispersion:
            See :meth:`build_from_df`.

        Raises
        ------
        ImportError
            If ``anndata`` is not installed.
        KeyError
            If the cell-type column is missing from ``adata.obs``.
        """
        try:
            import scipy.sparse as sp
        except ImportError as exc:
            raise ImportError(
                "scipy is required.  Install it with: pip install scipy"
            ) from exc

        validate_cell_type_column(adata.obs, self.cfg.celltype_col)

        X = adata.X
        if sp.issparse(X):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float64).T   # transpose: (G, C)

        gene_names = list(adata.var_names)
        validate_no_duplicate_genes(gene_names, source="adata.var_names")

        cell_types = adata.obs[self.cfg.celltype_col].astype(str).tolist()
        donors: Optional[list[str]] = None
        if self.cfg.donor_col and self.cfg.donor_col in adata.obs.columns:
            donors = adata.obs[self.cfg.donor_col].astype(str).tolist()

        return self._aggregate(
            X=X,
            gene_names=gene_names,
            cell_types=cell_types,
            donors=donors,
            estimate_overdispersion=estimate_overdispersion,
        )

    def build_from_h5ad(
        self,
        path: Union[Path, str],
        *,
        estimate_overdispersion: bool = False,
        backed: bool = False,
    ) -> ReferenceSignature:
        """Build from a ``.h5ad`` file.

        Parameters
        ----------
        path:
            Path to the ``.h5ad`` file.
        estimate_overdispersion:
            See :meth:`build_from_df`.
        backed:
            Open in backed (memory-mapped) mode.  Set ``True`` for large
            references that do not fit in RAM.

        Raises
        ------
        FileNotFoundError
        ImportError
        """
        from tissueresolve.io.reference import load_reference_h5ad
        adata = load_reference_h5ad(
            path,
            cell_type_col=self.cfg.celltype_col,
            backed=backed,
        )
        return self.build_from_adata(
            adata, estimate_overdispersion=estimate_overdispersion
        )

    # ------------------------------------------------------------------
    # Core aggregation
    # ------------------------------------------------------------------

    def _aggregate(
        self,
        X: np.ndarray,            # (G, C) raw counts
        gene_names: list[str],
        cell_types: list[str],
        donors: Optional[list[str]],
        *,
        estimate_overdispersion: bool,
    ) -> ReferenceSignature:
        G, C = X.shape
        ct_arr = np.asarray(cell_types)
        unique_cts = sorted(set(cell_types))

        cpm_profiles: dict[str, np.ndarray] = {}   # ct -> (G,) CPM
        cv_profiles: dict[str, np.ndarray] = {}    # ct -> (G,) CV
        dropped: list[str] = []

        for ct in unique_cts:
            mask = ct_arr == ct
            n_cells = int(mask.sum())

            if n_cells < self.cfg.min_cells:
                dropped.append(ct)
                warnings.warn(
                    f"Cell type '{ct}' has only {n_cells} cell(s) "
                    f"(min_cells={self.cfg.min_cells}); dropping.",
                    stacklevel=5,
                )
                continue

            ct_X = X[:, mask]   # (G, n_cells)

            if donors is not None:
                donor_arr = np.asarray(donors)[mask]
                unique_d = np.unique(donor_arr)

                if len(unique_d) >= 2:
                    # Donor-aware: mean of per-donor means
                    d_means = np.stack(
                        [ct_X[:, donor_arr == d].mean(axis=1) for d in unique_d],
                        axis=1,
                    )   # (G, D)
                    mean_expr = d_means.mean(axis=1)   # (G,)

                    # Cross-donor CV
                    mu_d = d_means.mean(axis=1)   # (G,)
                    sd_d = d_means.std(axis=1, ddof=1)   # (G,)
                    cv = np.where(mu_d > 1e-9, sd_d / mu_d, 0.0)
                    cv_profiles[ct] = cv.astype(np.float64)
                else:
                    mean_expr = ct_X.mean(axis=1)
            else:
                mean_expr = ct_X.mean(axis=1)

            # CPM-normalise
            total = float(mean_expr.sum())
            cpm_profiles[ct] = (mean_expr / total * 1e6) if total > 0 else mean_expr

        if dropped:
            logger.info("Dropped %d type(s): %s", len(dropped), ", ".join(dropped))

        if not cpm_profiles:
            raise ValueError(
                "No cell types survived the min_cells filter.  "
                f"Lower ReferenceConfig.min_cells (current: {self.cfg.min_cells})."
            )

        labels = sorted(cpm_profiles)
        K = len(labels)

        # R_cpm: (K, G) with pseudocount
        R_cpm_mat = np.stack(
            [(cpm_profiles[lbl] + self.pseudocount).astype(np.float32) for lbl in labels],
            axis=0,
        )

        R_log_mat = np.log1p(R_cpm_mat).astype(np.float32)

        # phi: (G, K) = L1-normalised R_cpm
        phi_raw = R_cpm_mat.T.copy()                          # (G, K)
        col_sums = phi_raw.sum(axis=0, keepdims=True)
        col_sums = np.where(col_sums == 0, 1.0, col_sums)
        phi_mat = (phi_raw / col_sums).astype(np.float64)

        # donor_cv: (G, K)
        donor_cv: Optional[np.ndarray] = None
        if cv_profiles:
            cv_mat = np.zeros((G, K), dtype=np.float64)
            for j, lbl in enumerate(labels):
                if lbl in cv_profiles:
                    cv_mat[:, j] = cv_profiles[lbl]
            donor_cv = cv_mat

        # NB overdispersion
        phi_g: Optional[np.ndarray] = None
        if estimate_overdispersion:
            phi_g = self._estimate_overdispersion(X)

        # n_cells_per_type (only for types that survived)
        n_cells_per_type = {ct: int((ct_arr == ct).sum()) for ct in labels}

        logger.info(
            "Reference built: %d genes × %d cell types (%s).",
            G, K, ", ".join(labels),
        )
        if dropped:
            logger.info("Dropped cell types (below min_cells): %s", ", ".join(dropped))

        return ReferenceSignature(
            gene_names=list(gene_names),
            cell_types=labels,
            phi=phi_mat,
            R_cpm=R_cpm_mat,
            R_log=R_log_mat,
            phi_g=phi_g,
            donor_cv=donor_cv,
            n_cells_per_type=n_cells_per_type,
            genome=self.cfg.genome,
        )

    # ------------------------------------------------------------------
    # Overdispersion estimation
    # ------------------------------------------------------------------

    def _estimate_overdispersion(
        self, X: np.ndarray
    ) -> np.ndarray:
        """Method-of-moments NB overdispersion estimator.

        Parameters
        ----------
        X:
            ``(G, C)`` raw count matrix.

        Returns
        -------
        np.ndarray
            Shape ``(G,)`` float32, clipped to [0.5, 100].

        Notes
        -----
        φ̂_g = μ̂_g² / max(σ̂_g² − μ̂_g, ε)

        When σ̂² ≤ μ (sub-Poisson), φ is set to ``_PHI_G_MAX`` (treated as
        Poisson for that gene).
        """
        XT = X.T.astype(np.float64)   # (C, G)
        n = XT.shape[0]
        if n < 2:
            logger.warning(
                "Only %d cell(s) — cannot estimate overdispersion reliably.  "
                "Returning phi_g = %.1f for all genes.",
                n, _PHI_G_MAX,
            )
            return np.full(X.shape[0], _PHI_G_MAX, dtype=np.float32)

        mean_g = XT.mean(axis=0)       # (G,)
        var_g = XT.var(axis=0, ddof=1) # (G,)
        excess_var = var_g - mean_g

        phi_g = np.where(
            excess_var > 0,
            mean_g ** 2 / np.maximum(excess_var, 1e-8),
            _PHI_G_MAX,
        )
        result = np.clip(phi_g, _PHI_G_MIN, _PHI_G_MAX).astype(np.float32)

        logger.debug(
            "Overdispersion: median=%.2f, mean=%.2f, range=[%.2f, %.2f].",
            float(np.median(result)), float(result.mean()),
            float(result.min()), float(result.max()),
        )
        return result
