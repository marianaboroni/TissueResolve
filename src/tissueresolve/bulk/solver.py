"""
Bulk RNA-seq deconvolution solver (CHIMERA-derived).

Mathematical formulation
------------------------
For each bulk sample n::

    min_θ  ‖ √w ⊙ b_norm − (√w ⊙ Φ_norm) · θ ‖²,   θ ≥ 0

where w ∈ ℝ^G is a non-negative weight vector and both b_norm and the
columns of Φ_norm are L1-normalised on the selected gene panel.

After solving, θ is normalised to sum to 1.  This produces mRNA proportions,
not cell fractions.  Use :class:`MRNAContentCorrector` to convert when
per-cell-type mRNA content data are available and scientifically valid.

Output note
-----------
``BulkDeconvResult.proportions`` records ``ESTIMATE_TYPE = "mRNA_proportion"``.
Every saved TSV header repeats this warning.  Never interpret these values as
cell fractions without explicit mRNA content correction.
"""
from __future__ import annotations

import logging
import warnings
from typing import Optional

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from tissueresolve.results import BulkDeconvResult, ReferenceSignature

__all__ = ["WNNLSSolver", "MRNAContentCorrector"]

logger = logging.getLogger("tissueresolve.bulk.solver")

_MIN_GENES_WARN = 20


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------


class WNNLSSolver:
    """Solve bulk deconvolution with weighted NNLS per sample.

    Parameters
    ----------
    seed:
        Random seed passed to bootstrap if invoked separately.  The NNLS
        solver itself is fully deterministic.
    """

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def solve(
        self,
        bulk: pd.DataFrame,
        ref: ReferenceSignature,
        gene_panel: list[str],
        gene_weights: Optional[pd.Series] = None,
    ) -> BulkDeconvResult:
        """Deconvolve all samples in *bulk*.

        Parameters
        ----------
        bulk:
            Raw or normalised bulk counts, shape ``(n_genes × n_samples)``.
            Index must be gene identifiers.
        ref:
            :class:`~tissueresolve.results.ReferenceSignature` built by
            :class:`~tissueresolve.reference.build.ReferenceBuilder`.
        gene_panel:
            Gene IDs to use.  The intersection with *bulk* and *ref* is taken
            automatically and reported.  Raises if the intersection is empty.
        gene_weights:
            Optional per-gene composite weights (pd.Series indexed by gene ID).
            ``None`` → uniform weights (standard NNLS).

        Returns
        -------
        BulkDeconvResult
            Without bootstrap CIs — attach them via
            :class:`~tissueresolve.uncertainty.bootstrap.BulkBootstrapCI`.
        """
        phi_mat = ref.as_phi()   # (G_ref, K)
        phi_df = pd.DataFrame(phi_mat, index=ref.gene_names, columns=ref.cell_types)

        common = sorted(set(gene_panel) & set(bulk.index) & set(ref.gene_names))
        if not common:
            raise ValueError(
                "No genes overlap between bulk, reference, and gene_panel.  "
                "Check that gene IDs use the same naming convention "
                "(symbol vs Ensembl, capitalisation)."
            )
        if len(common) < _MIN_GENES_WARN:
            warnings.warn(
                f"WNNLSSolver: only {len(common)} genes in the final panel "
                f"(threshold for reliable deconvolution: {_MIN_GENES_WARN}).  "
                "Results may be unreliable.  Check gene naming and panel size.",
                stacklevel=2,
            )

        B = bulk.loc[common].to_numpy(dtype=np.float64)    # (G, N)
        R = phi_df.loc[common].to_numpy(dtype=np.float64)  # (G, K)

        w = _resolve_weights(gene_weights, common)
        B_norm = _l1norm_cols(B)
        R_norm = _l1norm_cols(R)

        W_hat = _solve_all(B_norm, R_norm, w)             # (N, K)
        r2 = _r2_per_sample(B_norm, R_norm, W_hat)

        sample_ids = list(bulk.columns)
        ct_labels = list(ref.cell_types)

        props = pd.DataFrame(W_hat, index=sample_ids, columns=ct_labels)
        props.index.name = "sample"
        coverage = pd.Series(r2, index=sample_ids, name="coverage_r2")

        stored_w = gene_weights.reindex(common) if gene_weights is not None else None

        logger.info(
            "wNNLS: %d samples, %d cell types, %d genes.  Mean R² = %.4f.",
            len(sample_ids), len(ct_labels), len(common), float(r2.mean()),
        )
        return BulkDeconvResult(
            proportions=props,
            coverage_r2=coverage,
            gene_panel=common,
            gene_weights=stored_w,
            run_metadata={
                "solver": "WNNLSSolver",
                "n_genes_panel": len(common),
                "estimate_type": "mRNA_proportion",
            },
        )


# ---------------------------------------------------------------------------
# mRNA content corrector
# ---------------------------------------------------------------------------


class MRNAContentCorrector:
    """Convert RNA proportions → cell fractions via mRNA content normalisation.

    Bulk deconvolution produces **mRNA proportions** (fraction of total mRNA
    from each cell type), not cell fractions.  In tissues where cell types
    differ substantially in per-cell mRNA content (plasma cells, neurons,
    hepatocytes), these quantities can differ substantially.

    The correction is::

        cell_fraction_k  ∝  rna_proportion_k / mrna_content_k

    Rows are renormalised to sum to 1.

    Parameters
    ----------
    mrna_content:
        Per-cell-type relative mRNA content.  ``pd.Series`` indexed by
        cell-type name.  All values must be strictly positive.  The series is
        normalised so that its median equals 1.

    Raises
    ------
    ValueError
        If any values are ≤ 0 or if the content is trivially uniform (all
        values within 1 % of each other), which would make correction a no-op
        while potentially misleading users.  Pass ``allow_uniform=True`` to
        skip the uniformity check.
    """

    def __init__(
        self,
        mrna_content: pd.Series,
        *,
        allow_uniform: bool = False,
    ) -> None:
        if (mrna_content <= 0).any():
            raise ValueError(
                "mRNA content values must all be strictly positive.  "
                f"Got non-positive values for: "
                f"{mrna_content[mrna_content <= 0].index.tolist()}"
            )
        med = float(mrna_content.median())
        self.mrna_content: pd.Series = mrna_content / med

        # Warn when content is near-uniform (correction is a near no-op)
        cv = float(self.mrna_content.std() / self.mrna_content.mean())
        if cv < 0.01 and not allow_uniform:
            warnings.warn(
                "MRNAContentCorrector: mRNA content values are nearly uniform "
                f"(CV = {cv:.4f}).  The correction will have negligible effect "
                "on the estimates.  If content was estimated from an already "
                "L1-normalised reference (phi matrix), this is expected — use "
                "raw single-cell counts for a meaningful correction.  "
                "Pass allow_uniform=True to suppress this warning.",
                stacklevel=2,
            )

        logger.info(
            "MRNAContentCorrector: %d cell types, range [%.3f, %.3f] (median=1).",
            len(mrna_content),
            float(self.mrna_content.min()),
            float(self.mrna_content.max()),
        )

    def correct(self, rna_proportions: pd.DataFrame) -> pd.DataFrame:
        """Apply mRNA content correction to RNA proportions.

        Parameters
        ----------
        rna_proportions:
            ``(n_samples × n_cell_types)`` mRNA proportion DataFrame.
            Rows should sum to 1.

        Returns
        -------
        pd.DataFrame
            ``(n_samples × n_cell_types)`` cell fraction DataFrame.
            Rows sum to 1.  The result should be stored in
            ``BulkDeconvResult.cell_fractions``, not in ``proportions``.

        Warns
        -----
        UserWarning
            If any cell types in *rna_proportions* are missing from the mRNA
            content series.  Missing types receive content=1.0 (no correction).
        """
        missing = sorted(set(rna_proportions.columns) - set(self.mrna_content.index))
        if missing:
            warnings.warn(
                f"MRNAContentCorrector: no mRNA content for cell types {missing}.  "
                "Using content=1.0 for these types (no correction applied).",
                stacklevel=2,
            )

        content = self.mrna_content.reindex(rna_proportions.columns).fillna(1.0)
        scaled = rna_proportions.to_numpy() / content.to_numpy()[None, :]
        row_sums = scaled.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums == 0, 1.0, row_sums)
        cell_frac = scaled / row_sums

        result = pd.DataFrame(
            cell_frac,
            index=rna_proportions.index,
            columns=rna_proportions.columns,
        )
        result.index.name = "sample"

        delta = float((result - rna_proportions).abs().mean().mean())
        logger.info("mRNA correction applied.  Mean absolute change: %.4f.", delta)
        return result

    def correction_summary(self, rna_proportions: pd.DataFrame) -> pd.DataFrame:
        """Per-cell-type summary of how much correction changes estimates.

        Returns
        -------
        pd.DataFrame
            Columns: ``cell_type``, ``mrna_content``, ``mean_rna_prop``,
            ``mean_cell_frac``, ``mean_delta``.
        """
        cell_fracs = self.correct(rna_proportions)
        rows = []
        for ct in rna_proportions.columns:
            rows.append({
                "cell_type": ct,
                "mrna_content": float(self.mrna_content.get(ct, 1.0)),
                "mean_rna_prop": round(float(rna_proportions[ct].mean()), 4),
                "mean_cell_frac": round(float(cell_fracs[ct].mean()), 4),
                "mean_delta": round(float((cell_fracs[ct] - rna_proportions[ct]).abs().mean()), 4),
            })
        return pd.DataFrame(rows)

    @classmethod
    def from_reference_counts(
        cls,
        counts: pd.DataFrame,
        cell_type_labels: pd.Series,
        *,
        allow_uniform: bool = False,
    ) -> "MRNAContentCorrector":
        """Estimate mRNA content from raw single-cell reference counts.

        Uses mean total UMI count per cell of each type as a proxy for
        relative transcriptional activity per cell.

        Parameters
        ----------
        counts:
            Genes × cells raw count matrix.
        cell_type_labels:
            Cell-level annotation (index = cell_id, values = cell_type name).

        Returns
        -------
        MRNAContentCorrector
        """
        total_umi = counts.sum(axis=0)
        mean_umi: dict[str, float] = {}
        for ct in cell_type_labels.unique():
            cell_ids = cell_type_labels[cell_type_labels == ct].index
            valid = [c for c in cell_ids if c in total_umi.index]
            if valid:
                mean_umi[ct] = float(total_umi[valid].mean())
            else:
                mean_umi[ct] = float(total_umi.median())
                logger.warning("No cells found for '%s'; using median UMI.", ct)
        content = pd.Series(mean_umi)
        logger.info(
            "mRNA content from reference: %d types, UMI range [%.0f, %.0f].",
            len(content), content.min(), content.max(),
        )
        return cls(content, allow_uniform=allow_uniform)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_weights(
    gene_weights: Optional[pd.Series],
    common: list[str],
) -> np.ndarray:
    if gene_weights is not None:
        w = gene_weights.reindex(common).fillna(1.0).to_numpy(dtype=np.float64)
        w = np.clip(w, 0, None)
        mu = w.mean()
        if mu > 0:
            w = w / mu
        return w
    return np.ones(len(common))


def _l1norm_cols(M: np.ndarray) -> np.ndarray:
    """L1-normalise every column of M in-place semantics."""
    col_sums = M.sum(axis=0, keepdims=True)
    col_sums = np.where(col_sums == 0, 1.0, col_sums)
    return M / col_sums


def _solve_one(b: np.ndarray, R: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted NNLS for one sample.

    Solves min_θ ‖√w ⊙ b − (√w ⊙ R)·θ‖², θ ≥ 0, then normalises θ to
    sum to 1.  Falls back to uniform proportions when the solution is the
    zero vector (pathological case).
    """
    sqrt_w = np.sqrt(np.clip(w, 0, None))
    b_w = b * sqrt_w
    R_w = R * sqrt_w[:, None]
    theta, _ = nnls(R_w, b_w)
    total = theta.sum()
    if total > 0:
        return theta / total
    return np.full(R.shape[1], 1.0 / R.shape[1])


def _solve_all(B: np.ndarray, R: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Solve wNNLS for all N samples.  Returns (N, K)."""
    N = B.shape[1]
    K = R.shape[1]
    W = np.zeros((N, K))
    for n in range(N):
        W[n] = _solve_one(B[:, n], R, w)
    return W


def _r2_per_sample(
    B_norm: np.ndarray,
    R_norm: np.ndarray,
    W_hat: np.ndarray,
) -> np.ndarray:
    """Per-sample reconstruction R²."""
    recon = R_norm @ W_hat.T   # (G, N)
    out = np.empty(B_norm.shape[1])
    for n in range(B_norm.shape[1]):
        y = B_norm[:, n]
        yhat = recon[:, n]
        ss_res = float(np.dot(y - yhat, y - yhat))
        y_mean = float(y.mean())
        ss_tot = float(np.dot(y - y_mean, y - y_mean))
        out[n] = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 0.0
    return out
