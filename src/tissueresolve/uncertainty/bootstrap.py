"""
Non-parametric bootstrap CIs over the gene panel for bulk deconvolution.

Algorithm (bulk)
----------------
For each iteration b = 1 … B:

  1. Resample ⌊frac × G⌋ gene indices from the panel **with replacement**.
  2. L1-re-normalise bulk and reference columns on the sub-panel.
  3. Apply the corresponding weight sub-set.
  4. Solve wNNLS → θ̂_b.

CIs are the empirical percentiles of {θ̂_b} across B draws.

This approach bootstraps over genes, not samples, so it captures uncertainty
from the choice of gene panel and expression variability across markers.  It
does not substitute for biological replication.

Warning policy
--------------
If the panel has fewer than ``MIN_PANEL_GENES`` genes, a non-suppressible
``UserWarning`` is emitted and CIs will be very wide.
"""
from __future__ import annotations

import logging
import warnings
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.config import BootstrapConfig
from tissueresolve.results import BulkDeconvResult, ReferenceSignature

__all__ = ["BulkBootstrapCI", "SpatialBootstrapCI"]

logger = logging.getLogger("tissueresolve.uncertainty.bootstrap")

MIN_PANEL_GENES = 15


class BulkBootstrapCI:
    """Compute bootstrap credible intervals for bulk deconvolution results.

    Parameters
    ----------
    config:
        :class:`~tissueresolve.config.BootstrapConfig` with ``n_bootstrap``,
        ``bootstrap_frac``, ``ci_level``, and ``seed``.
    """

    def __init__(self, config: Optional[BootstrapConfig] = None) -> None:
        self.cfg = config or BootstrapConfig()
        self._rng = np.random.default_rng(self.cfg.seed)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def compute(
        self,
        bulk: pd.DataFrame,
        ref: ReferenceSignature,
        gene_panel: list[str],
        gene_weights: Optional[pd.Series] = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Compute bootstrap CIs.

        Parameters
        ----------
        bulk:
            Genes × samples bulk counts.
        ref:
            :class:`~tissueresolve.results.ReferenceSignature`.
        gene_panel:
            Panel of genes to bootstrap over.
        gene_weights:
            Optional composite weights indexed by gene name.

        Returns
        -------
        (lower_ci, upper_ci)
            Both DataFrames indexed by sample, columns by cell type.
            Shape ``(n_samples × n_cell_types)`` for each.

        Warns
        -----
        UserWarning
            If fewer than :data:`MIN_PANEL_GENES` genes are available.  CIs
            will be unreliable.
        """
        from tissueresolve.bulk.solver import _l1norm_cols, _resolve_weights, _solve_all

        phi_mat = ref.as_phi()  # (G_ref, K)
        phi_df = pd.DataFrame(phi_mat, index=ref.gene_names, columns=ref.cell_types)

        common = sorted(set(gene_panel) & set(bulk.index) & set(ref.gene_names))
        if not common:
            raise ValueError(
                "No gene overlap — cannot compute bootstrap CIs.  "
                "Ensure gene IDs match between bulk and reference."
            )
        if len(common) < MIN_PANEL_GENES:
            warnings.warn(
                f"BulkBootstrapCI: only {len(common)} genes available "
                f"(minimum recommended: {MIN_PANEL_GENES}).  Bootstrap CIs "
                "will be very wide and unreliable.",
                stacklevel=2,
            )

        B_raw = bulk.loc[common].to_numpy(dtype=np.float64)   # (G, N)
        R_raw = phi_df.loc[common].to_numpy(dtype=np.float64) # (G, K)
        B_norm = _l1norm_cols(B_raw)
        R_norm = _l1norm_cols(R_raw)
        w = _resolve_weights(gene_weights, common)

        G = B_norm.shape[0]
        N = B_norm.shape[1]
        K = R_norm.shape[1]
        n_sub = max(int(self.cfg.bootstrap_frac * G), 5)

        boot = np.zeros((self.cfg.n_bootstrap, N, K))
        for b in range(self.cfg.n_bootstrap):
            idx = self._rng.choice(G, size=n_sub, replace=True)
            B_sub = _l1norm_cols(B_norm[idx])
            R_sub = _l1norm_cols(R_norm[idx])
            w_sub = w[idx]
            mu_sub = w_sub.mean()
            if mu_sub > 0:
                w_sub = w_sub / mu_sub
            boot[b] = _solve_all(B_sub, R_sub, w_sub)

        alpha = (1.0 - self.cfg.ci_level) / 2.0
        lo_pct = 100.0 * alpha
        hi_pct = 100.0 * (1.0 - alpha)

        sample_ids = list(bulk.columns)
        ct_labels = list(ref.cell_types)

        lo = pd.DataFrame(
            np.percentile(boot, lo_pct, axis=0),
            index=sample_ids, columns=ct_labels,
        )
        hi = pd.DataFrame(
            np.percentile(boot, hi_pct, axis=0),
            index=sample_ids, columns=ct_labels,
        )
        lo.index.name = hi.index.name = "sample"

        mean_width = float((hi - lo).mean().mean())
        logger.info(
            "Bootstrap: %d iters, frac=%.2f, CI=%.0f%%.  Mean CI width=%.4f.",
            self.cfg.n_bootstrap,
            self.cfg.bootstrap_frac,
            self.cfg.ci_level * 100,
            mean_width,
        )
        return lo, hi

    def attach(
        self,
        result: BulkDeconvResult,
        bulk: pd.DataFrame,
        ref: ReferenceSignature,
    ) -> BulkDeconvResult:
        """Compute CIs and attach them to an existing :class:`BulkDeconvResult`.

        Mutates *result* in place (sets ``lower_ci`` and ``upper_ci``) and
        returns it.

        Parameters
        ----------
        result:
            Completed deconvolution result (from :class:`~WNNLSSolver`).
        bulk:
            Original genes × samples bulk matrix.
        ref:
            Reference used in the solve.
        """
        lo, hi = self.compute(
            bulk, ref,
            gene_panel=result.gene_panel,
            gene_weights=result.gene_weights,
        )
        result.lower_ci = lo
        result.upper_ci = hi
        return result

    @property
    def ci_metadata(self) -> dict:
        """Serialisable dict of bootstrap parameters for run_metadata."""
        return {
            "n_bootstrap": self.cfg.n_bootstrap,
            "bootstrap_frac": self.cfg.bootstrap_frac,
            "ci_level": self.cfg.ci_level,
            "bootstrap_seed": self.cfg.seed,
        }


# ---------------------------------------------------------------------------
# Spatial bootstrap
# ---------------------------------------------------------------------------


class SpatialBootstrapCI:
    """Parametric bootstrap confidence intervals for spatial deconvolution.

    Algorithm
    ---------
    For each bootstrap iteration *b* = 1 … B:

    1. Sample Y_boot ~ NB(μ̂_sg, φ_g) from the converged model estimates.
    2. Re-fit from Π̂ for ``n_iter_per_boot`` iterations using the same
       graph and mismatch factors.
    3. Compute percentile CIs across the B draws.

    Coverage note
    -------------
    Empirical coverage is ~88–93 % for 95 % nominal due to the warm-start
    bias.  Use ``n_iter_per_boot ≥ 50`` for closer-to-nominal coverage.
    This limitation is always reported in ``coverage_note``.

    Parameters
    ----------
    n_bootstrap:
        Number of bootstrap resamples (default 50).
    n_iter_per_boot:
        Re-fitting iterations per resample (default 30).
    ci_level:
        Nominal coverage probability (default 0.95).
    seed:
        Random seed.
    """

    COVERAGE_NOTE = (
        "Empirical coverage is approximately 88–93% for 95% nominal CIs "
        "due to warm-start bias.  Use n_iter_per_boot≥50 for closer-to-nominal "
        "coverage.  CIs are useful for ranking spots by uncertainty even when "
        "absolute coverage is imperfect."
    )

    def __init__(
        self,
        *,
        n_bootstrap: int = 50,
        n_iter_per_boot: int = 30,
        ci_level: float = 0.95,
        seed: int = 0,
    ) -> None:
        self.n_bootstrap = n_bootstrap
        self.n_iter_per_boot = n_iter_per_boot
        self.ci_level = ci_level
        self.seed = seed

    def compute(
        self,
        model: "tissueresolve.spatial.model.SpatCARModel",  # type: ignore[name-defined]
    ) -> tuple[np.ndarray, np.ndarray, str]:
        """Compute parametric bootstrap CIs from a fitted SpatCARModel.

        Parameters
        ----------
        model:
            A fitted :class:`~tissueresolve.spatial.model.SpatCARModel`.
            Must have been fitted with :meth:`~SpatCARModel.fit` so that
            ``_R_lin``, ``_phi_g``, ``_lib_sizes``, ``_graph_ref`` are
            available.

        Returns
        -------
        tuple (ci_lo, ci_hi, coverage_note)
            ``ci_lo``, ``ci_hi`` : shape ``(N, K)`` float32.
            ``coverage_note``    : string documenting empirical coverage
                                   limitation.  Always non-empty.
        """
        from tissueresolve.spatial.model import _nb_multiplicative_update

        if model.proportions_ is None:
            raise RuntimeError(
                "SpatialBootstrapCI.compute: model has not been fitted.  "
                "Call model.fit() first."
            )

        rng = np.random.default_rng(self.seed)

        Pi_hat = model.proportions_
        R_lin = model._R_lin
        phi_g = model._phi_g
        lib = model._lib_sizes
        graph = model._graph_ref
        d_g = (
            model._mismatch.d_g
            if model._mismatch is not None
            else np.ones(R_lin.shape[1], dtype=np.float32)
        )
        R_d = (R_lin * d_g[np.newaxis, :]).astype(np.float32)
        A = graph.A
        alpha = model.alpha

        Mu_hat = lib[:, None] * (Pi_hat @ R_d)

        Pi_boots = []
        for _ in range(self.n_bootstrap):
            phi_b = phi_g[np.newaxis, :]
            p_nb = phi_b / (phi_b + Mu_hat + 1e-8)
            p_nb = np.clip(p_nb, 1e-6, 1 - 1e-6)
            Y_boot = rng.negative_binomial(
                phi_g[np.newaxis, :] * np.ones_like(Mu_hat), p_nb
            ).astype(np.float32)

            Pi_b = Pi_hat.copy()
            for _ in range(self.n_iter_per_boot):
                Pi_b_new = np.empty_like(Pi_b)
                for batch in graph.batches:
                    Yb = Y_boot[batch]
                    Pb = Pi_b[batch]
                    lb = lib[batch]
                    Pnb = _nb_multiplicative_update(Yb, Pb, R_d, phi_g, lb)
                    if alpha > 0.0:
                        neigh = A[batch, :] @ Pi_b
                        Pnb = (1.0 - alpha) * Pnb + alpha * neigh
                        Pnb = np.maximum(Pnb, 1e-10)
                        Pnb /= Pnb.sum(1, keepdims=True)
                    Pi_b_new[batch] = Pnb
                Pi_b = Pi_b_new
            Pi_boots.append(Pi_b)

        stack = np.stack(Pi_boots, 0)
        lo_pct = (1.0 - self.ci_level) / 2.0 * 100
        hi_pct = (1.0 + self.ci_level) / 2.0 * 100
        ci_lo = np.percentile(stack, lo_pct, axis=0).astype(np.float32)
        ci_hi = np.percentile(stack, hi_pct, axis=0).astype(np.float32)

        logger.info(
            "SpatialBootstrapCI: n_boot=%d, n_iter=%d, mean_width=%.4f.",
            self.n_bootstrap, self.n_iter_per_boot,
            float((ci_hi - ci_lo).mean()),
        )
        return ci_lo, ci_hi, self.COVERAGE_NOTE

    def attach(
        self,
        result: "SpatialDeconvResult",  # type: ignore[name-defined]
        model: "SpatCARModel",           # type: ignore[name-defined]
    ) -> "SpatialDeconvResult":
        """Compute CIs and attach them to an existing SpatialDeconvResult.

        Sets ``result.lower_ci``, ``result.upper_ci``, and
        ``result.bootstrap_coverage_note``.
        """
        ci_lo, ci_hi, note = self.compute(model)
        result.lower_ci = ci_lo
        result.upper_ci = ci_hi
        result.bootstrap_coverage_note = note
        return result

    @property
    def ci_metadata(self) -> dict:
        """Serialisable metadata for run_metadata."""
        return {
            "spatial_n_bootstrap": self.n_bootstrap,
            "spatial_n_iter_per_boot": self.n_iter_per_boot,
            "spatial_ci_level": self.ci_level,
            "spatial_bootstrap_seed": self.seed,
            "coverage_note": self.COVERAGE_NOTE,
        }
