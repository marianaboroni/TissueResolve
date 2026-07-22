"""Count-likelihood (Poisson / negative-binomial) GLM deconvolution backbone.

A :class:`BaseSolver`-compatible wrapper so ``deconv_bulk(solver="poisson")`` (and
``"nb"``), ``AutoSolver``, and the benchmark can use a noise model matched to
RNA-seq counts instead of Gaussian/L2 on L1-normalised profiles.

It **reuses** the validated majorise–minimise fit in
``tissueresolve.experimental.nb_bulk_solver.fit_bulk_nb_glm`` (no re-derivation);
this module only adapts inputs/outputs to the ``SolverResult`` interface and adds
an explicit, recorded NNLS fallback on non-convergence or error.

Outputs are **RNA-derived (mRNA) proportions**, never cell fractions — identical
semantics to the other backbones. This backbone is **opt-in / non-default**; it does
not change the production ``BulkPipeline`` (weighted-NNLS) default. Benchmarks
(``benchmarks/rectangle_comparison.py``) show it is the strongest TissueResolve
bulk backbone on donor-held-out pseudobulk (breast + lung), but promotion to
default is gated on real-bulk validation.

Runtime: the MM update is vectorised over samples and genes; on the benchmark
references it converges in well under a second per call.
"""
from __future__ import annotations

import warnings

import numpy as np

from tissueresolve.solver.base import (
    BaseSolver, SolverResult, align_query_to_reference, _normalize)


class PoissonGLMSolver(BaseSolver):
    """Poisson-likelihood GLM backbone (``loss='poisson'``)."""

    name = "poisson"
    _loss = "poisson"

    def __init__(self, genes=None, max_iter: int = 500, tol: float = 1e-6,
                 fallback: bool = True):
        self.genes = genes
        self.max_iter = max_iter
        self.tol = tol
        self.fallback = fallback

    def solve(self, query, ref) -> SolverResult:
        sub, B = align_query_to_reference(query, ref, self.genes)      # B: genes×samples
        R = sub.as_R_cpm().T                                          # genes×K (CPM)
        Phi = R / 1e6                                                 # proportion scale
        cond = self.condition_number(R)

        phi_g = None
        if self._loss == "nb":
            pg = getattr(sub, "phi_g", None)
            if pg is not None:
                pg = np.asarray(pg, dtype=float)
                if pg.shape == (R.shape[0],):
                    phi_g = pg   # dispersion aligned to the panel genes

        diagnostics = {"solver": self.name, "loss": self._loss,
                       "n_genes": len(sub.gene_names), "condition_number": cond,
                       "feature_status": "experimental"}
        try:
            from tissueresolve.experimental.nb_bulk_solver import fit_bulk_nb_glm
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = fit_bulk_nb_glm(
                    B, Phi, gene_dispersion=phi_g, loss=self._loss,
                    max_iter=self.max_iter, tol=self.tol)
            theta = res.theta                                         # K×N (simplex)
            props = np.vstack([_normalize(theta[:, j]) for j in range(B.shape[1])])
            diagnostics.update(converged=bool(res.converged), n_iter=int(res.n_iter),
                               loss_monotonic=bool(res.loss_monotonic),
                               likelihood=res.metadata.get("likelihood", self._loss))
            if not res.converged and self.fallback:
                return self._nnls_fallback(query, ref, sub,
                                           "GLM did not converge", diagnostics)
            return SolverResult(self._frame(props, query.columns, sub.cell_types),
                                list(sub.gene_names), diagnostics)
        except Exception as exc:  # noqa: BLE001 — explicit, recorded fallback
            if not self.fallback:
                raise
            return self._nnls_fallback(query, ref, sub,
                                       f"{type(exc).__name__}: {exc}", diagnostics)

    def _nnls_fallback(self, query, ref, sub, reason, diagnostics):
        from tissueresolve.solver.nnls import NNLSSolver
        fb = NNLSSolver(genes=list(sub.gene_names)).solve(query, ref)
        diagnostics = {**diagnostics, "fallback": "nnls", "fallback_reason": reason,
                       **{f"nnls_{k}": v for k, v in fb.diagnostics.items()}}
        warnings.warn(
            f"{self.name} solver fell back to NNLS ({reason}); estimate is NNLS, "
            "recorded in diagnostics['fallback'].", stacklevel=2)
        return SolverResult(fb.proportions, fb.genes_used, diagnostics)


class NBGLMSolver(PoissonGLMSolver):
    """Negative-binomial GLM backbone (``loss='nb'``); uses reference ``phi_g``."""

    name = "nb"
    _loss = "nb"
