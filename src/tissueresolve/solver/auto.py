"""
``solver="auto"`` — pick the backbone that generalises best.

Candidates (NNLS, weighted NNLS, marker NNLS, ridge NNLS) are scored by
gene-masking cross-validation (reconstruction of held-out genes) and penalised
for ill-conditioning, so auto does **not** simply maximise in-sample
reconstruction (which can overfit / increase spillover).  The selected solver,
its score, and the comparison table are recorded.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.solver.base import BaseSolver, SolverResult
from tissueresolve.solver.nnls import NNLSSolver
from tissueresolve.solver.weighted_nnls import WeightedNNLSSolver
from tissueresolve.solver.marker_nnls import MarkerNNLSSolver
from tissueresolve.solver.ridge_nnls import RidgeNNLSSolver

__all__ = ["AutoSolver", "candidate_solvers"]


def candidate_solvers() -> list[BaseSolver]:
    return [NNLSSolver(), WeightedNNLSSolver(), MarkerNNLSSolver(), RidgeNNLSSolver()]


class AutoSolver(BaseSolver):
    name = "auto"

    def __init__(self, candidates=None, mask_fraction: float = 0.15,
                 n_splits: int = 3, seed: int = 0):
        self.candidates = candidates or candidate_solvers()
        self.mask_fraction = mask_fraction
        self.n_splits = n_splits
        self.seed = seed

    def _cv_table(self, query, ref) -> pd.DataFrame:
        from tissueresolve.validation.gene_masking import compare_solvers_by_gene_masking
        from tissueresolve.solver.ensemble import _solve_on
        fns = {sv.name: (lambda b, r, g, _sv=sv: _solve_on(_sv, b, r, g))
               for sv in self.candidates}
        return compare_solvers_by_gene_masking(
            query, ref, fns, mask_fraction=self.mask_fraction,
            n_splits=self.n_splits, seed=self.seed)

    def select(self, query, ref):
        """Return (best_solver, comparison_df, reason)."""
        cv = self._cv_table(query, ref)
        # multi-objective: CV reconstruction minus a small conditioning penalty
        rows = []
        for sv in self.candidates:
            full = sv.solve(query, ref)
            cond = full.diagnostics.get("condition_number", np.inf)
            cv_score = float(cv.loc[sv.name, "masked_gene_pearson"]) \
                if sv.name in cv.index else float("nan")
            cond_pen = 0.05 * np.tanh((np.log10(cond + 1)) / 3.0) if np.isfinite(cond) else 0.05
            score = (0.0 if np.isnan(cv_score) else cv_score) - cond_pen
            rows.append({"solver": sv.name, "cv_masked_gene_pearson": cv_score,
                         "condition_number": cond, "objective": score})
        comp = pd.DataFrame(rows).set_index("solver").sort_values(
            "objective", ascending=False)
        best_name = comp.index[0]
        best = next(sv for sv in self.candidates if sv.name == best_name)
        reason = (f"selected '{best_name}' by gene-masking CV "
                  f"(masked-gene Pearson={comp.loc[best_name, 'cv_masked_gene_pearson']:.3f}, "
                  f"condition≈{comp.loc[best_name, 'condition_number']:.1f}); "
                  "objective = CV reconstruction − conditioning penalty (not "
                  "in-sample reconstruction alone).")
        return best, comp, reason

    def solve(self, query, ref) -> SolverResult:
        best, comp, reason = self.select(query, ref)
        res = best.solve(query, ref)
        res.diagnostics.update({
            "solver": self.name, "selected_solver": best.name,
            "selection_reason": reason,
            "solver_comparison": comp.reset_index().to_dict(orient="records"),
        })
        return res
