"""
Ensemble reference-based solver.

Combines internal backbones (NNLS / weighted / marker / ridge) with weights
derived from gene-masking-CV reconstruction (better CV → more weight) and a
conditioning penalty.  Proportions are combined at the fine level and
renormalised; weights and reasons are recorded.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.solver.base import BaseSolver, SolverResult
from tissueresolve.solver.nnls import NNLSSolver
from tissueresolve.solver.weighted_nnls import WeightedNNLSSolver
from tissueresolve.solver.marker_nnls import MarkerNNLSSolver
from tissueresolve.solver.ridge_nnls import RidgeNNLSSolver

__all__ = ["EnsembleSolver"]


class EnsembleSolver(BaseSolver):
    name = "ensemble_nnls"

    def __init__(self, components=None, mask_fraction: float = 0.15,
                 n_splits: int = 2, seed: int = 0):
        self.components = components or [
            NNLSSolver(), WeightedNNLSSolver(), MarkerNNLSSolver(), RidgeNNLSSolver()]
        self.mask_fraction = mask_fraction
        self.n_splits = n_splits
        self.seed = seed

    def solve(self, query, ref) -> SolverResult:
        from tissueresolve.validation.gene_masking import compare_solvers_by_gene_masking

        def make_fn(sv):
            return lambda b, r, genes: sv.__class__(
                **{k: v for k, v in sv.__dict__.items() if k != "genes"}
            ).solve(b, r).proportions if "genes" not in sv.__dict__ else \
                _solve_on(sv, b, r, genes)
        # CV weights from masked-gene reconstruction
        solve_fns = {sv.name: (lambda b, r, g, _sv=sv: _solve_on(_sv, b, r, g))
                     for sv in self.components}
        cv = compare_solvers_by_gene_masking(query, ref, solve_fns,
                                             mask_fraction=self.mask_fraction,
                                             n_splits=self.n_splits, seed=self.seed)
        scores = cv["masked_gene_pearson"].fillna(0.0).clip(lower=0.0)
        if scores.sum() <= 0:
            scores = pd.Series(1.0, index=scores.index)
        weights = (scores / scores.sum()).to_dict()

        # full-data component predictions, aligned to a common cell-type set
        preds, cols = {}, None
        for sv in self.components:
            p = sv.solve(query, ref).proportions
            preds[sv.name] = p
            cols = list(p.columns) if cols is None else cols
        combined = sum(weights.get(n, 0.0) * preds[n][cols] for n in preds)
        combined = combined.div(combined.sum(axis=1), axis=0).fillna(0.0)
        return SolverResult(combined, cols,
                            {"solver": self.name, "weights": weights,
                             "cv_scores": scores.to_dict()})


def _solve_on(sv, bulk_train, ref, genes):
    """Re-run a solver restricted to *genes* (for CV)."""
    kind = sv.__class__
    if isinstance(sv, MarkerNNLSSolver):
        # marker solver picks its own genes; intersect with train genes
        res = sv.solve(bulk_train, ref.subset_genes(list(genes)))
        return res.proportions
    inst = kind(genes=list(genes)) if "genes" in sv.__dict__ else kind()
    return inst.solve(bulk_train, ref).proportions
