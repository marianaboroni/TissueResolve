"""Specificity-weighted NNLS (lightweight; not the full protocol pipeline)."""
from __future__ import annotations

import numpy as np

from tissueresolve.solver.base import BaseSolver, SolverResult, align_query_to_reference, _normalize


class WeightedNNLSSolver(BaseSolver):
    name = "weighted_nnls"

    def __init__(self, genes=None):
        self.genes = genes

    def solve(self, query, ref) -> SolverResult:
        from scipy.optimize import nnls
        sub, B = align_query_to_reference(query, ref, self.genes)
        R = sub.as_R_cpm().T
        row = R / (R.sum(axis=1, keepdims=True) + 1e-9)
        w = np.sqrt(np.clip(row.max(axis=1), 1e-9, None))[:, None]
        Rw, Bw = R * w, B * w[:, 0:1]
        props = np.vstack([_normalize(nnls(Rw, Bw[:, j])[0]) for j in range(B.shape[1])])
        return SolverResult(self._frame(props, query.columns, sub.cell_types),
                            list(sub.gene_names),
                            {"solver": self.name, "n_genes": len(sub.gene_names),
                             "condition_number": self.condition_number(R)})
