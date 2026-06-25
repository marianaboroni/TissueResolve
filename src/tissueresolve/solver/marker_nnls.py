"""NNLS on a per-cell-type top-marker panel (selected from the reference)."""
from __future__ import annotations

import numpy as np

from tissueresolve.solver.base import BaseSolver, SolverResult, align_query_to_reference, _normalize


class MarkerNNLSSolver(BaseSolver):
    name = "marker_nnls"

    def __init__(self, top_n: int = 50):
        self.top_n = top_n

    def solve(self, query, ref) -> SolverResult:
        from scipy.optimize import nnls
        logR = np.log1p(ref.as_R_cpm())
        K = logR.shape[0]
        idx = set()
        for k in range(K):
            others = np.delete(logR, k, axis=0).mean(axis=0)
            idx.update(np.argsort(logR[k] - others)[::-1][:self.top_n].tolist())
        genes = [ref.gene_names[i] for i in sorted(idx)]
        sub, B = align_query_to_reference(query, ref, genes)
        R = sub.as_R_cpm().T
        props = np.vstack([_normalize(nnls(R, B[:, j])[0]) for j in range(B.shape[1])])
        return SolverResult(self._frame(props, query.columns, sub.cell_types),
                            list(sub.gene_names),
                            {"solver": self.name, "n_genes": len(sub.gene_names),
                             "top_n": self.top_n,
                             "condition_number": self.condition_number(R)})
