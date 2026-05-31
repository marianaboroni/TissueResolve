"""Ridge-regularised non-negative deconvolution (closed form + clip)."""
from __future__ import annotations

import numpy as np

from tissueresolve.solver.base import BaseSolver, SolverResult, align_query_to_reference, _normalize


class RidgeNNLSSolver(BaseSolver):
    name = "ridge_nnls"

    def __init__(self, alpha: float = 0.1, genes=None):
        self.alpha = alpha
        self.genes = genes

    def solve(self, query, ref) -> SolverResult:
        sub, B = align_query_to_reference(query, ref, self.genes)
        R = sub.as_R_cpm().T
        K = R.shape[1]
        RtR = R.T @ R
        lam = self.alpha * (np.trace(RtR) / max(K, 1))
        A = np.linalg.pinv(RtR + lam * np.eye(K)) @ R.T
        props = np.vstack([_normalize(np.clip(A @ B[:, j], 0, None))
                           for j in range(B.shape[1])])
        return SolverResult(self._frame(props, query.columns, sub.cell_types),
                            list(sub.gene_names),
                            {"solver": self.name, "n_genes": len(sub.gene_names),
                             "alpha": self.alpha,
                             "condition_number": self.condition_number(R)})
