"""Plain non-negative least squares on all (or given) shared genes."""
from __future__ import annotations

import numpy as np

from tissueresolve.solver.base import BaseSolver, SolverResult, align_query_to_reference, _normalize


class NNLSSolver(BaseSolver):
    name = "nnls"

    def __init__(self, genes=None):
        self.genes = genes

    def solve(self, query, ref) -> SolverResult:
        from scipy.optimize import nnls
        sub, B = align_query_to_reference(query, ref, self.genes)
        R = sub.as_R_cpm().T  # genes × K
        props = np.vstack([_normalize(nnls(R, B[:, j])[0]) for j in range(B.shape[1])])
        return SolverResult(self._frame(props, query.columns, sub.cell_types),
                            list(sub.gene_names),
                            {"solver": self.name, "n_genes": len(sub.gene_names),
                             "condition_number": self.condition_number(R)})
