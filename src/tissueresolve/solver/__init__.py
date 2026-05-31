"""Reference-based deconvolution solver backbones and auto-selection."""
from tissueresolve.solver.base import BaseSolver, SolverResult, align_query_to_reference
from tissueresolve.solver.nnls import NNLSSolver
from tissueresolve.solver.weighted_nnls import WeightedNNLSSolver
from tissueresolve.solver.marker_nnls import MarkerNNLSSolver
from tissueresolve.solver.ridge_nnls import RidgeNNLSSolver
from tissueresolve.solver.ensemble import EnsembleSolver
from tissueresolve.solver.auto import AutoSolver, candidate_solvers

SOLVERS = {
    "nnls": NNLSSolver, "weighted_nnls": WeightedNNLSSolver,
    "marker_nnls": MarkerNNLSSolver, "ridge_nnls": RidgeNNLSSolver,
    "ensemble_nnls": EnsembleSolver, "auto": AutoSolver,
}

__all__ = ["BaseSolver", "SolverResult", "align_query_to_reference",
           "NNLSSolver", "WeightedNNLSSolver", "MarkerNNLSSolver",
           "RidgeNNLSSolver", "EnsembleSolver", "AutoSolver",
           "candidate_solvers", "SOLVERS", "get_solver"]


def get_solver(name: str, **kwargs) -> BaseSolver:
    if name not in SOLVERS:
        raise ValueError(f"unknown solver {name!r}; choose from {list(SOLVERS)}")
    return SOLVERS[name](**kwargs)
