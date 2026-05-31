"""
Benchmark utilities for TissueResolve.

``spillover``
    Spillover (cross-type leakage) simulation and estimation, shared by the
    bulk and spatial resolution-aware reporting layer.
"""
from tissueresolve.benchmark.spillover import (
    compute_spillover_risk,
    estimate_spillover_matrix,
    expression_spillover_proxy,
    simulate_pairwise_mixtures,
    simulate_pure_profiles,
    simulate_realistic_mixtures,
    summarize_spillover_partners,
)

__all__ = [
    "simulate_pure_profiles",
    "simulate_pairwise_mixtures",
    "simulate_realistic_mixtures",
    "estimate_spillover_matrix",
    "expression_spillover_proxy",
    "compute_spillover_risk",
    "summarize_spillover_partners",
]
