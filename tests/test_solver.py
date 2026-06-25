"""Offline tests for solver backbones, gene-masking CV, ensemble, and auto."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature


@pytest.fixture
def toy_bulk_ref():
    rng = np.random.default_rng(0)
    G = 80
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0); v[list(a)] = lvl; return v
    R = np.vstack([prof(range(0, 20)), prof(range(20, 40)),
                   prof(range(40, 60)), prof(range(60, 80))]).astype(np.float32)
    cts = ["A", "B", "C", "D"]
    ref = ReferenceSignature(gene_names=genes, cell_types=cts, R_cpm=R,
                             R_log=np.log1p(R).astype(np.float32),
                             n_cells_per_type={c: 100 for c in cts})
    true = rng.dirichlet(np.ones(4), size=6)
    bulk = pd.DataFrame((true @ R).T, index=genes,
                        columns=[f"s{i}" for i in range(6)])
    truth = pd.DataFrame(true, index=[f"s{i}" for i in range(6)], columns=cts)
    return bulk, ref, truth


def test_each_solver_runs_and_sums_to_one(toy_bulk_ref):
    from tissueresolve.solver import get_solver
    bulk, ref, _ = toy_bulk_ref
    for name in ("nnls", "weighted_nnls", "marker_nnls", "ridge_nnls"):
        res = get_solver(name).solve(bulk, ref)
        np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)
        assert res.genes_used


def test_auto_selects_and_records_reason(toy_bulk_ref):
    from tissueresolve.solver import AutoSolver
    bulk, ref, _ = toy_bulk_ref
    auto = AutoSolver(n_splits=2)
    res = auto.solve(bulk, ref)
    assert res.diagnostics["selected_solver"] in (
        "nnls", "weighted_nnls", "marker_nnls", "ridge_nnls")
    assert "selection_reason" in res.diagnostics
    assert "solver_comparison" in res.diagnostics


def test_auto_not_chosen_by_in_sample_reconstruction_alone(toy_bulk_ref):
    """The auto objective is CV reconstruction minus a conditioning penalty."""
    from tissueresolve.solver import AutoSolver
    bulk, ref, _ = toy_bulk_ref
    best, comp, reason = AutoSolver(n_splits=2).select(bulk, ref)
    assert "objective" in comp.columns
    assert "cv_masked_gene_pearson" in comp.columns
    assert "condition_number" in comp.columns


def test_gene_masking_reproducible(toy_bulk_ref):
    from tissueresolve.validation.gene_masking import split_genes_for_masking
    bulk, ref, _ = toy_bulk_ref
    genes = list(bulk.index)
    a = split_genes_for_masking(genes, n_splits=3, seed=0)
    b = split_genes_for_masking(genes, n_splits=3, seed=0)
    assert [m for _, m in a] == [m for _, m in b]


def test_gene_masking_scores_reconstruction(toy_bulk_ref):
    from tissueresolve.validation.gene_masking import run_gene_masking_cv_bulk
    from tissueresolve.solver import NNLSSolver
    bulk, ref, _ = toy_bulk_ref

    def solve_fn(b, r, genes):
        return NNLSSolver(genes=list(genes)).solve(b, r).proportions
    res = run_gene_masking_cv_bulk(bulk, ref, solve_fn, n_splits=2)
    assert -1.0 <= res.score <= 1.0
    assert res.per_fold


def test_ensemble_weights_sum_to_one_and_mass_preserved(toy_bulk_ref):
    from tissueresolve.solver import EnsembleSolver
    bulk, ref, _ = toy_bulk_ref
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = EnsembleSolver(n_splits=2).solve(bulk, ref)
    w = res.diagnostics["weights"]
    assert abs(sum(w.values()) - 1.0) < 1e-6
    np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)


def test_deconv_bulk_solver_path(toy_bulk_ref):
    import tissueresolve as tr
    bulk, ref, _ = toy_bulk_ref
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = tr.deconv_bulk(bulk, ref, solver="nnls")
    assert res.deconv.run_metadata["solver"] == "nnls"
    np.testing.assert_allclose(res.deconv.proportions.sum(axis=1), 1.0, atol=1e-6)
    assert res.deconv.coverage_r2 is not None
