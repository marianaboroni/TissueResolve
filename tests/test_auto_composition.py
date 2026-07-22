"""Offline tests for the composition-calibrated auto solver (experimental, opt-in)."""
from __future__ import annotations

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
    ref = ReferenceSignature(
        gene_names=genes, cell_types=cts, R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 10.0, dtype=np.float32),
        donor_cv=np.full((G, 4), 0.2, dtype=np.float32),
        n_cells_per_type={c: 100 for c in cts})
    true = rng.dirichlet(np.ones(4), size=6)
    bulk = pd.DataFrame(np.rint((true @ R).T * 5).astype(np.int64), index=genes,
                        columns=[f"s{i}" for i in range(6)])
    truth = pd.DataFrame(true, index=[f"s{i}" for i in range(6)], columns=cts)
    return bulk, ref, truth


def test_selects_valid_candidate_and_records(toy_bulk_ref):
    from tissueresolve.solver import AutoCompositionSolver
    bulk, ref, _ = toy_bulk_ref
    res = AutoCompositionSolver(n_sim=20, seed=0).solve(bulk, ref)
    np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)
    assert res.diagnostics["selected_solver"] in (
        "nnls", "weighted_nnls", "marker_nnls", "ridge_nnls", "poisson")
    assert "sim_composition_comparison" in res.diagnostics
    assert "Poisson" in res.diagnostics["selection_reason"]


def test_poisson_is_a_candidate(toy_bulk_ref):
    """Unlike gene-masking 'auto', the composition auto considers the Poisson GLM."""
    from tissueresolve.solver import AutoCompositionSolver
    bulk, ref, _ = toy_bulk_ref
    _, comp, _ = AutoCompositionSolver(n_sim=20, seed=0).select(bulk, ref)
    assert "poisson" in comp.index


def test_reproducible(toy_bulk_ref):
    from tissueresolve.solver import AutoCompositionSolver
    bulk, ref, _ = toy_bulk_ref
    a = AutoCompositionSolver(n_sim=20, seed=1).select(bulk, ref)[1]
    b = AutoCompositionSolver(n_sim=20, seed=1).select(bulk, ref)[1]
    pd.testing.assert_frame_equal(a, b)


def test_registered_and_dispatch():
    from tissueresolve.solver import SOLVERS, get_solver
    assert "auto_composition" in SOLVERS
    assert get_solver("auto_composition").name == "auto_composition"


def test_recovers_dominant_type(toy_bulk_ref):
    from tissueresolve.solver import AutoCompositionSolver
    bulk, ref, truth = toy_bulk_ref
    res = AutoCompositionSolver(n_sim=30, seed=0).solve(bulk, ref)
    est = res.proportions.reindex(index=truth.index, columns=truth.columns)
    # dominant cell type per sample recovered
    assert (est.idxmax(axis=1) == truth.idxmax(axis=1)).mean() >= 0.8


def test_deconv_bulk_integration(toy_bulk_ref):
    import tissueresolve as tr
    bulk, ref, _ = toy_bulk_ref
    res = tr.deconv_bulk(bulk, ref, solver="auto_composition", resolution_mode="flat")
    np.testing.assert_allclose(res.deconv.proportions.sum(axis=1), 1.0, atol=1e-6)
    assert res.deconv.run_metadata.get("selected_solver") in (
        "nnls", "weighted_nnls", "marker_nnls", "ridge_nnls", "poisson")


def test_handles_missing_donor_cv():
    """No donor_cv → falls back to a default CV, still selects a solver."""
    from tissueresolve.solver import AutoCompositionSolver
    rng = np.random.default_rng(0)
    G = 40; genes = [f"g{i}" for i in range(G)]
    R = np.vstack([np.where(np.arange(G) < 20, 100.0, 2.0),
                   np.where(np.arange(G) >= 20, 100.0, 2.0)]).astype("float32")
    ref = ReferenceSignature(gene_names=genes, cell_types=["A", "B"], R_cpm=R)
    true = rng.dirichlet(np.ones(2), size=4)
    bulk = pd.DataFrame(np.rint((true @ R).T * 5), index=genes,
                        columns=[f"s{i}" for i in range(4)])
    res = AutoCompositionSolver(n_sim=15, seed=0).solve(bulk, ref)
    np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)
