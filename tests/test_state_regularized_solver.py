"""Tests for the experimental in-solver state-regularized optimiser (Option A).

Numerical optimiser tests (loss decreases, valid compositions, no NaN, penalty
behaviour, rare protection, convergence metadata) plus config/preset/CLI/pipeline
integration (default off, production unchanged, bulk unchanged). The optimiser is a
SEPARATE experimental path; the production NB-CAR solver is never modified.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.state_regularized_solver import (
    fit_state_regularized_spatial, StateRegularizedSolverResult)
from tissueresolve.experimental.state_similarity_regularization import (
    compute_state_similarity_graph)


CELL_TYPES = ["A1", "A2", "B1", "B2"]
FAMILY = {"A1": "A", "A2": "A", "B1": "B", "B2": "B"}


def _toy(n_spots=40, seed=0):
    """Reference with collinear family A (A1≈A2), distinct family B; counts from truth."""
    rng = np.random.default_rng(seed)
    G = 30
    base = rng.random(G) * 5
    R = np.vstack([base + rng.normal(0, 0.05, G), base + rng.normal(0, 0.05, G),
                   np.r_[np.zeros(15), rng.random(15) * 5],
                   np.r_[rng.random(15) * 5, np.zeros(15)]])
    R = np.clip(R, 0, None)
    R = R / R.sum(1, keepdims=True)
    truth = rng.dirichlet(np.ones(4), size=n_spots)
    lib = rng.uniform(500, 1500, n_spots)
    Y = (truth @ R) * lib[:, None]
    sg = compute_state_similarity_graph(R, CELL_TYPES, FAMILY, within_family_only=True)
    return Y, R, lib, truth, sg


def _eff_n(theta):
    return float(np.mean([1.0 / np.sum(r ** 2) for r in theta]))


# ---- numerical optimiser ----

def test_loss_decreases_on_toy_data():
    Y, R, lib, _, sg = _toy()
    r = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, lambda_spatial=0.0, lambda_state=0.05,
                                      lambda_sparse=0.001, state_penalty="competition")
    assert r.final_loss <= r.loss_trace[0] + 1e-9
    assert r.loss_monotonic


def test_output_nonneg_valid_no_nan():
    Y, R, lib, _, sg = _toy()
    r = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, lambda_state=0.05)
    th = r.theta
    assert (th >= -1e-9).all()
    assert np.allclose(th.sum(1), 1.0, atol=1e-5)
    assert not np.any(np.isnan(th)) and np.all(np.isfinite(th))


def test_lambda_state_zero_matches_no_state_penalty():
    Y, R, lib, _, sg = _toy()
    a = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, lambda_state=0.0, lambda_sparse=0.0,
                                      lambda_spatial=0.0)
    b = fit_state_regularized_spatial(Y, R, state_graph=None, cell_types=CELL_TYPES,
                                      lib_sizes=lib, lambda_state=0.5, lambda_sparse=0.0,
                                      lambda_spatial=0.0)  # no graph → state term inactive
    assert np.allclose(a.theta, b.theta, atol=1e-4)


def test_lambda_sparse_reduces_effective_n():
    Y, R, lib, _, sg = _toy()
    base = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                         lib_sizes=lib, lambda_state=0.0, lambda_sparse=0.0,
                                         lambda_spatial=0.0)
    sparse = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                           lib_sizes=lib, lambda_state=0.0, lambda_sparse=0.05,
                                           lambda_spatial=0.0)
    assert _eff_n(sparse.theta) <= _eff_n(base.theta) + 1e-9


def test_competition_discourages_co_assignment_of_similar_states():
    # competition penalty should reduce simultaneous A1 & A2 mass vs no penalty
    Y, R, lib, _, sg = _toy()
    base = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                         lib_sizes=lib, lambda_state=0.0, lambda_spatial=0.0)
    comp = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                         lib_sizes=lib, lambda_state=0.3, lambda_spatial=0.0,
                                         state_penalty="competition")
    co_base = np.mean(np.minimum(base.theta[:, 0], base.theta[:, 1]))
    co_comp = np.mean(np.minimum(comp.theta[:, 0], comp.theta[:, 1]))
    assert co_comp <= co_base + 1e-9


def test_laplacian_smooths_similar_state_abundances():
    # From an asymmetric A1/A2 split (recon is flat along the collinear direction),
    # the laplacian penalty should pull the A1 and A2 abundances together.
    Y, R, lib, _, sg = _toy()
    rng = np.random.default_rng(1)
    init = rng.dirichlet(np.ones(4), size=Y.shape[0])
    init[:, 0] += 0.3 * init[:, 1]  # push A1 above A2
    init = init / init.sum(1, keepdims=True)
    base = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                         lib_sizes=lib, init_theta=init, lambda_state=0.0,
                                         lambda_spatial=0.0)
    lap = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                        lib_sizes=lib, init_theta=init, lambda_state=0.5,
                                        lambda_spatial=0.0, state_penalty="laplacian")
    gap_base = np.mean(np.abs(base.theta[:, 0] - base.theta[:, 1]))
    gap_lap = np.mean(np.abs(lap.theta[:, 0] - lap.theta[:, 1]))
    assert gap_lap <= gap_base + 1e-9


def test_within_family_only_prevents_cross_family_interaction():
    # cross-family similarity must not be coupled: a graph restricted to families
    # has no A-B edges, so the penalty matrix is block-diagonal by family.
    Y, R, lib, _, _ = _toy()
    sg = compute_state_similarity_graph(R, CELL_TYPES, FAMILY, within_family_only=True)
    for e in range(sg.n_edges):
        a, b = sg.cell_types[sg.state_i[e]], sg.cell_types[sg.state_j[e]]
        assert FAMILY[a] == FAMILY[b]


def test_rare_protection_prevents_zeroing():
    Y, R, lib, _, sg = _toy()
    r = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, lambda_state=0.3, lambda_sparse=0.2,
                                      state_penalty="competition",
                                      rare_protection={"B2": 0.05})
    assert (r.theta[:, CELL_TYPES.index("B2")] >= 0.05 - 1e-6).all()


def test_convergence_metadata_recorded():
    Y, R, lib, _, sg = _toy()
    r = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, lambda_state=0.05)
    assert isinstance(r, StateRegularizedSolverResult)
    for k in ("optimizer", "state_penalty", "n_iter", "converged", "loss_monotonic"):
        assert k in r.metadata
    assert r.n_iter >= 1


def test_preserve_broad_mass_conserves_family_mass():
    Y, R, lib, _, sg = _toy()
    init = np.full((Y.shape[0], 4), 0.25)
    r = fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, init_theta=init, lambda_state=0.3,
                                      lambda_sparse=0.05, preserve_broad_mass=True,
                                      family_map=FAMILY)
    famA = r.theta[:, [0, 1]].sum(1)
    famB = r.theta[:, [2, 3]].sum(1)
    assert np.allclose(famA, init[:, [0, 1]].sum(1), atol=1e-4)
    assert np.allclose(famB, init[:, [2, 3]].sum(1), atol=1e-4)


def test_unknown_optimizer_raises():
    Y, R, lib, _, sg = _toy()
    with pytest.raises(ValueError):
        fit_state_regularized_spatial(Y, R, state_graph=sg, cell_types=CELL_TYPES,
                                      lib_sizes=lib, optimizer="lbfgs")


# ---- config / preset / CLI ----

def test_default_solver_config_disabled():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    assert cfg.state_regularized_solver.enabled is False
    assert cfg.spatial_solver.lambda_spatial == 0.1


def test_solver_presets_enable_and_record_penalty():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_presets import apply_spatial_preset
    expected = {
        "state_regularized_solver_experimental": "competition",
        "state_regularized_solver_competition": "competition",
        "state_regularized_solver_laplacian": "laplacian",
        "state_regularized_solver_weak": "competition",
    }
    for preset, pen in expected.items():
        cfg = TissueResolveConfig()
        info = apply_spatial_preset(cfg, preset)
        assert cfg.state_regularized_solver.enabled is True
        assert cfg.state_regularized_solver.state_penalty == pen
        assert info.experimental and info.state_regularized_solver
        md = info.to_metadata()
        assert md["state_regularized_solver_used"] is True
        assert md["state_solver_penalty"] == pen
        # production spatial lambda untouched
        assert cfg.spatial_solver.lambda_spatial == 0.1


def test_existing_presets_do_not_enable_solver():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_presets import apply_spatial_preset
    for p in ("default", "weak_smoothing", "edge_aware_smoothing",
              "combined_weak_edge_smoothing", "state_regularized_experimental"):
        cfg = TissueResolveConfig()
        apply_spatial_preset(cfg, p)
        assert cfg.state_regularized_solver.enabled is False


def test_unknown_preset_raises():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_presets import apply_spatial_preset
    with pytest.raises(ValueError):
        apply_spatial_preset(TissueResolveConfig(), "no_such_preset")


def test_cli_accepts_solver_presets():
    from tissueresolve.cli import spatial_run
    for p in spatial_run.params:
        if p.name == "spatial_preset":
            for preset in ("state_regularized_solver_experimental",
                           "state_regularized_solver_competition",
                           "state_regularized_solver_laplacian",
                           "state_regularized_solver_weak"):
                assert preset in p.type.choices
            break
    else:
        raise AssertionError("spatial_preset option not found")


# ---- pipeline integration ----

def test_pipeline_default_solver_off_unchanged():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    from tissueresolve.config import TissueResolveConfig
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    res = SpatialPipeline(TissueResolveConfig()).run(Y, ref, r, c, lib, genes, marker_genes=genes)
    assert res.deconv.run_metadata.get("state_regularized_solver_used") is False


def test_pipeline_solver_valid_and_records_diagnostics():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    from tissueresolve.config import TissueResolveConfig
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    cts = list(ref.cell_types)
    family_map = {ct: ("famA" if i % 2 == 0 else "famB") for i, ct in enumerate(cts)}
    cfg = TissueResolveConfig()
    cfg.state_regularized_solver.enabled = True
    cfg.state_regularized_solver.state_penalty = "competition"
    cfg.state_regularized_solver.lambda_state = 0.05
    cfg.state_regularized_solver.max_iter = 100
    res = SpatialPipeline(cfg).run(Y, ref, r, c, lib, genes, marker_genes=genes,
                                   family_map=family_map)
    P = res.deconv.proportions.to_numpy()
    assert not np.any(np.isnan(P)) and (P >= -1e-9).all()
    assert np.allclose(P.sum(1), 1.0, atol=1e-4)
    md = res.deconv.run_metadata
    assert md.get("state_regularized_solver_used") is True
    assert "state_solver_final_loss" in md and "state_solver_n_iter" in md
    assert md.get("state_solver_converged") in (True, False)


def test_bulk_behaviour_unchanged_by_solver_config():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    before = cfg.bulk_solver.__dict__.copy()
    cfg.state_regularized_solver.enabled = True
    assert cfg.bulk_solver.__dict__ == before
