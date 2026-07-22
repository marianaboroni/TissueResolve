"""Offline tests for the Poisson/NB GLM solver backbone (opt-in, non-default).

Covers: solve + simplex, dispersion path (NB), recovery on clean toy data,
registration/dispatch, deconv_bulk(solver=...) integration, mass conservation,
and the explicit NNLS fallback on failure (never a silent wrong answer).
"""
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
    ref = ReferenceSignature(
        gene_names=genes, cell_types=cts, R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 10.0, dtype=np.float32),      # finite dispersion for NB path
        n_cells_per_type={c: 100 for c in cts})
    true = rng.dirichlet(np.ones(4), size=6)
    # integer-ish counts (GLM expects counts): scale then round
    bulk = pd.DataFrame(np.rint((true @ R).T * 5).astype(np.int64), index=genes,
                        columns=[f"s{i}" for i in range(6)])
    truth = pd.DataFrame(true, index=[f"s{i}" for i in range(6)], columns=cts)
    return bulk, ref, truth


def test_poisson_and_nb_run_and_sum_to_one(toy_bulk_ref):
    from tissueresolve.solver import get_solver
    bulk, ref, _ = toy_bulk_ref
    for name in ("poisson", "nb"):
        res = get_solver(name).solve(bulk, ref)
        np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)
        assert (res.proportions.to_numpy() >= -1e-9).all()
        assert res.diagnostics["solver"] == name
        assert res.diagnostics["feature_status"] == "experimental"


def test_registered_in_dispatch():
    from tissueresolve.solver import SOLVERS, get_solver
    assert {"poisson", "nb"} <= set(SOLVERS)
    assert get_solver("poisson").name == "poisson"
    assert get_solver("nb").name == "nb"


def test_recovers_clean_proportions(toy_bulk_ref):
    from tissueresolve.solver import get_solver
    bulk, ref, truth = toy_bulk_ref
    res = get_solver("poisson").solve(bulk, ref)
    est = res.proportions.reindex(index=truth.index, columns=truth.columns)
    rmse = float(np.sqrt(np.mean((est.to_numpy() - truth.to_numpy()) ** 2)))
    assert rmse < 0.05, f"clean-data RMSE too high: {rmse}"


def test_nb_uses_dispersion_when_shape_matches(toy_bulk_ref):
    from tissueresolve.solver import get_solver
    bulk, ref, _ = toy_bulk_ref
    res = get_solver("nb").solve(bulk, ref)
    assert res.diagnostics["loss"] == "nb"
    assert res.diagnostics.get("fallback") is None


def test_explicit_fallback_on_error(monkeypatch, toy_bulk_ref):
    """A solver failure must fall back to NNLS, warn, and record it — never silent."""
    from tissueresolve.solver import get_solver
    import tissueresolve.experimental.nb_bulk_solver as nb
    bulk, ref, _ = toy_bulk_ref

    def _boom(*a, **k):
        raise RuntimeError("forced")
    monkeypatch.setattr(nb, "fit_bulk_nb_glm", _boom)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = get_solver("poisson").solve(bulk, ref)
    np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)
    assert res.diagnostics["fallback"] == "nnls"
    assert "forced" in res.diagnostics["fallback_reason"]
    assert any("fell back to NNLS" in str(x.message) for x in w)


def test_fallback_can_be_disabled(monkeypatch, toy_bulk_ref):
    from tissueresolve.solver import PoissonGLMSolver
    import tissueresolve.experimental.nb_bulk_solver as nb
    bulk, ref, _ = toy_bulk_ref
    monkeypatch.setattr(nb, "fit_bulk_nb_glm",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
    with pytest.raises(RuntimeError):
        PoissonGLMSolver(fallback=False).solve(bulk, ref)


def test_deconv_bulk_solver_poisson_integration(toy_bulk_ref):
    """deconv_bulk(solver='poisson') routes through the additive solver path."""
    import tissueresolve as tr
    bulk, ref, _ = toy_bulk_ref
    res = tr.deconv_bulk(bulk, ref, solver="poisson", resolution_mode="flat")
    props = res.deconv.proportions
    np.testing.assert_allclose(props.sum(axis=1), 1.0, atol=1e-6)
    assert res.deconv.run_metadata.get("solver") == "poisson"


def test_deconv_bulk_poisson_hierarchical(toy_bulk_ref):
    """solver='poisson' composes with hierarchical: the GLM drives BOTH the
    family and fine solves (via cfg.bulk_solver.method), not just a flat estimate."""
    import tissueresolve as tr
    bulk, ref, _ = toy_bulk_ref
    mapping = {"A": "fam1", "B": "fam1", "C": "fam2", "D": "fam2"}
    res = tr.deconv_bulk(bulk, ref, solver="poisson", resolution_mode="hierarchical",
                         hierarchy_mapping=mapping, n_bootstrap=0)
    props = res.deconv.proportions
    np.testing.assert_allclose(props.sum(axis=1), 1.0, atol=1e-6)      # mass conserved
    assert res.deconv.run_metadata.get("resolution_mode") == "hierarchical"
    assert res.deconv.run_metadata.get("solver") == "poisson_glm_experimental"


def test_deconv_bulk_poisson_hierarchical_does_not_change_wnnls_default(toy_bulk_ref):
    """The default (no solver) hierarchical run stays wNNLS — Poisson is opt-in."""
    import tissueresolve as tr
    bulk, ref, _ = toy_bulk_ref
    mapping = {"A": "fam1", "B": "fam1", "C": "fam2", "D": "fam2"}
    res = tr.deconv_bulk(bulk, ref, resolution_mode="hierarchical",
                         hierarchy_mapping=mapping, n_bootstrap=0)
    assert res.deconv.run_metadata.get("solver") != "poisson_glm_experimental"
