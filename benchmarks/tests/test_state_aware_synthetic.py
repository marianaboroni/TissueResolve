"""Offline tests for the synthetic three-level state benchmark."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("anndata")

from benchmarks.synthetic.state_hierarchy import (
    SCENARIOS, simulate_scenario, build_state_reference, three_level_hierarchy)
from benchmarks.synthetic import run_state_aware_synthetic as RUN

_SMALL = dict(n_cells_per_state=30, n_samples=6)


def test_simulate_three_levels_and_truth():
    sc = simulate_scenario("separable", **_SMALL)
    # truth sums to 1 per sample at each level
    for t in (sc.truth_state, sc.truth_celltype, sc.truth_broad):
        assert np.allclose(t.sum(axis=1).to_numpy(), 1.0, atol=1e-6)
    # state→cell_type→broad column counts
    assert sc.truth_state.shape[1] == 7
    assert sc.truth_celltype.shape[1] == 4
    assert sc.truth_broad.shape[1] == 2


def test_missing_state_dropped_from_reference_only():
    sc = simulate_scenario("missing_state", **_SMALL)
    assert "D2" in sc.missing_states
    assert "D2" not in sc.state_cell_types          # absent from reference
    assert "D2" in sc.truth_state.columns           # still in the bulk truth


def test_scenarios_differ_in_state_separability():
    sep = simulate_scenario("separable", **_SMALL)
    non = simulate_scenario("non_separable", **_SMALL)
    assert sep.resolvable_states and not non.resolvable_states


def test_run_methods_returns_all_five():
    sc = simulate_scenario("separable", **_SMALL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = RUN.run_methods(sc)
    assert set(res) == {
        "NNLS_baseline", "TissueResolve_auto",
        "TissueResolve_hierarchical_standard",
        "TissueResolve_hierarchical_granular", "TissueResolve_state_aware"}


def test_aggregate_levels_routes_unresolved_mass():
    sc = simulate_scenario("separable", **_SMALL)
    # a synthetic prediction: half resolved to A1, half unresolved at broad X
    pred = pd.DataFrame({"A1": [0.5], "unresolved_X": [0.5]}, index=["s0"])
    bp, cp, sp = RUN.aggregate_levels(pred, sc)
    assert bp.loc["s0", "X"] == pytest.approx(1.0)   # both map to broad X
    assert sp.loc["s0", "A1"] == pytest.approx(0.5)
    assert cp.loc["s0", "A"] == pytest.approx(0.5)   # only resolved A1 → cell type A


def test_metrics_columns_present():
    sc = simulate_scenario("separable", **_SMALL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = RUN.compute_metrics(sc, RUN.run_methods(sc))
    for col in ("broad_pearson", "celltype_pearson", "state_pearson",
                "resolvable_state_pearson", "false_resolution_rate",
                "unresolved_mass_fraction", "unresolved_precision",
                "unresolved_recall", "runtime_seconds"):
        assert col in m.columns


def test_state_aware_avoids_false_resolution_on_non_separable():
    """The headline honesty check: on non-separable states, state-aware leaves
    mass UNRESOLVED (low false-resolution) while flat NNLS resolves everything."""
    sc = simulate_scenario("non_separable", **_SMALL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = RUN.compute_metrics(sc, RUN.run_methods(sc)).set_index("method")
    sa_fr = m.loc["TissueResolve_state_aware", "false_resolution_rate"]
    nn_fr = m.loc["NNLS_baseline", "false_resolution_rate"]
    assert sa_fr < nn_fr            # state-aware abstains, NNLS over-resolves
    # state-aware keeps the non-separable mass unresolved
    assert m.loc["TissueResolve_state_aware", "unresolved_mass_fraction"] > 0.3
    assert m.loc["NNLS_baseline", "unresolved_mass_fraction"] == pytest.approx(0.0)


def test_broad_accuracy_preserved_by_state_aware():
    """State-aware must not harm broad-level accuracy vs the flat baseline."""
    sc = simulate_scenario("separable", **_SMALL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = RUN.compute_metrics(sc, RUN.run_methods(sc)).set_index("method")
    assert m.loc["TissueResolve_state_aware", "broad_pearson"] >= \
        m.loc["NNLS_baseline", "broad_pearson"] - 0.2


def test_write_outputs(tmp_path):
    metrics = pd.concat(
        [RUN.compute_metrics(simulate_scenario(n, **_SMALL),
                             RUN.run_methods(simulate_scenario(n, **_SMALL)))
         for n in ("separable", "non_separable")], ignore_index=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        paths = RUN.write_outputs(metrics, tmp_path)
    assert paths["benchmark"].exists() and paths["family_summary"].exists()
    assert paths["report"].exists()
    html = paths["report"].read_text()
    assert "synthetic" in html.lower() and "Verdict" in html
