"""
Tests for tissueresolve.spatial.benchmark — synthetic Visium benchmark.

These are fast, deterministic, offline unit/integration tests using small
synthetic sections.  They check that:

* simulation produces self-consistent, well-shaped data;
* dataset transforms (mismatch, shuffle) preserve invariants;
* every baseline returns valid simplex proportions;
* metrics behave sensibly (perfect prediction → 0 RMSE / 0 JSD / PCC≈1);
* the full suite recovers structure better than the shuffle negative control;
* the suite saves its underlying data and returns unified BenchmarkResults.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import BenchmarkResult
from tissueresolve.spatial.benchmark import (
    METHODS,
    SCENARIOS,
    BenchmarkDataset,
    apply_mismatch,
    benchmark_results_to_frame,
    compute_benchmark_metrics,
    run_baseline_dwls,
    run_baseline_nnls,
    run_baseline_spatial_nnls,
    run_benchmark_suite,
    run_spatcar,
    shuffle_coordinates,
    simulate_visium,
)

SMALL = dict(n_spots=40, n_types=3, n_genes=30)


@pytest.fixture
def dataset() -> BenchmarkDataset:
    return simulate_visium(**SMALL, seed=0)


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


class TestSimulateVisium:
    def test_shapes_consistent(self, dataset: BenchmarkDataset):
        N, G, K = dataset.n_spots, len(dataset.gene_names), dataset.n_types
        assert dataset.Y_sparse.shape == (N, G)
        assert dataset.Pi_true.shape == (N, K)
        assert dataset.R_cpm.shape == (K, G)
        assert dataset.phi_g.shape == (G,)
        assert dataset.lib_sizes.shape == (N,)

    def test_proportions_on_simplex(self, dataset: BenchmarkDataset):
        assert np.all(dataset.Pi_true >= 0)
        np.testing.assert_allclose(dataset.Pi_true.sum(1), 1.0, atol=1e-5)

    def test_counts_non_negative_integers(self, dataset: BenchmarkDataset):
        Y = dataset.dense_counts()
        assert np.all(Y >= 0)
        np.testing.assert_array_equal(Y, np.round(Y))

    def test_deterministic_given_seed(self):
        a = simulate_visium(**SMALL, seed=3)
        b = simulate_visium(**SMALL, seed=3)
        np.testing.assert_array_equal(a.dense_counts(), b.dense_counts())
        np.testing.assert_array_equal(a.Pi_true, b.Pi_true)

    def test_different_seeds_differ(self):
        a = simulate_visium(**SMALL, seed=1)
        b = simulate_visium(**SMALL, seed=2)
        assert not np.array_equal(a.Pi_true, b.Pi_true)

    def test_to_reference_roundtrip(self, dataset: BenchmarkDataset):
        ref = dataset.to_reference()
        ref.validate()
        assert ref.n_cell_types == dataset.n_types
        assert ref.n_genes == len(dataset.gene_names)
        assert ref.phi_g is not None

    def test_to_graph_matches_spots(self, dataset: BenchmarkDataset):
        graph = dataset.to_graph()
        assert graph.n_spots == dataset.n_spots


# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------


class TestTransforms:
    def test_mismatch_changes_reference_not_counts(self, dataset: BenchmarkDataset):
        mis, true_d = apply_mismatch(dataset, seed=5)
        assert mis.scenario == "mismatch"
        # Counts shared (same underlying matrix).
        np.testing.assert_array_equal(mis.dense_counts(), dataset.dense_counts())
        # Reference altered.
        assert not np.allclose(mis.R_cpm, dataset.R_cpm)
        assert true_d.shape == (len(dataset.gene_names),)
        assert np.any(true_d != 1.0)

    def test_shuffle_preserves_counts_and_truth(self, dataset: BenchmarkDataset):
        sh = shuffle_coordinates(dataset, seed=2)
        assert sh.scenario == "shuffle"
        np.testing.assert_array_equal(sh.dense_counts(), dataset.dense_counts())
        np.testing.assert_array_equal(sh.Pi_true, dataset.Pi_true)
        # Coordinate multiset preserved, ordering changed.
        assert sorted(sh.array_row.tolist()) == sorted(dataset.array_row.tolist())
        assert not np.array_equal(sh.array_row, dataset.array_row) or \
               not np.array_equal(sh.array_col, dataset.array_col)


# ---------------------------------------------------------------------------
# Baseline methods
# ---------------------------------------------------------------------------


def _assert_valid_proportions(Pi: np.ndarray, N: int, K: int):
    assert Pi.shape == (N, K)
    assert np.all(Pi >= -1e-6)
    np.testing.assert_allclose(Pi.sum(1), 1.0, atol=1e-4)


class TestBaselines:
    def test_nnls(self, dataset: BenchmarkDataset):
        Pi = run_baseline_nnls(dataset)
        _assert_valid_proportions(Pi, dataset.n_spots, dataset.n_types)

    def test_dwls(self, dataset: BenchmarkDataset):
        Pi = run_baseline_dwls(dataset, n_iter=5)
        _assert_valid_proportions(Pi, dataset.n_spots, dataset.n_types)

    def test_spatial_nnls(self, dataset: BenchmarkDataset):
        graph = dataset.to_graph()
        Pi = run_baseline_spatial_nnls(dataset, graph, n_smooth=2)
        _assert_valid_proportions(Pi, dataset.n_spots, dataset.n_types)

    def test_spatcar(self, dataset: BenchmarkDataset):
        graph = dataset.to_graph()
        Pi = run_spatcar(dataset, graph, max_iter=20)
        _assert_valid_proportions(Pi, dataset.n_spots, dataset.n_types)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


class TestMetrics:
    def test_perfect_prediction(self, dataset: BenchmarkDataset):
        m = compute_benchmark_metrics(dataset.Pi_true, dataset.Pi_true)
        assert m["rmse"] == pytest.approx(0.0, abs=1e-6)
        assert m["jsd"] == pytest.approx(0.0, abs=1e-6)
        assert m["pcc"] == pytest.approx(1.0, abs=1e-6)

    def test_worse_prediction_has_higher_rmse(self, dataset: BenchmarkDataset):
        rng = np.random.default_rng(0)
        noisy = dataset.Pi_true + rng.normal(0, 0.2, dataset.Pi_true.shape)
        noisy = np.clip(noisy, 1e-6, None)
        noisy /= noisy.sum(1, keepdims=True)
        good = compute_benchmark_metrics(dataset.Pi_true, dataset.Pi_true)
        bad = compute_benchmark_metrics(noisy, dataset.Pi_true)
        assert bad["rmse"] > good["rmse"]

    def test_metric_keys(self, dataset: BenchmarkDataset):
        m = compute_benchmark_metrics(dataset.Pi_true, dataset.Pi_true)
        assert set(m) == {"rmse", "jsd", "pcc", "pcc_std"}


# ---------------------------------------------------------------------------
# Full suite
# ---------------------------------------------------------------------------


class TestRunBenchmarkSuite:
    def test_returns_benchmark_results(self):
        results = run_benchmark_suite(
            **SMALL, scenarios=("basic",), spatcar_max_iter=20, seed=0
        )
        assert len(results) == len(METHODS)
        assert all(isinstance(r, BenchmarkResult) for r in results)
        for r in results:
            assert r.scenario == "basic"
            assert "rmse" in r.metrics
            assert "morans_i_mean" in r.metrics
            assert r.per_celltype_metrics is not None

    def test_all_scenarios_and_methods(self):
        results = run_benchmark_suite(
            **SMALL, scenarios="all", spatcar_max_iter=15, seed=0
        )
        assert len(results) == len(SCENARIOS) * len(METHODS)
        got = {(r.scenario, r.method) for r in results}
        assert got == {(s, m) for s in SCENARIOS for m in METHODS}

    def test_saves_underlying_data(self, tmp_path):
        run_benchmark_suite(
            **SMALL, scenarios=("basic",), output_dir=tmp_path,
            spatcar_max_iter=15, seed=0,
        )
        csv = tmp_path / "benchmark_results.csv"
        assert csv.exists()
        frame = pd.read_csv(csv)
        assert len(frame) == len(METHODS)
        assert {"scenario", "method", "rmse"}.issubset(frame.columns)

    def test_to_frame(self):
        results = run_benchmark_suite(
            **SMALL, scenarios=("basic",), spatcar_max_iter=15, seed=0
        )
        frame = benchmark_results_to_frame(results)
        assert len(frame) == len(METHODS)
        assert set(frame["method"]) == set(METHODS)

    def test_shuffle_negative_control_drops_moran(self):
        # On a spatially-structured section, the fitted Moran's I should track
        # the structure; once coordinates are shuffled, the *true* spatial
        # autocorrelation collapses — so true Moran's I must be far lower in
        # the shuffle scenario than in basic.  This guards the negative control.
        results = run_benchmark_suite(
            n_spots=120, n_types=3, n_genes=40,
            scenarios=("basic", "shuffle"), methods=("NNLS",),
            spatcar_max_iter=15, seed=0,
        )
        by = {r.scenario: r.metrics for r in results}
        assert by["basic"]["morans_i_true_mean"] > by["shuffle"]["morans_i_true_mean"]

    def test_rejects_unknown_scenario(self):
        with pytest.raises(ValueError, match="scenario"):
            run_benchmark_suite(**SMALL, scenarios=("nonsense",))

    def test_rejects_unknown_method(self):
        with pytest.raises(ValueError, match="method"):
            run_benchmark_suite(**SMALL, methods=("FancyNet",))
