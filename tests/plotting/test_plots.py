"""
Tests for tissueresolve.plotting (Stage 5).

Verifies that every plotting function:
* returns a PlotResult carrying the figure object,
* writes figure files (png/pdf/svg) that exist on disk,
* writes the underlying source data as readable TSV,
* (where relevant) produces a caption naming the estimate type.

Guarded by importorskip so the suite still runs without matplotlib.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

mpl = pytest.importorskip("matplotlib")

from tissueresolve.plotting import (  # noqa: E402
    benchmark_plots,
    bulk_plots,
    qc_plots,
    spatial_plots,
)
from tissueresolve.plotting.style import PlotResult  # noqa: E402

N_SPOTS = 20


def _assert_plot_result(pr: PlotResult):
    assert isinstance(pr, PlotResult)
    assert pr.figure is not None
    assert pr.figure_paths, "no figure files saved"
    for fmt, path in pr.figure_paths.items():
        assert path.exists(), f"missing {fmt} figure {path}"
    assert pr.data_paths, "no source data saved (DESIGN_SPEC violation)"
    for path in pr.data_paths.values():
        assert path.exists(), f"missing data file {path}"
        # Readable (skip comment lines)
        pd.read_csv(path, sep="\t", comment="#")


def _coords(n=N_SPOTS):
    r = np.repeat(np.arange((n + 4) // 5), 5)[:n]
    c = np.tile(np.arange(5) * 2, (n + 4) // 5)[:n]
    return r.astype(int), c.astype(int)


# ---------------------------------------------------------------------------
# Bulk
# ---------------------------------------------------------------------------


class TestBulkPlots:
    def test_composition_barplot(self, bulk_result, tmp_path):
        pr = bulk_plots.composition_barplot(bulk_result, tmp_path)
        _assert_plot_result(pr)
        assert "mRNA proportion" in pr.caption

    def test_proportion_heatmap(self, bulk_result, tmp_path):
        _assert_plot_result(bulk_plots.proportion_heatmap(bulk_result, tmp_path))

    def test_qc_barplot(self, bulk_qc_report, tmp_path):
        _assert_plot_result(bulk_plots.qc_barplot(bulk_qc_report, tmp_path, r2_warn=0.5))

    def test_marker_recall_plot(self, bulk_qc_report, tmp_path):
        _assert_plot_result(bulk_plots.marker_recall_plot(bulk_qc_report, tmp_path))

    def test_bootstrap_ci_plot(self, bulk_result_with_ci, tmp_path):
        pr = bulk_plots.bootstrap_ci_plot(bulk_result_with_ci, tmp_path)
        _assert_plot_result(pr)

    def test_bootstrap_ci_plot_requires_ci(self, bulk_result, tmp_path):
        with pytest.raises(ValueError, match="bootstrap"):
            bulk_plots.bootstrap_ci_plot(bulk_result, tmp_path)

    def test_three_formats_saved(self, bulk_result, tmp_path):
        pr = bulk_plots.composition_barplot(bulk_result, tmp_path)
        assert set(pr.figure_paths) == {"png", "pdf", "svg"}


# ---------------------------------------------------------------------------
# Spatial
# ---------------------------------------------------------------------------


class TestSpatialPlots:
    def test_abundance_map(self, spatial_result, tmp_path):
        r, c = _coords()
        ct = spatial_result.proportions.columns[0]
        pr = spatial_plots.abundance_map(spatial_result, r, c, ct, tmp_path)
        _assert_plot_result(pr)
        # coordinates recorded in source data
        df = pd.read_csv(next(iter(pr.data_paths.values())), sep="\t", comment="#")
        assert {"array_row", "array_col"}.issubset(df.columns)

    def test_abundance_map_multi(self, spatial_result, tmp_path):
        r, c = _coords()
        cts = list(spatial_result.proportions.columns[:2])
        _assert_plot_result(
            spatial_plots.abundance_map_multi(spatial_result, r, c, cts, tmp_path)
        )

    def test_dominant_type_map(self, spatial_result, tmp_path):
        r, c = _coords()
        pr = spatial_plots.dominant_type_map(spatial_result, r, c, tmp_path)
        _assert_plot_result(pr)
        assert "single-cell counts" in pr.caption

    def test_spatial_qc_map(self, spatial_qc_report, tmp_path):
        r, c = _coords()
        pr = spatial_plots.spatial_qc_map(
            spatial_qc_report.spot_qc, r, c, "nb_loglik", tmp_path, lambda_spatial=0.1
        )
        _assert_plot_result(pr)

    def test_morans_i_barplot(self, spatial_qc_report, tmp_path):
        _assert_plot_result(
            spatial_plots.morans_i_barplot(spatial_qc_report.morans_i, tmp_path)
        )

    def test_abundance_map_unknown_type_errors(self, spatial_result, tmp_path):
        r, c = _coords()
        with pytest.raises(ValueError):
            spatial_plots.abundance_map(spatial_result, r, c, "NoSuchType", tmp_path)


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------


class TestQCPlots:
    def test_separability_heatmap(self, separability_report, tmp_path):
        cts = ["CellTypeA", "CellTypeB", "CellTypeC", "CellTypeD"]
        pr = qc_plots.separability_heatmap(separability_report, cts, tmp_path)
        _assert_plot_result(pr)

    def test_gene_overlap_plot(self, tmp_path):
        q = [f"G{i}" for i in range(20)]
        r = [f"G{i}" for i in range(10, 30)]
        _assert_plot_result(qc_plots.gene_overlap_plot(q, r, tmp_path))

    def test_warning_summary_plot(self, tmp_path):
        _assert_plot_result(
            qc_plots.warning_summary_plot(["w1", "w2 long warning"], tmp_path)
        )

    def test_warning_summary_plot_empty(self, tmp_path):
        _assert_plot_result(qc_plots.warning_summary_plot([], tmp_path))


# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------


class TestBenchmarkPlots:
    def test_estimated_vs_true(self, tmp_path):
        rng = np.random.default_rng(0)
        true = rng.dirichlet(np.ones(3), size=30)
        pred = np.clip(true + rng.normal(0, 0.05, true.shape), 0, None)
        pred /= pred.sum(1, keepdims=True)
        pr = benchmark_plots.estimated_vs_true_scatter(
            pred, true, ["A", "B", "C"], tmp_path
        )
        _assert_plot_result(pr)

    def test_per_celltype_rmse(self, benchmark_result, tmp_path):
        _assert_plot_result(
            benchmark_plots.per_celltype_rmse_barplot(
                benchmark_result.per_celltype_metrics, tmp_path
            )
        )

    def test_method_comparison(self, benchmark_result, tmp_path):
        import dataclasses

        other = dataclasses.replace(benchmark_result, method="NNLS")
        _assert_plot_result(
            benchmark_plots.method_comparison_plot([benchmark_result, other], tmp_path)
        )
