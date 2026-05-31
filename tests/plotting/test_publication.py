"""
Tests for the Plotly publication layer (figures + export contract).

Offline; skipped entirely if plotly is unavailable.  kaleido is intentionally
absent in CI, so these also exercise the "static export skipped + warning" path.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting import (  # noqa: E402
    bulk_plots, separability_plots, spatial_plots, spillover_plots,
)
from tissueresolve.plotting.export import FigureResult, kaleido_available  # noqa: E402


def _bulk_props(n=6, k=4, seed=0):
    rng = np.random.default_rng(seed)
    raw = rng.dirichlet(np.ones(k), size=n)
    return pd.DataFrame(raw, index=[f"s{i}" for i in range(n)],
                        columns=[f"CellType_{j}" for j in range(k)])


def _spatial_props(n=12, k=4, seed=1):
    rng = np.random.default_rng(seed)
    raw = rng.dirichlet(np.ones(k), size=n)
    return pd.DataFrame(raw, index=[f"spot{i}" for i in range(n)],
                        columns=[f"CT{j}" for j in range(k)])


def _coords(n=12):
    r = np.repeat(np.arange((n + 3) // 4), 4)[:n]
    c = np.tile(np.arange(4) * 2, (n + 3) // 4)[:n]
    return r, c


def _assert_figure_result(res: FigureResult):
    assert isinstance(res, FigureResult)
    assert res.figure is not None
    assert res.html_path is not None and res.html_path.exists()
    assert res.data_paths, "figure produced no source data (contract violation)"
    for p in res.data_paths.values():
        assert p.exists()
        pd.read_csv(p, sep="\t", comment="#")


# ---------------------------------------------------------------------------
# Bulk
# ---------------------------------------------------------------------------


def test_bulk_clustered_barplot_html_and_data(tmp_path):
    res = bulk_plots.plot_bulk_composition_clustered_barplot(_bulk_props(), tmp_path)
    _assert_figure_result(res)
    assert "data" in res.data_paths
    assert "sample_order" in res.data_paths
    assert "cell_type_order" in res.data_paths
    order = pd.read_csv(res.data_paths["sample_order"], sep="\t", comment="#")
    assert len(order) == 6  # non-empty sample order


def test_bulk_clustered_barplot_top_n_other(tmp_path):
    res = bulk_plots.plot_bulk_composition_clustered_barplot(
        _bulk_props(k=6), tmp_path, top_n=3)
    data = pd.read_csv(res.data_paths["data"], sep="\t", comment="#", index_col=0)
    assert "Other" in data.columns
    assert data.shape[1] == 4  # 3 kept + Other


def test_bulk_composition_heatmap_saves_data(tmp_path):
    res = bulk_plots.plot_bulk_composition_heatmap(_bulk_props(), tmp_path)
    _assert_figure_result(res)


def test_bulk_qc_summary_saves_data(tmp_path):
    qc = pd.DataFrame({"recon_r2": [0.9, 0.8, 0.7],
                       "profile_corr": [0.95, 0.85, 0.75]},
                      index=["s0", "s1", "s2"])
    res = bulk_plots.plot_bulk_qc_summary(qc, tmp_path, r2_warn=0.5)
    _assert_figure_result(res)


def test_bulk_uncertainty_without_ci_warns(tmp_path):
    res = bulk_plots.plot_bulk_uncertainty(_bulk_props(), None, None, tmp_path)
    _assert_figure_result(res)
    assert any("CI" in w or "uncertainty" in w for w in res.warnings)


def test_kaleido_missing_records_warning(tmp_path):
    # kaleido is not installed in CI → static export skipped, HTML + data kept.
    res = bulk_plots.plot_bulk_composition_heatmap(_bulk_props(), tmp_path)
    if not kaleido_available():
        assert res.html_path.exists()
        assert res.data_paths
        assert any("kaleido" in w.lower() for w in res.warnings)
        assert not res.static_paths


# ---------------------------------------------------------------------------
# Spatial
# ---------------------------------------------------------------------------


def test_spatial_mean_composition_saves_data(tmp_path):
    res = spatial_plots.plot_spatial_mean_composition_barplot(_spatial_props(), tmp_path)
    _assert_figure_result(res)


def test_spatial_pie_charts_small_coords_and_other(tmp_path):
    props = _spatial_props(k=6)
    r, c = _coords(len(props))
    res = spatial_plots.plot_spatial_spot_pie_charts(props, r, c, tmp_path, top_n=3)
    _assert_figure_result(res)
    data = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "Other" in set(data["cell_type"])  # low-abundance collapsed
    assert {"array_row", "array_col", "fraction"}.issubset(data.columns)


def test_spatial_abundance_map_saves_coords(tmp_path):
    props = _spatial_props()
    r, c = _coords(len(props))
    res = spatial_plots.plot_spatial_abundance_maps(props, r, c, tmp_path)
    _assert_figure_result(res)
    data = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert {"array_row", "array_col"}.issubset(data.columns)


def test_spatial_dominant_map_saves_data(tmp_path):
    props = _spatial_props()
    r, c = _coords(len(props))
    res = spatial_plots.plot_spatial_dominant_cell_type_map(props, r, c, tmp_path)
    _assert_figure_result(res)
    data = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "dominant_type" in data.columns and "dominant_fraction" in data.columns


def test_spatial_qc_maps_saves_data(tmp_path):
    n = 12
    r, c = _coords(n)
    qc = pd.DataFrame({"prop_entropy": np.linspace(0, 2, n),
                       "nb_loglik": np.linspace(-100, -10, n)},
                      index=[f"spot{i}" for i in range(n)])
    res = spatial_plots.plot_spatial_qc_maps(qc, r, c, tmp_path)
    _assert_figure_result(res)


def test_spatial_morans_barplot(tmp_path):
    s = pd.Series([0.7, 0.2, -0.1], index=["A", "B", "C"])
    res = spatial_plots.plot_spatial_morans_i_barplot(s, tmp_path)
    _assert_figure_result(res)


# ---------------------------------------------------------------------------
# Separability & spillover
# ---------------------------------------------------------------------------


def test_separability_heatmap_saves_data(tmp_path):
    M = pd.DataFrame([[1.0, 0.3], [0.3, 1.0]], index=["A", "B"], columns=["A", "B"])
    res = separability_plots.plot_separability_heatmap(M, tmp_path)
    _assert_figure_result(res)


def test_spillover_heatmap_and_network_save_data(tmp_path):
    M = pd.DataFrame([[0.8, 0.2], [0.15, 0.85]], index=["A", "B"], columns=["A", "B"])
    _assert_figure_result(spillover_plots.plot_spillover_heatmap(M, tmp_path))
    _assert_figure_result(spillover_plots.plot_spillover_network(M, tmp_path, threshold=0.1))
