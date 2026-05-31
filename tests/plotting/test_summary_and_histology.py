"""
Tests for the main summary figures and H&E histology overlays (offline).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting.histology import (  # noqa: E402
    plot_abundance_on_he, plot_dominant_cell_type_on_he, plot_he_with_spots,
)
from tissueresolve.plotting.palette import assign_family_palette  # noqa: E402
from tissueresolve.plotting.summary_figures import (  # noqa: E402
    bulk_main_summary_figure, spatial_main_summary_figure,
)


def _props(n=8, cols=None, seed=0):
    cols = cols or ["CD4 T cell", "CD8 T cell", "macrophage", "monocyte",
                    "luminal epithelial cell", "fibroblast"]
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.dirichlet(np.ones(len(cols)), size=n),
                        index=[f"s{i}" for i in range(n)], columns=cols)


def _toy_coords(props):
    rng = np.random.default_rng(1)
    n = len(props)
    return pd.DataFrame({
        "spot": props.index, "x": rng.random(n), "y": rng.random(n),
        "image_x": rng.random(n) * 50, "image_y": rng.random(n) * 50,
        "array_row": rng.integers(0, 10, n), "array_col": rng.integers(0, 10, n),
        "tissue_status": 1,
    })


# ---------------------------------------------------------------------------
# Summary figures
# ---------------------------------------------------------------------------


def test_bulk_main_summary_figure(tmp_path):
    props = _props(cols=[f"ct{i}" for i in range(14)])  # > top_n → "Other"
    res = bulk_main_summary_figure(props, tmp_path, top_n=10)
    assert res.html_path.exists()
    assert res.data_paths  # source data saved
    assert (tmp_path / "bulk_main_summary_figure.caption.txt").exists()
    data = pd.read_csv(res.data_paths["data"], sep="\t", comment="#", index_col=0)
    assert "Other" in data.columns


def test_spatial_main_summary_figure(tmp_path):
    props = _props()
    coords = _toy_coords(props)[["array_row", "array_col"]]
    coords.index = props.index
    morans = pd.Series(np.linspace(0.1, 0.8, props.shape[1]), index=props.columns)
    res = spatial_main_summary_figure(props, tmp_path, coords=coords, morans=morans)
    assert res.html_path.exists()
    assert (tmp_path / "spatial_main_summary_figure.caption.txt").exists()
    assert res.data_paths


# ---------------------------------------------------------------------------
# H&E overlays
# ---------------------------------------------------------------------------


def test_he_with_spots_toy_image(tmp_path):
    props = _props()
    coords = _toy_coords(props)
    img = (np.random.default_rng(0).random((30, 30, 3)) * 255).astype("uint8")
    res = plot_he_with_spots(img, coords, tmp_path)
    assert res.html_path.exists()
    # coordinate-transform metadata recorded in the source-data comment header
    header = res.data_paths["data"].read_text().splitlines()[:6]
    joined = "\n".join(header)
    assert "y_axis_inverted" in joined
    assert "has_he_image: True" in joined
    assert "image_x" in joined  # column manifest recorded


def test_he_without_image_warns(tmp_path):
    props = _props()
    coords = _toy_coords(props)
    res = plot_he_with_spots(None, coords, tmp_path, name="he_noimg")
    assert res.html_path.exists()
    assert any("No H&E image" in w for w in res.warnings)
    header = "\n".join(res.data_paths["data"].read_text().splitlines()[:6])
    assert "has_he_image: False" in header


def test_dominant_and_abundance_on_he(tmp_path):
    props = _props()
    coords = _toy_coords(props)
    img = (np.random.default_rng(0).random((30, 30, 3)) * 255).astype("uint8")
    cmap = assign_family_palette(list(props.columns))
    rd = plot_dominant_cell_type_on_he(img, coords, props, cmap, tmp_path)
    assert rd.html_path.exists()
    d = pd.read_csv(rd.data_paths["data"], sep="\t", comment="#")
    assert "dominant_type" in d.columns and "dominant_fraction" in d.columns
    ra = plot_abundance_on_he(img, coords, props.iloc[:, 0], tmp_path,
                              cell_type=str(props.columns[0]))
    assert ra.html_path.exists()
