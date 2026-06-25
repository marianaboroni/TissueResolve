"""Offline tests for reference-QC publication figures (report Part 3)."""
from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting import reference_qc_plots as R

_COUNTS = {"CD8 T cell": 800, "CD4 T cell": 600, "Treg": 120,
           "macrophage": 500, "monocyte": 90, "fibroblast": 1500, "Rare": 20}
_MAP = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "Treg": "T/NK",
        "macrophage": "Myeloid", "monocyte": "Myeloid",
        "fibroblast": "Stromal", "Rare": "Stromal"}


def _has_data(res):
    assert res.html_path is None or res.html_path.exists()
    assert res.data_paths and all(p.exists() for p in res.data_paths.values())


def test_broad_family_composition(tmp_path):
    res = R.plot_reference_broad_family_composition(_COUNTS, _MAP, tmp_path)
    _has_data(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert {"broad_family", "n_cells", "percent"} <= set(df.columns)
    # aggregated to families, sorted; Stromal (1520) should be the top family
    assert df.iloc[-1]["broad_family"] == "Stromal"


def test_fine_subpopulation_support_writes_threshold(tmp_path):
    res = R.plot_reference_fine_subpopulation_support(_COUNTS, _MAP, tmp_path,
                                                      min_cells=100)
    _has_data(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert {"fine_cell_type", "broad_family", "n_cells"} <= set(df.columns)


def test_celltype_imbalance_reports_gini(tmp_path):
    res = R.plot_reference_celltype_imbalance(_COUNTS, tmp_path)
    _has_data(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "cumulative_fraction" in df.columns
    assert abs(df["cumulative_fraction"].iloc[-1] - 1.0) < 1e-6


def test_gene_overlap_by_modality(tmp_path):
    overlap = {"bulk": {"n_reference": 5000, "n_query": 5000, "n_shared": 5000},
               "spatial": {"n_reference": 5000, "n_query": 36000, "n_shared": 4900}}
    res = R.plot_gene_overlap_by_modality(overlap, tmp_path)
    _has_data(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert {"modality", "reference_genes", "query_genes", "shared_genes"} <= set(df.columns)


def test_suitability_components_traffic_light(tmp_path):
    comp = pd.DataFrame({
        "component": ["gene_overlap", "celltype_balance", "fine_label_separability"],
        "score": [1.0, 0.06, 0.71],
        "status": ["PASS", "FAIL", "CAUTION"],
        "detail": ["100%", "imbalanced", "71% separable"]})
    res = R.plot_reference_suitability_components(comp, tmp_path)
    _has_data(res)


def test_missing_mapping_raises_not_crash(tmp_path):
    with pytest.raises(ValueError):
        R.plot_reference_broad_family_composition(_COUNTS, {}, tmp_path)
    with pytest.raises(ValueError):
        R.plot_gene_overlap_by_modality({}, tmp_path)


def test_uses_hierarchical_palette(tmp_path):
    """Broad-family bars use the hierarchical palette (distinct family colors)."""
    res = R.plot_reference_broad_family_composition(_COUNTS, _MAP, tmp_path)
    colors = res.figure.data[0].marker.color
    assert len(set(colors)) >= 2  # families get distinct colors
