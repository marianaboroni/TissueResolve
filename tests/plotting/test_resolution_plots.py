"""Offline tests for resolution/spillover/unresolved figures (report Part 8)."""
from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting import resolution_plots as RP

_MAP = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "Treg": "T/NK",
        "macrophage": "Myeloid", "monocyte": "Myeloid", "fibroblast": "Stromal"}

_SEP = pd.DataFrame({
    "type_a": ["CD8 T cell", "CD4 T cell", "macrophage", "CD8 T cell"],
    "type_b": ["CD4 T cell", "Treg", "monocyte", "fibroblast"],
    "separability_score": [0.02, 0.06, 0.55, 0.95],
})


def _ok(res):
    assert res.data_paths and all(p.exists() for p in res.data_paths.values())


def test_separability_distribution(tmp_path):
    res = RP.plot_separability_distribution(_SEP, tmp_path, high_risk_threshold=0.25)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "separability_score" in df.columns


def test_unresolved_mass_by_family(tmp_path):
    um = pd.DataFrame({
        "sample": ["s0", "s1"],
        "unresolved_T/NK": [0.24, 0.37],
        "unresolved_Myeloid": [0.08, 0.12],
        "unresolved_Epithelial": [0.0, 0.24],
    }).set_index("sample")
    res = RP.plot_unresolved_mass_by_family(um, tmp_path)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "mean_unresolved_mass" in df.columns
    # prefix stripped
    assert "T/NK" in set(df["broad_family"])


def test_trusted_resolution_summary_levels(tmp_path):
    res = RP.plot_trusted_resolution_summary(
        _MAP, _SEP, tmp_path, unresolved_families=["T/NK"])
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert {"broad_family", "recommended_level", "separability_status"} <= set(df.columns)
    # T/NK marked unresolved -> family-level
    tnk = df[df["broad_family"] == "T/NK"].iloc[0]
    assert tnk["recommended_level"] == "broad / family-level"


def test_missing_columns_raise(tmp_path):
    with pytest.raises(ValueError):
        RP.plot_separability_distribution(pd.DataFrame({"x": [1]}), tmp_path)
    with pytest.raises(ValueError):
        RP.plot_trusted_resolution_summary({}, _SEP, tmp_path)
