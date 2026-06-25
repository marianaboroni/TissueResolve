"""Offline tests for signature-quality & hierarchy figures (report Part 4)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting import signature_qc_plots as S

_MAP = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "Treg": "T/NK",
        "macrophage": "Myeloid", "monocyte": "Myeloid"}

_SEP = pd.DataFrame({
    "type_a": ["CD8 T cell", "CD4 T cell", "macrophage", "CD8 T cell", "Treg"],
    "type_b": ["CD4 T cell", "Treg", "monocyte", "macrophage", "CD4 T cell"],
    "bhattacharyya": [0.99, 0.95, 0.90, 0.30, 0.92],
    "separability_score": [0.01, 0.05, 0.10, 0.70, 0.08],
    "n_discriminating_genes": [12, 20, 30, 120, 18],
    "resolvability": ["unresolved", "unresolved", "poorly_resolved",
                      "resolved", "unresolved"],
})


def _ok(res):
    assert res.data_paths and all(p.exists() for p in res.data_paths.values())


def test_signature_matrix_heatmap_caps_top_genes(tmp_path):
    rng = np.random.default_rng(0)
    genes = [f"GENE_{i}" for i in range(120)]
    cts = list(_MAP)
    M = rng.gamma(1.0, 1.0, (120, len(cts)))
    res = S.plot_signature_matrix_heatmap(M, genes, cts, _MAP, tmp_path,
                                          top_n_genes=30)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#", index_col=0)
    assert df.shape[0] == 30  # only top genes shown, not all 120
    assert df.shape[1] == len(cts)


def test_top_confusable_pairs_specific(tmp_path):
    res = S.plot_top_confusable_pairs(_SEP, tmp_path, top_n=3)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    # most confusable (lowest separability) included
    assert "CD8 T cell" in set(df["type_a"]) | set(df["type_b"])
    assert "Visual summary" not in (res.caption or "")


def test_within_vs_between_family(tmp_path):
    res = S.plot_within_vs_between_family_separability(_SEP, _MAP, tmp_path)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert set(df["group"]) <= {"within-family", "between-family"}
    # CD8↔CD4 is within T/NK; macrophage↔monocyte within Myeloid
    assert "within-family" in set(df["group"])


def test_hierarchy_map_has_broad_and_fine(tmp_path):
    counts = {"CD8 T cell": 800, "CD4 T cell": 600, "Treg": 100,
              "macrophage": 500, "monocyte": 90}
    res = S.plot_hierarchy_map(counts, _MAP, tmp_path,
                               unresolved_families=["T/NK"])
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert set(df["level"]) == {"broad", "fine"}
    assert df["unresolved"].any()


def test_marker_support_by_family(tmp_path):
    res = S.plot_marker_support_by_family(_SEP, _MAP, tmp_path)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "mean_discriminating_genes" in df.columns


def test_missing_columns_raise(tmp_path):
    bad = pd.DataFrame({"x": [1, 2]})
    with pytest.raises(ValueError):
        S.plot_top_confusable_pairs(bad, tmp_path)
    with pytest.raises(ValueError):
        S.plot_marker_support_by_family(bad, _MAP, tmp_path)
