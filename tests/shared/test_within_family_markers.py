"""
Offline tests for within-family HVG / DE / pairwise marker selection.

A toy AnnData encodes the key structure: genes that separate broad *families*
are flat *within* a family, while distinct genes separate subtypes *within* each
family.  Within-family selection must recover the latter, not the former.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("anndata")
import anndata as ad

from tissueresolve.reference import within_family_markers as W


@pytest.fixture
def toy_adata():
    rng = np.random.default_rng(0)
    # cells: family F1 {A,B}, family F2 {C,D}, plus singleton family F3 {E}
    spec = [("F1", "A", 60), ("F1", "B", 60), ("F2", "C", 60),
            ("F2", "D", 60), ("F3", "E", 5)]
    broad, fine = [], []
    for b, f, n in spec:
        broad += [b] * n
        fine += [f] * n
    n_cells = len(fine)
    gene_names = (["BROAD_F1", "BROAD_F2"]                 # separate families
                  + ["A_MARK", "B_MARK", "C_MARK", "D_MARK"]  # within-family
                  + ["MT-CO1", "RPS12", "FOS"]              # to be filtered
                  + [f"BG_{i}" for i in range(40)])          # background
    gi = {g: i for i, g in enumerate(gene_names)}
    X = rng.poisson(0.4, (n_cells, len(gene_names))).astype(np.float64)
    fine_arr = np.array(fine)
    broad_arr = np.array(broad)
    # broad markers (flat within a family, differ between families)
    X[broad_arr == "F1", gi["BROAD_F1"]] += rng.poisson(15, (broad_arr == "F1").sum())
    X[broad_arr == "F2", gi["BROAD_F2"]] += rng.poisson(15, (broad_arr == "F2").sum())
    # within-family subtype markers
    for sub, g in [("A", "A_MARK"), ("B", "B_MARK"),
                   ("C", "C_MARK"), ("D", "D_MARK")]:
        X[fine_arr == sub, gi[g]] += rng.poisson(12, (fine_arr == sub).sum())
    # mito/ribo/stress: subtype-A-variable (so they ARE candidates) but must be
    # filtered out of the panel for being mito/ribo/stress.
    for g in ("MT-CO1", "RPS12", "FOS"):
        X[fine_arr == "A", gi[g]] += rng.poisson(12, (fine_arr == "A").sum())
    obs = pd.DataFrame({"broad": broad, "fine": fine},
                       index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=gene_names)
    return ad.AnnData(X=X, obs=obs, var=var)


def test_hvgs_selected_per_family_and_differ(toy_adata):
    h1 = W.select_within_family_hvgs(toy_adata, "broad", "fine", "F1",
                                     n_hvgs=10, min_cells_per_subtype=20,
                                     min_detection=0.0)
    h2 = W.select_within_family_hvgs(toy_adata, "broad", "fine", "F2",
                                     n_hvgs=10, min_cells_per_subtype=20,
                                     min_detection=0.0)
    # within-family subtype markers are the within-family HVGs
    assert "A_MARK" in h1.index and "B_MARK" in h1.index
    assert "C_MARK" in h2.index and "D_MARK" in h2.index
    # the two families pick different HVGs (not just one global set)
    assert set(h1.index) != set(h2.index)
    # the cross-family broad marker is NOT a within-F1 HVG (flat within F1)
    assert "BROAD_F2" not in h1.index


def test_singleton_family_raises(toy_adata):
    with pytest.raises(ValueError):
        W.select_within_family_hvgs(toy_adata, "broad", "fine", "F3")


def test_pairwise_markers_have_family_and_pair_columns(toy_adata):
    pw = W.select_pairwise_discriminative_genes_within_family(
        toy_adata, "broad", "fine", "F1", top_n_per_pair=20,
        min_logfc=0.25, min_cells_per_subtype=20)
    assert {"gene", "type_a", "type_b", "log2fc", "pair_score"} <= set(pw.columns)
    # A_MARK / B_MARK discriminate A vs B
    assert {"A_MARK", "B_MARK"} & set(pw["gene"])


def test_de_genes_one_vs_rest(toy_adata):
    de = W.select_within_family_de_genes(toy_adata, "broad", "fine", "F1",
                                         top_n_per_pair=20, min_logfc=0.25)
    assert {"gene", "subtype", "log2fc"} <= set(de.columns)
    a_markers = set(de.loc[de["subtype"] == "A", "gene"])
    assert "A_MARK" in a_markers


def test_build_panels_filters_and_scores(toy_adata, tmp_path):
    res = W.build_family_specific_gene_panels(
        toy_adata, "broad", "fine", n_hvgs_per_family=20,
        top_n_de_per_pair=20, min_cells_per_subtype=20, min_logfc=0.25)
    # F1 and F2 get panels; F3 (singleton) is skipped
    assert "F1" in res.family_panels and "F2" in res.family_panels
    assert "F3" in res.families_skipped
    # panels differ between families (family-specific, not global)
    assert set(res.family_panels["F1"]) != set(res.family_panels["F2"])
    # within-family subtype markers are in the panel
    assert "A_MARK" in res.family_panels["F1"]
    # mito/ribo/stress are excluded with reasons
    excl = set(res.excluded_genes["gene"])
    assert {"MT-CO1", "RPS12", "FOS"} <= excl
    assert "MT-CO1" not in res.family_panels["F1"]
    # weights normalized 0..1 within family
    w = res.gene_weights
    assert (w["final_within_family_weight"] >= 0).all()
    assert (w["final_within_family_weight"] <= 1.0 + 1e-9).all()
    # support tables present
    assert {"family", "n_panel_genes"} <= set(res.marker_support_by_family.columns)
    # save writes the Part-3 TSVs
    paths = res.save(tmp_path)
    for key in ("within_family_marker_scores", "within_family_selected_genes",
                "within_family_pairwise_markers", "within_family_gene_weights",
                "within_family_excluded_genes", "marker_support_by_family",
                "marker_support_by_pair"):
        assert paths[key].exists()


def test_query_detection_filter(toy_adata):
    # query missing the within-family markers entirely -> excluded by query detection
    qgenes = [g for g in toy_adata.var_names if g not in ("A_MARK", "B_MARK")]
    q = toy_adata[:, qgenes].copy()
    res = W.build_family_specific_gene_panels(
        toy_adata, "broad", "fine", query_adata=q, min_query_detection=0.01,
        n_hvgs_per_family=20, top_n_de_per_pair=20, min_cells_per_subtype=20)
    reasons = res.excluded_genes
    # A_MARK absent in query -> excluded for query detection (in F1)
    amark = reasons[reasons["gene"] == "A_MARK"]
    assert not amark.empty and amark["reason"].str.contains("query_detection").any()
