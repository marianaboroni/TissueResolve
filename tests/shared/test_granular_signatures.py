"""Offline tests for multi-granularity (broad / cell type / state) gene panels."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("anndata")
import anndata as ad

from tissueresolve.reference import granular_signatures as G


@pytest.fixture
def toy_three_level():
    """Broad {T/NK, Myeloid}; cell types within; states within cell types.

    Gene blocks: broad markers (flat within a family), cell-type markers (flat
    within a cell type), state markers (vary within a cell type).
    """
    rng = np.random.default_rng(0)
    # (broad, cell_type, state, n)
    spec = [
        ("T/NK", "CD8 T cell", "CD8_eff", 50), ("T/NK", "CD8 T cell", "CD8_mem", 50),
        ("T/NK", "CD4 T cell", "CD4_naive", 50), ("T/NK", "CD4 T cell", "Treg", 50),
        ("Myeloid", "macrophage", "macro_M1", 50), ("Myeloid", "macrophage", "macro_M2", 50),
        ("Myeloid", "monocyte", "mono_classical", 50),
    ]
    broad, ctype, state = [], [], []
    for b, c, s, n in spec:
        broad += [b] * n; ctype += [c] * n; state += [s] * n
    nc = len(state)
    genes = (["BROAD_TNK", "BROAD_MYE"]                       # broad
             + ["CT_CD8", "CT_CD4", "CT_MAC", "CT_MONO"]      # cell type
             + ["ST_CD8eff", "ST_CD8mem", "ST_CD4naive", "ST_Treg",
                "ST_M1", "ST_M2"]                              # state
             + [f"BG_{i}" for i in range(30)])
    gi = {g: i for i, g in enumerate(genes)}
    X = rng.poisson(0.4, (nc, len(genes))).astype(np.float64)
    b_arr, c_arr, s_arr = np.array(broad), np.array(ctype), np.array(state)
    X[b_arr == "T/NK", gi["BROAD_TNK"]] += rng.poisson(15, (b_arr == "T/NK").sum())
    X[b_arr == "Myeloid", gi["BROAD_MYE"]] += rng.poisson(15, (b_arr == "Myeloid").sum())
    for c, g in [("CD8 T cell", "CT_CD8"), ("CD4 T cell", "CT_CD4"),
                 ("macrophage", "CT_MAC"), ("monocyte", "CT_MONO")]:
        X[c_arr == c, gi[g]] += rng.poisson(12, (c_arr == c).sum())
    for st, g in [("CD8_eff", "ST_CD8eff"), ("CD8_mem", "ST_CD8mem"),
                  ("CD4_naive", "ST_CD4naive"), ("Treg", "ST_Treg"),
                  ("macro_M1", "ST_M1"), ("macro_M2", "ST_M2")]:
        X[s_arr == st, gi[g]] += rng.poisson(12, (s_arr == st).sum())
    obs = pd.DataFrame({"broad": broad, "cell_type": ctype, "state": state},
                       index=[f"c{i}" for i in range(nc)])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))


def test_broad_panel_distinguishes_families(toy_three_level):
    panel = G.build_broad_gene_panel(toy_three_level, "broad",
                                     n_hvgs_per_family=20, min_cells_per_subtype=20)
    assert "BROAD_TNK" in panel or "BROAD_MYE" in panel
    # temp root column is cleaned up
    assert "__tr_root__" not in toy_three_level.obs.columns


def _disc_genes(fp, group):
    """Discriminative (pairwise) genes selected for a group — the meaningful set
    (HVGs may also pull in flat noise genes, which carry no discriminatory
    weight; the pairwise set is what distinguishes the children)."""
    pw = fp.pairwise_markers
    if pw.empty:
        return set()
    return set(pw.loc[pw["family"] == group, "gene"])


def test_celltype_panels_within_family(toy_three_level):
    fp = G.build_celltype_within_family_gene_panels(
        toy_three_level, "broad", "cell_type", n_hvgs_per_family=20,
        min_cells_per_subtype=20, min_logfc=1.0)
    assert "T/NK" in fp.family_panels and "Myeloid" in fp.family_panels
    # cell-type markers are DISCRIMINATIVE within their family
    assert "CT_CD8" in _disc_genes(fp, "T/NK")
    assert "CT_MAC" in _disc_genes(fp, "Myeloid")
    # a broad marker active-but-flat within the family (BROAD_TNK is high in all
    # T/NK cells) is NOT a discriminative within-family gene
    assert "BROAD_TNK" not in _disc_genes(fp, "T/NK")


def test_state_panels_within_celltype_not_global(toy_three_level):
    fp = G.build_state_within_celltype_gene_panels(
        toy_three_level, "cell_type", "state", n_hvgs_per_family=20,
        min_cells_per_subtype=20, min_logfc=1.0)
    # panels keyed by cell type; states resolved within their cell type
    assert "CD8 T cell" in fp.family_panels and "macrophage" in fp.family_panels
    assert "ST_CD8eff" in _disc_genes(fp, "CD8 T cell")
    assert "ST_M1" in _disc_genes(fp, "macrophage")
    # the cell-type marker (high in ALL CD8 cells, hence flat across CD8 states)
    # is NOT a within-cell-type state discriminator — i.e. selection is within
    # the cell type, not reusing the cell-type-level genes
    assert "CT_CD8" not in _disc_genes(fp, "CD8 T cell")
    # monocyte has a single state -> skipped (cannot sub-resolve)
    assert "monocyte" in fp.families_skipped


def test_panels_differ_by_granularity(toy_three_level):
    res = G.build_multigranularity_panels(
        toy_three_level, "broad", "cell_type", "state",
        n_hvgs_per_family=20, min_cells_per_subtype=20)
    broad = set(res["broad_panel"])
    ct = set(g for genes in res["celltype_panels"].family_panels.values() for g in genes)
    st = set(g for genes in res["state_panels"].family_panels.values() for g in genes)
    # the three levels are not identical sets
    assert broad != ct and ct != st


def test_granularity_weights_and_outputs(tmp_path, toy_three_level):
    res = G.build_multigranularity_panels(
        toy_three_level, "broad", "cell_type", "state",
        n_hvgs_per_family=20, min_cells_per_subtype=20)
    w = res["gene_weights"]
    assert set(w["granularity"]) <= {"broad", "cell_type", "state"}
    assert {"broad", "cell_type", "state"} <= set(w["granularity"])
    paths = G.write_granular_signature_outputs(res, tmp_path)
    for key in ("broad_gene_panel", "celltype_within_family_gene_panels",
                "state_within_celltype_gene_panels", "granularity_gene_weights",
                "pairwise_state_marker_support"):
        assert paths[key].exists()
    # pairwise state markers carry cell_type + state pair columns
    pw = pd.read_csv(paths["pairwise_state_marker_support"], sep="\t")
    if not pw.empty:
        assert {"cell_type", "state_a", "state_b"} <= set(pw.columns)


def test_two_level_skips_state(toy_three_level):
    res = G.build_multigranularity_panels(
        toy_three_level, "broad", "cell_type", state_col=None,
        n_hvgs_per_family=20, min_cells_per_subtype=20)
    assert res["state_panels"] is None
    assert (res["gene_weights"]["granularity"] == "state").sum() == 0
