"""Offline tests for donor-aware gene selection (experimental, opt-in).

Toy per-cell AnnData with explicit donor structure and per-type marker blocks;
verifies DE recovers the right markers, ML-minimal returns a small set, hybrid
composes, single-donor genes are penalised (not silently removed), and small/edge
inputs degrade gracefully.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


@pytest.fixture
def toy_adata():
    """3 types (A/B/C), 6 donors, per-type marker blocks + shared background."""
    rng = np.random.default_rng(0)
    G, types, donors = 60, ["A", "B", "C"], [f"d{i}" for i in range(6)]
    genes = [f"g{i}" for i in range(G)]
    blocks = {"A": range(0, 10), "B": range(10, 20), "C": range(20, 30)}
    X, obs_ct, obs_dn = [], [], []
    for d in donors:
        for t in types:
            for _ in range(25):                       # 25 cells per (donor, type)
                base = rng.poisson(2.0, size=G).astype(float)
                base[list(blocks[t])] += rng.poisson(40.0, size=10)
                X.append(base); obs_ct.append(t); obs_dn.append(d)
    ad = anndata.AnnData(np.vstack(X).astype("float32"))
    ad.var_names = genes
    ad.obs["cell_type"] = pd.Categorical(obs_ct)
    ad.obs["donor_id"] = pd.Categorical(obs_dn)
    return ad, blocks


def test_donor_aware_de_recovers_markers(toy_adata):
    from tissueresolve.reference.gene_selection import donor_aware_de
    ad, blocks = toy_adata
    res = donor_aware_de(ad, "cell_type", "donor_id", "A", top_n=10)
    assert res.genes, "no genes returned"
    expected = {f"g{i}" for i in blocks["A"]}
    hits = sum(g in expected for g in res.genes[:10])
    assert hits >= 7, f"DE recovered only {hits}/10 true A markers"
    # stats present and log2fc positive for top markers
    assert (res.stats.loc[res.genes[0], "log2fc"] > 0)


def test_de_penalises_single_donor_genes(toy_adata):
    """A gene expressed in A but only in ONE donor must be penalised, not dropped."""
    from tissueresolve.reference.gene_selection import donor_aware_de
    ad, _ = toy_adata
    X = ad.X.copy()
    # inject a spurious "marker" high only in donor d0's A cells
    mask = (ad.obs["cell_type"].astype(str).to_numpy() == "A") & \
           (ad.obs["donor_id"].astype(str).to_numpy() == "d0")
    X[mask, 55] += 500.0
    ad2 = ad.copy(); ad2.X = X
    res = donor_aware_de(ad2, "cell_type", "donor_id", "A", top_n=60)
    assert "g55" in res.stats.index
    assert bool(res.stats.loc["g55", "single_donor_penalised"])
    # penalised gene should not outrank the genuine multi-donor markers
    assert res.genes[0] != "g55"


def test_ml_minimal_returns_small_set(toy_adata):
    from tissueresolve.reference.gene_selection import ml_minimal_genes
    ad, _ = toy_adata
    res = ml_minimal_genes(ad, "cell_type", "donor_id",
                           sizes=(5, 10, 15, 20, 30), tol=0.98, seed=0)
    assert 5 <= len(res.genes) <= 30
    assert res.metadata.get("best_balanced_acc", 0) > 0.8


def test_hybrid_composes(toy_adata):
    from tissueresolve.reference.gene_selection import hybrid_gene_set
    ad, _ = toy_adata
    res = hybrid_gene_set(ad, "cell_type", "donor_id", ["A", "B", "C"],
                          de_top_n=30, ml_sizes=(10, 20, 30), seed=0)
    assert res.genes
    assert res.metadata.get("n_de_candidates", 0) >= 5


def test_donor_pseudobulk_shapes(toy_adata):
    from tissueresolve.reference.gene_selection import donor_pseudobulk
    ad, _ = toy_adata
    X, labels, donors, genes = donor_pseudobulk(ad, "cell_type", "donor_id", min_cells=10)
    assert X.shape == (len(labels), len(genes))
    assert set(labels) == {"A", "B", "C"} and len(set(donors)) == 6


def test_sibling_mode_runs(toy_adata):
    from tissueresolve.reference.gene_selection import donor_aware_de
    ad, blocks = toy_adata
    res = donor_aware_de(ad, "cell_type", "donor_id", "A", mode="sibling",
                         siblings=["A", "B"], top_n=10)
    assert res.metadata["mode"] == "sibling"
    expected = {f"g{i}" for i in blocks["A"]}
    assert sum(g in expected for g in res.genes[:10]) >= 6


def test_select_donor_aware_genes_union(toy_adata):
    from tissueresolve.reference.gene_selection import select_donor_aware_genes
    ad, blocks = toy_adata
    genes = select_donor_aware_genes(ad, "cell_type", "donor_id", top_n_per_type=8)
    assert genes and len(genes) == len(set(genes))          # deduped
    all_markers = {f"g{i}" for b in blocks.values() for i in b}
    assert sum(g in all_markers for g in genes) >= 15       # recovers markers of all types


# ---- reference-build integration (opt-in, non-default) --------------------
def _cfg_with_donor():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    cfg.reference.celltype_col = "cell_type"
    cfg.reference.donor_col = "donor_id"
    return cfg


def test_build_reference_stores_donor_de_panel(toy_adata):
    import tissueresolve as tr
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor(),
                             gene_selection="donor_de")
    assert ref.selected_genes, "donor_de panel not stored"
    assert set(ref.selected_genes) <= set(ref.gene_names)
    assert len(ref.selected_genes) < len(ref.gene_names)


def test_default_build_has_no_stored_panel(toy_adata):
    import tissueresolve as tr
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor())
    assert ref.selected_genes is None                       # opt-in only


def test_pipeline_uses_stored_panel(toy_adata):
    """A reference built with a stored panel drives the pipeline's default panel."""
    import numpy as np, tissueresolve as tr
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor(),
                             gene_selection="donor_de")
    rng = np.random.default_rng(1)
    true = rng.dirichlet(np.ones(len(ref.cell_types)), size=4)
    bulk = pd.DataFrame((true @ ref.as_R_cpm()).T, index=list(ref.gene_names),
                        columns=[f"s{i}" for i in range(4)])
    res = tr.deconv_bulk(bulk, ref, resolution_mode="flat", n_bootstrap=0)
    n_used = res.deconv.run_metadata.get("n_genes_panel")
    assert n_used == len([g for g in ref.selected_genes if g in set(bulk.index)])
    assert n_used < len(ref.gene_names)


def test_selected_genes_roundtrip(tmp_path, toy_adata):
    import tissueresolve as tr
    from tissueresolve.results import ReferenceSignature
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor(),
                             gene_selection="donor_de")
    ref.save(tmp_path / "ref")
    loaded = ReferenceSignature.load(tmp_path / "ref")
    assert loaded.selected_genes == ref.selected_genes


def test_solver_backbone_honors_stored_panel(toy_adata):
    """The solver= backbone (flat) uses ref.selected_genes, like the pipeline does."""
    import numpy as np, tissueresolve as tr
    from tissueresolve.solver import get_solver
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor(),
                             gene_selection="donor_de")
    ref_plain = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor())
    rng = np.random.default_rng(2)
    true = rng.dirichlet(np.ones(len(ref.cell_types)), size=4)
    bulk = pd.DataFrame((true @ ref.as_R_cpm()).T, index=list(ref.gene_names),
                        columns=[f"s{i}" for i in range(4)])
    res = get_solver("poisson").solve(bulk, ref)
    res_plain = get_solver("poisson").solve(bulk, ref_plain)
    n_panel = len([g for g in ref.selected_genes if g in set(bulk.index)])
    assert res.diagnostics["n_genes"] == n_panel               # used the stored panel
    assert res_plain.diagnostics["n_genes"] == len(ref.gene_names)  # no panel -> all genes
    assert res.diagnostics["n_genes"] < res_plain.diagnostics["n_genes"]


def test_explicit_genes_override_stored_panel(toy_adata):
    import numpy as np, tissueresolve as tr
    from tissueresolve.solver import PoissonGLMSolver
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type", config=_cfg_with_donor(),
                             gene_selection="donor_de")
    genes = list(ref.gene_names[:12])
    true = np.random.default_rng(3).dirichlet(np.ones(len(ref.cell_types)), size=3)
    bulk = pd.DataFrame((true @ ref.as_R_cpm()).T, index=list(ref.gene_names),
                        columns=[f"s{i}" for i in range(3)])
    res = PoissonGLMSolver(genes=genes).solve(bulk, ref)
    assert res.diagnostics["n_genes"] == len(genes)            # explicit genes win over panel


def test_gene_selection_warns_without_donor(toy_adata):
    import warnings, tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    ad, _ = toy_adata
    cfg = TissueResolveConfig()
    cfg.reference.celltype_col = "cell_type"
    cfg.reference.donor_col = None                           # no donor column
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ref = tr.build_reference(ad, cell_type_col="cell_type", config=cfg,
                                 gene_selection="donor_de")
    assert ref.selected_genes is None
    assert any("donor" in str(x.message).lower() for x in w)
