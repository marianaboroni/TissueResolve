"""
Offline synthetic tests for hierarchical broad-to-fine deconvolution.

All tests are deterministic and offline (no network, no dataset downloads).
They cover the hierarchy arithmetic, annotation validation, within-family
resolvability gating, the bulk and spatial hierarchical workflows, the
reproducible colour map, and the CLI surface.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature


# ---------------------------------------------------------------------------
# Fixtures: a reference with one separable and one non-separable family
# ---------------------------------------------------------------------------

FAMILY_MAP = {
    "X_sub1": "FamX", "X_sub2": "FamX",   # separable subtypes
    "Y_sub1": "FamY", "Y_sub2": "FamY",   # near-identical subtypes
    "Z_solo": "FamZ",                      # singleton family
}


@pytest.fixture
def hier_ref() -> ReferenceSignature:
    rng = np.random.default_rng(0)
    G = 120
    genes = [f"g{i}" for i in range(G)]

    def prof(active, lvl=200.0):
        v = np.full(G, 2.0)
        v[list(active)] = lvl
        return v

    X1 = prof(range(0, 20))
    X2 = prof(range(20, 40))
    Y1 = prof(range(40, 60))
    Y2 = Y1 * (1.0 + rng.normal(0, 0.004, G))   # nearly identical to Y1
    Z = prof(range(60, 80))
    R = np.vstack([X1, X2, Y1, Y2, Z]).astype(np.float32)
    cts = ["X_sub1", "X_sub2", "Y_sub1", "Y_sub2", "Z_solo"]
    return ReferenceSignature(
        gene_names=genes, cell_types=cts, R_cpm=R,
        R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={c: 100 for c in cts},
    )


# ---------------------------------------------------------------------------
# Hierarchy inference and mapping I/O
# ---------------------------------------------------------------------------


def test_infer_broad_cell_type_family_keywords():
    from tissueresolve.reference.hierarchy import infer_broad_cell_type_family as f
    assert f("CD8-positive memory T cell") == "T/NK"
    assert f("natural killer cell") == "T/NK"
    assert f("IgA plasma cell") == "B/plasma"
    assert f("classical monocyte") == "Myeloid"
    assert f("capillary endothelial cell") == "Endothelial"
    assert f("luminal epithelial cell of mammary gland") == "Epithelial"
    assert f("myofibroblast") == "Stromal/fibroblast"
    assert f("pericyte") == "Mural"
    assert f("adipocyte") == "Adipocyte"
    assert f("some unknown cell") == "Other"


def test_build_cell_type_hierarchy_with_mapping():
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    cts = ["a", "b", "c"]
    mp = {"a": "F1", "b": "F1", "c": "F2"}
    out = build_cell_type_hierarchy(cts, mp)
    assert out == mp


def test_build_cell_type_hierarchy_missing_raises():
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    with pytest.raises(ValueError):
        build_cell_type_hierarchy(["a", "b"], {"a": "F1"})


def test_build_cell_type_hierarchy_infers_with_warning():
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    with pytest.warns(UserWarning):
        out = build_cell_type_hierarchy(["CD8 T cell", "macrophage"], None)
    assert out["CD8 T cell"] == "T/NK"
    assert out["macrophage"] == "Myeloid"


def test_save_load_hierarchy_mapping_roundtrip(tmp_path):
    from tissueresolve.reference.hierarchy import (
        save_hierarchy_mapping, load_hierarchy_mapping,
    )
    mp = {"a": "F1", "b": "F1", "c": "F2"}
    path = save_hierarchy_mapping(mp, tmp_path / "h.tsv")
    assert load_hierarchy_mapping(path) == mp


def test_load_hierarchy_mapping_alt_headers(tmp_path):
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping
    p = tmp_path / "m.tsv"
    p.write_text("fine_cell_type\tbroad_cell_type\nCD8 T cell\tT/NK\nmacrophage\tMyeloid\n")
    mp = load_hierarchy_mapping(p)
    assert mp == {"CD8 T cell": "T/NK", "macrophage": "Myeloid"}


def test_load_hierarchy_mapping_conflict_raises(tmp_path):
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping
    p = tmp_path / "m.tsv"
    p.write_text("cell_type\tfamily\nx\tA\nx\tB\n")
    with pytest.raises(ValueError):
        load_hierarchy_mapping(p)


# ---------------------------------------------------------------------------
# Reference aggregation: preserves gene order & total mass
# ---------------------------------------------------------------------------


def test_aggregate_reference_by_family(hier_ref):
    from tissueresolve.reference.hierarchy import aggregate_reference_by_family
    fam_ref = aggregate_reference_by_family(hier_ref, FAMILY_MAP)
    assert set(fam_ref.cell_types) == {"FamX", "FamY", "FamZ"}
    # gene order preserved
    assert list(fam_ref.gene_names) == list(hier_ref.gene_names)
    # cell counts preserved (sum over members)
    assert fam_ref.n_cells_per_type["FamX"] == 200


def test_aggregate_predictions_by_family_preserves_mass():
    from tissueresolve.reference.hierarchy import aggregate_predictions_by_family
    df = pd.DataFrame({"a": [0.2, 0.5], "b": [0.3, 0.1], "c": [0.5, 0.4]})
    fam = aggregate_predictions_by_family(df, {"a": "F", "b": "F", "c": "G"})
    assert list(fam.columns) == ["F", "G"]
    np.testing.assert_allclose(fam.sum(axis=1), df.sum(axis=1))


# ---------------------------------------------------------------------------
# Conditional / combine / unresolved arithmetic (mass preservation)
# ---------------------------------------------------------------------------


def test_conditional_sums_to_one_within_family():
    from tissueresolve.reference.hierarchy import (
        compute_conditional_subtype_proportions,
    )
    fam = pd.DataFrame({"F": [0.6], "G": [0.4]})
    fine = pd.DataFrame({"a": [0.3], "b": [0.1], "c": [0.6]})
    mp = {"a": "F", "b": "F", "c": "G"}
    cond = compute_conditional_subtype_proportions(fine, fam, mp)
    np.testing.assert_allclose(cond[["a", "b"]].sum(axis=1), 1.0)
    np.testing.assert_allclose(cond[["c"]].sum(axis=1), 1.0)


def test_combine_and_unresolved_preserve_total_mass():
    from tissueresolve.reference.hierarchy import (
        compute_conditional_subtype_proportions,
        combine_family_and_conditional_estimates,
        add_unresolved_family_mass,
    )
    fam = pd.DataFrame({"F": [0.6, 0.5], "G": [0.4, 0.5]})
    fine = pd.DataFrame({"a": [0.3, 0.2], "b": [0.1, 0.3], "c": [0.6, 0.5]})
    mp = {"a": "F", "b": "F", "c": "G"}
    cond = compute_conditional_subtype_proportions(fine, fam, mp)
    absolute = combine_family_and_conditional_estimates(fam, cond, mp)
    np.testing.assert_allclose(absolute.sum(axis=1), 1.0)
    resolved, unresolved = add_unresolved_family_mass(absolute, fam, mp, ["G"])
    # G's subtype (c) is zeroed; its mass moves to unresolved_G
    assert (resolved["c"] == 0).all()
    np.testing.assert_allclose(unresolved["unresolved_G"], fam["G"])
    total = resolved.sum(axis=1) + unresolved.sum(axis=1)
    np.testing.assert_allclose(total, 1.0)


# ---------------------------------------------------------------------------
# Within-family resolvability gating
# ---------------------------------------------------------------------------


def test_resolvability_flags_nonseparable_family(hier_ref):
    from tissueresolve.reference.hierarchy import (
        evaluate_within_family_resolvability, decide_unresolved_families,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = evaluate_within_family_resolvability(
            hier_ref, FAMILY_MAP, min_discriminating_genes=5)
    assert bool(res.loc["FamX", "resolvable"]) is True
    assert bool(res.loc["FamY", "resolvable"]) is False
    # singleton family is trivially resolvable
    assert bool(res.loc["FamZ", "resolvable"]) is True
    unresolved = decide_unresolved_families(res, allow_unresolved=True)
    assert unresolved == ["FamY"]


def test_no_unresolved_when_disallowed(hier_ref):
    from tissueresolve.reference.hierarchy import (
        evaluate_within_family_resolvability, decide_unresolved_families,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = evaluate_within_family_resolvability(hier_ref, FAMILY_MAP)
    assert decide_unresolved_families(res, allow_unresolved=False) == []


# ---------------------------------------------------------------------------
# Within-family gene panels (additive; default behaviour unchanged)
# ---------------------------------------------------------------------------

def test_family_panel_changes_separability(hier_ref):
    """The resolvability call must actually USE the family panel: restricting
    FamX to flat (non-discriminating) genes makes it unresolvable, while the
    discriminating panel keeps it resolvable."""
    from tissueresolve.reference.hierarchy import evaluate_within_family_resolvability
    genes = [str(g) for g in hier_ref.gene_names]
    flat_panel = {"FamX": genes[80:120]}        # all flat (value 2.0) → not separable
    disc_panel = {"FamX": genes[0:40]}          # the X_sub1/X_sub2 discriminating block
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_flat = evaluate_within_family_resolvability(
            hier_ref, FAMILY_MAP, min_discriminating_genes=5,
            family_gene_panels=flat_panel)
        res_disc = evaluate_within_family_resolvability(
            hier_ref, FAMILY_MAP, min_discriminating_genes=5,
            family_gene_panels=disc_panel)
    assert bool(res_flat.loc["FamX", "resolvable"]) is False   # flat panel kills it
    assert bool(res_disc.loc["FamX", "resolvable"]) is True    # discriminating panel keeps it
    assert int(res_disc.loc["FamX", "n_panel_genes"]) == 40


def test_default_behaviour_unchanged_without_panels(hier_ref):
    """No panels → identical verdict to the original global-gene path."""
    from tissueresolve.reference.hierarchy import evaluate_within_family_resolvability
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = evaluate_within_family_resolvability(hier_ref, FAMILY_MAP,
                                                 min_discriminating_genes=5)
        b = evaluate_within_family_resolvability(hier_ref, FAMILY_MAP,
                                                 min_discriminating_genes=5,
                                                 family_gene_panels=None)
    assert a.drop(columns=["n_panel_genes"]).equals(b.drop(columns=["n_panel_genes"]))


def test_run_hierarchical_bulk_accepts_family_panels(hier_ref):
    """The bulk hierarchical entry point threads family_gene_panels without error
    and records it in metadata."""
    from tissueresolve.bulk.hierarchical import run_hierarchical_bulk
    # tiny pseudobulk from the reference profiles (samples × genes)
    R = hier_ref.as_R_cpm()
    samples = pd.DataFrame(
        {f"s{i}": R[i % R.shape[0]] for i in range(3)},
        index=hier_ref.gene_names).T
    genes = [str(g) for g in hier_ref.gene_names]
    panels = {"FamX": genes[0:40], "FamY": genes[40:60]}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_hierarchical_bulk(samples.T, hier_ref, FAMILY_MAP,
                                    family_gene_panels=panels,
                                    min_discriminating_genes=5)
    assert res.estimates.metadata.get("within_family_panels") is True
    # combined fine proportions still sum to ~1 per sample (mass preserved)
    s = res.estimates.combined_fine.sum(axis=1)
    assert np.allclose(s.to_numpy(), 1.0, atol=1e-6)


def test_within_family_resolution_summary(hier_ref):
    from tissueresolve.reference.within_family_markers import (
        within_family_resolution_summary)
    genes = [str(g) for g in hier_ref.gene_names]
    panels = {"FamX": genes[0:40], "FamY": genes[40:60]}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        glob, within, summary = within_family_resolution_summary(
            hier_ref, FAMILY_MAP, panels, min_discriminating_genes=5)
    assert {"family", "global_mean_separability",
            "within_family_mean_separability", "global_resolvable",
            "within_family_resolvable", "separability_improved",
            "verdict_changed"} <= set(summary.columns)


# ---------------------------------------------------------------------------
# Annotation detection & validation
# ---------------------------------------------------------------------------


def _toy_obs():
    return pd.DataFrame({
        "broad_cell_type": ["T/NK", "T/NK", "Myeloid", "Epithelial", "Epithelial"],
        "cell_type": ["CD8 T cell", "CD4 T cell", "macrophage", "luminal", "basal"],
    })


def test_detect_columns():
    from tissueresolve.io import validation as v
    obs = _toy_obs()
    assert v.detect_broad_cell_type_col(obs) == "broad_cell_type"
    assert v.detect_fine_cell_type_col(obs, exclude="broad_cell_type") == "cell_type"


def test_validate_hierarchical_annotations_ok():
    from tissueresolve.io import validation as v
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = v.validate_hierarchical_annotations(_toy_obs(), "broad_cell_type", "cell_type")
    assert r["n_broad"] == 3 and r["n_fine"] == 5
    assert r["mapping"]["CD8 T cell"] == "T/NK"


def test_validate_missing_column_raises():
    from tissueresolve.io import validation as v
    with pytest.raises(KeyError):
        v.validate_hierarchical_annotations(_toy_obs(), "nope", "cell_type")


def test_validate_ambiguous_fine_raises():
    from tissueresolve.io import validation as v
    obs = pd.DataFrame({"broad_cell_type": ["A", "B"], "cell_type": ["x", "x"]})
    with pytest.raises(ValueError):
        v.validate_hierarchical_annotations(obs, "broad_cell_type", "cell_type")


def test_summarize_hierarchical_annotations():
    from tissueresolve.io import validation as v
    s = v.summarize_hierarchical_annotations(_toy_obs(), "broad_cell_type", "cell_type")
    assert set(s.columns) == {"broad_cell_type", "fine_cell_type", "n_cells",
                              "n_fine_in_family"}
    assert int(s.loc[s["broad_cell_type"] == "T/NK", "n_fine_in_family"].iloc[0]) == 2


# ---------------------------------------------------------------------------
# Bulk hierarchical workflow
# ---------------------------------------------------------------------------


def test_hierarchical_bulk_workflow(hier_ref):
    from tissueresolve.bulk.hierarchical import (
        run_hierarchical_bulk, save_hierarchical_bulk_outputs,
    )
    rng = np.random.default_rng(7)
    R = hier_ref.as_R_cpm()
    true = rng.dirichlet(np.ones(5), size=4)
    bulk_mat = true @ R
    bulk = pd.DataFrame(bulk_mat.T, index=hier_ref.gene_names,
                        columns=[f"s{i}" for i in range(4)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_hierarchical_bulk(bulk, hier_ref, FAMILY_MAP, n_bootstrap=0,
                                    min_discriminating_genes=5)
    # mass preserved including unresolved column
    np.testing.assert_allclose(res.deconv.proportions.sum(axis=1), 1.0, atol=1e-9)
    np.testing.assert_allclose(res.estimates.family_proportions.sum(axis=1), 1.0,
                               atol=1e-9)
    assert "FamY" in res.estimates.unresolved_families
    assert "unresolved_FamY" in res.deconv.proportions.columns


def test_hierarchical_bulk_outputs_written(hier_ref, tmp_path):
    from tissueresolve.bulk.hierarchical import (
        run_hierarchical_bulk, save_hierarchical_bulk_outputs,
    )
    R = hier_ref.as_R_cpm()
    bulk = pd.DataFrame(R.T, index=hier_ref.gene_names,
                        columns=[f"c{i}" for i in range(R.shape[0])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_hierarchical_bulk(bulk, hier_ref, FAMILY_MAP, n_bootstrap=0,
                                    min_discriminating_genes=5)
        written = save_hierarchical_bulk_outputs(res, tmp_path)
    for name in ("bulk_family_proportions.tsv",
                 "bulk_conditional_fine_proportions.tsv",
                 "bulk_hierarchical_fine_proportions.tsv",
                 "bulk_unresolved_family_mass.tsv",
                 "bulk_hierarchical_qc.tsv",
                 "cell_type_hierarchy.tsv"):
        assert (tmp_path / name).exists()


# ---------------------------------------------------------------------------
# Spatial hierarchical workflow (toy coordinates)
# ---------------------------------------------------------------------------


def test_hierarchical_spatial_workflow(hier_ref):
    from tissueresolve.spatial.hierarchical import run_hierarchical_spatial
    from tissueresolve.config import TissueResolveConfig

    rng = np.random.default_rng(11)
    n_spots = 24
    R = hier_ref.as_R_cpm()
    true = rng.dirichlet(np.ones(5), size=n_spots)
    lib = rng.integers(800, 1500, size=n_spots).astype("float32")
    lam = (true @ R) / R.sum(axis=1, keepdims=True).mean()
    Y = rng.poisson(np.clip(lam * lib[:, None] / lam.sum(axis=1, keepdims=True), 0, None))
    Y = Y.astype("float32")
    # simple hex-like grid coordinates
    rows = np.repeat(np.arange(4), 6)
    cols = np.tile(np.arange(6), 4)
    cfg = TissueResolveConfig()
    cfg.spatial_solver.max_iter = 10
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_hierarchical_spatial(
            Y, hier_ref, rows, cols, lib, list(hier_ref.gene_names), FAMILY_MAP,
            spot_ids=[f"spot{i}" for i in range(n_spots)],
            config=cfg, min_discriminating_genes=5)
    np.testing.assert_allclose(res.deconv.proportions.sum(axis=1), 1.0, atol=1e-6)
    assert "FamY" in res.estimates.unresolved_families
    # spatial smoothing parameter is recorded for both stages
    assert "lambda_spatial_family" in res.run_metadata
    assert "lambda_spatial_fine" in res.run_metadata


# ---------------------------------------------------------------------------
# Reproducible family-aware colour map
# ---------------------------------------------------------------------------


def test_color_map_family_aware_and_neutral():
    from tissueresolve.plotting import palette as p
    fine = ["CD8 T cell", "CD4 T cell", "macrophage", "unresolved_T/NK", "Other"]
    mp = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "macrophage": "Myeloid"}
    cm = p.build_hierarchical_color_map(fine, mp)
    by_label = dict(zip(cm["fine_cell_type"], cm["color_hex"]))
    assert by_label["Other"] == p.OTHER_COLOR
    assert by_label["unresolved_T/NK"] == p.UNRESOLVED_COLOR
    # CD4 and CD8 (same family) get different but related shades
    assert by_label["CD4 T cell"] != by_label["CD8 T cell"]


def test_color_map_reuse_preserves_existing():
    from tissueresolve.plotting import palette as p
    fine = ["CD8 T cell", "CD4 T cell"]
    mp = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK"}
    cm = p.build_hierarchical_color_map(fine, mp)
    old = dict(zip(cm["fine_cell_type"], cm["color_hex"]))
    cm2 = p.build_hierarchical_color_map(
        fine + ["NK cell"], {**mp, "NK cell": "T/NK"}, existing=cm)
    new = dict(zip(cm2["fine_cell_type"], cm2["color_hex"]))
    assert new["CD8 T cell"] == old["CD8 T cell"]


def test_color_map_save_json_tsv(tmp_path):
    from tissueresolve.plotting import palette as p
    cm = p.build_hierarchical_color_map(["macrophage"], {"macrophage": "Myeloid"})
    paths = p.save_hierarchical_color_map(cm, tmp_path)
    assert paths["tsv"].exists() and paths["json"].exists()


# ---------------------------------------------------------------------------
# Within-family marker refinement
# ---------------------------------------------------------------------------


def test_within_family_marker_panels(hier_ref):
    from tissueresolve.reference.pairwise_markers import (
        build_family_specific_gene_panels,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        panels, disc = build_family_specific_gene_panels(hier_ref, FAMILY_MAP, top_n=5)
    # multi-member families get panels; singleton FamZ does not
    assert "FamX" in panels and "FamY" in panels
    assert "FamZ" not in panels
    assert set(disc.columns) == {"family", "gene", "resolution_score",
                                 "best_pair", "abs_log2fc", "leakage"}


# ---------------------------------------------------------------------------
# Partial hierarchical resolution (Part 3)
# ---------------------------------------------------------------------------


def _partial_ref():
    rng = np.random.default_rng(0)
    G = 120
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0); v[list(a)] = lvl; return v
    X1 = prof(range(0, 20)); X2 = prof(range(20, 40))          # separable family
    P1 = prof(range(40, 60))                                    # distinct in FamP
    P2 = prof(range(60, 80)); P3 = P2 * (1 + rng.normal(0, 0.004, G))  # near-identical
    R = np.vstack([X1, X2, P1, P2, P3]).astype(np.float32)
    cts = ["X1", "X2", "P1", "P2", "P3"]
    ref = ReferenceSignature(gene_names=genes, cell_types=cts, R_cpm=R,
                             R_log=np.log1p(R).astype(np.float32),
                             n_cells_per_type={c: 100 for c in cts})
    mapping = {"X1": "FamX", "X2": "FamX", "P1": "FamP", "P2": "FamP", "P3": "FamP"}
    return ref, mapping


def test_partial_resolution_splits_confident_keeps_residual():
    from tissueresolve.reference.hierarchy import assemble_hierarchical_estimates
    ref, mapping = _partial_ref()
    fam = pd.DataFrame({"FamX": [0.5], "FamP": [0.5]}, index=["s0"])
    fine = pd.DataFrame({"X1": [0.3], "X2": [0.2], "P1": [0.2], "P2": [0.15],
                         "P3": [0.15]}, index=["s0"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = assemble_hierarchical_estimates(
            fam, fine, ref, mapping, min_discriminating_genes=5,
            allow_partial_resolution=True)
    c = est.combined_fine
    # separable family fully resolved
    assert c["X1"].iloc[0] > 0 and c["X2"].iloc[0] > 0
    # distinct subtype in the hard family keeps confident mass...
    assert c["P1"].iloc[0] > 0
    # ...while the near-identical pair is abstained (residual to unresolved)
    assert c["P2"].iloc[0] == pytest.approx(0.0)
    assert c["P3"].iloc[0] == pytest.approx(0.0)
    assert "unresolved_FamP" in c.columns and c["unresolved_FamP"].iloc[0] > 0
    # mass preserved
    np.testing.assert_allclose(c.sum(axis=1), 1.0, atol=1e-9)


def test_partial_resolution_mass_preserved_separable_family_no_unresolved():
    from tissueresolve.reference.hierarchy import assemble_hierarchical_estimates
    ref, mapping = _partial_ref()
    fam = pd.DataFrame({"FamX": [1.0], "FamP": [0.0]}, index=["s0"])
    fine = pd.DataFrame({"X1": [0.6], "X2": [0.4], "P1": [0.0], "P2": [0.0],
                         "P3": [0.0]}, index=["s0"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = assemble_hierarchical_estimates(
            fam, fine, ref, mapping, min_discriminating_genes=5,
            allow_partial_resolution=True)
    # FamX is fully separable → its unresolved residual is ~0
    fam_x_unresolved = est.combined_fine.get("unresolved_FamX")
    if fam_x_unresolved is not None:
        assert fam_x_unresolved.iloc[0] == pytest.approx(0.0, abs=1e-9)
    np.testing.assert_allclose(est.combined_fine.sum(axis=1), 1.0, atol=1e-9)
