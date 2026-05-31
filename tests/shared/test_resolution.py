"""
Tests for the resolution-aware layer (reference.resolution + hierarchy).

Synthetic, offline.  Covers classification (resolved vs unresolved), family
inference, the abstain/unresolved mode (mass preservation), hierarchical
aggregation (mass preservation), and the per-cell-type annotation.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.reference.hierarchy import (
    aggregate_reference_by_group,
    build_hierarchical_reference,
    decompose_within_family,
    infer_broad_groups_from_labels,
)
from tissueresolve.reference.resolution import (
    ResolutionConfig,
    annotate_resolution,
    apply_unresolved_mode,
    build_resolution_report,
    classify_resolvability,
    infer_cell_type_families,
)
from tissueresolve.reference.separability import compute_separability
from tissueresolve.results import ReferenceSignature


def _make_ref(R_cpm: np.ndarray, cell_types: list[str]) -> ReferenceSignature:
    K, G = R_cpm.shape
    R_cpm = R_cpm.astype(np.float32)
    return ReferenceSignature(
        gene_names=[f"G{i:03d}" for i in range(G)],
        cell_types=list(cell_types),
        R_cpm=R_cpm,
        R_log=np.log1p(R_cpm).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={ct: 100 for ct in cell_types},
    )


def _separable_ref():
    # Two disjoint blocks → highly separable.
    G = 20
    R = np.full((2, G), 5.0)
    R[0, :10] = 1000.0
    R[1, 10:] = 1000.0
    return _make_ref(R, ["TypeA", "TypeB"])


def _identical_ref():
    G = 20
    base = np.linspace(50, 500, G)
    R = np.vstack([base, base * 1.001])  # nearly identical
    return _make_ref(R, ["TypeA", "TypeB"])


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def test_classify_thresholds():
    cfg = ResolutionConfig()
    assert classify_resolvability(0.9, cfg) == "resolved"
    assert classify_resolvability(0.15, cfg) == "partially_resolved"
    assert classify_resolvability(0.05, cfg) == "poorly_resolved"
    assert classify_resolvability(0.0, cfg) == "unresolved"


def test_separable_pair_is_resolved():
    ref = _separable_ref()
    rep = compute_separability(ref)
    res = build_resolution_report(rep, ref.cell_types)
    assert res.pairwise.iloc[0]["resolvability"] == "resolved"


def test_identical_pair_is_unresolved():
    ref = _identical_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep = compute_separability(ref)
    res = build_resolution_report(rep, ref.cell_types)
    row = res.pairwise.iloc[0]
    assert row["resolvability"] == "unresolved"
    # And they cluster into one (unresolvable) family.
    fams = [f for f in res.families if f.size > 1]
    assert len(fams) == 1
    assert fams[0].resolvability == "unresolved"
    assert set(fams[0].members) == {"TypeA", "TypeB"}


def test_resolution_report_save(tmp_path):
    ref = _separable_ref()
    res = build_resolution_report(compute_separability(ref), ref.cell_types)
    res.save(tmp_path)
    assert (tmp_path / "pairwise_resolvability.tsv").exists()
    assert (tmp_path / "cell_type_families.tsv").exists()
    assert (tmp_path / "resolution_config.json").exists()


# ---------------------------------------------------------------------------
# Abstain / unresolved mode
# ---------------------------------------------------------------------------


def test_unresolved_mode_collapses_and_preserves_mass():
    ref = _identical_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = build_resolution_report(compute_separability(ref), ref.cell_types)
    props = pd.DataFrame({"TypeA": [0.6, 0.3], "TypeB": [0.4, 0.7]},
                         index=["s0", "s1"])
    out, info = apply_unresolved_mode(props, res)
    assert info["collapsed_families"]
    # collapsed into a single unresolved_* column; mass preserved per sample.
    assert any(c.startswith("unresolved_") for c in out.columns)
    np.testing.assert_allclose(out.sum(axis=1).to_numpy(),
                               props.sum(axis=1).to_numpy(), atol=1e-9)


def test_unresolved_mode_disabled_returns_input():
    ref = _identical_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = build_resolution_report(compute_separability(ref), ref.cell_types,
                                      ResolutionConfig(allow_unresolved=False))
    props = pd.DataFrame({"TypeA": [0.6], "TypeB": [0.4]}, index=["s0"])
    out, info = apply_unresolved_mode(props, res, cfg=res.config)
    assert info["allow_unresolved"] is False
    assert list(out.columns) == ["TypeA", "TypeB"]


def test_resolved_family_not_collapsed():
    ref = _separable_ref()
    res = build_resolution_report(compute_separability(ref), ref.cell_types)
    props = pd.DataFrame({"TypeA": [0.6], "TypeB": [0.4]}, index=["s0"])
    out, info = apply_unresolved_mode(props, res)
    assert not info["collapsed_families"]
    assert set(out.columns) == {"TypeA", "TypeB"}


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------


def test_annotate_resolution_recommends_level():
    ref = _identical_ref()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = build_resolution_report(compute_separability(ref), ref.cell_types)
    mean_est = pd.Series({"TypeA": 0.5, "TypeB": 0.5})
    ann = annotate_resolution(mean_est, res, modality="bulk")
    assert set(ann.index) == {"TypeA", "TypeB"}
    assert (ann["recommended_level"] == "family/broad").all()
    assert "family_resolvability" in ann.columns


# ---------------------------------------------------------------------------
# Hierarchy
# ---------------------------------------------------------------------------


def test_infer_broad_groups_keywords():
    m = infer_broad_groups_from_labels(
        ["CD8-positive, alpha-beta T cell", "naive B cell", "macrophage"])
    assert m["CD8-positive, alpha-beta T cell"] == "T cell"
    assert m["naive B cell"] == "B cell"
    assert m["macrophage"] == "Myeloid"


def test_aggregate_reference_by_group_preserves_genes():
    R = np.vstack([np.full(12, 100.0), np.full(12, 120.0), np.full(12, 900.0)])
    ref = _make_ref(R, ["CD8 T cell", "CD4 T cell", "macrophage"])
    mapping = infer_broad_groups_from_labels(ref.cell_types)
    broad = aggregate_reference_by_group(ref, mapping)
    assert broad.n_genes == ref.n_genes
    assert set(broad.cell_types) == {"T cell", "Myeloid"}
    broad.validate()


def test_build_hierarchical_reference_structure():
    R = np.vstack([np.full(10, 100.0), np.full(10, 120.0), np.full(10, 900.0)])
    ref = _make_ref(R, ["CD8 T cell", "CD4 T cell", "macrophage"])
    h = build_hierarchical_reference(ref)
    assert "broad" in h and "groups" in h
    assert sorted(h["groups"]["T cell"]) == ["CD4 T cell", "CD8 T cell"]


def test_decompose_within_family_preserves_total_mass():
    broad = pd.DataFrame({"T cell": [0.6, 0.5], "Myeloid": [0.4, 0.5]},
                         index=["s0", "s1"])
    cond = pd.DataFrame({
        "CD8 T cell": [0.5, 0.25], "CD4 T cell": [0.5, 0.75],
        "macrophage": [1.0, 1.0],
    }, index=["s0", "s1"])
    mapping = {"CD8 T cell": "T cell", "CD4 T cell": "T cell", "macrophage": "Myeloid"}
    out = decompose_within_family(broad, cond, mapping)
    abs_props = out["absolute_subtype_proportions"]
    # absolute subtype mass per sample equals broad total (no unresolved here).
    np.testing.assert_allclose(abs_props.sum(axis=1).to_numpy(), [1.0, 1.0], atol=1e-9)


def test_decompose_with_unresolved_family_holds_mass():
    broad = pd.DataFrame({"T cell": [0.6], "Myeloid": [0.4]}, index=["s0"])
    cond = pd.DataFrame({"CD8 T cell": [0.5], "CD4 T cell": [0.5],
                         "macrophage": [1.0]}, index=["s0"])
    mapping = {"CD8 T cell": "T cell", "CD4 T cell": "T cell", "macrophage": "Myeloid"}
    out = decompose_within_family(broad, cond, mapping, unresolved_groups=["T cell"])
    # T cell kept at broad level as unresolved mass; total preserved.
    total = (out["absolute_subtype_proportions"].sum(axis=1)
             + out["unresolved_family_mass"].sum(axis=1))
    np.testing.assert_allclose(total.to_numpy(), [1.0], atol=1e-9)
    assert out["unresolved_family_mass"].loc["s0", "unresolved_T cell"] == pytest.approx(0.6)
