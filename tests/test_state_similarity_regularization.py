"""Tests for experimental state-similarity / Redeconve-inspired regularization.

Covers the state-similarity graph (finite/bounded weights, no self-loops,
within-family-only, similarity ordering, isolated states), the within-family
state-regularized + sparsity refinement (valid compositions, broad-mass
conservation, lambda=0 baseline, effective-N reduction, rare protection, no
cross-family redistribution), and the diagnostic state grouping. Everything is
opt-in; bulk/default behaviour is untouched (covered in the config/CLI tests).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.state_similarity_regularization import (
    compute_state_similarity_graph,
    apply_state_regularized_refinement,
    apply_sparsity_aware_refinement,
    recommend_state_groups,
    StateSimilarityGraph,
)


# --------------------------------------------------------------------------
# Fixtures: 4 states in 2 families. Within family A, states are near-collinear;
# within family B, states are distinct. Families are biologically unrelated.
# --------------------------------------------------------------------------

CELL_TYPES = ["A1", "A2", "B1", "B2"]
FAMILY = {"A1": "A", "A2": "A", "B1": "B", "B2": "B"}


def _profiles():
    rng = np.random.default_rng(0)
    base_a = rng.random(20) * 10
    # A1, A2 nearly collinear (similar); B1, B2 distinct from each other and from A
    A1 = base_a + rng.normal(0, 0.05, 20)
    A2 = base_a + rng.normal(0, 0.05, 20)
    B1 = np.r_[np.zeros(10), rng.random(10) * 10]
    B2 = np.r_[rng.random(10) * 10, np.zeros(10)]
    return np.clip(np.vstack([A1, A2, B1, B2]), 0, None)


# ---- Phase 1: state-similarity graph ----

def test_graph_finite_and_bounded_weights():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY)
    assert isinstance(g, StateSimilarityGraph)
    assert np.all(np.isfinite(g.weights))
    assert (g.weights >= 0).all() and (g.weights <= 1.0 + 1e-9).all()


def test_graph_no_self_loops():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY)
    assert np.all(g.state_i != g.state_j)


def test_within_family_only_blocks_cross_family_edges():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY, within_family_only=True)
    for e in range(g.n_edges):
        a, b = g.cell_types[g.state_i[e]], g.cell_types[g.state_j[e]]
        assert FAMILY[a] == FAMILY[b]


def test_similar_states_connected_with_high_weight():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY, within_family_only=False)
    # A1 and A2 are near-collinear → should be an edge with high weight
    edge_w = {}
    for e in range(g.n_edges):
        a, b = g.cell_types[g.state_i[e]], g.cell_types[g.state_j[e]]
        edge_w[frozenset((a, b))] = g.weights[e]
    assert frozenset(("A1", "A2")) in edge_w
    assert edge_w[frozenset(("A1", "A2"))] > 0.9


def test_dissimilar_states_not_over_connected():
    # B1 vs B2 are anti-structured; with min_similarity they should not connect
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY,
                                       within_family_only=False, min_similarity=0.3)
    for e in range(g.n_edges):
        a, b = g.cell_types[g.state_i[e]], g.cell_types[g.state_j[e]]
        assert frozenset((a, b)) != frozenset(("B1", "B2"))


def test_isolated_states_allowed_and_recorded():
    # singleton families → no within-family edges → all isolated
    fam = {c: c for c in CELL_TYPES}
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, fam, within_family_only=True)
    assert g.n_edges == 0
    assert set(g.metadata["isolated_states"]) == set(CELL_TYPES)


def test_within_family_only_no_family_info_is_noop():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, None, within_family_only=True)
    assert g.n_edges == 0
    assert g.metadata["no_family_info"] is True


def test_graph_metadata_complete():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY)
    md = g.metadata
    for k in ("n_states", "n_edges", "within_family_only", "similarity_method",
              "weight_min", "weight_max", "weight_mean", "weight_median",
              "n_connected_components", "families", "highly_connected_states",
              "isolated_states"):
        assert k in md
    assert md["n_states"] == 4 and set(md["families"]) == {"A", "B"}


# ---- helpers for refinement: proportions frame with two families + unresolved ----

def _props(n_spots=8, seed=0):
    rng = np.random.default_rng(seed)
    # family A diffuse across A1,A2; family B concentrated; plus unresolved_A column
    A1 = rng.uniform(0.10, 0.20, n_spots)
    A2 = rng.uniform(0.10, 0.20, n_spots)
    B1 = rng.uniform(0.20, 0.30, n_spots)
    B2 = rng.uniform(0.02, 0.05, n_spots)
    unres_A = rng.uniform(0.05, 0.10, n_spots)
    df = pd.DataFrame({"A1": A1, "A2": A2, "B1": B1, "B2": B2, "unresolved_A": unres_A})
    df = df.div(df.sum(axis=1), axis=0)  # valid compositions
    return df


def _graph_for_props():
    return compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY, within_family_only=False)


# ---- Phase 2: state-regularized refinement ----

def test_refinement_valid_compositions_and_nonneg():
    df = _props()
    g = _graph_for_props()
    out, md = apply_state_regularized_refinement(df, g, FAMILY, lambda_state=0.2, lambda_sparse=0.05)
    P = out.to_numpy(float)
    assert (P >= -1e-12).all()
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-6)
    assert not np.any(np.isnan(P))


def test_refinement_lambda_zero_is_baseline():
    df = _props()
    g = _graph_for_props()
    out, _ = apply_state_regularized_refinement(df, g, FAMILY, lambda_state=0.0, lambda_sparse=0.0)
    pd.testing.assert_frame_equal(out, df)


def test_refinement_conserves_broad_family_mass():
    df = _props()
    g = _graph_for_props()
    out, md = apply_state_regularized_refinement(df, g, FAMILY, lambda_state=0.2, lambda_sparse=0.1)
    # family A mass = A1+A2+unresolved_A ; family B mass = B1+B2
    famA_before = df[["A1", "A2", "unresolved_A"]].sum(axis=1).to_numpy()
    famA_after = out[["A1", "A2", "unresolved_A"]].sum(axis=1).to_numpy()
    famB_before = df[["B1", "B2"]].sum(axis=1).to_numpy()
    famB_after = out[["B1", "B2"]].sum(axis=1).to_numpy()
    assert np.allclose(famA_before, famA_after, atol=1e-9)
    assert np.allclose(famB_before, famB_after, atol=1e-9)
    assert md["broad_mass_max_deviation"] < 1e-6


def test_refinement_preserves_unresolved_mass():
    df = _props()
    g = _graph_for_props()
    out, _ = apply_state_regularized_refinement(df, g, FAMILY, lambda_state=0.2, lambda_sparse=0.1)
    assert np.allclose(out["unresolved_A"].to_numpy(), df["unresolved_A"].to_numpy(), atol=1e-9)


def test_refinement_no_cross_family_redistribution():
    df = _props()
    g = _graph_for_props()
    out, _ = apply_state_regularized_refinement(df, g, FAMILY, lambda_state=0.3, lambda_sparse=0.2)
    # total family A mass and family B mass each unchanged → no cross-family flow
    assert np.allclose(df[["A1", "A2", "unresolved_A"]].sum(axis=1),
                       out[["A1", "A2", "unresolved_A"]].sum(axis=1), atol=1e-9)


def test_sparsity_reduces_effective_n_in_diffuse_family():
    # family A is diffuse (A1~A2). sparsity should concentrate → lower effective-N.
    df = _props()
    g = _graph_for_props()
    out, md = apply_state_regularized_refinement(df, g, FAMILY, lambda_state=0.0, lambda_sparse=0.5)

    def eff_n(sub):
        s = sub.div(sub.sum(axis=1), axis=0).to_numpy()
        return float(np.mean([1.0 / np.sum(r ** 2) for r in s]))

    assert eff_n(out[["A1", "A2"]]) <= eff_n(df[["A1", "A2"]]) + 1e-9
    assert md["mass_shifted_sparse"] > 0


# ---- Phase 3: sparsity-aware refinement + rare protection ----

def test_sparsity_unsupported_small_state_shrinks():
    df = _props()
    out, _ = apply_sparsity_aware_refinement(df, FAMILY, lambda_sparse=0.5)
    # B2 is the small/unsupported state in family B → should shrink relative to B1
    assert out["B2"].mean() <= df["B2"].mean() + 1e-9


def test_sparsity_protects_rare_state():
    df = _props()
    before = df["B2"].to_numpy().copy()
    out, md = apply_sparsity_aware_refinement(df, FAMILY, lambda_sparse=0.9,
                                              rare_protection=["B2"])
    # protected B2 must not be shrunk to a smaller within-family share
    shareB_before = df["B2"] / df[["B1", "B2"]].sum(axis=1)
    shareB_after = out["B2"] / out[["B1", "B2"]].sum(axis=1)
    assert (shareB_after.to_numpy() >= shareB_before.to_numpy() - 1e-9).all()
    assert "B2" in md["rare_protected"]


def test_sparsity_conserves_broad_mass_no_nan():
    df = _props()
    out, md = apply_sparsity_aware_refinement(df, FAMILY, lambda_sparse=0.5)
    assert not np.any(np.isnan(out.to_numpy()))
    assert md["broad_mass_max_deviation"] < 1e-6
    assert np.allclose(out.to_numpy().sum(axis=1), 1.0, atol=1e-6)


# ---- Phase 4: adaptive resolution / state grouping ----

def test_grouping_merges_similar_states():
    g = _graph_for_props()
    rec = recommend_state_groups(g, threshold=0.9, family_map=FAMILY)
    assert ["A1", "A2"] in rec["recommended_groups"]


def test_grouping_does_not_merge_dissimilar_states():
    g = _graph_for_props()
    rec = recommend_state_groups(g, threshold=0.9, family_map=FAMILY)
    flat = {s for grp in rec["recommended_groups"] for s in grp}
    # B1/B2 dissimilar → not grouped together
    for grp in rec["recommended_groups"]:
        assert set(grp) != {"B1", "B2"}


def test_grouping_respects_family_boundaries():
    g = compute_state_similarity_graph(_profiles(), CELL_TYPES, FAMILY, within_family_only=False)
    rec = recommend_state_groups(g, threshold=0.5, family_map=FAMILY)
    for grp in rec["recommended_groups"]:
        fams = {FAMILY[s] for s in grp}
        assert len(fams) == 1


def test_grouping_protects_rare_state():
    g = _graph_for_props()
    rec = recommend_state_groups(g, threshold=0.9, family_map=FAMILY, rare_protection=["A2"])
    for grp in rec["recommended_groups"]:
        assert "A2" not in grp
    assert rec["used_for_estimation"] is False


def test_grouping_metadata_recorded():
    g = _graph_for_props()
    rec = recommend_state_groups(g, threshold=0.9, family_map=FAMILY)
    for k in ("recommended_groups", "n_groups", "families_with_unresolved_fine_structure",
              "used_for", "reasons"):
        assert k in rec


# ---- config / preset / CLI (default behaviour unchanged) ----

def test_default_config_state_regularization_disabled():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    assert cfg.state_regularization.enabled is False
    assert cfg.spatial_solver.lambda_spatial == 0.1  # default untouched


def test_existing_presets_do_not_enable_state_regularization():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_presets import apply_spatial_preset
    for p in ("default", "weak_smoothing", "no_smoothing", "edge_aware_smoothing",
              "combined_weak_edge_smoothing"):
        cfg = TissueResolveConfig()
        apply_spatial_preset(cfg, p)
        assert cfg.state_regularization.enabled is False


def test_state_presets_enable_and_record_mode():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_presets import apply_spatial_preset
    expected = {
        "state_regularized_experimental": "state_regularized",
        "sparsity_state_regularized_experimental": "sparsity",
        "adaptive_resolution_experimental": "adaptive_resolution",
    }
    for preset, mode in expected.items():
        cfg = TissueResolveConfig()
        info = apply_spatial_preset(cfg, preset)
        assert cfg.state_regularization.enabled is True
        assert cfg.state_regularization.mode == mode
        assert info.experimental and info.state_regularization
        md = info.to_metadata()
        assert md["state_regularization_used"] is True
        assert md["state_regularization_mode"] == mode
        # spatial lambda stays at the package default for state presets
        assert cfg.spatial_solver.lambda_spatial == 0.1


def test_unknown_preset_raises():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.experimental.spatial_presets import apply_spatial_preset
    with pytest.raises(ValueError):
        apply_spatial_preset(TissueResolveConfig(), "no_such_preset")


def test_cli_accepts_state_presets():
    from tissueresolve.cli import spatial_run
    for p in spatial_run.params:
        if p.name == "spatial_preset":
            for preset in ("state_regularized_experimental",
                           "sparsity_state_regularized_experimental",
                           "adaptive_resolution_experimental"):
                assert preset in p.type.choices
            break
    else:
        raise AssertionError("spatial_preset option not found")


def test_config_yaml_roundtrip_state_regularization(tmp_path):
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    cfg.state_regularization.enabled = True
    cfg.state_regularization.mode = "sparsity"
    cfg.state_regularization.lambda_sparse = 0.02
    p = tmp_path / "cfg.yaml"
    cfg.to_yaml(p)
    loaded = TissueResolveConfig.from_yaml(p)
    assert loaded.state_regularization.enabled is True
    assert loaded.state_regularization.mode == "sparsity"
    assert loaded.state_regularization.lambda_sparse == 0.02


# ---- end-to-end pipeline: opt-in, default off, valid compositions ----

def test_pipeline_default_off_unchanged():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    from tissueresolve.config import TissueResolveConfig
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    res = SpatialPipeline(TissueResolveConfig()).run(Y, ref, r, c, lib, genes, marker_genes=genes)
    assert res.deconv.run_metadata.get("state_regularization_used") is False


def test_pipeline_state_regularized_valid_and_conserves_mass():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    from tissueresolve.config import TissueResolveConfig
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    cell_types = list(ref.cell_types)
    # toy fine->broad family map: split cell types into two families
    family_map = {ct: ("famA" if i % 2 == 0 else "famB") for i, ct in enumerate(cell_types)}
    cfg = TissueResolveConfig()
    cfg.state_regularization.enabled = True
    cfg.state_regularization.mode = "state_regularized"
    cfg.state_regularization.lambda_state = 0.2
    cfg.state_regularization.lambda_sparse = 0.05
    cfg.state_regularization.within_family_only = True
    res = SpatialPipeline(cfg).run(Y, ref, r, c, lib, genes, marker_genes=genes,
                                   family_map=family_map)
    P = res.deconv.proportions.to_numpy()
    assert not np.any(np.isnan(P)) and (P >= -1e-9).all()
    assert np.allclose(P.sum(1), 1.0, atol=1e-4)
    md = res.deconv.run_metadata
    assert md.get("state_regularization_used") is True
    assert md.get("state_reg_broad_mass_max_deviation", 0.0) < 1e-3


def test_pipeline_state_regularized_no_family_map_is_noop():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    from tissueresolve.config import TissueResolveConfig
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    cfg_off = TissueResolveConfig()
    base = SpatialPipeline(cfg_off).run(Y, ref, r, c, lib, genes, marker_genes=genes)
    cfg = TissueResolveConfig()
    cfg.state_regularization.enabled = True
    cfg.state_regularization.within_family_only = True
    with pytest.warns(UserWarning):
        res = SpatialPipeline(cfg).run(Y, ref, r, c, lib, genes, marker_genes=genes)
    # within_family_only + no family_map → no edges → estimates unchanged
    assert np.allclose(base.deconv.proportions.to_numpy(),
                       res.deconv.proportions.to_numpy(), atol=1e-6)


def test_bulk_behaviour_unchanged_by_state_config():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    before = cfg.bulk_solver.__dict__.copy()
    cfg.state_regularization.enabled = True
    assert cfg.bulk_solver.__dict__ == before
