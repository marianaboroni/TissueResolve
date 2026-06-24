"""Tests for in-solver edge-aware spatial smoothing (experimental, opt-in).

Covers the weighted-graph builder (finite/bounded weights, boundary attenuation,
zero-distance → high weight, fallback), the preset (default unchanged, marked
experimental, metadata), and an end-to-end run (non-negative, valid compositions,
no NaN). Bulk behaviour is unaffected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.config import TissueResolveConfig
from tissueresolve.spatial.graph import (
    compute_edge_aware_graph, build_edge_aware_spatial_graph,
    build_hex_graph_from_arrays, WeightedSpatialGraph,
)
from tissueresolve.experimental.spatial_presets import apply_spatial_preset


def _hex_two_domains(n_rows=6, n_cols=6):
    """Hex grid + an expression matrix split into two domains at the col midline."""
    rows, cols = [], []
    for r in range(n_rows):
        for c in range(n_cols):
            rows.append(r); cols.append(c * 2 + (r % 2))
    rows = np.array(rows, np.int32); cols = np.array(cols, np.int32)
    rng = np.random.default_rng(0)
    mid = cols.mean()
    expr = []
    for cc in cols:
        left = cc < mid
        base = np.r_[np.ones(6), np.zeros(6)] if left else np.r_[np.zeros(6), np.ones(6)]
        expr.append(base * 50 + rng.normal(0, 0.5, 12))
    return rows, cols, np.clip(np.array(expr), 0, None), (cols < mid)


def test_weighted_graph_finite_and_bounded():
    r, c, X, _ = _hex_two_domains()
    wg = compute_edge_aware_graph(X, array_row=r, array_col=c,
                                  min_edge_weight=0.05, max_edge_weight=1.0)
    assert isinstance(wg, WeightedSpatialGraph)
    assert np.all(np.isfinite(wg.weights))
    assert wg.weights.min() >= 0.05 - 1e-9 and wg.weights.max() <= 1.0 + 1e-9


def test_cross_boundary_edges_lower_weight():
    r, c, X, left = _hex_two_domains()
    wg = compute_edge_aware_graph(X, array_row=r, array_col=c)
    cross, within = [], []
    for e in range(wg.weights.size):
        i, j = wg.indices_i[e], wg.indices_j[e]
        (cross if left[i] != left[j] else within).append(wg.weights[e])
    assert np.mean(cross) < np.mean(within)


def test_zero_distance_neighbours_get_high_weight():
    # identical expression everywhere => all edges at max weight
    r, c, _, _ = _hex_two_domains()
    X = np.ones((len(r), 8), dtype=float) * 7.0
    wg = compute_edge_aware_graph(X, array_row=r, array_col=c, max_edge_weight=1.0)
    assert np.allclose(wg.weights, 1.0)


def test_build_weighted_spatial_graph_valid_AL():
    r, c, X, _ = _hex_two_domains()
    g = build_edge_aware_spatial_graph(r, c, X, spot_ids=[f"s{i}" for i in range(len(r))])
    n = len(r)
    assert g.A.shape == (n, n) and g.L.shape == (n, n)
    # row-normalized: rows with neighbours sum to 1
    rs = np.asarray(g.A.sum(1)).ravel()
    assert np.allclose(rs[g.degree > 0], 1.0, atol=1e-5)
    assert g.metadata.get("edge_aware") is True
    for k in ("edge_weight_min", "edge_weight_max", "edge_weight_mean", "edge_weight_median"):
        assert k in g.metadata and np.isfinite(g.metadata[k])


def test_fallback_to_unweighted_when_no_edges():
    # single isolated spot => no hex edges => fallback to standard graph
    r = np.array([0], np.int32); c = np.array([0], np.int32)
    X = np.ones((1, 5))
    g = build_edge_aware_spatial_graph(r, c, X, spot_ids=["s0"])
    assert g.metadata.get("fallback") is True


def test_default_config_and_preset_unchanged():
    assert TissueResolveConfig().spatial_solver.lambda_spatial == 0.1
    assert TissueResolveConfig().spatial_solver.edge_aware is False
    cfg = TissueResolveConfig()
    info = apply_spatial_preset(cfg, "default")
    assert cfg.spatial_solver.edge_aware is False and info.lambda_spatial == 0.1


def test_edge_aware_preset_experimental_and_metadata():
    cfg = TissueResolveConfig()
    info = apply_spatial_preset(cfg, "edge_aware_smoothing")
    assert cfg.spatial_solver.edge_aware is True
    assert info.experimental and info.edge_aware
    md = info.to_metadata()
    assert md["edge_aware_smoothing_used"] is True and md["spatial_preset"] == "edge_aware_smoothing"


# ---- end-to-end: predictions stay valid; bulk unaffected ----

def test_end_to_end_edge_aware_valid_compositions():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    cfg = TissueResolveConfig()
    cfg.spatial_solver.edge_aware = True
    res = SpatialPipeline(cfg).run(Y, ref, r, c, lib, genes, marker_genes=genes)
    P = res.deconv.proportions.to_numpy()
    assert not np.any(np.isnan(P))
    assert (P >= -1e-9).all()
    assert np.allclose(P.sum(1), 1.0, atol=1e-4)
    assert res.deconv.run_metadata.get("edge_aware_smoothing_used") is True


def test_bulk_behaviour_unchanged_by_edge_aware_config():
    # toggling the spatial edge_aware flag must not touch bulk config
    cfg = TissueResolveConfig()
    before = cfg.bulk_solver.__dict__.copy()
    cfg.spatial_solver.edge_aware = True
    assert cfg.bulk_solver.__dict__ == before


# ---- combined weak + edge-aware preset (experimental) ----

def test_combined_preset_sets_low_lambda_and_edge_aware():
    cfg = TissueResolveConfig()
    info = apply_spatial_preset(cfg, "combined_weak_edge_smoothing")
    assert cfg.spatial_solver.lambda_spatial == 0.02
    assert cfg.spatial_solver.edge_aware is True
    assert info.experimental and info.combined and info.edge_aware
    assert info.min_edge_weight == 0.05


def test_combined_preset_metadata_records_combined():
    md = apply_spatial_preset(TissueResolveConfig(), "combined_weak_edge_smoothing").to_metadata()
    assert md["combined_preset"] is True
    assert md["spatial_preset"] == "combined_weak_edge_smoothing"
    assert md["lambda_spatial"] == 0.02 and md["edge_aware_smoothing_used"] is True


def test_existing_presets_unchanged_by_combined_addition():
    c1 = TissueResolveConfig(); apply_spatial_preset(c1, "weak_smoothing")
    assert c1.spatial_solver.lambda_spatial == 0.02 and c1.spatial_solver.edge_aware is False
    c2 = TissueResolveConfig(); apply_spatial_preset(c2, "edge_aware_smoothing")
    assert c2.spatial_solver.lambda_spatial == 0.05 and c2.spatial_solver.edge_aware is True
    c3 = TissueResolveConfig(); i3 = apply_spatial_preset(c3, "default")
    assert c3.spatial_solver.lambda_spatial == 0.1 and c3.spatial_solver.edge_aware is False
    assert i3.combined is False


def test_unknown_preset_raises():
    import pytest as _pt
    with _pt.raises(ValueError):
        apply_spatial_preset(TissueResolveConfig(), "no_such_preset")


def test_cli_accepts_combined_preset():
    from tissueresolve.cli import spatial_run
    for p in spatial_run.params:
        if p.name == "spatial_preset":
            assert "combined_weak_edge_smoothing" in p.type.choices
            break
    else:
        raise AssertionError("spatial_preset option not found")


def test_end_to_end_combined_valid_compositions():
    from tests.spatial.test_pipeline import _make_synthetic_dataset
    from tissueresolve.spatial.pipeline import SpatialPipeline
    Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
    cfg = TissueResolveConfig()
    apply_spatial_preset(cfg, "combined_weak_edge_smoothing")
    res = SpatialPipeline(cfg).run(Y, ref, r, c, lib, genes, marker_genes=genes)
    P = res.deconv.proportions.to_numpy()
    assert not np.any(np.isnan(P)) and (P >= -1e-9).all()
    assert np.allclose(P.sum(1), 1.0, atol=1e-4)
    assert res.deconv.run_metadata.get("edge_aware_smoothing_used") is True
    assert abs(res.deconv.lambda_spatial - 0.02) < 1e-9
