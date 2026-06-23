"""Tests for experimental edge-aware + family-specific spatial smoothing.

Fast, offline, deterministic. Verifies edge weights are finite/bounded, boundary
edges get lower weights, family lambdas are bounded with fine ≤ broad, rare/immune
families are smoothed less than structural ones, predictions stay valid
compositions, metadata is recorded, and the package default lambda is unchanged.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.spatial_adaptive_smoothing import (
    compute_edge_aware_spatial_weights, apply_edge_aware_smoothing,
    compute_family_lambda, family_adaptive_smoothing,
)
from tissueresolve.config import TissueResolveConfig


def _toy_two_domain(n_side=6, k=3, seed=0):
    """Two domains split at the vertical midline with distinct expression/composition."""
    rng = np.random.default_rng(seed)
    rows, cols, expr, prop, dom = [], [], [], [], []
    for r in range(n_side):
        for c in range(n_side):
            left = c < n_side // 2
            base_expr = np.r_[np.ones(5), np.zeros(5)] if left else np.r_[np.zeros(5), np.ones(5)]
            expr.append(base_expr * 10 + rng.normal(0, 0.1, 10))
            prop.append([0.8, 0.2] if left else [0.2, 0.8])
            rows.append(r); cols.append(c); dom.append("L" if left else "R")
    coords = np.c_[rows, cols].astype(float)
    spots = [f"s{i}" for i in range(len(rows))]
    X = np.array(expr)
    P = pd.DataFrame(np.array(prop), index=spots, columns=["A", "B"])
    return coords, X, P, np.array(dom)


def test_edge_weights_finite_bounded():
    coords, X, P, dom = _toy_two_domain()
    ew = compute_edge_aware_spatial_weights(X, coords, initial_proportions=P,
                                            min_edge_weight=0.1, max_edge_weight=1.0)
    assert np.all(np.isfinite(ew.edge_weights))
    assert ew.edge_weights.min() >= 0.1 - 1e-9 and ew.edge_weights.max() <= 1.0 + 1e-9
    assert np.all(np.isfinite(ew.spot_smoothness))


def test_cross_boundary_edges_lower_weight():
    coords, X, P, dom = _toy_two_domain()
    ew = compute_edge_aware_spatial_weights(X, coords, initial_proportions=P)
    cross, within = [], []
    for i in range(coords.shape[0]):
        for jj, j in enumerate(ew.neighbors[i]):
            (cross if dom[i] != dom[j] else within).append(ew.edge_weights[i, jj])
    assert np.mean(cross) < np.mean(within)   # boundaries smoothed less


def test_edge_aware_smoothing_valid_compositions_and_metadata():
    coords, X, P, dom = _toy_two_domain()
    out, meta = apply_edge_aware_smoothing(P, coords, expression=X, base_lambda=0.1)
    assert not np.any(np.isnan(out.to_numpy()))
    assert np.allclose(out.sum(axis=1), 1.0)
    assert (out.to_numpy() >= -1e-9).all()
    assert meta["adaptive_smoothing_used"] and meta["edge_aware"]
    assert meta["base_lambda"] == 0.1 and meta["default_lambda_unchanged"]


def test_default_lambda_unchanged_by_module_import():
    assert TissueResolveConfig().spatial_solver.lambda_spatial == 0.1


def test_family_lambda_bounded_and_classes():
    lam_struct, _ = compute_family_lambda("Endothelial", base_lambda=0.1,
                                          min_lambda=0.01, max_lambda=0.1)
    lam_immune, _ = compute_family_lambda("Lymphoid", base_lambda=0.1,
                                          min_lambda=0.01, max_lambda=0.1)
    assert 0.01 <= lam_immune <= lam_struct <= 0.1
    assert lam_immune < lam_struct          # immune/rare smoothed less than structural


def test_family_lambda_rare_collinear_reduces():
    lam, rule = compute_family_lambda("Epithelial",
                                      family_metrics={"signature_collinearity": 0.95},
                                      base_lambda=0.1, min_lambda=0.01, max_lambda=0.1)
    assert lam < 0.1 and "reduce" in rule


def test_family_adaptive_smoothing_fine_le_broad_and_valid():
    coords, X, P, dom = _toy_two_domain()
    P3 = P.copy()
    P3.columns = ["Endothelial_1", "Lymphoid_1"]
    mapping = {"Endothelial_1": "Endothelial", "Lymphoid_1": "Lymphoid"}
    out, meta = family_adaptive_smoothing(P3, coords, mapping, base_lambda=0.1)
    assert np.allclose(out.sum(axis=1), 1.0) and not np.any(np.isnan(out.to_numpy()))
    fl = meta["family_lambdas"]
    assert fl["Lymphoid"] <= fl["Endothelial"]   # immune ≤ structural
    assert all(0.01 <= v <= 0.1 for v in fl.values())
    assert meta["default_lambda_unchanged"]
