"""Tests for the experimental discriminative-marker within-family split."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.discriminative_within_family import (
    discriminative_gene_weights, solve_conditional_split, refine_family_split)


def _family():
    # 3 states, 20 genes: genes 0-2 discriminate (differ across states), rest shared.
    rng = np.random.default_rng(0)
    shared = rng.random(17) * 10 + 5
    P = np.vstack([
        np.r_[10, 0, 0, shared],
        np.r_[0, 10, 0, shared],
        np.r_[0, 0, 10, shared],
    ]).astype(float)
    P = P / P.sum(1, keepdims=True)
    return P


def test_weights_upweight_discriminative_genes():
    P = _family()
    w = discriminative_gene_weights(P, power=1.0)
    assert w.shape == (P.shape[1],)
    # discriminating genes 0-2 should weigh more than shared genes 3+
    assert w[:3].mean() > w[3:].mean()
    assert np.isclose(w.mean(), 1.0, atol=1e-6)


def test_power_zero_is_uniform():
    P = _family()
    w = discriminative_gene_weights(P, power=0.0)
    assert np.allclose(w, 1.0)


def test_min_weight_floor_compresses_range():
    P = _family()
    w_lo = discriminative_gene_weights(P, power=4.0, min_weight=0.02)
    w_hi = discriminative_gene_weights(P, power=4.0, min_weight=0.5)
    # a higher floor raises the minimum weight and compresses the dynamic range
    assert (w_hi > 0).all() and (w_lo > 0).all()
    assert w_hi.min() > w_lo.min()
    assert (w_hi.max() / w_hi.min()) < (w_lo.max() / w_lo.min())


def test_solve_nonneg_and_conserves_mass():
    P = _family()
    true = np.array([0.5, 0.3, 0.2])
    b = true @ P
    w = discriminative_gene_weights(P, power=2.0)
    theta = solve_conditional_split(b, P.T, w, total_mass=1.0)
    assert (theta >= -1e-12).all()
    assert np.isclose(theta.sum(), 1.0, atol=1e-9)


def test_discriminative_recovers_better_on_diluted_signal():
    # build a near-collinear family: states differ only in 2 of many genes (diluted).
    rng = np.random.default_rng(1)
    G = 200
    shared = rng.random(G) * 20 + 5
    a = shared.copy(); a[0] += 8
    b_ = shared.copy(); b_[1] += 8
    P = np.vstack([a, b_]); P = P / P.sum(1, keepdims=True)
    true = np.array([0.7, 0.3])
    obs = true @ P
    uni = solve_conditional_split(obs, P.T, None, 1.0)
    disc = solve_conditional_split(obs, P.T, discriminative_gene_weights(P, power=2.0), 1.0)
    # discriminative should be at least as close to truth as uniform (usually closer)
    assert np.sum((disc - true) ** 2) <= np.sum((uni - true) ** 2) + 1e-9


def test_refine_family_split_metadata():
    P = _family()
    b = np.array([0.4, 0.4, 0.2]) @ P
    theta, meta = refine_family_split(b, P, family_mass=0.5, power=2.0)
    assert np.isclose(theta.sum(), 0.5, atol=1e-9)
    assert meta["n_states"] == 3 and meta["power"] == 2.0 and "weight_max" in meta
