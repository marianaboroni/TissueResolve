"""Deterministic unit tests for the TCGA sparsity-audit complexity metrics.

Covers PART 18 items: richness / entropy / effective-N / Gini / complexity-row
behaviour on known compositions. Offline, no real data, no network.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

# Load the audit module by path (benchmarks/audit is not an importable package).
_AUDIT = Path(__file__).resolve().parents[1] / "audit" / "tcga_sparsity_audit.py"
_spec = importlib.util.spec_from_file_location("tcga_sparsity_audit", _AUDIT)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)


def test_uniform_composition_is_maximally_complex():
    k = 8
    v = np.full(k, 1.0 / k)
    r = mod.complexity_row(v, n_total=k)
    assert r["n_gt_0"] == k
    assert r["effective_n_populations"] == pytest.approx(k, rel=1e-6)  # exp(H)=K
    assert r["gini"] == pytest.approx(0.0, abs=1e-9)                   # perfect equality
    assert r["dominant_fraction"] == pytest.approx(1.0 / k, rel=1e-6)
    assert r["shannon_entropy_nats"] == pytest.approx(np.log(k), rel=1e-6)


def test_one_hot_composition_is_maximally_sparse():
    k = 8
    v = np.zeros(k)
    v[0] = 1.0
    r = mod.complexity_row(v, n_total=k)
    assert r["n_gt_0"] == 1
    assert r["effective_n_populations"] == pytest.approx(1.0, rel=1e-6)
    assert r["dominant_fraction"] == pytest.approx(1.0, rel=1e-6)
    assert r["top1_cumulative"] == pytest.approx(1.0, rel=1e-6)
    assert r["shannon_entropy_nats"] == pytest.approx(0.0, abs=1e-9)


def test_threshold_counts_are_monotone_nonincreasing():
    v = np.array([0.5, 0.3, 0.05, 0.004, 0.0005, 0.00005, 0.0])
    r = mod.complexity_row(v, n_total=v.size)
    counts = [r["n_gt_0"], r["n_gt_0.0001"], r["n_gt_0.001"],
              r["n_gt_0.005"], r["n_gt_0.01"]]
    assert counts == sorted(counts, reverse=True)
    assert r["n_gt_0"] == 6          # six strictly-positive entries
    assert r["n_gt_0.01"] == 3       # 0.5, 0.3, 0.05


def test_effective_n_between_one_and_k_and_normalization_invariant():
    rng = np.random.default_rng(0)
    v = rng.random(12)
    r1 = mod.complexity_row(v, n_total=12)
    r2 = mod.complexity_row(v * 7.0, n_total=12)   # scale-invariant after L1 norm
    assert 1.0 <= r1["effective_n_populations"] <= 12.0
    assert r1["effective_n_populations"] == pytest.approx(
        r2["effective_n_populations"], rel=1e-9)


def test_gini_unequal_greater_than_equal():
    eq = mod._gini(np.array([1.0, 1.0, 1.0, 1.0]))
    uneq = mod._gini(np.array([0.97, 0.01, 0.01, 0.01]))
    assert eq == pytest.approx(0.0, abs=1e-9)
    assert uneq > 0.6


def test_entropy_handles_zeros_and_empty():
    assert mod._entropy_nats(np.array([0.0, 0.0])) != mod._entropy_nats(np.array([0.0, 0.0])) \
        or np.isnan(mod._entropy_nats(np.array([0.0, 0.0])))  # nan on all-zero
    # zeros are ignored, not -inf
    h = mod._entropy_nats(np.array([0.5, 0.5, 0.0]))
    assert h == pytest.approx(np.log(2), rel=1e-9)
