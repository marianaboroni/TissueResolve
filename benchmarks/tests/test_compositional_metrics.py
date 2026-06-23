"""Tests for compositional distances + complexity panel (PART 5/6/10/12/18).

Offline, deterministic.  Explicitly covers the PART-18 edge cases: JSD with
zero entries, Aitchison distance with zero replacement, and the richness /
entropy / effective-N metrics.
"""
import numpy as np
import pandas as pd
import pytest

from benchmarks.shared import metrics as M


# --- Jensen-Shannon divergence ----------------------------------------------
def test_jsd_identical_is_zero():
    p = np.array([0.5, 0.3, 0.2])
    assert M.jensen_shannon_divergence(p, p) == pytest.approx(0.0, abs=1e-12)


def test_jsd_disjoint_support_is_one_base2():
    p = np.array([1.0, 0.0, 0.0])
    q = np.array([0.0, 1.0, 0.0])
    assert M.jensen_shannon_divergence(p, q) == pytest.approx(1.0, rel=1e-9)


def test_jsd_handles_zeros_without_nan_or_inf():
    p = np.array([0.0, 0.7, 0.3, 0.0])
    q = np.array([0.1, 0.6, 0.0, 0.3])
    v = M.jensen_shannon_divergence(p, q)
    assert np.isfinite(v) and 0.0 <= v <= 1.0


def test_jsd_symmetric():
    p = np.array([0.2, 0.5, 0.3])
    q = np.array([0.4, 0.1, 0.5])
    assert M.jensen_shannon_divergence(p, q) == pytest.approx(
        M.jensen_shannon_divergence(q, p), rel=1e-12)


# --- total variation ---------------------------------------------------------
def test_total_variation_bounds():
    assert M.total_variation_distance([1, 0, 0], [1, 0, 0]) == pytest.approx(0.0)
    assert M.total_variation_distance([1, 0, 0], [0, 1, 0]) == pytest.approx(1.0)


# --- Aitchison distance ------------------------------------------------------
def test_aitchison_zero_replacement_finite():
    # zeros would make clr undefined; multiplicative replacement keeps it finite
    p = np.array([0.0, 0.5, 0.5])
    q = np.array([0.5, 0.5, 0.0])
    d = M.aitchison_distance(p, q)
    assert np.isfinite(d) and d > 0.0


def test_aitchison_identical_is_zero():
    p = np.array([0.2, 0.3, 0.5])
    assert M.aitchison_distance(p, p) == pytest.approx(0.0, abs=1e-9)


def test_aitchison_scale_invariant():
    p = np.array([0.2, 0.3, 0.5])
    q = np.array([0.1, 0.6, 0.3])
    assert M.aitchison_distance(p, q) == pytest.approx(
        M.aitchison_distance(p * 10, q * 3), rel=1e-6)


# --- concordance correlation -------------------------------------------------
def test_ccc_perfect_agreement():
    t = np.array([0.1, 0.4, 0.2, 0.3])
    assert M.concordance_correlation_coefficient(t, t) == pytest.approx(1.0, rel=1e-9)


def test_ccc_penalises_scale_shift_more_than_pearson():
    t = np.array([0.1, 0.4, 0.2, 0.3])
    e = t * 0.5  # perfectly correlated but scaled -> CCC < 1
    assert M.concordance_correlation_coefficient(t, e) < 0.999


# --- complexity panel --------------------------------------------------------
def test_effective_n_uniform_and_onehot():
    assert M.effective_n_populations([0.25] * 4) == pytest.approx(4.0, rel=1e-6)
    assert M.effective_n_populations([1, 0, 0, 0]) == pytest.approx(1.0, rel=1e-6)


def test_gini_equal_vs_concentrated():
    assert M.gini_coefficient([1, 1, 1, 1]) == pytest.approx(0.0, abs=1e-9)
    assert M.gini_coefficient([0.97, 0.01, 0.01, 0.01]) > 0.6


def test_composition_complexity_dataframe_shape_and_values():
    df = pd.DataFrame(
        {"a": [0.7, 0.25], "b": [0.2, 0.25], "c": [0.1, 0.25], "d": [0.0, 0.25]},
        index=["s1", "s2"])
    out = M.composition_complexity(df)
    assert list(out.index) == ["s1", "s2"]
    assert out.loc["s1", "n_gt_0"] == 3
    assert out.loc["s2", "n_gt_0"] == 4
    assert out.loc["s2", "effective_n_populations"] == pytest.approx(4.0, rel=1e-6)
    assert out.loc["s1", "dominant_fraction"] == pytest.approx(0.7, rel=1e-6)


def test_compositional_metrics_bundle():
    true = pd.DataFrame({"a": [0.6, 0.2], "b": [0.4, 0.8]}, index=["s1", "s2"])
    est = pd.DataFrame({"a": [0.5, 0.3], "b": [0.5, 0.7]}, index=["s1", "s2"])
    m = M.compositional_metrics(true, est)
    for k in ("jsd_mean", "tv_mean", "aitchison_mean", "ccc"):
        assert k in m and np.isfinite(m[k])
    assert 0.0 <= m["jsd_mean"] <= 1.0
