"""Tests for benchmark statistical helpers (PART 14 / PART 18). Offline, deterministic."""
import numpy as np
import pandas as pd
import pytest

from benchmarks.shared import stats as S


def test_bootstrap_ci_brackets_point_and_is_reproducible():
    rng = np.random.default_rng(0)
    vals = rng.normal(0.8, 0.05, 40)
    c1 = S.bootstrap_ci(vals, seed=1)
    c2 = S.bootstrap_ci(vals, seed=1)
    assert c1 == c2                          # deterministic
    assert c1["lo"] <= c1["point"] <= c1["hi"]
    assert c1["n"] == 40


def test_bootstrap_ci_handles_nan_and_empty():
    assert S.bootstrap_ci([np.nan, np.nan])["n"] == 0
    c = S.bootstrap_ci([1.0, np.nan, 1.0, 1.0])
    assert c["n"] == 3 and c["point"] == pytest.approx(1.0)


def test_paired_wilcoxon_detects_consistent_difference():
    a = np.array([0.9, 0.85, 0.92, 0.88, 0.91, 0.87])
    b = a - 0.1                              # a consistently larger
    r = S.paired_wilcoxon(a, b)
    assert r["p_value"] < 0.05
    assert r["effect_size"] == pytest.approx(1.0)   # a always > b


def test_paired_wilcoxon_identical_arms():
    a = [0.5, 0.6, 0.7]
    r = S.paired_wilcoxon(a, a)
    assert r["effect_size"] == 0.0


def test_effect_size_sign():
    a = np.array([1.0, 2.0, 3.0, 4.0])
    b = np.array([0.5, 2.5, 2.0, 3.0])       # a > b in 3 of 4
    assert S.matched_pairs_effect_size(a, b) > 0
    assert S.matched_pairs_effect_size(b, a) < 0


def test_rank_stability_picks_clear_winner():
    # method "good" beats "bad" on every replicate
    df = pd.DataFrame({"good": [0.9, 0.88, 0.91, 0.92],
                       "bad": [0.4, 0.45, 0.42, 0.41]})
    rs = S.rank_stability(df, higher_is_better=True, seed=0)
    assert rs.loc["good", "prob_best"] == pytest.approx(1.0)
    assert rs.loc["good", "mean_rank"] < rs.loc["bad", "mean_rank"]
