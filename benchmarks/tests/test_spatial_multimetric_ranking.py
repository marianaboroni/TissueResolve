"""Offline tests for the scenario-separated spatial multi-metric ranking."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from benchmarks.shared import spatial_multimetric_ranking as R


# --- primitives -------------------------------------------------------------

def test_jsd_identical_is_zero_and_symmetric():
    p = [0.2, 0.3, 0.5]
    assert R.jensen_shannon_divergence(p, p) == pytest.approx(0.0, abs=1e-9)
    q = [0.5, 0.3, 0.2]
    assert R.jensen_shannon_divergence(p, q) == pytest.approx(
        R.jensen_shannon_divergence(q, p))


def test_jsd_handles_zeros():
    # zeros must not produce nan/inf via log(0)
    v = R.jensen_shannon_divergence([0.0, 1.0, 0.0], [0.5, 0.5, 0.0])
    assert np.isfinite(v) and 0.0 <= v <= 1.0


def test_jsd_bounded_01():
    assert R.jensen_shannon_divergence([1, 0], [0, 1]) == pytest.approx(1.0, abs=1e-9)


def test_normalize_direction_inverts():
    s = pd.Series({"a": 1.0, "b": 2.0, "c": 3.0})
    hi = R.normalize_metric_direction(s, higher_is_better=True)
    lo = R.normalize_metric_direction(s, higher_is_better=False)
    assert hi["c"] == pytest.approx(1.0) and hi["a"] == pytest.approx(0.0)
    assert lo["a"] == pytest.approx(1.0) and lo["c"] == pytest.approx(0.0)


def test_normalize_constant_column():
    s = pd.Series({"a": 2.0, "b": 2.0})
    out = R.normalize_metric_direction(s)
    assert (out == 1.0).all()


def test_marker_recovery_on_toy_data():
    spots = [f"s{i}" for i in range(20)]
    rng = np.random.default_rng(0)
    a = pd.Series(rng.random(20), index=spots)
    expr = pd.DataFrame({"GENE_A": a * 5 + 0.1, "GENE_B": rng.random(20)},
                        index=spots)
    pred = pd.DataFrame({"TypeA": a, "TypeB": rng.random(20)}, index=spots)
    rec = R.marker_recovery_score(pred, expr, {"TypeA": ["GENE_A"], "TypeB": ["GENE_B"]})
    # TypeA abundance tracks its marker → high recovery
    assert rec["TypeA"] > 0.9


# --- scoring / ranking ------------------------------------------------------

_REAL = pd.DataFrame({
    "marker_recovery": [0.8, 0.6, 0.4, np.nan],
    "method_concordance": [0.7, 0.65, 0.5, np.nan],
    "spatial_structure": [0.5, 0.4, 0.3, np.nan],
    "stability": [0.9, 0.8, 0.7, np.nan],
    "runtime_completeness": [0.6, 0.9, 0.5, np.nan],
    "interpretability": [0.8, 0.6, 0.7, np.nan],
    "status": ["executed", "imported", "executed", "skipped"],
}, index=["CARD", "cell2location", "SPOTlight", "RCTD"])

_SYNTH = pd.DataFrame({
    "rmse": [0.04, 0.06, 0.10],
    "jsd": [0.10, 0.15, 0.30],
    "pearson": [0.9, 0.8, 0.6],
    "dominant_accuracy": [0.85, 0.75, 0.6],
    "status": ["executed", "executed", "skipped"],
}, index=["CARD", "cell2location", "SPOTlight"])


def test_real_score_excludes_accuracy_and_skipped():
    out = R.compute_real_spatial_score(_REAL)
    assert "real_score" in out.columns
    # skipped method gets no score
    assert np.isnan(out.loc["RCTD", "real_score"])
    # no accuracy dimension leaks in
    assert not any("rmse" in c or "pearson" in c for c in out.columns)
    assert out.loc["CARD", "real_score"] >= out.loc["SPOTlight", "real_score"]


def test_synthetic_score_uses_accuracy():
    out = R.compute_synthetic_spatial_accuracy_score(_SYNTH)
    assert "synthetic_accuracy_score" in out.columns
    assert "inverse_error_score" in out.columns and "correlation_score" in out.columns
    assert np.isnan(out.loc["SPOTlight", "synthetic_accuracy_score"])  # skipped
    # CARD (low error, high corr) outranks cell2location
    assert out.loc["CARD", "synthetic_accuracy_score"] >= out.loc["cell2location", "synthetic_accuracy_score"]


def test_rank_excludes_skipped_marks_imported():
    out = R.compute_real_spatial_score(_REAL)
    ranked = R.rank_spatial_methods(out, score_col="real_score")
    assert np.isnan(ranked.loc["RCTD", "rank"])          # skipped not ranked
    assert bool(ranked.loc["cell2location", "imported"]) is True
    assert (ranked["rank"].dropna() >= 1).all()


def test_one_method_no_ranking_claim():
    one = _REAL.iloc[[0]]
    out = R.compute_real_spatial_score(one)
    ranked = R.rank_spatial_methods(out, score_col="real_score")
    assert "no ranking claimed" in ranked["ranking_claim"].iloc[0]


def test_overall_scorecard_labeled():
    real = R.compute_real_spatial_score(_REAL)
    synth = R.compute_synthetic_spatial_accuracy_score(_SYNTH)
    card = R.compute_spatial_overall_scorecard(real, synth)
    assert "scorecard" in card.columns
    assert "scorecard_note" in card.columns
    assert "NOT objective truth" in card["scorecard_note"].iloc[0]


def test_explain_ranking_states_no_accuracy_for_real():
    out = R.compute_real_spatial_score(_REAL)
    ranked = R.rank_spatial_methods(out, score_col="real_score")
    md = R.explain_spatial_ranking(ranked, score_col="real_score",
                                   scenario="real spatial (no ground truth)",
                                   weights=R.DEFAULT_REAL_WEIGHTS)
    assert "not accuracy" in md.lower()
    assert "| rank | method |" in md
