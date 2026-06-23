"""Tests for experimental diagnostic-driven auto-mode recommendation (advisory)."""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.auto_mode import recommend_mode, ModeRecommendation


def test_missing_evidence_is_conservative_bulk():
    rec = recommend_mode("bulk")
    assert rec.recommended_mode == "broad_only"
    assert rec.is_experimental and rec.reasons


def test_missing_evidence_is_conservative_spatial():
    rec = recommend_mode("spatial", reference_qc={"ok": 1}, query_qc={"ok": 1})
    assert rec.recommended_preset == "default"
    assert rec.reasons


def test_high_boundary_recommends_edge_aware():
    rec = recommend_mode(
        "spatial", reference_qc={"ok": 1}, query_qc={"ok": 1},
        spatial_diagnostics={"boundary_score": 0.8, "graph_quality": 1.0})
    assert rec.recommended_preset == "edge_aware_smoothing"
    assert rec.is_experimental


def test_strong_rare_niche_avoids_aggressive_smoothing():
    rec = recommend_mode(
        "spatial", reference_qc={"ok": 1}, query_qc={"ok": 1},
        spatial_diagnostics={"boundary_score": 0.1, "rare_niche_signal": 0.9,
                             "graph_quality": 1.0})
    assert rec.recommended_preset == "weak_smoothing"


def test_poor_fine_reliability_recommends_broad():
    rec = recommend_mode(
        "bulk", reference_qc={"spillover_risk": 0.6}, query_qc={"ok": 1},
        resolution_decisions={"F1": "broad_only", "F2": "broad_only", "F3": "selected_fine"})
    assert rec.recommended_mode == "broad_only"
    assert any("spillover" in w for w in rec.warnings)


def test_user_mode_not_overridden():
    rec = recommend_mode(
        "spatial", reference_qc={"ok": 1}, query_qc={"ok": 1},
        spatial_diagnostics={"boundary_score": 0.9, "graph_quality": 1.0},
        config={"user_preset": "default", "user_mode": "none"})
    assert rec.recommended_preset == "default"   # user choice kept
    assert rec.recommended_mode == "none"
    assert any("kept" in w for w in rec.warnings)


def test_reasons_recorded():
    rec = recommend_mode("spatial", reference_qc={"ok": 1}, query_qc={"ok": 1},
                         spatial_diagnostics={"oversmoothing_risk": 0.7, "graph_quality": 1.0})
    assert isinstance(rec, ModeRecommendation) and len(rec.reasons) >= 1
