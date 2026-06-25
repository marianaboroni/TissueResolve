"""Tests for adaptive resolution reporting (P4). Reporting/QC only — never alters
estimates. Collinear subtypes are grouped/diagnostic, separable ones stay fine,
broad-only families stay broad-only, and rare protected states are not collapsed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.results import ReferenceSignature
from tissueresolve.experimental.adaptive_resolution import (
    build_adaptive_resolution_report, classify_family_resolution, AdaptiveResolutionReport)


def _ref():
    """4 types in 2 families: famA = {A1≈A2 collinear}, famB = {B1,B2 distinct}."""
    rng = np.random.default_rng(0)
    G = 200
    base = rng.random(G) * 50
    A1 = base + rng.normal(0, 0.02, G)
    A2 = base + rng.normal(0, 0.02, G)          # near-identical to A1
    B1 = np.r_[np.zeros(100), rng.random(100) * 50]
    B2 = np.r_[rng.random(100) * 50, np.zeros(100)]   # disjoint from B1
    R = np.clip(np.vstack([A1, A2, B1, B2]), 1e-3, None)
    R_cpm = (R / R.sum(1, keepdims=True) * 1e6).astype(np.float32)
    genes = [f"g{i}" for i in range(G)]
    cts = ["A1", "A2", "B1", "B2"]
    ref = ReferenceSignature(
        gene_names=genes, cell_types=cts,
        phi=(R_cpm.T / R_cpm.T.sum(0, keepdims=True)).astype(float),
        R_cpm=R_cpm, R_log=np.log1p(R_cpm).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={c: 100 for c in cts}, genome="hg38")
    mapping = {"A1": "famA", "A2": "famA", "B1": "famB", "B2": "famB"}
    return ref, mapping


def test_report_table_produced_with_schema():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping)
    assert isinstance(rep, AdaptiveResolutionReport)
    for col in ("family", "recommended_resolution", "supported_states", "grouped_states",
                "unresolved_mass", "reason", "confidence"):
        assert col in rep.table.columns
    assert set(rep.table["family"]) == {"famA", "famB"}


def test_collinear_family_grouped_or_broad_or_diagnostic():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping)
    famA = next(f for f in rep.families if f.family == "famA")
    assert famA.recommended_resolution in (
        "broad_only", "partially_resolved_group", "diagnostic_only")
    # A1/A2 must NOT be reported as independently resolved fine states
    assert famA.recommended_resolution != "resolved_fine"


def test_separable_family_stays_fine():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping)
    famB = next(f for f in rep.families if f.family == "famB")
    assert famB.recommended_resolution == "resolved_fine"
    assert set(famB.supported_states) == {"B1", "B2"}


def test_broad_only_for_fully_collinear_family():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping)
    famA = next(f for f in rep.families if f.family == "famA")
    # A1≈A2 are mutually collinear (BC>0.97) → broad_only
    assert famA.recommended_resolution == "broad_only"
    assert famA.max_within_bc > 0.97


def test_rare_protected_state_not_collapsed():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping, rare_protection=["A2"])
    famA = next(f for f in rep.families if f.family == "famA")
    # protected A2 must not be swept into a collapsed/grouped bucket
    grouped_flat = {s for g in famA.grouped_states for s in g}
    assert "A2" not in grouped_flat


def test_unresolved_family_when_mass_high():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping,
                                           unresolved_mass={"famB": 0.8})
    famB = next(f for f in rep.families if f.family == "famB")
    assert famB.recommended_resolution == "unresolved_family"


def test_reasons_recorded():
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping)
    assert all(isinstance(f.reason, str) and f.reason for f in rep.families)
    assert rep.metadata["estimates_modified"] is False


def test_write_outputs(tmp_path):
    ref, mapping = _ref()
    rep = build_adaptive_resolution_report(ref, mapping)
    paths = rep.write(tmp_path)
    assert paths["adaptive_resolution"].exists()
