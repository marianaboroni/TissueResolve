"""Tests for the experimental spatial weak-smoothing preset.

Fast, offline, deterministic. Verifies the preset layer does not change default
spatial behaviour, sets the documented lambdas, marks experimental presets, and
records provenance metadata — and that bulk config is untouched.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.config import TissueResolveConfig
from tissueresolve.experimental.spatial_presets import (
    DEFAULT_LAMBDA_SPATIAL, apply_spatial_preset, available_presets,
)


def test_default_lambda_remains_0_1():
    assert TissueResolveConfig().spatial_solver.lambda_spatial == 0.1
    assert DEFAULT_LAMBDA_SPATIAL == 0.1


def test_default_preset_is_noop():
    cfg = TissueResolveConfig()
    info = apply_spatial_preset(cfg, "default")
    assert cfg.spatial_solver.lambda_spatial == 0.1
    assert info.lambda_spatial == 0.1
    assert info.experimental is False
    assert info.is_default is True


def test_weak_smoothing_sets_lower_lambda_and_is_experimental():
    cfg = TissueResolveConfig()
    info = apply_spatial_preset(cfg, "weak_smoothing")
    assert cfg.spatial_solver.lambda_spatial == 0.02
    assert info.lambda_spatial == 0.02
    assert info.experimental is True
    assert info.is_default is False
    assert info.smoothing_used is True


def test_no_smoothing_disables_smoothing():
    cfg = TissueResolveConfig()
    info = apply_spatial_preset(cfg, "no_smoothing")
    assert cfg.spatial_solver.lambda_spatial == 0.0
    assert info.smoothing_used is False
    assert info.experimental is True


def test_metadata_records_lambda_and_preset_and_smoothing_status():
    cfg = TissueResolveConfig()
    md = apply_spatial_preset(cfg, "weak_smoothing").to_metadata()
    assert md["spatial_preset"] == "weak_smoothing"
    assert md["lambda_spatial"] == 0.02
    assert md["smoothing_used"] is True
    assert md["spatial_preset_experimental"] is True
    assert md["is_default_spatial_behaviour"] is False
    assert "spatial_preset_description" in md


def test_bulk_behaviour_unchanged_by_spatial_preset():
    cfg = TissueResolveConfig()
    before = (cfg.bulk_solver.__dict__.copy(), cfg.genes.__dict__.copy())
    apply_spatial_preset(cfg, "weak_smoothing")
    after = (cfg.bulk_solver.__dict__.copy(), cfg.genes.__dict__.copy())
    assert before == after  # only spatial_solver.lambda_spatial changed


def test_unknown_preset_raises():
    with pytest.raises(ValueError):
        apply_spatial_preset(TissueResolveConfig(), "ultra_smoothing")


def test_available_presets():
    presets = available_presets()
    assert {"default", "weak_smoothing", "no_smoothing"} <= set(presets)


def test_report_documents_lambda_and_smoothing_status():
    report = (REPO / "docs" / "SPATIAL_WEAK_SMOOTHING_BENCHMARK_REPORT.md").read_text()
    assert "0.02" in report and "0.1" in report
    assert "oversmoothing" in report.lower()
    assert "weak_smoothing" in report


# --- generalized synthetic spatial generator (domain_families) ---------------
def test_build_spatial_targets_respects_domain_families():
    import numpy as np
    from benchmarks.shared.synthetic_spatial import build_spatial_targets, BREAST_DOMAIN_FAMILIES
    ref_types = ["epi1", "epi2", "imm1", "imm2", "str1"]
    mapping = {"epi1": "Epi", "epi2": "Epi", "imm1": "Imm", "imm2": "Imm", "str1": "Str"}
    fams = {"left": ["Epi"], "right": ["Imm"], "gradient": ["Str"]}
    targets, coords, domains, meta = build_spatial_targets(
        ref_types, mapping, n_side=6, seed=0, domain_families=fams)
    # left half should be Epi-enriched, right half Imm-enriched
    mid = 3
    left = targets[coords[:, 1] < mid]
    right = targets[coords[:, 1] >= mid]
    assert left[["epi1", "epi2"]].sum(1).mean() > left[["imm1", "imm2"]].sum(1).mean()
    assert right[["imm1", "imm2"]].sum(1).mean() > right[["epi1", "epi2"]].sum(1).mean()
    assert set(domains) <= {"L", "R", "niche"}
    # backward-compat default exists and is the breast layout
    assert BREAST_DOMAIN_FAMILIES["left"] == ["Epithelial"]


def test_build_spatial_targets_default_is_breast_layout():
    from benchmarks.shared.synthetic_spatial import build_spatial_targets
    ref_types = ["T", "B"]
    mapping = {"T": "Epithelial", "B": "Myeloid"}
    # no domain_families -> uses breast defaults without error
    targets, coords, domains, meta = build_spatial_targets(ref_types, mapping, n_side=4, seed=1)
    assert targets.shape[0] == 16 and abs(targets.sum(1).mean() - 1.0) < 1e-6
