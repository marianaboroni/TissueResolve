"""Tests for experimental family-specific unresolved calibration.

Verifies mass conservation (overall + per broad family), non-negative unresolved
mass, broad-only stays broad-only, supported subtypes retain more mass than
unsupported, rare marker-supported subtypes are not zeroed, and calibration can
be disabled.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.unresolved_calibration import (
    calibrate_unresolved_by_family, CalibrationParams,
)


def _toy():
    # family F: A1 (supported), A2 (unsupported); family G: B1 (single)
    df = pd.DataFrame({
        "A1": [0.4, 0.5], "A2": [0.2, 0.1], "B1": [0.3, 0.3],
        "unresolved_F": [0.0, 0.0], "unresolved_G": [0.1, 0.1],
    }, index=["s0", "s1"])
    mapping = {"A1": "F", "A2": "F", "B1": "G"}
    return df, mapping


def test_mode_none_is_identity():
    df, mapping = _toy()
    res = calibrate_unresolved_by_family(df, mode="none", mapping=mapping)
    pd.testing.assert_frame_equal(res.predictions, df)
    assert res.metadata["unresolved_calibration"] == "none"


def test_mass_conservation_total_and_per_family():
    df, mapping = _toy()
    res = calibrate_unresolved_by_family(
        df, confidence_scores={"A1": 0.9, "A2": 0.1, "B1": 0.9},
        resolution_decisions={"F": "selected_fine", "G": "full_fine"},
        mapping=mapping)
    out = res.predictions
    # overall mass conserved
    assert np.allclose(out.sum(axis=1).to_numpy(), df.sum(axis=1).to_numpy())
    # per-family mass conserved (members + unresolved_<fam>)
    f_before = df[["A1", "A2"]].sum(1) + df["unresolved_F"]
    f_after = out[["A1", "A2"]].sum(1) + out["unresolved_F"]
    assert np.allclose(f_before.to_numpy(), f_after.to_numpy())


def test_unresolved_non_negative():
    df, mapping = _toy()
    res = calibrate_unresolved_by_family(
        df, confidence_scores={"A1": 0.9, "A2": 0.1, "B1": 0.9},
        resolution_decisions={"F": "selected_fine", "G": "full_fine"}, mapping=mapping)
    assert (res.predictions["unresolved_F"] >= -1e-12).all()
    assert (res.predictions["unresolved_G"] >= -1e-12).all()


def test_supported_subtype_retains_more_than_unsupported():
    df, mapping = _toy()
    res = calibrate_unresolved_by_family(
        df, confidence_scores={"A1": 0.9, "A2": 0.1, "B1": 0.9},
        resolution_decisions={"F": "selected_fine", "G": "full_fine"}, mapping=mapping)
    out = res.predictions
    assert out["A1"].sum() > out["A2"].sum()
    assert np.allclose(out["A2"].to_numpy(), 0.0)   # unsupported dropped to unresolved


def test_broad_only_family_stays_broad_only():
    df, mapping = _toy()
    res = calibrate_unresolved_by_family(
        df, resolution_decisions={"F": "broad_only", "G": "full_fine"}, mapping=mapping)
    out = res.predictions
    assert np.allclose(out["A1"].to_numpy(), 0.0) and np.allclose(out["A2"].to_numpy(), 0.0)
    assert (out["unresolved_F"] > 0).all()


def test_rare_marker_supported_subtype_not_zeroed():
    df = pd.DataFrame({"A1": [0.7, 0.7], "A2": [0.005, 0.004],
                       "unresolved_F": [0.295, 0.296]}, index=["s0", "s1"])
    mapping = {"A1": "F", "A2": "F"}
    res = calibrate_unresolved_by_family(
        df, confidence_scores={"A1": 0.9, "A2": 0.1},   # A2 low confidence but rare+marker
        resolution_decisions={"F": "selected_fine"},
        family_metrics={"F": {"marker_supported": {"A2": True}}},
        calibration_params=CalibrationParams(rare_abundance=0.01), mapping=mapping)
    assert res.predictions["A2"].sum() > 0   # not zeroed
