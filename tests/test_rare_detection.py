"""Tests for the rare-subtype detection / calibration layer (P2). Decision layer
on top of any solver — never modifies the solver; the gate conserves mass.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.results import ReferenceSignature
from tissueresolve.experimental.rare_detection import (
    compute_marker_support, rare_detection_probability, calibrate_detection,
    precision_recall_curve, apply_rare_detection_gate)


def _ref():
    rng = np.random.default_rng(0)
    G, K = 60, 3
    genes = [f"g{i}" for i in range(G)]
    cts = ["common", "other", "RARE"]
    R = rng.random((K, G)) * 5 + 0.5
    # give RARE a private marker block (genes 0..9 high, others low)
    R[2, :10] += 50
    R[:2, :10] *= 0.05
    R_cpm = (R / R.sum(1, keepdims=True) * 1e6).astype(np.float32)
    ref = ReferenceSignature(gene_names=genes, cell_types=cts,
                             phi=(R_cpm.T / R_cpm.T.sum(0, keepdims=True)).astype(float),
                             R_cpm=R_cpm, R_log=np.log1p(R_cpm).astype(np.float32),
                             phi_g=np.full(G, 8.0, dtype=np.float32),
                             n_cells_per_type={c: 100 for c in cts}, genome="hg38")
    return ref, genes


def _bulk_with_without_rare(ref, genes, n=6):
    """Samples 0..n-1 contain RARE (marker block expressed), n..2n-1 do not."""
    rng = np.random.default_rng(1)
    R = ref.as_R_cpm() / 1e6
    cols, truth_present = [], []
    for present in (True, False):
        for _ in range(n):
            theta = np.array([0.5, 0.5, 0.0]) if not present else np.array([0.45, 0.40, 0.15])
            mu = (theta @ R) * 5000
            cols.append(rng.poisson(mu))
            truth_present.append(present)
    B = pd.DataFrame(np.array(cols).T, index=genes,
                     columns=[f"s{i}" for i in range(2 * n)])
    return B, np.array(truth_present)


def test_marker_support_higher_when_rare_present():
    ref, genes = _ref()
    bulk, present = _bulk_with_without_rare(ref, genes)
    sup = compute_marker_support(bulk, ref, ["RARE"], n_markers=10)
    s = sup["RARE"].to_numpy()
    assert np.nanmean(s[present]) > np.nanmean(s[~present])


def test_detection_probability_monotonic_in_mass_and_support():
    ref, genes = _ref()
    props = pd.DataFrame({"common": [0.5, 0.5], "other": [0.4, 0.49], "RARE": [0.1, 0.01]},
                         index=["a", "b"])
    sup = pd.DataFrame({"RARE": [1.0, 1.0]}, index=["a", "b"])
    score = rare_detection_probability(props, sup, ["RARE"])
    assert score.loc["a", "RARE"] > score.loc["b", "RARE"]   # more mass → higher
    sup2 = pd.DataFrame({"RARE": [1.0, 0.0]}, index=["a", "b"])
    score2 = rare_detection_probability(props, sup2, ["RARE"])
    assert score2.loc["b", "RARE"] == 0.0                    # no marker support → 0


def test_precision_recall_curve_computed():
    scores = np.array([0.9, 0.8, 0.2, 0.1, 0.6, 0.05])
    labels = np.array([1, 1, 0, 0, 1, 0])
    pr = precision_recall_curve(scores, labels)
    assert {"threshold", "precision", "recall", "fpr"} <= set(pr.columns)
    assert pr["recall"].iloc[0] == 1.0          # threshold 0 → all called → recall 1
    assert (pr["fpr"].dropna() <= 1.0).all()


def test_calibration_returns_valid_probabilities():
    rng = np.random.default_rng(0)
    s = rng.random(200)
    y = (s + rng.normal(0, 0.1, 200) > 0.5).astype(int)
    cal = calibrate_detection(s, y, kind="isotonic")
    p = cal(np.array([0.0, 0.5, 1.0]))
    assert np.all((p >= 0) & (p <= 1))
    assert cal(np.array([0.9]))[0] >= cal(np.array([0.1]))[0]  # monotone increasing


def test_gate_reduces_false_positives_and_conserves_mass():
    ref, genes = _ref()
    bulk, present = _bulk_with_without_rare(ref, genes)
    # solver over-calls RARE everywhere (false positives in the absent samples)
    props = pd.DataFrame(0.0, index=bulk.columns, columns=["common", "other", "RARE"])
    props["RARE"] = 0.12
    props["common"] = 0.48
    props["other"] = 0.40
    sup = compute_marker_support(bulk, ref, ["RARE"], n_markers=10)
    prob = rare_detection_probability(props, sup, ["RARE"])
    gated, meta = apply_rare_detection_gate(props, prob, sup, ["RARE"],
                                            family_map={"RARE": "famR"},
                                            threshold=0.3, min_marker_support=0.5)
    # mass conserved per row
    assert np.allclose(gated.sum(axis=1).to_numpy(), props.sum(axis=1).to_numpy(), atol=1e-9)
    # RARE dropped (→ unresolved) more often in the truly-absent samples
    rare_after = gated["RARE"].to_numpy()
    assert (rare_after[~present] == 0).mean() > (rare_after[present] == 0).mean()
    assert meta["mass_conserving"] and meta["n_calls_dropped"] >= 1


def test_gate_protects_listed_state():
    ref, genes = _ref()
    bulk, _ = _bulk_with_without_rare(ref, genes)
    props = pd.DataFrame(0.0, index=bulk.columns, columns=["common", "other", "RARE"])
    props["RARE"] = 0.1; props["common"] = 0.5; props["other"] = 0.4
    sup = pd.DataFrame(0.0, index=bulk.columns, columns=["RARE"])  # zero support everywhere
    prob = rare_detection_probability(props, sup, ["RARE"])
    gated, meta = apply_rare_detection_gate(props, prob, sup, ["RARE"],
                                            threshold=0.9, min_marker_support=0.9,
                                            protected=["RARE"])
    assert np.allclose(gated["RARE"].to_numpy(), props["RARE"].to_numpy())
    assert meta["n_calls_dropped"] == 0


def test_gate_does_not_touch_nonrare_columns():
    ref, genes = _ref()
    bulk, _ = _bulk_with_without_rare(ref, genes)
    props = pd.DataFrame(0.0, index=bulk.columns, columns=["common", "other", "RARE"])
    props["RARE"] = 0.5; props["common"] = 0.3; props["other"] = 0.2  # force all dropped
    sup = pd.DataFrame(0.0, index=bulk.columns, columns=["RARE"])
    prob = rare_detection_probability(props, sup, ["RARE"])
    common_before = props["common"].copy()
    gated, _ = apply_rare_detection_gate(props, prob, sup, ["RARE"], threshold=0.9)
    assert np.allclose(gated["common"].to_numpy(), common_before.to_numpy())
