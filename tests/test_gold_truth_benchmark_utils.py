from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from benchmarks.diagnostics import gold_truth_performance_benchmark as G


def test_split_donors_is_disjoint_and_complete():
    donors = [f"D{i}" for i in range(12)]
    split = G.split_donors(donors, seed=5)
    G.assert_donor_disjoint(split)
    recovered = set().union(*(set(v) for v in split.values()))
    assert recovered == set(donors)
    assert split["reference_train"]
    assert split["final_test"]


def test_validate_truth_tables_accepts_consistent_truth():
    mapping = {"a": "F1", "b": "F1", "c": "F2"}
    fine = pd.DataFrame(
        {"a": [0.2, 0.0], "b": [0.3, 0.5], "c": [0.5, 0.5]},
        index=["s1", "s2"],
    )
    broad = G.aggregate_truth_to_broad(fine, mapping)
    cond = G.conditional_truth_from_fine(fine, mapping)
    G.validate_truth_tables(fine, broad, cond, mapping)
    assert np.allclose(broad.sum(axis=1), 1.0)
    assert np.allclose(cond.loc[["s1"], ["a", "b"]].sum(axis=1), 1.0)


def test_validate_truth_tables_rejects_leaky_broad_truth():
    mapping = {"a": "F1", "b": "F1", "c": "F2"}
    fine = pd.DataFrame({"a": [0.2], "b": [0.3], "c": [0.5]}, index=["s1"])
    broad = pd.DataFrame({"F1": [0.1], "F2": [0.9]}, index=["s1"])
    cond = G.conditional_truth_from_fine(fine, mapping)
    with pytest.raises(ValueError, match="broad truth"):
        G.validate_truth_tables(fine, broad, cond, mapping)


def test_conditional_truth_handles_zero_family_safely():
    mapping = {"a": "F1", "b": "F1", "c": "F2"}
    fine = pd.DataFrame({"a": [0.0], "b": [0.0], "c": [1.0]}, index=["s1"])
    cond = G.conditional_truth_from_fine(fine, mapping)
    assert cond.loc["s1", ["a", "b"]].sum() == 0.0
    assert cond.loc["s1", "c"] == 1.0


def test_composition_metrics_handles_zero_truth():
    true = pd.DataFrame({"a": [0.0, 1.0], "b": [1.0, 0.0]}, index=["s1", "s2"])
    pred = pd.DataFrame({"a": [0.1, 0.8], "b": [0.9, 0.2]}, index=["s1", "s2"])
    m = G.composition_metrics(true, pred)
    assert m["rmse"] >= 0.0
    assert np.isfinite(m["jensen_shannon"])
    assert np.isfinite(m["aitchison"])


def test_report_level_donor_leak_detection():
    split = {
        "reference_train": ["D1", "D2"],
        "calibration_unused": ["D3"],
        "validation_unused": ["D2"],
        "final_test": ["D4"],
    }
    with pytest.raises(ValueError, match="donor split leak"):
        G.assert_donor_disjoint(split)
