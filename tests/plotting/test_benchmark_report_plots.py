"""Offline tests for benchmark report figures (report Part 9)."""
from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting import benchmark_report_plots as B

_STATUS = pd.DataFrame({
    "method": ["TissueResolve_auto", "NNLS_baseline", "MuSiC", "CARD",
               "cell2location", "BayesPrism", "SPOTlight"],
    "modality": ["bulk", "bulk", "bulk", "spatial", "spatial", "bulk", "spatial"],
    "status": ["executed", "executed", "executed", "executed", "executed",
               "failed", "skipped"],
    "runtime_seconds": [13.9, 1.7, 40.0, 120.0, 90.0, 0.0, 0.0],
    "accuracy": [0.77, 0.77, 0.65, float("nan"), float("nan"),
                 float("nan"), float("nan")],
}).set_index("method")

_COMPOSITE = pd.DataFrame({
    "method": ["TissueResolve_auto", "NNLS_baseline"],
    "modality": ["bulk", "bulk"],
    "accuracy_score": [0.77, 0.77],
    "robustness_score": [0.85, 0.5],
    "usability_score": [0.9, 0.7],
    "final_score": [0.80, 0.54],
}).set_index("method")


def _ok(res):
    assert res.data_paths and all(p.exists() for p in res.data_paths.values())


def test_method_status_summary(tmp_path):
    res = B.plot_benchmark_method_status(_STATUS, tmp_path)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert "n_methods" in df.columns


def test_bulk_accuracy_leaderboard_excludes_spatial(tmp_path):
    res = B.plot_bulk_accuracy_leaderboard(_STATUS, tmp_path)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    methods = set(df["method"])
    # only bulk executed/imported with accuracy
    assert "CARD" not in methods and "cell2location" not in methods
    assert "BayesPrism" not in methods  # failed
    assert "TissueResolve_auto" in methods


def test_runtime_comparison(tmp_path):
    res = B.plot_runtime_comparison(_STATUS, tmp_path)
    _ok(res)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert (df["runtime_seconds"] > 0).all()


def test_composite_scorecard_labeled(tmp_path):
    res = B.plot_composite_scorecard(_COMPOSITE, tmp_path)
    _ok(res)
    assert "scorecard" in (res.caption or "").lower()
    assert "not an objective accuracy" in (res.caption or "").lower() \
        or "not objective accuracy" in (res.caption or "").lower()


def test_missing_columns_raise(tmp_path):
    with pytest.raises(ValueError):
        B.plot_benchmark_method_status(pd.DataFrame({"x": [1]}), tmp_path)
    with pytest.raises(ValueError):
        B.plot_bulk_accuracy_leaderboard(pd.DataFrame({"method": ["a"]}), tmp_path)
