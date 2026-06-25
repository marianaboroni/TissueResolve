"""Offline tests for the split bulk/spatial benchmark plotting modules."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.plotting import bulk_benchmark_plots as BB
from tissueresolve.plotting import spatial_benchmark_plots as SB

_STATUS = pd.DataFrame({
    "method": ["TissueResolve_auto", "NNLS_baseline", "MuSiC",
               "CARD", "cell2location", "SPOTlight"],
    "modality": ["bulk", "bulk", "bulk", "spatial", "spatial", "spatial"],
    "status": ["executed", "executed", "executed", "executed", "executed", "skipped"],
    "runtime_seconds": [13.9, 1.7, 40.0, 120.0, 90.0, 0.0],
    "accuracy": [0.77, 0.77, 0.65, np.nan, np.nan, np.nan],
}).set_index("method")


def _ok(res):
    assert res.data_paths and all(p.exists() for p in res.data_paths.values())


# --- bulk benchmark ---------------------------------------------------------

def test_bulk_status_summary_bulk_only(tmp_path):
    res = BB.plot_bulk_method_status_summary(_STATUS, tmp_path)
    _ok(res)


def test_bulk_fine_leaderboard_excludes_spatial(tmp_path):
    res = BB.plot_bulk_fine_accuracy_leaderboard(_STATUS, tmp_path)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    methods = set(df["method"])
    assert "CARD" not in methods and "cell2location" not in methods
    assert "TissueResolve_auto" in methods
    assert "ground truth" in (res.caption or "").lower()


def test_bulk_runtime_comparison(tmp_path):
    res = BB.plot_bulk_runtime_comparison(_STATUS, tmp_path)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert set(df["method"]) <= {"TissueResolve_auto", "NNLS_baseline", "MuSiC"}


def test_bulk_rmse_mae(tmp_path):
    metrics = pd.DataFrame({"method": ["A", "B"], "rmse": [0.04, 0.06],
                            "mae": [0.02, 0.03]})
    res = BB.plot_bulk_rmse_mae_comparison(metrics, tmp_path)
    _ok(res)


# --- spatial benchmark ------------------------------------------------------

def test_spatial_status_summary_spatial_only(tmp_path):
    res = SB.plot_spatial_method_status_summary(_STATUS, tmp_path)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    # spatial-only: total counted methods == spatial methods (3)
    assert int(df["n_methods"].sum()) == 3
    assert "no spot-level ground truth" in (res.caption or "").lower()


def test_spatial_concordance_heatmap(tmp_path):
    conc = pd.DataFrame([[1.0, 0.6], [0.6, 1.0]],
                        index=["CARD", "cell2location"],
                        columns=["CARD", "cell2location"])
    res = SB.plot_spatial_concordance_heatmap(conc, tmp_path)
    assert "not accuracy" in (res.caption or "").lower()


def test_spatial_concordance_needs_two_methods(tmp_path):
    conc = pd.DataFrame([[1.0]], index=["CARD"], columns=["CARD"])
    with pytest.raises(ValueError):
        SB.plot_spatial_concordance_heatmap(conc, tmp_path)


def test_spatial_structure_metrics(tmp_path):
    metrics = pd.DataFrame({
        "population": ["T/NK", "Myeloid", "T/NK", "Myeloid"],
        "metric": ["morans_i", "morans_i", "entropy", "entropy"],
        "value": [0.4, 0.2, 1.1, 1.3]})
    res = SB.plot_spatial_structure_metrics_summary(metrics, tmp_path)
    _ok(res)
    assert "not accuracy" in (res.caption or "").lower()


def test_spatial_runtime(tmp_path):
    res = SB.plot_spatial_runtime_comparison(_STATUS, tmp_path)
    df = pd.read_csv(res.data_paths["data"], sep="\t", comment="#")
    assert set(df["method"]) <= {"CARD", "cell2location"}


def test_spatial_no_accuracy_functions(tmp_path):
    """Spatial benchmark module exposes no accuracy/RMSE/Pearson-vs-truth plot."""
    for fn in SB.__all__:
        low = fn.lower()
        assert "rmse" not in low and "accuracy" not in low and "pearson" not in low
    # the spatial caption carries the no-ground-truth caveat
    res = SB.plot_spatial_method_status_summary(_STATUS, tmp_path)
    assert "not accuracy" in (res.caption or "").lower()
