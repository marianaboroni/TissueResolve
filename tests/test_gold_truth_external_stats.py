"""Phase 1 tests for gold-truth external benchmark statistics layer.

Fast, offline, deterministic. Uses tiny synthetic fixtures only.
Covers: per-sample prediction saving, internal/external sample-ID alignment,
paired-bootstrap refusal when predictions missing, skipped tools not ranked,
and estimate-type recording.
"""
from __future__ import annotations

import json
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import benchmarks.diagnostics.gold_truth_performance_benchmark as gtp
import benchmarks.diagnostics.gold_truth_paired_stats as ps


def _toy_ctx():
    mapping = {"A1": "A", "A2": "A", "B1": "B"}
    spec = types.SimpleNamespace(name="toy")
    split = {"reference_train": ["d1", "d2"], "final_test": ["d3"]}
    return types.SimpleNamespace(spec=spec, mapping=mapping, split=split)


def _toy_bundle():
    idx = ["s0", "s1", "s2"]
    fine = pd.DataFrame({"A1": [0.5, 0.2, 0.1], "A2": [0.3, 0.3, 0.4],
                         "B1": [0.2, 0.5, 0.5]}, index=idx)
    return types.SimpleNamespace(fine_truth=fine)


def test_save_mode_predictions_structure_and_estimate_type(tmp_path):
    ctx, bundle = _toy_ctx(), _toy_bundle()
    pred = bundle.fine_truth.copy()  # well-formed RNA-derived proportions
    meta = {"runtime_seconds": 0.1, "trusted_resolution": {"A": "selected_fine"}}
    rows = gtp._save_mode_predictions(tmp_path, ctx, bundle, pred, meta,
                                      mode="flat", group_name="full_overlap", seed=2026)
    fine_file = tmp_path / "predictions" / "TissueResolve_flat__toy.tsv"
    broad_file = tmp_path / "predictions_broad" / "TissueResolve_flat__toy.tsv"
    assert fine_file.exists() and broad_file.exists()
    saved = pd.read_csv(fine_file, sep="\t", index_col=0)
    # rows = sample IDs, columns = fine labels
    assert list(saved.index) == ["s0", "s1", "s2"]
    assert set(saved.columns) == {"A1", "A2", "B1"}
    # estimate type recorded in sidecar
    sidecar = json.loads((fine_file.with_suffix(".meta.json")).read_text())
    assert sidecar["estimate_type"] == "RNA_derived"
    assert sidecar["donor_split"]["final_test"] == ["d3"]
    assert any(r["estimate_type"] == "RNA_derived" for r in rows)


def test_broad_aggregation_matches_mapping():
    fine = pd.DataFrame({"A1": [0.5], "A2": [0.3], "B1": [0.2]}, index=["s0"])
    truth_broad = pd.DataFrame({"A": [0.8], "B": [0.2]}, index=["s0"])
    mapping = {"A1": "A", "A2": "A", "B1": "B"}
    broad = ps._broad_from_fine(fine, mapping, truth_broad)
    assert broad.loc["s0", "A"] == pytest.approx(0.8)
    assert broad.loc["s0", "B"] == pytest.approx(0.2)


def test_fine_aligned_drops_unresolved_and_reindexes():
    truth = pd.DataFrame({"A1": [0.5], "A2": [0.5]}, index=["s0"])
    pred = pd.DataFrame({"A1": [0.4], "unresolved_A": [0.6]}, index=["s0"])
    aligned = ps._fine_aligned(pred, truth)
    assert list(aligned.columns) == ["A1", "A2"]  # unresolved dropped, A2 filled
    assert aligned.loc["s0", "A1"] == pytest.approx(0.4)
    assert aligned.loc["s0", "A2"] == pytest.approx(0.0)


def test_per_sample_rmse_correct():
    pred = pd.DataFrame({"A": [0.0, 1.0], "B": [1.0, 0.0]}, index=["s0", "s1"])
    truth = pd.DataFrame({"A": [0.0, 0.0], "B": [1.0, 1.0]}, index=["s0", "s1"])
    rmse = ps._per_sample_rmse(pred, truth)
    assert rmse[0] == pytest.approx(0.0)
    assert rmse[1] == pytest.approx(1.0)  # sqrt(mean(1,1))


def test_internal_external_sample_ids_align():
    """If predictions exist, internal and external full_overlap IDs must match."""
    internal = ps._read_internal("TissueResolve_flat", "breast")
    external = ps._read_external("MuSiC", "breast")
    if internal is None or external is None:
        pytest.skip("benchmark predictions not generated in this environment")
    assert list(internal.index) == list(external.index)


def test_paired_bootstrap_refuses_when_predictions_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(ps, "PRED_DIR", tmp_path / "no_such_predictions")
    rc = ps.run(n_boot=10, seed=0)
    assert rc == 2  # refuses, does not fabricate CIs


def test_skipped_tool_not_ranked(tmp_path, monkeypatch):
    """A tool with no prediction file must not appear in rank_by_task."""
    # point all dirs at empty tmp except provide one internal + truth pair
    monkeypatch.setattr(ps, "PRED_DIR", tmp_path / "predictions")
    monkeypatch.setattr(ps, "PERF_DIR", tmp_path / "perf")
    monkeypatch.setattr(ps, "EXTERNAL_DIR", tmp_path / "ext")
    monkeypatch.setattr(ps, "INPUT_DIR", tmp_path / "inp")
    monkeypatch.setattr(ps, "DATASETS", ["toy"])

    (tmp_path / "predictions").mkdir(parents=True)
    (tmp_path / "perf").mkdir(parents=True)
    (tmp_path / "ext").mkdir(parents=True)

    idx = ["s0", "s1", "s2"]
    fine = pd.DataFrame({"A1": [0.5, 0.2, 0.1], "A2": [0.3, 0.3, 0.4],
                         "B1": [0.2, 0.5, 0.5], "sample": idx,
                         "dataset": "toy", "group": "full_overlap"})
    broad = pd.DataFrame({"A": [0.8, 0.5, 0.5], "B": [0.2, 0.5, 0.5], "sample": idx,
                          "dataset": "toy", "group": "full_overlap"})
    fine.to_csv(tmp_path / "perf" / "fine_truth.tsv", sep="\t", index=False)
    broad.to_csv(tmp_path / "perf" / "broad_truth.tsv", sep="\t", index=False)
    fine.to_csv(tmp_path / "perf" / "conditional_truth.tsv", sep="\t", index=False)
    # only TissueResolve_flat has predictions; MuSiC/BisqueRNA/etc do not
    pd.DataFrame({"A1": [0.5, 0.2, 0.1], "A2": [0.3, 0.3, 0.4], "B1": [0.2, 0.5, 0.5]},
                 index=idx).to_csv(tmp_path / "predictions" / "TissueResolve_flat__toy.tsv", sep="\t")

    rc = ps.run(n_boot=50, seed=0)
    assert rc == 0
    rank = pd.read_csv(tmp_path / "ext" / "rank_by_task.tsv", sep="\t")
    present = set(rank["method"].unique())
    assert "TissueResolve_flat" in present
    assert "MuSiC" not in present  # skipped (no predictions) -> not ranked
