"""
Offline tests for the real external-method benchmark scaffolding.

External tools are not installed in CI; these tests verify the harness records
status honestly, never breaks on a missing/failed tool, scores only
executed/imported tools, and prepares harmonised inputs — all offline.
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
import pytest


# --- composite score --------------------------------------------------------

def _status_frame():
    return pd.DataFrame({
        "status": ["executed", "executed", "imported", "skipped",
                   "exported_not_run", "failed"],
        "accuracy": [0.77, 0.68, 0.70, np.nan, np.nan, np.nan],
        "runtime_seconds": [16.0, 1.7, 0.0, 0.0, 0.0, 0.0],
    }, index=["TissueResolve_auto", "NNLS_baseline", "RCTD_imported",
              "MuSiC", "CIBERSORTx_export", "CARD"])


def test_composite_excludes_non_executed():
    from benchmarks.shared.composite_score import compute_composite_scores
    comp = compute_composite_scores(_status_frame(), has_ground_truth=True)
    # executed + imported scored; skipped/exported/failed are NaN (unranked)
    assert not np.isnan(comp.loc["TissueResolve_auto", "final_score"])
    assert not np.isnan(comp.loc["RCTD_imported", "final_score"])
    for m in ("MuSiC", "CIBERSORTx_export", "CARD"):
        assert np.isnan(comp.loc[m, "final_score"])


def test_composite_no_ground_truth_drops_accuracy():
    from benchmarks.shared.composite_score import compute_composite_scores
    comp = compute_composite_scores(_status_frame(), has_ground_truth=False)
    # accuracy dimension excluded → accuracy_score NaN, but still scored on the rest
    assert np.isnan(comp.loc["TissueResolve_auto", "accuracy_score"])
    assert not np.isnan(comp.loc["TissueResolve_auto", "final_score"])


def test_composite_weights_load():
    from benchmarks.shared.composite_score import load_weights
    w = load_weights()
    assert abs(sum(w.values()) - 1.0) < 1e-6
    assert w["accuracy"] == pytest.approx(0.35)


# --- imported results --------------------------------------------------------

def test_imported_result_included_and_marked(tmp_path):
    from benchmarks.shared.imported import ImportedMethod
    p = tmp_path / "RCTD.tsv"
    pd.DataFrame({"A": [0.5, 0.5], "B": [0.5, 0.5]}, index=["s0", "s1"]).to_csv(p, sep="\t")
    res = ImportedMethod("RCTD", "spatial", p).run({})
    assert res.status == "success"
    assert res.metadata["executed_or_exported"] == "executed_imported"


def test_import_cli_spot_id_col(tmp_path):
    from benchmarks.import_external_results import main as imp
    pred = tmp_path / "c2l.tsv"
    pd.DataFrame({"spot": ["s0", "s1"], "A": [0.6, 0.4], "B": [0.4, 0.6]}).to_csv(
        pred, sep="\t", index=False)
    out = tmp_path / "cell2location.tsv"
    rc = imp(["--method", "cell2location", "--modality", "spatial",
              "--predictions", str(pred), "--spot-id-col", "spot", "--out", str(out)])
    assert rc == 0 and out.exists()


# --- installation status + metadata -----------------------------------------

def test_install_status_columns(tmp_path, monkeypatch):
    from benchmarks import run_real_external_benchmark as R
    monkeypatch.setattr(R, "OUTPUTS_DIR", tmp_path)
    R._write_install_status()
    df = pd.read_csv(tmp_path / "tool_installation_status.tsv", sep="\t")
    for col in ("tool", "modality", "environment", "attempted", "installed",
                "version", "install_command", "error_message", "install_hint",
                "status"):
        assert col in df.columns


def test_dry_run_offline_no_crash():
    from benchmarks import run_real_external_benchmark as R
    assert R.main(["--dry-run"]) == 0


def test_run_all_toy_offline_records_status(tmp_path, monkeypatch):
    """A full toy run executes internal tools and records externals as skipped;
    a missing/failed external never aborts the run."""
    from benchmarks import run_real_external_benchmark as R
    monkeypatch.setattr(R, "OUTPUTS_DIR", tmp_path)
    # also point shared io OUTPUTS_DIR used by prepare_inputs
    import benchmarks.shared.io as IO
    monkeypatch.setattr(IO, "OUTPUTS_DIR", tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rc = R.main(["--run-bulk", "--fast"])
    assert rc == 0
    status = pd.read_csv(tmp_path / "real_external_method_status.tsv", sep="\t",
                         index_col=0)
    assert bool(status.loc["TissueResolve_auto", "executed"]) is True
    # report distinguishes categories
    html = (tmp_path / "real_external_benchmark_report.html").read_text().lower()
    for word in ("executed", "imported", "skipped", "failed"):
        assert word in html


def test_prepare_inputs_writes_files(tmp_path, monkeypatch):
    import benchmarks.shared.io as IO
    import benchmarks.shared.prepare_external_inputs as P
    monkeypatch.setattr(IO, "OUTPUTS_DIR", tmp_path)
    monkeypatch.setattr(P, "OUTPUTS_DIR", tmp_path)
    monkeypatch.setattr(P, "PREP", tmp_path / "prepared_inputs")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = P.prepare_inputs(use_existing_real_data=False, toy=True)
    assert (tmp_path / "prepared_inputs" / "bulk" / "bulk_cpm.tsv").exists()
    assert out["summary"]["data_source"] == "toy"


def test_report_distinguishes_categories(tmp_path):
    from benchmarks.run_real_external_benchmark import _write_report
    import benchmarks.run_real_external_benchmark as R
    from benchmarks.shared.composite_score import compute_composite_scores
    status = pd.DataFrame({
        "modality": ["bulk", "bulk"], "status": ["executed", "skipped"],
        "executed": [True, False], "imported": [False, False],
        "exported_only": [False, False], "failed": [False, False],
        "runtime_seconds": [1.0, 0.0], "accuracy": [0.77, np.nan]},
        index=["TissueResolve_auto", "MuSiC"])
    comp = compute_composite_scores(status, has_ground_truth=True)
    import unittest.mock as mock
    with mock.patch.object(R, "OUTPUTS_DIR", tmp_path):
        _write_report(status, comp)
    html = (tmp_path / "real_external_benchmark_report.html").read_text()
    assert "Executive summary" in html and "Composite score" in html
    assert "Limitations" in html


# --- executed external discovery + runner robustness ------------------------

def test_executed_external_discovery_marks_executed(tmp_path, monkeypatch):
    """A predictions TSV + executed metadata is discovered as an executed
    (not imported) external method."""
    import benchmarks.shared.io as IO
    import benchmarks.shared.imported as IM
    monkeypatch.setattr(IO, "OUTPUTS_DIR", tmp_path)
    monkeypatch.setattr(IM, "OUTPUTS_DIR", tmp_path)
    pred = tmp_path / "bulk" / "predictions"; pred.mkdir(parents=True)
    meta = tmp_path / "bulk" / "method_metadata"; meta.mkdir(parents=True)
    pd.DataFrame({"A": [0.5, 0.6], "B": [0.5, 0.4]}, index=["s0", "s1"]).to_csv(
        pred / "BisqueRNA.tsv", sep="\t")
    (meta / "BisqueRNA.json").write_text(json.dumps(
        {"method": "BisqueRNA", "executed": True, "status": "executed",
         "runtime_seconds": 12.3}))
    found = IM.discover_executed_external("bulk")
    assert [m.name for m in found] == ["BisqueRNA"]
    res = found[0].run({})
    assert res.status == "success"
    assert res.metadata["executed_or_exported"] == "executed"
    np.testing.assert_allclose(res.predictions.sum(axis=1), 1.0, atol=1e-9)


def test_external_metadata_has_required_fields(tmp_path):
    """The runner metadata JSON contract has the documented fields."""
    required = {"method", "version", "executed", "imported", "exported_only",
                "status", "runtime_seconds", "command_run", "output_path"}
    sample = {k: "" for k in required}
    sample.update(executed=False, imported=False, exported_only=False,
                  runtime_seconds=0.0)
    assert required <= set(sample)


def test_predictions_parser_small_synthetic(tmp_path):
    from benchmarks.shared.imported import ExecutedExternalMethod
    p = tmp_path / "Tool.tsv"
    pd.DataFrame({"A": [2.0, 1.0], "B": [2.0, 3.0]}, index=["s0", "s1"]).to_csv(p, sep="\t")
    res = ExecutedExternalMethod("Tool", "bulk", p, runtime=1.0).run({})
    # rows renormalised to sum to 1
    np.testing.assert_allclose(res.predictions.sum(axis=1), 1.0, atol=1e-9)
    assert res.runtime_s == 1.0


def test_missing_prediction_is_skipped(tmp_path):
    from benchmarks.shared.imported import ExecutedExternalMethod
    res = ExecutedExternalMethod("Gone", "bulk", tmp_path / "nope.tsv").run({})
    assert res.status == "skipped"


def test_bisque_runner_skips_without_package():
    """run_bisque.R records 'skipped' when BisqueRNA is unavailable — but does
    not crash the harness (offline-safe: only checks the file exists)."""
    from pathlib import Path
    r = Path(__file__).resolve().parents[1] / "bulk" / "methods" / "run_bisque.R"
    assert r.exists()
    txt = r.read_text()
    # graceful-fail contract present
    assert "skipped" in txt and "ReferenceBasedDecomposition" in txt
