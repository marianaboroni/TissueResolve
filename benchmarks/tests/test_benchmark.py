"""
Offline tests for the benchmark framework.

Deterministic and offline: external tools must be absent-safe (skipped, never
failing the suite), internal baselines must run on toy data, and the
normalization/protocol/library/batch diagnostics must produce sane outputs.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from benchmarks.shared import metrics as M
from benchmarks.shared import normalization as norm
from benchmarks.shared import library_type as lt
from benchmarks.shared import batch_effects as be
from benchmarks.shared.method_registry import bulk_methods, spatial_methods
from benchmarks.shared.method_selection import build_compatibility_table
from benchmarks.shared.synthetic import toy_reference, toy_bulk, toy_spatial


# --- metrics ----------------------------------------------------------------

def test_accuracy_metrics_perfect():
    df = pd.DataFrame({"a": [0.6, 0.4], "b": [0.4, 0.6]}, index=["s0", "s1"])
    m = M.accuracy_metrics(df, df.copy())
    assert m["rmse"] == pytest.approx(0.0)
    assert m["pearson"] == pytest.approx(1.0)


def test_pairwise_method_correlation():
    a = pd.DataFrame({"x": [0.5, 0.5], "y": [0.5, 0.5]}, index=["s0", "s1"])
    out = M.pairwise_method_correlation({"m1": a, "m2": a.copy()})
    assert "m1__vs__m2" in out


def test_proportion_entropy_and_near_zero():
    df = pd.DataFrame({"a": [1.0, 0.5], "b": [0.0, 0.5]}, index=["s0", "s1"])
    ent = M.proportion_entropy(df)
    assert ent.iloc[0] < ent.iloc[1]
    assert 0.0 <= M.near_zero_fraction(df) <= 1.0


# --- normalization detection ------------------------------------------------

def test_detect_counts_vs_cpm_vs_log():
    rng = np.random.default_rng(0)
    counts = rng.integers(0, 500, size=(50, 8)).astype(float)
    assert norm.detect_normalization_status(counts) == "counts"
    cpm = counts / counts.sum(axis=1, keepdims=True) * 1e6
    assert norm.detect_normalization_status(cpm) in ("cpm", "counts")  # row-sum heuristic
    logn = np.log1p(cpm / 100.0)
    assert norm.detect_normalization_status(logn) in ("log_normalized", "unknown")


def test_convert_counts_to_cpm():
    df = pd.DataFrame({"s0": [1.0, 1.0], "s1": [2.0, 2.0]}, index=["g0", "g1"])
    cpm = norm.convert_counts_to_cpm(df)
    assert cpm["s0"].sum() == pytest.approx(1e6)


# --- library-type detection -------------------------------------------------

def test_library_type_detection_scrna_and_sn_and_mixed():
    sc = pd.DataFrame({"assay": ["10x 3' v3"] * 4, "cell_type": list("abcd")})
    assert lt.detect_library_type(sc)["overall"] == "scrna_10x_3p"
    sn = pd.DataFrame({"suspension_type": ["nucleus"] * 4, "cell_type": list("abcd")})
    assert lt.detect_library_type(sn)["overall"] == "single_nucleus"
    mixed = pd.DataFrame({"assay": ["10x 3'", "10x 3'", "single-nucleus", "nucleus"],
                          "cell_type": list("abcd")})
    out = lt.detect_library_type(mixed)
    assert out["is_mixed"] and out["overall"] == "mixed"


def test_cell_type_by_library_type():
    obs = pd.DataFrame({"lib": ["sc", "sc", "sn", "sn"],
                        "ct": ["T", "T", "T", "B"]})
    frac = lt.cell_type_by_library_type(obs, "lib", "ct")
    assert "dominant_library_fraction" in frac.columns
    assert frac.loc["B", "dominant_library_fraction"] == pytest.approx(1.0)


def test_simulate_library_shift_preserves_cpm():
    ref, _ = toy_reference()
    shifted = lt.simulate_library_shift(ref.as_R_cpm(), list(ref.gene_names),
                                        kind="single_nucleus", seed=1)
    np.testing.assert_allclose(shifted.sum(axis=1), 1e6, rtol=1e-4)


# --- batch diagnostics ------------------------------------------------------

def test_detect_batch_columns_and_confounding():
    obs = pd.DataFrame({"donor": ["d1", "d1", "d2", "d2"],
                        "cell_type": ["T", "T", "T", "B"]})
    assert "donor" in be.detect_batch_columns(obs)
    conf = be.compute_celltype_batch_confounding(obs, "cell_type", "donor")
    assert bool(conf.loc["B", "single_batch_only"]) is True
    assert 0.0 <= be.compute_batch_mixing_score(obs, "donor") <= 1.0


def test_marker_batch_stability_handles_missing_donor_cv():
    ref, _ = toy_reference()
    stab = be.compute_marker_batch_stability(ref)
    assert "mean_marker_donor_cv" in stab.columns


# --- method registry / availability ----------------------------------------

def test_external_methods_absent_are_skipped_not_failed():
    ref, mapping = toy_reference()
    bulk, _ = toy_bulk(ref)
    scenario = {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
                "normalization_status": "counts"}
    for m in bulk_methods(include_external=True):
        if m.external and not m.is_available():
            res = m.run(scenario)
            assert res.status in ("skipped",)
            assert res.install_hint  # actionable hint recorded


def test_internal_bulk_baselines_run_on_toy():
    ref, mapping = toy_reference()
    bulk, truth = toy_bulk(ref)
    scenario = {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
                "normalization_status": "counts"}
    from benchmarks.bulk.methods.nnls_baseline import NNLSBaseline
    from benchmarks.bulk.methods.tissueresolve import (
        TissueResolveBulkFlat, TissueResolveBulkHierarchical,
    )
    for m in (NNLSBaseline(), TissueResolveBulkFlat(), TissueResolveBulkHierarchical()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = m.run(scenario)
        assert res.status == "success"
        np.testing.assert_allclose(res.predictions.sum(axis=1), 1.0, atol=1e-6)


def test_internal_spatial_baselines_run_on_toy():
    ref, mapping = toy_reference()
    sc, truth = toy_spatial(ref)
    sc["reference"] = ref
    sc["hierarchy_mapping"] = mapping
    sc["normalization_status"] = "counts"
    from benchmarks.spatial.methods.nnls_spot_baseline import NNLSSpotBaseline
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = NNLSSpotBaseline().run(sc)
    assert res.status == "success"
    np.testing.assert_allclose(res.predictions.sum(axis=1), 1.0, atol=1e-6)


def test_compatibility_table():
    ref, mapping = toy_reference()
    bulk, _ = toy_bulk(ref)
    scenario = {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
                "normalization_status": "counts", "reference_library_type": "scrna_10x_3p"}
    tab = build_compatibility_table(bulk_methods(), scenario)
    assert "available" in tab.columns
    assert "TissueResolve_hierarchical" in tab.index


def test_cibersortx_export_status(tmp_path, monkeypatch):
    ref, mapping = toy_reference()
    bulk, _ = toy_bulk(ref)
    from benchmarks.bulk.methods.cibersortx_export import CIBERSORTxExport
    res = CIBERSORTxExport().run({"reference": ref, "bulk": bulk})
    assert res.status == "exported_not_run"


def test_dry_run_offline():
    """Dry-run main() returns 0 without touching the network."""
    from benchmarks.bulk import run_bulk_benchmark as B
    from benchmarks.spatial import run_spatial_benchmark as S
    assert B.main(["--dry-run"]) == 0
    assert S.main(["--dry-run"]) == 0


def test_readme_documents_new_capabilities():
    from pathlib import Path
    readme = (Path(__file__).resolve().parents[2] / "README.md").read_text().lower()
    for kw in ("benchmark", "normalization", "protocol", "library", "batch",
               "hierarchical broad→fine", "resolution-mode auto"):
        assert kw.lower() in readme, f"README missing mention of {kw!r}"


def test_docs_pages_exist():
    from pathlib import Path
    docs = Path(__file__).resolve().parents[2] / "docs"
    for name in ("benchmarking.md", "normalization_and_protocols.md",
                 "batch_effects.md", "library_type_references.md"):
        assert (docs / name).exists(), f"missing docs/{name}"


# --- fair hierarchical / family-level / unresolved-aware metrics ------------

def test_family_level_metrics_aggregates_unresolved_columns():
    mapping = {"X1": "FamX", "X2": "FamX", "Y1": "FamY", "Y2": "FamY"}
    true = pd.DataFrame({"X1": [0.3, 0.1], "X2": [0.2, 0.2],
                         "Y1": [0.3, 0.4], "Y2": [0.2, 0.3]}, index=["s0", "s1"])
    est = pd.DataFrame({"X1": [0.25, 0.15], "X2": [0.25, 0.15],
                        "Y1": [0, 0], "Y2": [0, 0],
                        "unresolved_FamY": [0.5, 0.7]}, index=["s0", "s1"])
    fam = M.family_level_metrics(true, est, mapping)
    assert fam["pearson"] == pytest.approx(1.0)
    assert fam["rmse"] == pytest.approx(0.0, abs=1e-9)


def test_unresolved_aware_metrics_rewards_appropriate_abstention():
    mapping = {"X1": "FamX", "X2": "FamX", "Y1": "FamY", "Y2": "FamY"}
    est = pd.DataFrame({"X1": [0.25], "X2": [0.25], "Y1": [0], "Y2": [0],
                        "unresolved_FamY": [0.5]})
    u = M.unresolved_aware_metrics(est, mapping, high_risk_families={"FamY"})
    assert u["abstains"] is True
    assert u["unresolved_precision"] == pytest.approx(1.0)
    assert u["unresolved_recall"] == pytest.approx(1.0)
    assert u["unresolved_mass_fraction"] == pytest.approx(0.5)


def test_hierarchical_not_judged_only_by_fine_level():
    """Hierarchical fair metrics must include family-level + resolvable-only fine."""
    mapping = {"X1": "FamX", "X2": "FamX", "Y1": "FamY", "Y2": "FamY"}
    true = pd.DataFrame({"X1": [0.3], "X2": [0.2], "Y1": [0.3], "Y2": [0.2]})
    est = pd.DataFrame({"X1": [0.3], "X2": [0.2], "Y1": [0], "Y2": [0],
                        "unresolved_FamY": [0.5]})
    fair = M.hierarchical_fair_metrics(true, est, mapping, {"FamY"})
    assert "family_pearson" in fair and "fine_resolvable_pearson" in fair
    # resolvable-only fine excludes the abstained FamY subtypes
    assert fair["abstains"] is True


# --- capability matrix / analysis / import ----------------------------------

def test_capability_matrix_generated():
    from benchmarks.shared.analysis import capability_matrix
    cap = capability_matrix(bulk_methods())
    assert "hierarchical" in cap.columns
    assert "open_local_execution" in cap.columns
    assert "TissueResolve_hierarchical" in cap.index


def test_categorize_separates_executed_exported_skipped():
    from benchmarks.shared.analysis import categorize
    ref, mapping = toy_reference()
    bulk, _ = toy_bulk(ref)
    scenario = {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
                "normalization_status": "counts"}
    results = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for m in bulk_methods(include_external=True):
            results.append(m.run(scenario))
    cats = categorize(results)
    assert cats["executed"], "expected executed methods"
    # CIBERSORTx export → exported (NOT counted as benchmarked)
    assert any(r.method == "CIBERSORTx_export" for r in cats["exported"])


def test_at_least_three_non_tissueresolve_executed():
    ref, mapping = toy_reference()
    bulk, _ = toy_bulk(ref)
    scenario = {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
                "normalization_status": "counts"}
    executed_nontr = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for m in bulk_methods(include_external=False):
            if not m.name.startswith("TissueResolve") and m.run(scenario).status == "success":
                executed_nontr += 1
    assert executed_nontr >= 3  # NNLS, WNNLS, MarkerOnly


def test_imported_results_loadable(tmp_path):
    from benchmarks.shared.imported import ImportedMethod
    p = tmp_path / "ext.tsv"
    pd.DataFrame({"X1": [0.5, 0.5], "X2": [0.5, 0.5]},
                 index=["s0", "s1"]).to_csv(p, sep="\t")
    m = ImportedMethod("ExtTool", "bulk", p)
    res = m.run({})
    assert res.status == "success"
    assert res.metadata["executed_or_exported"] == "executed_imported"
    np.testing.assert_allclose(res.predictions.sum(axis=1), 1.0)


def test_import_external_results_cli(tmp_path):
    from benchmarks.import_external_results import main as import_main
    pred = tmp_path / "p.tsv"
    pd.DataFrame({"A": [0.6], "B": [0.4]}, index=["s0"]).to_csv(pred, sep="\t")
    out = tmp_path / "ExtTool.tsv"
    rc = import_main(["--method", "ExtTool", "--modality", "bulk",
                      "--predictions", str(pred), "--out", str(out)])
    assert rc == 0 and out.exists()


# --- unified report / audit / README ----------------------------------------

def test_unified_report_builder_links_sections(tmp_path):
    from tissueresolve.report.unified import Section, build_unified_report
    secs = [Section("bulk", "Bulk results", "<p>x</p>",
                    links=[("detail", "bulk/report.html")]),
            Section("benchmark", "Benchmark comparison", "<p>y</p>")]
    out = build_unified_report(tmp_path / "report.html", secs)
    html = out.read_text()
    assert "href='#bulk'" in html and "href='#benchmark'" in html
    assert "bulk/report.html" in html


def test_audit_doc_exists():
    from pathlib import Path
    assert (Path(__file__).resolve().parents[2] / "docs" /
            "project_structure_audit.md").exists()


def test_readme_documents_results_location_and_simplified_workflow():
    from pathlib import Path
    readme = (Path(__file__).resolve().parents[2] / "README.md").read_text().lower()
    assert "where are my results" in readme
    assert "report.html" in readme
    assert "not intended to replace every specialized method" in readme
