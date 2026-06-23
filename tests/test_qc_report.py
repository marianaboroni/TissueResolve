"""Tests for the QC-first report logic + structure (Task A).

Pure-logic tests run offline; structural tests validate the built report if
present (skip otherwise, keeping the default suite offline).
"""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "build_qc_first_report", REPO / "benchmarks" / "build_qc_first_report.py")
qc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qc)

QCDIR = REPO / "benchmarks" / "outputs" / "qc_report"


# ---- pure logic ----
def test_trusted_resolution_uses_full_panel_not_cell_auroc():
    mix = pd.DataFrame({
        "family": ["B/Plasma", "T/NK", "Epithelial"],
        "median_shared_lineage_fraction": [0.86, 0.99, 0.94],
        "median_condition_number": [2.8, 8.0, 4.6],
    })
    cond = pd.DataFrame({"family": ["B/Plasma", "T/NK", "Epithelial"],
                         "conditional_rmse_fullpanel": [0.10, 0.40, 0.12]})
    t = qc.trusted_resolution(mix, cond).set_index("family")
    # high shared-lineage (T/NK 0.99) → broad_only regardless of cell AUROC
    assert t.loc["T/NK", "trusted_resolution"] == "broad_only"
    # low shared-lineage + low conditional RMSE → selected_fine
    assert t.loc["B/Plasma", "trusted_resolution"] == "selected_fine"
    # basis explicitly excludes cell-AUROC
    assert "NOT cell-AUROC" in t.loc["B/Plasma", "basis"]


def test_no_family_marked_full_fine():
    mix = pd.DataFrame({"family": ["X"], "median_shared_lineage_fraction": [0.5],
                        "median_condition_number": [2.0]})
    cond = pd.DataFrame({"family": ["X"], "conditional_rmse_fullpanel": [0.05]})
    t = qc.trusted_resolution(mix, cond)
    assert set(t["trusted_resolution"]) <= {"selected_fine", "broad_only"}  # never 'fine'


def test_status_card_has_metric_threshold_and_affected():
    html = qc.card("CAUTION", "warn", "Fine reliability", "RMSE=0.31",
                   "fragile", "fine predictions")
    assert "RMSE=0.31" in html and "fragile" in html and "fine predictions" in html


# ---- structural (skip if report not built) ----
def _need_report():
    if not (QCDIR / "report.html").exists():
        pytest.skip("QC report not built (run build_qc_first_report.py --run-real-data)")


def test_qc_sections_before_predictions():
    _need_report()
    html = (QCDIR / "report.html").read_text()
    i_refq = html.index("2. Reference quality")
    i_predq = html.index("5. Prediction quality")
    i_bulk = html.index("6. Final bulk predictions")
    assert i_refq < i_predq < i_bulk    # quality precedes predictions


def test_three_identifiability_levels_distinguished():
    _need_report()
    html = (QCDIR / "report.html").read_text()
    assert "three distinct" in html.lower()      # 3 evidence levels, in section 3
    assert "do not merge" in html.lower()


def test_fine_predictions_gated_and_no_spatial_accuracy_claim():
    _need_report()
    html = (QCDIR / "report.html").read_text()
    assert "GATED OUT" in html or "gated" in html.lower()
    assert "not accuracy" in html.lower()


def test_every_main_figure_has_caption_and_source_data():
    _need_report()
    man = pd.read_csv(QCDIR / "figures" / "figure_manifest.tsv", sep="\t")
    assert len(man) >= 8
    assert man["caption"].astype(str).str.len().min() > 10
    assert man["source_data"].astype(str).str.len().min() > 3
    assert man["purpose"].astype(str).str.len().min() > 5


def test_donor_support_present():
    _need_report()
    man = pd.read_csv(QCDIR / "figures" / "figure_manifest.tsv", sep="\t")
    figs = set(man["figure"])
    assert "figR1_family_support" in figs and "figR3_donor_coverage" in figs


def test_technical_appendix_holds_heavy_material():
    _need_report()
    appx = (QCDIR / "technical_appendix.html").read_text()
    main = (QCDIR / "report.html").read_text()
    # full pairwise table referenced in appendix, not embedded as a main figure
    assert "pairs.tsv" in appx
    man = pd.read_csv(QCDIR / "figures" / "figure_manifest.tsv", sep="\t")
    assert not any("pairwise" in str(p).lower() and "heatmap" in str(p).lower()
                   for p in man["figure"])
