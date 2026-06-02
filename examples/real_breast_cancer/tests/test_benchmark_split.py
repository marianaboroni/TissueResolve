"""
Offline tests for the benchmark split + main-summary-figure relocation.

Contract-level (caption/constant/figure-card) tests plus checks against the
generated report/manifest when they are present (skipped otherwise, so CI stays
offline).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("plotly")

from tissueresolve.report.figures import FigureManifest

REPORT = Path("examples/real_breast_cancer/outputs/report.html")
MANIFEST = Path("examples/real_breast_cancer/outputs/figures/figure_manifest.tsv")

_BULK_BENCH = ["bulk_method_status_summary", "bulk_fine_accuracy_leaderboard",
               "bulk_family_accuracy_leaderboard", "bulk_rmse_mae_comparison",
               "bulk_runtime_comparison"]
_SPATIAL_BENCH = ["spatial_method_status_summary", "spatial_concordance_heatmap",
                  "spatial_structure_metrics_summary", "spatial_runtime_comparison",
                  "spatial_output_completeness_summary"]


@pytest.fixture
def rep(load_script):
    return load_script("07_generate_reports.py")


# --- constants / captions ---------------------------------------------------

def test_summary_figs_excluded_constant(rep):
    assert rep._SUMMARY_FIGS == {"bulk_main_summary_figure",
                                 "spatial_main_summary_figure"}


def test_bulk_benchmark_captions_mention_ground_truth(rep):
    for stem in _BULK_BENCH:
        cap = rep._CAPTIONS[stem]["caption"].lower()
        # accuracy figures must justify themselves via ground truth;
        # runtime/status are cost/robustness and need not claim accuracy
        if "accuracy" in stem or "rmse" in stem:
            assert "ground truth" in cap


def test_spatial_benchmark_captions_disclaim_accuracy(rep):
    for stem in _SPATIAL_BENCH:
        cap = rep._CAPTIONS[stem]["caption"].lower()
        assert ("not accuracy" in cap or "no spot-level ground truth" in cap
                or "no ground truth" in cap)


def test_no_spatial_accuracy_vs_truth_caption(rep):
    """No spatial benchmark caption claims Pearson/RMSE vs truth."""
    for stem in _SPATIAL_BENCH:
        cap = rep._CAPTIONS[stem]["caption"].lower()
        assert "rmse" not in cap and "pearson" not in cap


# --- figure-card exclude / only ---------------------------------------------

def _toy_fig(figdir: Path, stem: str):
    figdir.mkdir(parents=True, exist_ok=True)
    (figdir / f"{stem}.html").write_text("<html><body>fig</body></html>")
    (figdir / f"{stem}.data.tsv").write_text("a\tb\n1\t2\n")


def test_main_summary_excluded_from_primary(tmp_path, rep):
    figdir = tmp_path / "figures"
    _toy_fig(figdir, "bulk_composition_clustered_barplot")
    _toy_fig(figdir, "bulk_main_summary_figure")
    m = FigureManifest()
    primary = rep._figure_cards(figdir, tmp_path, m, "bulk", rep._CAPTIONS,
                                exclude=rep._SUMMARY_FIGS)
    assert "bulk_composition_clustered_barplot" in primary
    assert "bulk_main_summary_figure" not in primary   # excluded from primary
    appendix = rep._figure_cards(figdir, tmp_path, FigureManifest(), "bulk",
                                 rep._CAPTIONS, only={"bulk_main_summary_figure"})
    assert "bulk_main_summary_figure" in appendix       # shown only in appendix


# --- generated report / manifest (skip if not present) ----------------------

@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_report_has_separate_benchmark_sections():
    html = REPORT.read_text()
    assert "9. Bulk benchmark" in html
    assert "10. Spatial benchmark" in html
    assert "Benchmark comparison" not in html          # old combined section gone
    # main summary figures relocated to collapsible technical appendix
    assert "technical_appendix.html" in html   # heavy/composite figs relocated to the appendix file


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_spatial_benchmark_states_no_ground_truth():
    html = REPORT.read_text().lower()
    i = html.find("id='spatial_benchmark'")      # the section body, not the nav link
    assert i > 0
    seg = html[i:i + 4000]
    assert "no spot-level ground truth" in seg or "not accuracy" in seg


@pytest.mark.skipif(not MANIFEST.exists(), reason="manifest not generated")
def test_manifest_distinguishes_benchmark_sections():
    m = pd.read_csv(MANIFEST, sep="\t")
    secs = set(m["section"])
    assert "bulk_benchmark" in secs and "spatial_benchmark" in secs
    assert "benchmark_scorecard" in secs
    # bulk and spatial benchmark figures never share a section
    bulk = set(m[m.section == "bulk_benchmark"]["figure_id"])
    spat = set(m[m.section == "spatial_benchmark"]["figure_id"])
    assert not (bulk & spat)
