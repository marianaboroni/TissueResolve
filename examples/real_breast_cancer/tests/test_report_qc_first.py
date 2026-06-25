"""
Tests for the QC-first report restructuring (REPORT_RESTRUCTURING_PLAN.md):
10-section layout, QC before predictions, ~18 reliability/interpretation figures
in the main report, everything technical/exploratory in the appendix.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPORT = Path("examples/real_breast_cancer/outputs/report.html")
APPENDIX = Path("examples/real_breast_cancer/outputs/technical_appendix.html")

# Figures that MUST be in the main report (18) and figures that MUST NOT be.
_MAIN = [
    "reference_suitability_components", "reference_broad_family_composition",
    "gene_overlap_by_modality", "top_confusable_pairs",
    "within_vs_between_family_separability", "trusted_resolution_summary",
    "unresolved_mass_by_family", "bulk_composition_clustered_barplot",
    "bulk_composition_heatmap", "bulk_qc_summary", "bulk_uncertainty_plot",
    "spatial_mean_composition_barplot", "he_dominant_cell_type",
    "spatial_dominant_cell_type_map", "spatial_abundance_maps",
    "spatial_morans_i_barplot", "bulk_fine_accuracy_leaderboard",
    "spatial_structure_metrics_summary",
]
_NOT_MAIN = [
    "hierarchy_map", "separability_distribution", "signature_matrix_heatmap",
    "marker_support_by_family", "reference_celltype_imbalance",
    "reference_fine_subpopulation_support", "bulk_main_summary_figure",
    "spatial_main_summary_figure", "spatial_spot_pie_charts",
    "composite_scorecard", "bulk_runtime_comparison", "he_spots_check",
]
_EXPECTED_SECTIONS = [
    "1. Executive decision summary",
    "2. Reference and signature quality",
    "3. Input compatibility",
    "4. Trusted resolution and uncertainty",
    "5. Final bulk predictions",
    "6. Final spatial predictions",
    "7. Bulk benchmark summary",
    "8. Spatial benchmark summary",
    "9. Warnings and limitations",
    "10. Methods and source data",
]


@pytest.fixture
def rep(load_script):
    return load_script("07_generate_reports.py")


# --- routing constants (always runnable) ------------------------------------

def test_routing_sets_keep_tweaks(rep):
    # user tweaks: these two stay in MAIN (not in any appendix set)
    appendix_all = (rep._BULK_APPENDIX_FIGS | rep._SPATIAL_APPENDIX_FIGS
                    | rep._SIGNATURE_APPENDIX_FIGS | rep._RESOLUTION_APPENDIX_FIGS
                    | rep._REFERENCE_APPENDIX_FIGS)
    assert "bulk_composition_heatmap" not in appendix_all
    assert "spatial_dominant_cell_type_map" not in appendix_all
    # and these move to APPENDIX
    assert "hierarchy_map" in rep._SIGNATURE_APPENDIX_FIGS
    assert "separability_distribution" in rep._RESOLUTION_APPENDIX_FIGS
    # per-type H&E panels routed by prefix
    assert "he_abundance_*" in rep._SPATIAL_APPENDIX_FIGS


def test_fig_match_prefix(rep):
    assert rep._fig_match("he_abundance_macrophage", {"he_abundance_*"})
    assert not rep._fig_match("spatial_abundance_maps", {"he_abundance_*"})
    assert rep._fig_match("x", {"x"})


# --- generated main report (skip if not generated) --------------------------

@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_exactly_ten_sections_in_order():
    html = REPORT.read_text()
    titles = re.findall(r"<h2>([^<]+)</h2>", html)
    # strip HTML entities for arrow etc.; compare against expected prefixes
    assert len(titles) == 10, f"expected 10 sections, got {len(titles)}: {titles}"
    for got, exp in zip(titles, _EXPECTED_SECTIONS):
        assert got.strip().startswith(exp), f"{got!r} != {exp!r}"


@pytest.mark.skipif(not (REPORT.exists() and APPENDIX.exists()),
                    reason="report/appendix not generated")
def test_main_figures_present_appendix_figures_absent():
    rep_html = REPORT.read_text()
    appx = APPENDIX.read_text()
    for stem in _MAIN:
        assert stem in rep_html, f"{stem} must be in the main report"
    for stem in _NOT_MAIN:
        assert stem not in rep_html, f"{stem} must NOT be in the main report"
        assert stem in appx, f"{stem} must be in the technical appendix"


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_main_figure_count_in_target_range():
    n = REPORT.read_text().count("class='figure-card'")
    assert 15 <= n <= 22, f"main report has {n} figure cards (target 15-22)"


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_no_generic_main_summary_in_main():
    html = REPORT.read_text()
    assert "bulk_main_summary_figure" not in html
    assert "spatial_main_summary_figure" not in html
