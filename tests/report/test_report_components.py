"""
Offline tests for the redesigned report design system (components, glossary,
figure manifest, unified report).  These exercise the component contracts so the
generated reports are navigable, labelled, and honest.
"""
from __future__ import annotations

import pandas as pd
import pytest

from tissueresolve.report import components as C
from tissueresolve.report import glossary as G
from tissueresolve.report.figures import FigureManifest, FigureRecord
from tissueresolve.report.unified import Section, build_unified_report


# --- design system / components ---------------------------------------------

def test_style_has_required_css_classes():
    from tissueresolve.report.style import REPORT_CSS
    for cls in ("report-container", "sidebar-nav", "section-card", "metric-card",
                "figure-card", "methods-card", "variable-card", "warning-card",
                "limitation-card", "source-data-link", "collapsible-table",
                "interpretation-guide", "figure-caption", "figure-legend",
                "status-pass", "status-caution", "status-warning", "status-critical",
                "how-to-read", "checklist", "status-card"):
        assert f".{cls}" in REPORT_CSS, f"missing CSS class .{cls}"


def test_metric_card_and_grid():
    html = C.metric_grid({"Samples": 6, "Genes": (5000, "shared")})
    assert "metric-card" in html and ">6<" in html and "5000" in html


def test_figure_card_has_title_caption_and_source_link():
    html = C.figure_card(title="Composition", subtitle="by family",
                         caption="Stacked barplot of RNA-derived proportions.",
                         how_to_read="Check clustering.",
                         source_links=[("source data", "x.data.tsv")])
    assert "fig-title" in html and "Composition" in html
    assert "Figure." in html and "Stacked barplot" in html
    assert "source-data-link" in html and "x.data.tsv" in html
    assert "how-to-read" in html


def test_warning_box_empty_when_no_warnings():
    assert C.warning_box([]) == ""
    assert C.warning_box(["", "  "]) == ""
    box = C.warning_box(["overlap is low"])
    assert "overlap is low" in box and "No warnings" not in box


def test_collapsible_table_renders_details():
    html = C.collapsible_table("Full matrix", "<table></table>")
    assert html.startswith("<details") and "Full matrix" in html


def test_status_badge_classes():
    assert "status-pass" in C.status_badge("PASS")
    assert "status-critical" in C.status_badge("FAIL")


def test_decision_status_grid_renders_badges():
    html = C.decision_status_grid([
        ("Reference suitability", "WARNING", "is the reference good enough?"),
        ("Bulk results", "yes"),
    ])
    assert "status-warning" in html and "Reference suitability" in html
    assert "is the reference good enough?" in html
    # a two-element item (no sub) still renders
    assert "Bulk results" in html


def test_checklist_marks_states():
    html = C.checklist([("done", True), ("not done", False), ("n/a", None)])
    assert "check-ok" in html and "check-bad" in html and "check-na" in html
    assert "done" in html and "not done" in html and "n/a" in html


def test_variable_dictionary_renders_terms():
    html = C.variable_dictionary(G.subset(["Pearson correlation", "RMSE"]))
    assert "variable-card" in html and "Pearson correlation" in html


# --- glossary ----------------------------------------------------------------

def test_glossary_contains_required_terms():
    for term in ("Pearson correlation", "RMSE", "entropy", "Moran's I",
                 "unresolved mass", "spillover", "separability", "gene overlap"):
        assert term in G.GLOSSARY and G.GLOSSARY[term]


# --- figure manifest ---------------------------------------------------------

def test_figure_manifest_columns_and_save(tmp_path):
    m = FigureManifest()
    m.add(FigureRecord("f1", "bulk", "Composition barplot",
                       source_data_path="figures/f1.data.tsv"))
    p = m.save(tmp_path / "figures")
    assert p.exists()
    df = pd.read_csv(p, sep="\t")
    for col in ("figure_id", "section", "title", "html_path", "png_path", "svg_path",
                "pdf_path", "source_data_path", "caption_path", "caption",
                "methodology", "variables_defined", "status", "reason_if_missing"):
        assert col in df.columns
    assert df.loc[0, "figure_id"] == "f1"


def test_figure_manifest_records_missing_data(tmp_path):
    m = FigureManifest()
    m.add(FigureRecord("missing_fig", "bulk", "Solver CV",
                       status="missing_data",
                       reason_if_missing="solver=auto CV scores not saved"))
    df = pd.read_csv(m.save(tmp_path / "figures"), sep="\t")
    row = df[df["figure_id"] == "missing_fig"].iloc[0]
    assert row["status"] == "missing_data"
    assert "not saved" in str(row["reason_if_missing"])


# --- unified report ----------------------------------------------------------

def _demo_report(tmp_path, warnings=False):
    secs = [
        Section("summary", "1. Executive summary",
                C.metric_grid({"Modality": "bulk + spatial"})
                + C.estimate_note("RNA-derived estimates, not cell fractions.")),
        Section("reference", "2. Reference quality",
                C.methodology_summary(["built from h5ad"])
                + C.figure_card(title="Reference composition",
                                caption="Reference composition barplot.",
                                source_links=[("source data", "r.data.tsv")])
                + C.collapsible_table("Full table", "<table></table>")
                + C.variable_dictionary(G.subset(["separability"]))),
        Section("spatial", "5. Spatial deconvolution",
                C.methodology_summary([
                    "Real Visium has no ground truth — concordance/structure, "
                    "not accuracy."])),
        Section("benchmark", "8. Benchmark comparison",
                C.methodology_summary([
                    "Only executed or imported methods are ranked; exported-only "
                    "and skipped tools are listed but not scored."])),
        Section("warnings", "9. Warnings & limitations",
                (C.warning_box(["low overlap"]) if warnings
                 else "<p class='muted'>No warnings recorded for this run.</p>")
                + C.limitation_box(["not absolute cell fractions"])),
    ]
    return build_unified_report(tmp_path / "report.html", secs, subtitle="t")


def test_unified_report_has_navigation(tmp_path):
    html = _demo_report(tmp_path).read_text()
    assert "sidebar-nav" in html
    for a in ("#summary", "#reference", "#spatial", "#benchmark", "#warnings"):
        assert f"href='{a}'" in html


def test_unified_sections_have_methodology_boxes(tmp_path):
    html = _demo_report(tmp_path).read_text()
    # reference / spatial / benchmark sections each carry a methodology box
    assert html.count("methods-card") >= 3


def test_every_figure_has_title_caption_sourcelink(tmp_path):
    html = _demo_report(tmp_path).read_text()
    # one figure card present, fully labelled
    assert "fig-title" in html and "Figure." in html and "source-data-link" in html


def test_visual_summary_before_raw_table(tmp_path):
    """A section must present visual summaries (metric/figure cards) before any
    full raw table (collapsible)."""
    html = _demo_report(tmp_path).read_text()
    assert "figure-card" in html and "collapsible-table" in html
    assert html.index("figure-card") < html.index("collapsible-table")


def test_no_false_no_warnings_when_warnings_exist(tmp_path):
    html = _demo_report(tmp_path, warnings=True).read_text()
    assert "low overlap" in html
    assert "No warnings recorded" not in html


def test_spatial_section_states_no_ground_truth(tmp_path):
    html = _demo_report(tmp_path).read_text()
    assert "no ground truth" in html.lower()
    assert "not accuracy" in html.lower()


def test_benchmark_section_distinguishes_tool_categories(tmp_path):
    html = _demo_report(tmp_path).read_text().lower()
    for word in ("executed", "imported", "exported", "skipped"):
        assert word in html


def test_color_map_deterministic():
    from tissueresolve.plotting import palette as p
    fine = ["CD8 T cell", "CD4 T cell", "macrophage", "Other"]
    mp = {"CD8 T cell": "T/NK", "CD4 T cell": "T/NK", "macrophage": "Myeloid"}
    a = p.build_hierarchical_color_map(fine, mp)
    b = p.build_hierarchical_color_map(fine, mp)
    assert dict(zip(a["fine_cell_type"], a["color_hex"])) == \
        dict(zip(b["fine_cell_type"], b["color_hex"]))
