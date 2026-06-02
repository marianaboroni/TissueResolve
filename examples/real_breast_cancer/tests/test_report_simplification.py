"""
Offline tests for report simplification: heavy/technical figures are routed out
of the primary grid into collapsible technical appendices (report-simplification
audit, Part 5), and the simplified structure holds.
"""
from __future__ import annotations

from pathlib import Path

import pytest

REPORT = Path("examples/real_breast_cancer/outputs/report.html")


@pytest.fixture
def rep(load_script):
    return load_script("07_generate_reports.py")


def test_appendix_figure_sets_defined(rep):
    assert {"bulk_separability_heatmap", "bulk_spillover_heatmap",
            "spillover_network"} <= rep._BULK_APPENDIX_FIGS
    assert {"spatial_separability_heatmap", "spatial_spillover_heatmap",
            "spatial_spot_pie_charts"} <= rep._SPATIAL_APPENDIX_FIGS
    assert "signature_matrix_heatmap" in rep._SIGNATURE_APPENDIX_FIGS


def test_figure_cards_excludes_appendix_from_primary(tmp_path, rep):
    from tissueresolve.report.figures import FigureManifest
    figdir = tmp_path / "figures"
    figdir.mkdir(parents=True)
    for stem in ("bulk_composition_clustered_barplot", "bulk_separability_heatmap",
                 "spillover_network"):
        (figdir / f"{stem}.html").write_text("<html>f</html>")
        (figdir / f"{stem}.data.tsv").write_text("a\tb\n1\t2\n")
    primary = rep._figure_cards(figdir, tmp_path, FigureManifest(), "bulk",
                                rep._CAPTIONS, exclude=rep._BULK_APPENDIX_FIGS)
    assert "bulk_composition_clustered_barplot" in primary
    assert "bulk_separability_heatmap" not in primary    # routed to appendix
    assert "spillover_network" not in primary
    appendix = rep._figure_cards(figdir, tmp_path, FigureManifest(), "bulk",
                                 rep._CAPTIONS, only=rep._BULK_APPENDIX_FIGS)
    assert "bulk_separability_heatmap" in appendix


# --- generated report (skip when not present) -------------------------------

def _section_pos(html):
    import re
    return {m.group(1): m.start() for m in re.finditer(r"id='([a-z_]+)'", html)}


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_heavy_figures_in_technical_appendix_not_primary():
    html = REPORT.read_text()
    assert "Technical appendix" in html          # collapsible appendix present
    pos = _section_pos(html)
    # each heavy figure sits within its section but AFTER the section's
    # "Technical appendix"/"exploratory" collapsible marker
    def in_appendix(stem, sec_a, sec_b, marker):
        i = html.find(stem)
        if i < 0:
            return None
        # the appendix collapsible marker for this section precedes the figure
        seg = html[pos.get(sec_a, 0):pos.get(sec_b, len(html))]
        return (marker in seg) and (seg.find(marker) < seg.find(stem))
    assert in_appendix("bulk_separability_heatmap", "bulk", "spatial",
                       "Technical appendix (bulk)")
    assert in_appendix("spatial_spot_pie_charts", "spatial", "hierarchical",
                       "Technical / exploratory figures (spatial)")


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_no_empty_figure_cards_and_benchmarks_separated():
    html = REPORT.read_text()
    assert "fig-body'></div>" not in html
    assert "9. Bulk benchmark" in html and "10. Spatial benchmark" in html
    # spatial benchmark still disclaims accuracy
    i = html.find("id='spatial_benchmark'")
    seg = html[i:i + 4000].lower()
    assert "not accuracy" in seg or "no spot-level ground truth" in seg
