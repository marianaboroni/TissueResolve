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


# --- separate technical appendix file (skip when not generated) -------------

APPENDIX = Path("examples/real_breast_cancer/outputs/technical_appendix.html")
_HEAVY = ("bulk_separability_heatmap", "bulk_spillover_heatmap",
          "spatial_separability_heatmap", "spatial_spillover_heatmap",
          "signature_matrix_heatmap", "spatial_spot_pie_charts")


@pytest.mark.skipif(not (REPORT.exists() and APPENDIX.exists()),
                    reason="report/appendix not generated")
def test_separate_appendix_exists_and_cross_links():
    rep = REPORT.read_text()
    appx = APPENDIX.read_text()
    assert "technical_appendix.html" in rep        # main → appendix
    assert "report.html" in appx                   # appendix → main


@pytest.mark.skipif(not (REPORT.exists() and APPENDIX.exists()),
                    reason="report/appendix not generated")
def test_heavy_figures_in_appendix_not_in_main():
    rep = REPORT.read_text()
    appx = APPENDIX.read_text()
    for stem in _HEAVY:
        assert stem not in rep, f"{stem} should not be in the main report"
        assert stem in appx, f"{stem} should be in the technical appendix"


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_main_report_clean_and_benchmarks_separated():
    html = REPORT.read_text()
    assert "fig-body'></div>" not in html                  # no empty cards
    assert "All source-data tables" not in html            # full dump moved out
    # QC-first 10-section layout: bulk benchmark = §7, spatial benchmark = §8
    assert "7. Bulk benchmark" in html and "8. Spatial benchmark" in html
    i = html.find("id='spatial_benchmark'")
    seg = html[i:i + 4000].lower()
    assert "not accuracy" in seg or "no spot-level ground truth" in seg


@pytest.mark.skipif(not REPORT.exists(), reason="report not generated")
def test_qc_before_predictions_in_main():
    import re
    html = REPORT.read_text()
    pos = {m.group(1): m.start() for m in re.finditer(r"id='([a-z_]+)'", html)}
    # reference+signature (merged → 'reference'), input, and trusted resolution
    # ('resolution') are all QC and must precede the prediction sections.
    for qc in ("reference", "input", "resolution"):
        assert pos[qc] < pos["bulk"] and pos[qc] < pos["spatial"]
    # predictions precede benchmarks
    assert pos["bulk"] < pos["bulk_benchmark"]
    assert pos["spatial"] < pos["spatial_benchmark"]
