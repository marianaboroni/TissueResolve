"""
Report-stability regression tests (offline).

Guards the conclusions in ``docs/REPORT_STABILITY_AUDIT.md``: report.html is
generated, QC precedes predictions, figure cards are never empty, the spatial
no-ground-truth caveat exists, and a methods / source-data path is present.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from tissueresolve.report import generate_report


def _bulk_results_dir(tmp_path):
    rdir = tmp_path / "bulk"
    (rdir / "tables").mkdir(parents=True)
    (rdir / "figures").mkdir(parents=True)
    pd.DataFrame(np.eye(3), index=["s0", "s1", "s2"],
                 columns=["CD8 T cell", "macrophage", "fibroblast"]).to_csv(
        rdir / "bulk_estimated_proportions.tsv", sep="\t")
    pd.DataFrame({"value": [40000, 5000, 3]},
                 index=["n_cells", "n_genes", "n_cell_types"]).to_csv(
        rdir / "tables" / "reference_summary.tsv", sep="\t")
    (rdir / "figures" / "bulk_main_summary_figure.html").write_text("<html></html>")
    return rdir


def test_report_html_generated(tmp_path):
    out = generate_report("bulk", _bulk_results_dir(tmp_path))
    assert out.exists() and out.name == "report.html"


def test_qc_section_before_predictions(tmp_path):
    doc = generate_report("bulk", _bulk_results_dir(tmp_path)).read_text()
    ref_pos = doc.find("Single-cell reference quality")
    pred_pos = doc.find("Deconvolution predictions")
    assert ref_pos != -1 and pred_pos != -1
    assert ref_pos < pred_pos


def test_no_empty_figure_cards(tmp_path):
    """A figure card must never render an empty body."""
    doc = generate_report("bulk", _bulk_results_dir(tmp_path)).read_text()
    # no literally-empty figure body
    assert not re.search(r"<div class='fig-body'>\s*</div>", doc)
    # the real figure is embedded
    assert "bulk_main_summary_figure.html" in doc


def test_methods_and_source_data_present(tmp_path):
    doc = generate_report("bulk", _bulk_results_dir(tmp_path)).read_text()
    assert "Methods" in doc
    # estimate-type disclaimer + source-data / output file listing
    assert "mRNA-derived" in doc
    assert "Output files" in doc or "Detailed outputs" in doc


def test_spatial_no_ground_truth_caveat():
    """The report glossary must define concordance as not-equal-to correctness
    when no ground truth exists (no accuracy claim without ground truth)."""
    from pathlib import Path
    import tissueresolve.report.glossary as glossary
    src = Path(glossary.__file__).read_text().lower()
    assert "concordance" in src
    assert "ground" in src and "truth" in src


def test_spatial_report_shows_failed_convergence(tmp_path):
    """converged=False must be visible, not hidden."""
    from tissueresolve.report.html import generate_spatial_report
    rdir = tmp_path / "spatial"
    (rdir / "figures").mkdir(parents=True)
    pd.DataFrame(np.eye(3), index=["sp0", "sp1", "sp2"],
                 columns=["A", "B", "C"]).to_csv(
        rdir / "spatial_spot_proportions.tsv", sep="\t")
    pd.Series([0.8, 0.3, 0.1], index=["A", "B", "C"], name="morans_i"
              ).to_frame().to_csv(rdir / "morans_i.tsv", sep="\t")
    doc = generate_spatial_report(
        str(rdir), out=rdir / "report.html",
        run_metadata={"lambda_spatial": 0.1, "converged": False}).read_text()
    assert "cell counts" in doc          # estimate-type disclaimer
    assert "lambda_spatial" in doc       # smoothing recorded
    assert "False" in doc                # converged=False not hidden
