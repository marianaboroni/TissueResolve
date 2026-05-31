"""
Tests for the redesigned report layout (executive summary, key findings,
interpretation, severity warnings, collapsible details).  Offline.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from tissueresolve.report import generate_report


def _bulk_results(tmp_path, *, with_separability=True):
    rdir = tmp_path / "bulk"
    tables = rdir / "tables"
    figs = rdir / "figures"
    tables.mkdir(parents=True)
    figs.mkdir(parents=True)
    pd.DataFrame(np.eye(3), index=["s0", "s1", "s2"],
                 columns=["CD8 T cell", "macrophage", "fibroblast"]).to_csv(
        rdir / "bulk_estimated_proportions.tsv", sep="\t")
    pd.DataFrame({"value": [40000, 5000, 3]},
                 index=["n_cells", "n_genes", "n_cell_types"]).to_csv(
        tables / "reference_summary.tsv", sep="\t")
    if with_separability:
        pd.DataFrame({
            "type_a": ["CD8 T cell"], "type_b": ["macrophage"],
            "separability_score": [0.04], "resolvability": ["unresolved"],
        }).to_csv(tables / "pairwise_resolvability.tsv", sep="\t", index=False)
    (figs / "bulk_main_summary_figure.html").write_text("<html></html>")
    return rdir


def _titles(doc: str, headings):
    return {h: doc.find(h) for h in headings}


def test_executive_summary_and_key_findings(tmp_path):
    doc = generate_report("bulk", _bulk_results(tmp_path)).read_text()
    assert "Executive summary" in doc
    assert "Key findings" in doc
    # summary cards present
    assert "class='cards'" in doc or "class=\"cards\"" in doc


def test_main_figure_before_detailed_tables(tmp_path):
    doc = generate_report("bulk", _bulk_results(tmp_path)).read_text()
    pos = _titles(doc, ["Main publication figure", "Detailed outputs"])
    assert pos["Main publication figure"] != -1
    assert pos["Detailed outputs"] != -1
    assert pos["Main publication figure"] < pos["Detailed outputs"]
    # main summary figure embedded
    assert "bulk_main_summary_figure.html" in doc


def test_interpretation_paragraphs_present(tmp_path):
    doc = generate_report("bulk", _bulk_results(tmp_path)).read_text()
    assert "Main results interpretation" in doc
    assert "mRNA-derived" in doc  # estimate-type-aware interpretation


def test_not_no_warnings_when_separability_exists(tmp_path):
    doc = generate_report("bulk", _bulk_results(tmp_path, with_separability=True)).read_text()
    assert "No warnings." not in doc
    assert "unresolved" in doc.lower()
    # severity chips rendered
    assert "sev-" in doc


def test_detailed_outputs_collapsible(tmp_path):
    doc = generate_report("bulk", _bulk_results(tmp_path)).read_text()
    assert "<details" in doc
    assert "Full prediction table" in doc


def test_missing_bootstrap_shows_card_not_plot(tmp_path):
    doc = generate_report("bulk", _bulk_results(tmp_path)).read_text()
    assert "Bootstrap uncertainty was not computed" in doc
    # no uncertainty figure iframe (no bootstrap figure was generated)
    assert "bulk_uncertainty_plot.html" not in doc


def test_spatial_report_redesign(tmp_path):
    rdir = tmp_path / "spatial"
    (rdir / "figures").mkdir(parents=True)
    pd.DataFrame(np.eye(3), index=["sp0", "sp1", "sp2"],
                 columns=["A", "B", "C"]).to_csv(
        rdir / "spatial_spot_proportions.tsv", sep="\t")
    pd.Series([0.8, 0.3, 0.1], index=["A", "B", "C"], name="morans_i").to_frame(
        ).to_csv(rdir / "morans_i.tsv", sep="\t")
    doc = generate_report("spatial", rdir,
                          run_metadata={"lambda_spatial": 0.1, "converged": False,
                                        "has_he_image": False}).read_text()
    assert "Executive summary" in doc
    assert "cell counts" in doc  # spatial estimate-type disclaimer
    assert "No warnings." not in doc  # missing H&E + no bootstrap → warnings exist
