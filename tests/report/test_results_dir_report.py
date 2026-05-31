"""
Tests for the results-directory-driven HTML report (publication layer).

Offline.  Builds a tiny results directory and checks the generated report
includes the required sections, estimate-type disclaimers, warnings and
methods — and that failed checks are not hidden.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from tissueresolve.report import generate_report
from tissueresolve.report.html import generate_bulk_report, generate_spatial_report


def _make_bulk_results(tmp_path):
    rdir = tmp_path / "bulk"
    tables = rdir / "tables"
    figs = rdir / "figures"
    tables.mkdir(parents=True)
    figs.mkdir(parents=True)
    props = pd.DataFrame(np.eye(3), index=["s0", "s1", "s2"],
                         columns=["Tcell", "Bcell", "Myeloid"])
    props.to_csv(rdir / "bulk_estimated_proportions.tsv", sep="\t")
    pd.DataFrame({"recon_r2": [0.9, 0.8, 0.7]},
                 index=["s0", "s1", "s2"]).to_csv(rdir / "bulk_qc.tsv", sep="\t")
    pd.DataFrame({"value": [40, 8, 3]},
                 index=["n_cells", "n_genes", "n_cell_types"]).to_csv(
        tables / "reference_summary.tsv", sep="\t")
    (tables / "warnings.json").write_text(json.dumps(
        {"qc_recommendations": ["sample s2 has low R²"]}))
    # a placeholder figure so the figure block links it
    (figs / "bulk_composition_clustered_barplot.html").write_text("<html></html>")
    return rdir


def test_bulk_results_dir_report_sections(tmp_path):
    rdir = _make_bulk_results(tmp_path)
    out = generate_report("bulk", rdir, rdir / "report.html")
    assert out.exists()
    doc = out.read_text()
    for heading in ("Single-cell reference quality", "Input data summary",
                    "Deconvolution predictions", "Prediction QC",
                    "Warnings and limitations", "Methods"):
        assert heading in doc, heading
    # estimate-type disclaimer present (bulk)
    assert "not" in doc and "cell fractions" in doc
    # surfaced warning
    assert "low R" in doc
    # figure embedded by iframe
    assert "bulk_composition_clustered_barplot.html" in doc


def test_bulk_report_dispatch_via_path(tmp_path):
    rdir = _make_bulk_results(tmp_path)
    out = generate_bulk_report(str(rdir), out=rdir / "r.html")
    assert out.exists()
    assert "mRNA-derived" in out.read_text()


def test_spatial_results_dir_report(tmp_path):
    rdir = tmp_path / "spatial"
    (rdir / "figures").mkdir(parents=True)
    props = pd.DataFrame(np.eye(3), index=["sp0", "sp1", "sp2"],
                         columns=["A", "B", "C"])
    props.to_csv(rdir / "spatial_spot_proportions.tsv", sep="\t")
    pd.Series([0.8, 0.3, 0.1], index=["A", "B", "C"], name="morans_i").to_frame(
        ).to_csv(rdir / "morans_i.tsv", sep="\t")
    out = generate_spatial_report(
        str(rdir), out=rdir / "report.html",
        run_metadata={"lambda_spatial": 0.1, "converged": False})
    doc = out.read_text()
    for heading in ("Single-cell reference quality", "Spatial predictions",
                    "Spatial structure", "Methods"):
        assert heading in doc, heading
    # spatial estimate-type disclaimer + recorded smoothing + failed convergence shown
    assert "cell counts" in doc
    assert "lambda_spatial" in doc
    assert "False" in doc  # converged=False not hidden


def test_report_offline(monkeypatch, tmp_path):
    import socket

    def _blocked(*a, **k):
        raise OSError("network blocked")

    monkeypatch.setattr(socket.socket, "connect", _blocked)
    rdir = _make_bulk_results(tmp_path)
    out = generate_report("bulk", rdir)
    assert out.exists()


def test_generate_report_bad_modality(tmp_path):
    with pytest.raises(ValueError):
        generate_report("nonsense", str(tmp_path))


def test_report_surfaces_recommended_merges_and_families(tmp_path):
    rdir = _make_bulk_results(tmp_path)
    tables = rdir / "tables"
    pd.DataFrame({
        "family_name": ["T/NK lymphocytes"], "members": ["CD4 T cell; CD8 T cell"],
        "n_members": [2], "resolvability": ["unresolved"],
        "mean_separability": [0.05], "recommended_merge": [True],
    }).to_csv(tables / "recommended_merges.tsv", sep="\t", index=False)
    pd.DataFrame(np.eye(2), index=["s0", "s1"],
                 columns=["T/NK lymphocytes", "Myeloid"]).to_csv(
        rdir / "bulk_family_proportions.tsv", sep="\t")
    doc = generate_report("bulk", rdir, rdir / "report.html").read_text()
    assert "Recommended merge families" in doc
    assert "Family-level estimates" in doc
    assert "T/NK lymphocytes" in doc
    # banner warning about confusable types
    assert "recommended merge family" in doc.lower()
