"""
Tests for the CLI report commands (offline).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from click.testing import CliRunner

from tissueresolve.cli import cli


def test_report_help():
    r = CliRunner().invoke(cli, ["report", "--help"])
    assert r.exit_code == 0
    assert "--modality" in r.output and "--results-dir" in r.output


def test_bulk_report_help():
    r = CliRunner().invoke(cli, ["bulk", "report", "--help"])
    assert r.exit_code == 0
    assert "--results-dir" in r.output


def test_spatial_report_help():
    r = CliRunner().invoke(cli, ["spatial", "report", "--help"])
    assert r.exit_code == 0
    assert "--results-dir" in r.output


def test_bulk_report_runs(tmp_path):
    rdir = tmp_path / "bulk"
    rdir.mkdir()
    pd.DataFrame(np.eye(2), index=["s0", "s1"], columns=["A", "B"]).to_csv(
        rdir / "bulk_estimated_proportions.tsv", sep="\t")
    out = rdir / "report.html"
    r = CliRunner().invoke(cli, ["bulk", "report", "--results-dir", str(rdir),
                                 "--out", str(out)])
    assert r.exit_code == 0, r.output
    assert out.exists()
    assert "Wrote bulk report" in r.output


def test_top_level_report_runs(tmp_path):
    rdir = tmp_path / "spatial"
    rdir.mkdir()
    pd.DataFrame(np.eye(2), index=["sp0", "sp1"], columns=["A", "B"]).to_csv(
        rdir / "spatial_spot_proportions.tsv", sep="\t")
    r = CliRunner().invoke(cli, ["report", "--modality", "spatial",
                                 "--results-dir", str(rdir)])
    assert r.exit_code == 0, r.output
    assert (rdir / "report.html").exists()
