"""
Tests for the output-modality guard on `tissueresolve run` and the explicit
combined bulk+spatial report command.  Offline + fast (the guard tests use
--dry-run; the combine tests fabricate run-shaped directories).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tissueresolve import cli


# ---------------------------------------------------------------------------
# Output modality guard
# ---------------------------------------------------------------------------


def _tiny_ref(tmp_path):
    import anndata as ad
    ref = tmp_path / "ref.h5ad"
    adata = ad.AnnData(X=np.ones((3, 2)))
    adata.obs["cell_type"] = ["A"] * 3
    adata.write(ref)
    return ref


def _tiny_query(tmp_path):
    q = tmp_path / "counts.tsv"
    q.write_text("gene\ts1\nG1\t1\n")
    return q


def _seed_prior_run(out: Path, mode: str):
    out.mkdir(parents=True, exist_ok=True)
    (out / "analysis_plan.json").write_text(json.dumps({"mode": mode}))


def _run(tmp_path, out, mode, *, force=False):
    args = ["--reference", str(_tiny_ref(tmp_path)), "--query", str(_tiny_query(tmp_path)),
            "--out", str(out), "--mode", mode, "--dry-run"]
    if force:
        args.append("--force")
    return cli.run(args)


def test_guard_blocks_spatial_over_bulk(tmp_path):
    out = tmp_path / "out"
    _seed_prior_run(out, "bulk")
    with pytest.raises(Exception) as e:
        _run(tmp_path, out, "spatial")
    msg = str(e.value)
    assert "single modality" in msg.lower()
    assert "--force" in msg
    assert "results/bulk" in msg and "results/spatial" in msg


def test_guard_blocks_bulk_over_spatial(tmp_path):
    out = tmp_path / "out"
    _seed_prior_run(out, "spatial")
    with pytest.raises(Exception) as e:
        _run(tmp_path, out, "bulk")
    assert "spatial" in str(e.value)


def test_guard_allows_same_modality(tmp_path):
    out = tmp_path / "out"
    _seed_prior_run(out, "bulk")
    assert _run(tmp_path, out, "bulk") == 0  # same modality refreshes, no error


def test_guard_allows_empty_dir(tmp_path):
    assert _run(tmp_path, tmp_path / "fresh", "bulk") == 0


def test_force_allows_cross_modality_overwrite(tmp_path):
    out = tmp_path / "out"
    _seed_prior_run(out, "bulk")
    assert _run(tmp_path, out, "spatial", force=True) == 0


def test_run_help_documents_force():
    from click.testing import CliRunner
    r = CliRunner().invoke(cli.cli, ["run", "--help"])
    assert "--force" in r.output


# ---------------------------------------------------------------------------
# Combined report
# ---------------------------------------------------------------------------


def _fake_bulk_run(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    (d / "run_metadata.json").write_text(json.dumps({"analysis_plan": {
        "mode": "bulk", "preset": "quick", "resolution_mode": "none",
        "detected_reference": "single_cell_h5ad_reference"}}))
    (d / "analysis_plan.json").write_text(json.dumps({"mode": "bulk"}))
    (d / "deconv").mkdir(exist_ok=True)
    pd.DataFrame(np.eye(2), index=["s0", "s1"], columns=["A", "B"]).to_csv(
        d / "deconv" / "proportions.tsv", sep="\t")
    pd.DataFrame({"coverage_r2": [0.9, 0.8]}, index=["s0", "s1"]).to_csv(
        d / "deconv" / "coverage_r2.tsv", sep="\t")
    (d / "qc").mkdir(exist_ok=True)
    (d / "qc" / "recommendations.txt").write_text("bulk sample s1 low R2\n")
    figs = d / "figures"; figs.mkdir(exist_ok=True)
    (figs / "bulk_composition.png").write_bytes(b"\x89PNG\r\n\x1a\n")  # stub png
    (figs / "bulk_composition.data.tsv").write_text("a\tb\n1\t2\n")
    (d / "methods.txt").write_text("Bulk methods text.")
    (d / "warnings.json").write_text(json.dumps([
        {"severity": "info", "category": "estimate_type",
         "message": "Estimates are RNA-derived mRNA proportions, not absolute cell fractions."},
        {"severity": "warning", "category": "qc", "message": "BULK_WARN_X"}]))
    return d


def _fake_spatial_run(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    (d / "run_metadata.json").write_text(json.dumps({"analysis_plan": {
        "mode": "spatial", "preset": "quick", "resolution_mode": "none",
        "detected_reference": "single_cell_h5ad_reference"}}))
    (d / "analysis_plan.json").write_text(json.dumps({"mode": "spatial"}))
    (d / "deconv").mkdir(exist_ok=True)
    pd.DataFrame(np.eye(2), index=["sp0", "sp1"], columns=["A", "B"]).to_csv(
        d / "deconv" / "proportions.tsv", sep="\t")
    (d / "qc").mkdir(exist_ok=True)
    pd.Series([0.8, 0.3], index=["A", "B"], name="morans_i").to_frame().to_csv(
        d / "qc" / "morans_i.tsv", sep="\t")
    (d / "qc" / "recommendations.txt").write_text("spatial low-quality spots\n")
    (d / "methods.txt").write_text("Spatial methods text.")
    (d / "warnings.json").write_text(json.dumps([
        {"severity": "info", "category": "estimate_type",
         "message": "Estimates are spot-level RNA-derived composition, not single-cell counts."},
        {"severity": "warning", "category": "qc", "message": "SPATIAL_WARN_Y"}]))
    return d


def _combined(tmp_path):
    from tissueresolve.report.combined import generate_combined_report
    b = _fake_bulk_run(tmp_path / "bulk")
    s = _fake_spatial_run(tmp_path / "spatial")
    out = tmp_path / "combined"
    return generate_combined_report(b, s, out), out


def test_combine_writes_bundle(tmp_path):
    rpt, out = _combined(tmp_path)
    assert rpt.exists() and rpt.name == "report.html"
    for f in ("report.html", "methods.txt", "warnings.json", "run_metadata.json"):
        assert (out / f).exists(), f


def test_combine_embeds_run_figures(tmp_path):
    rpt, _ = _combined(tmp_path)
    doc = rpt.read_text()
    # bulk figure from the run dir is embedded (figure-driven, not table-only)
    assert "<img" in doc
    assert "bulk_composition.png" in doc


def test_combine_has_bulk_and_spatial_sections(tmp_path):
    rpt, _ = _combined(tmp_path)
    doc = rpt.read_text()
    assert "Bulk QC and predictions" in doc
    assert "Spatial QC and predictions" in doc
    # estimate-type caveats for both modalities present
    assert "mRNA proportions" in doc and "single-cell counts" in doc


def test_combine_keeps_benchmarks_separate(tmp_path):
    rpt, _ = _combined(tmp_path)
    doc = rpt.read_text()
    assert "Bulk benchmark summary" in doc
    assert "Spatial benchmark summary" in doc
    # spatial benchmark must not claim accuracy without ground truth
    assert "not accuracy" in doc or "no ground truth" in doc.lower()


def test_combine_warnings_include_both(tmp_path):
    _, out = _combined(tmp_path)
    warns = json.loads((out / "warnings.json").read_text())
    modalities = {w.get("modality") for w in warns}
    assert {"bulk", "spatial"} <= modalities
    messages = " ".join(w.get("message", "") for w in warns)
    assert "BULK_WARN_X" in messages and "SPATIAL_WARN_Y" in messages


def test_combine_states_not_joint_model(tmp_path):
    rpt, out = _combined(tmp_path)
    doc = rpt.read_text().lower()
    assert "not" in doc and "joint" in doc  # "not a single joint ... model"
    md = (out / "run_metadata.json").read_text().lower()
    assert "not a joint" in md or "not a single joint" in md


def test_combine_validates_modality(tmp_path):
    """Passing a spatial run as --bulk-dir is rejected."""
    from tissueresolve.report.combined import generate_combined_report
    b = _fake_bulk_run(tmp_path / "bulk")
    s = _fake_spatial_run(tmp_path / "spatial")
    with pytest.raises(ValueError):
        generate_combined_report(s, b, tmp_path / "bad")  # swapped


def test_combine_cli(tmp_path):
    from click.testing import CliRunner
    b = _fake_bulk_run(tmp_path / "bulk")
    s = _fake_spatial_run(tmp_path / "spatial")
    out = tmp_path / "combined"
    r = CliRunner().invoke(cli.cli, ["combine-report", "--bulk-dir", str(b),
                                     "--spatial-dir", str(s), "--out", str(out)])
    assert r.exit_code == 0, r.output
    assert (out / "report.html").exists()
