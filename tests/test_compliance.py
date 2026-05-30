"""
Cross-cutting compliance tests enforcing the non-negotiable scientific rules
from CLAUDE.md / DESIGN_SPEC.md.

These guard the *contract*, independent of any single module:

1.  No silent gene removal (subset emits a warning).
2.  Bulk outputs are mRNA/RNA-derived (estimate type recorded everywhere).
3.  Spatial outputs record smoothing parameters.
4.  Plots save underlying source data.
5.  Reports surface warnings.
6.  Reports include methods text.
7.  Captions mention the estimate type.
8.  Separability problems are surfaced.
9.  The default suite runs offline (no network).
"""
from __future__ import annotations

import socket
from types import SimpleNamespace

import numpy as np
import pytest

from tissueresolve.config import TissueResolveConfig
from tissueresolve.plotting import captions
from tissueresolve.report import html, methods_text
from tissueresolve.results import (
    BulkDeconvResult,
    SpatialDeconvResult,
)


# ---------------------------------------------------------------------------
# 1. No silent gene removal
# ---------------------------------------------------------------------------


def test_subset_genes_warns_on_missing(ref_sig_full):
    with pytest.warns(UserWarning, match="not in reference"):
        ref_sig_full.subset_genes(list(ref_sig_full.gene_names[:3]) + ["FAKE_GENE_X"])


# ---------------------------------------------------------------------------
# 2. Bulk estimate type is mRNA proportion, recorded in saved output
# ---------------------------------------------------------------------------


def test_bulk_estimate_type_constant(bulk_result):
    assert bulk_result.ESTIMATE_TYPE == "mRNA_proportion"


def test_bulk_saved_output_records_estimate_type(bulk_result, tmp_path):
    bulk_result.save(tmp_path / "deconv")
    header = (tmp_path / "deconv" / "proportions.tsv").read_text()
    assert "estimate_type: mRNA_proportion" in header
    assert "NOT cell fractions" in header


def test_spatial_estimate_type_constant(spatial_result):
    assert spatial_result.ESTIMATE_TYPE == "spot_rna_composition"


# ---------------------------------------------------------------------------
# 3. Spatial smoothing parameters are always recorded
# ---------------------------------------------------------------------------


def test_spatial_saved_output_records_smoothing(spatial_result, tmp_path):
    spatial_result.save(tmp_path / "deconv")
    header = (tmp_path / "deconv" / "proportions.tsv").read_text()
    assert "lambda_spatial" in header
    import json

    meta = json.loads((tmp_path / "deconv" / "metadata.json").read_text())
    assert "lambda_spatial" in meta
    assert "n_smooth" in meta


def test_spatial_pipeline_records_lambda_end_to_end():
    """A real (tiny) spatial run must record λ in the result + metadata."""
    from tissueresolve.spatial.benchmark import simulate_visium
    import tissueresolve as tr

    ds = simulate_visium(n_spots=40, n_types=3, n_genes=30, seed=0)
    cfg = TissueResolveConfig()
    cfg.spatial_solver.max_iter = 15
    result = tr.deconv_spatial(
        ds.dense_counts(), ds.to_reference(), ds.array_row, ds.array_col,
        ds.lib_sizes, list(ds.gene_names),
        marker_genes=list(ds.gene_names), config=cfg,
    )
    assert result.deconv.lambda_spatial == cfg.spatial_solver.lambda_spatial
    assert "lambda_spatial" in result.run_metadata


# ---------------------------------------------------------------------------
# 4. Plots save underlying source data
# ---------------------------------------------------------------------------


def test_plots_save_source_data(bulk_result, tmp_path):
    pytest.importorskip("matplotlib")
    from tissueresolve.plotting import bulk_plots

    pr = bulk_plots.composition_barplot(bulk_result, tmp_path)
    assert pr.data_paths, "plot produced no source data"
    for p in pr.data_paths.values():
        assert p.exists()


# ---------------------------------------------------------------------------
# 5 & 6. Reports surface warnings and include methods text
# ---------------------------------------------------------------------------


def test_report_surfaces_warnings_and_methods(bulk_result, bulk_qc_report, tmp_path):
    res = SimpleNamespace(
        deconv=bulk_result, qc=bulk_qc_report, protocol_risk=None,
        run_metadata={"genome": "hg38"},
    )
    doc = html.generate_bulk_report(res, tmp_path / "r.html").read_text()
    assert "Warnings" in doc
    assert "Methods" in doc
    # estimate-type caveat present
    assert "mRNA" in doc and "cell fractions" in doc


def test_spatial_report_shows_failed_convergence(spatial_result, spatial_qc_report, tmp_path):
    import dataclasses

    nc = dataclasses.replace(spatial_result, converged=False)
    res = SimpleNamespace(
        deconv=nc, qc=spatial_qc_report, morans_i=spatial_qc_report.morans_i,
        spot_qc=spatial_qc_report.spot_qc, run_metadata={"alpha": 0.1},
    )
    doc = html.generate_spatial_report(res, tmp_path / "s.html").read_text()
    assert "did NOT converge" in doc  # failed check not hidden


# ---------------------------------------------------------------------------
# 7. Captions mention the estimate type
# ---------------------------------------------------------------------------


def test_bulk_caption_mentions_mrna():
    cap = captions.bulk_composition_caption(3, 4)
    assert "mRNA proportions" in cap
    assert "NOT absolute cell fractions" in cap


def test_spatial_caption_mentions_composition():
    cap = captions.spatial_abundance_caption("Tcell", lambda_spatial=0.1, n_spots=100)
    assert "composition estimates" in cap
    assert "NOT" in cap and "single-cell counts" in cap

    dom = captions.dominant_type_caption(lambda_spatial=0.1, n_spots=100)
    assert "single-cell counts" in dom


# ---------------------------------------------------------------------------
# 8. Separability problems are surfaced
# ---------------------------------------------------------------------------


def test_separability_warning_surfaced(separability_report):
    assert separability_report.has_problems
    cap = captions.separability_heatmap_caption(
        n_critical=separability_report.n_critical,
        n_high=separability_report.n_high,
    )
    assert "WARNING" in cap
    # And in the methods text
    txt = methods_text.separability_methods(
        n_critical=separability_report.n_critical,
        n_high=separability_report.n_high,
    )
    assert "poorly-separable" in txt


# ---------------------------------------------------------------------------
# 9. Default suite runs offline
# ---------------------------------------------------------------------------


def test_core_pipeline_runs_without_network(monkeypatch):
    """Block outbound sockets, then run a tiny spatial deconvolution."""
    def _blocked(*args, **kwargs):
        raise OSError("network access is blocked in the default test suite")

    monkeypatch.setattr(socket.socket, "connect", _blocked)

    from tissueresolve.spatial.benchmark import simulate_visium
    import tissueresolve as tr

    ds = simulate_visium(n_spots=30, n_types=3, n_genes=24, seed=1)
    cfg = TissueResolveConfig()
    cfg.spatial_solver.max_iter = 10
    result = tr.deconv_spatial(
        ds.dense_counts(), ds.to_reference(), ds.array_row, ds.array_col,
        ds.lib_sizes, list(ds.gene_names),
        marker_genes=list(ds.gene_names), config=cfg,
    )
    assert result.deconv.proportions.shape[0] == ds.n_spots
