"""
Tests for tissueresolve.report (methods_text + html) — Stage 5.

Reports must surface warnings, include methods text, and never hide failed
checks.  Methods text must distinguish RNA-derived proportions / cell fractions
/ spot composition / single-cell counts and must not overclaim.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from tissueresolve.report import html, methods_text


# ---------------------------------------------------------------------------
# methods_text
# ---------------------------------------------------------------------------


class TestMethodsText:
    def test_bulk_estimate_statement(self):
        s = methods_text.estimate_type_statement("bulk")
        assert "mRNA proportions" in s
        assert "not" in s.lower() and "cell fractions" in s

    def test_spatial_estimate_statement(self):
        s = methods_text.estimate_type_statement("spatial")
        assert "composition" in s
        assert "single-cell counts" in s

    def test_spatial_deconv_records_lambda(self):
        s = methods_text.spatial_deconvolution_methods(
            lambda_spatial=0.25, n_marker_genes=100, converged=True, n_iter=42
        )
        assert "0.25" in s and "recorded" in s

    def test_non_convergence_flagged(self):
        s = methods_text.spatial_deconvolution_methods(
            lambda_spatial=0.1, n_marker_genes=50, converged=False, n_iter=200
        )
        assert "did NOT converge" in s

    def test_protocol_risk_mentions_no_silent_removal(self):
        s = methods_text.protocol_risk_methods(
            risk_level="high", n_excluded=12, mismatch_types=["polyA_vs_snRNA"]
        )
        assert "silently" in s

    def test_compose_bulk(self, bulk_result, bulk_qc_report):
        res = SimpleNamespace(
            deconv=bulk_result, qc=bulk_qc_report, protocol_risk=None,
            run_metadata={"genome": "hg38"},
        )
        text = methods_text.compose_bulk_methods(res)
        assert "mRNA proportions" in text
        assert "wNNLS" in text or "non-negative least squares" in text

    def test_compose_spatial(self, spatial_result, spatial_qc_report):
        res = SimpleNamespace(
            deconv=spatial_result, morans_i=spatial_qc_report.morans_i,
            run_metadata={"genome": "hg38"},
        )
        text = methods_text.compose_spatial_methods(res)
        assert "composition estimates" in text
        assert "λ_spatial" in text or "lambda" in text.lower()


# ---------------------------------------------------------------------------
# HTML reports
# ---------------------------------------------------------------------------


def _bulk_pipeline_result(deconv, qc):
    return SimpleNamespace(
        deconv=deconv, qc=qc, protocol_risk=None,
        run_metadata={"genome": "hg38", "n_reference_genes": 100},
    )


def _spatial_pipeline_result(deconv, qc):
    return SimpleNamespace(
        deconv=deconv, qc=qc, morans_i=qc.morans_i, spot_qc=qc.spot_qc,
        run_metadata={"genome": "hg38", "alpha": 0.1, "lambda_spatial": deconv.lambda_spatial},
    )


class TestBulkReport:
    def test_writes_html(self, bulk_result, bulk_qc_report, tmp_path):
        res = _bulk_pipeline_result(bulk_result, bulk_qc_report)
        out = html.generate_bulk_report(res, tmp_path / "bulk.html")
        assert out.exists()
        doc = out.read_text()
        assert "<html" in doc.lower()

    def test_surfaces_estimate_warning(self, bulk_result, bulk_qc_report, tmp_path):
        res = _bulk_pipeline_result(bulk_result, bulk_qc_report)
        doc = html.generate_bulk_report(res, tmp_path / "bulk.html").read_text()
        assert "Warnings" in doc
        assert "mRNA" in doc and "cell fractions" in doc

    def test_includes_methods(self, bulk_result, bulk_qc_report, tmp_path):
        res = _bulk_pipeline_result(bulk_result, bulk_qc_report)
        doc = html.generate_bulk_report(res, tmp_path / "bulk.html").read_text()
        assert "Methods" in doc

    def test_surfaces_qc_recommendation(self, bulk_result, bulk_qc_report, tmp_path):
        # bulk_qc_report fixture carries a recommendation + a 'high' mismatch flag.
        res = _bulk_pipeline_result(bulk_result, bulk_qc_report)
        doc = html.generate_bulk_report(res, tmp_path / "bulk.html").read_text()
        assert "high" in doc.lower()


class TestSpatialReport:
    def test_writes_html(self, spatial_result, spatial_qc_report, tmp_path):
        res = _spatial_pipeline_result(spatial_result, spatial_qc_report)
        out = html.generate_spatial_report(res, tmp_path / "spatial.html")
        assert out.exists()

    def test_records_smoothing(self, spatial_result, spatial_qc_report, tmp_path):
        res = _spatial_pipeline_result(spatial_result, spatial_qc_report)
        doc = html.generate_spatial_report(res, tmp_path / "spatial.html").read_text()
        assert "lambda_spatial" in doc
        assert "Smoothing parameters" in doc

    def test_shows_convergence(self, spatial_result, spatial_qc_report, tmp_path):
        res = _spatial_pipeline_result(spatial_result, spatial_qc_report)
        doc = html.generate_spatial_report(res, tmp_path / "spatial.html").read_text()
        assert "Convergence" in doc

    def test_non_convergence_shown(self, spatial_result, spatial_qc_report, tmp_path):
        import dataclasses

        nc = dataclasses.replace(spatial_result, converged=False)
        res = _spatial_pipeline_result(nc, spatial_qc_report)
        doc = html.generate_spatial_report(res, tmp_path / "spatial.html").read_text()
        assert "did NOT converge" in doc

    def test_surfaces_separability_warning(self, spatial_result, spatial_qc_report,
                                           separability_report, tmp_path):
        res = _spatial_pipeline_result(spatial_result, spatial_qc_report)
        doc = html.generate_spatial_report(
            res, tmp_path / "spatial.html", separability=separability_report
        ).read_text()
        assert "poorly-separable" in doc or "CRITICAL" in doc
