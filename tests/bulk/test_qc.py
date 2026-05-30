"""
Tests for tissueresolve.bulk.qc (BulkQC).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.bulk.qc import BulkQC
from tissueresolve.bulk.solver import WNNLSSolver
from tissueresolve.config import BulkQCConfig
from tissueresolve.results import BulkDeconvResult, QCReport, ReferenceSignature


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ref(n_genes: int = 80, n_types: int = 4, seed: int = 3) -> ReferenceSignature:
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    phi_raw = np.zeros((n_genes, n_types))
    for k in range(n_types):
        phi_raw[k * block:(k + 1) * block, k] = rng.uniform(2.0, 5.0, block)
    phi_raw += rng.uniform(0.0, 0.05, phi_raw.shape)
    col_sums = phi_raw.sum(axis=0, keepdims=True)
    phi = (phi_raw / col_sums).astype(np.float64)
    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_types)]
    return ReferenceSignature(gene_names=gene_names, cell_types=cell_types, phi=phi)


def _make_bulk_and_result(
    ref: ReferenceSignature,
    n_samples: int = 8,
    seed: int = 17,
) -> tuple[pd.DataFrame, BulkDeconvResult]:
    rng = np.random.default_rng(seed)
    K = ref.n_cell_types
    props = rng.dirichlet(np.ones(K), size=n_samples)
    counts = (ref.as_phi() @ props.T) * 1e4
    bulk = pd.DataFrame(counts, index=ref.gene_names,
                        columns=[f"S{i}" for i in range(n_samples)])
    solver = WNNLSSolver()
    result = solver.solve(bulk, ref, ref.gene_names)
    return bulk, result


# ---------------------------------------------------------------------------
# Basic output structure
# ---------------------------------------------------------------------------


class TestBulkQCOutput:
    def test_returns_qc_report(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        qc = BulkQC()
        report = qc.compute(result, ref, bulk)
        assert isinstance(report, QCReport)

    def test_modality_is_bulk(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.modality == "bulk"

    def test_recon_r2_present(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.recon_r2 is not None
        assert len(report.recon_r2) == bulk.shape[1]

    def test_profile_corr_present(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.profile_corr is not None
        assert len(report.profile_corr) == bulk.shape[1]

    def test_mismatch_flag_present(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.mismatch_flag is not None
        assert set(report.mismatch_flag.unique()).issubset({"low", "medium", "high"})

    def test_marker_recall_present(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.marker_recall is not None
        assert len(report.marker_recall) == ref.n_cell_types

    def test_spillover_risk_present(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.spillover_risk is not None
        assert len(report.spillover_risk) == ref.n_cell_types

    def test_condition_number_positive(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.condition_number is not None
        assert report.condition_number > 0

    def test_recommendations_is_list(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert isinstance(report.recommendations, list)


# ---------------------------------------------------------------------------
# Metric value ranges
# ---------------------------------------------------------------------------


class TestBulkQCMetricRanges:
    def test_recon_r2_in_minus_one_to_one(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert (report.recon_r2 <= 1.0 + 1e-9).all()

    def test_profile_corr_in_minus_one_to_one(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        finite = report.profile_corr.dropna()
        assert (finite >= -1.0 - 1e-9).all()
        assert (finite <= 1.0 + 1e-9).all()

    def test_marker_recall_in_zero_one(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert (report.marker_recall >= 0.0).all()
        assert (report.marker_recall <= 1.0 + 1e-9).all()

    def test_spillover_risk_bounded(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert (report.spillover_risk >= -1.0 - 1e-9).all()
        assert (report.spillover_risk <= 1.0 + 1e-9).all()

    def test_high_r2_for_noiseless_data(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.recon_r2.mean() > 0.85, "Expected high R² for noiseless bulk data"

    def test_metadata_is_heuristic_flagged(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.metadata.get("is_heuristic") is True


# ---------------------------------------------------------------------------
# RMSE
# ---------------------------------------------------------------------------


class TestBulkQCRMSE:
    def test_rmse_stored_in_metadata(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert "recon_rmse_mean" in report.metadata

    def test_rmse_non_negative(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.metadata["recon_rmse_mean"] >= 0.0

    def test_rmse_small_for_good_fit(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.metadata["recon_rmse_mean"] < 0.01, (
            "RMSE should be near-zero for noiseless well-separated data"
        )


# ---------------------------------------------------------------------------
# Bootstrap CI width (optional field)
# ---------------------------------------------------------------------------


class TestBulkQCWithBootstrap:
    def test_mean_ci_width_populated_when_ci_present(self):
        from tissueresolve.config import BootstrapConfig
        from tissueresolve.uncertainty.bootstrap import BulkBootstrapCI
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        cfg = BootstrapConfig(n_bootstrap=20, seed=0)
        BulkBootstrapCI(cfg).attach(result, bulk, ref)

        report = BulkQC().compute(result, ref, bulk)
        assert report.mean_ci_width is not None
        assert len(report.mean_ci_width) == ref.n_cell_types

    def test_mean_ci_width_positive(self):
        from tissueresolve.config import BootstrapConfig
        from tissueresolve.uncertainty.bootstrap import BulkBootstrapCI
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        cfg = BootstrapConfig(n_bootstrap=20, seed=0)
        BulkBootstrapCI(cfg).attach(result, bulk, ref)
        report = BulkQC().compute(result, ref, bulk)
        assert (report.mean_ci_width >= 0.0).all()


# ---------------------------------------------------------------------------
# Protocol risk integration
# ---------------------------------------------------------------------------


class TestBulkQCProtocolRisk:
    def test_protocol_risk_level_from_report(self):
        from tissueresolve.protocol.risk import ProtocolRiskReport

        fake_risk = ProtocolRiskReport(
            risk_level="medium",
            mismatch_types=["length_bias"],
            affected_gene_fraction=0.1,
            n_genes_excluded=5,
            excluded_genes=["G0000", "G0001", "G0002", "G0003", "G0004"],
            downweighted_genes=[],
            rationale="Test rationale",
            recommended_actions=[],
            literature_refs=[],
        )
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk, protocol_risk_report=fake_risk)
        assert report.protocol_risk_level == "medium"
        assert report.n_genes_excluded_protocol == 5

    def test_no_protocol_report_gives_unknown(self):
        ref = _make_ref()
        bulk, result = _make_bulk_and_result(ref)
        report = BulkQC().compute(result, ref, bulk)
        assert report.protocol_risk_level == "unknown"


# ---------------------------------------------------------------------------
# Recommendations
# ---------------------------------------------------------------------------


class TestBulkQCRecommendations:
    def test_recommendations_non_empty_for_poor_fit(self):
        # Force poor fit by using mismatched bulk
        rng = np.random.default_rng(99)
        ref = _make_ref()
        n_genes = ref.n_genes
        n_samples = 4
        # Random bulk (not matched to reference)
        bulk = pd.DataFrame(
            rng.exponential(100, (n_genes, n_samples)),
            index=ref.gene_names,
            columns=[f"S{i}" for i in range(n_samples)],
        )
        # Solve anyway
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)

        cfg = BulkQCConfig(r2_fail=0.99, profile_corr_fail=0.99)  # very strict
        report = BulkQC(config=cfg).compute(result, ref, bulk)
        assert len(report.recommendations) > 0
