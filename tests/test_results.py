"""
Tests: result class instantiation, save/load round-trips, and
scientific contract enforcement.
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import (
    BenchmarkResult,
    BulkDeconvResult,
    PairSeparability,
    QCReport,
    ReferenceSignature,
    SeparabilityReport,
    SpatialDeconvResult,
)

# ---------------------------------------------------------------------------
# ReferenceSignature
# ---------------------------------------------------------------------------

class TestReferenceSignature:

    def test_basic_construction(self, ref_sig_full):
        assert ref_sig_full.n_genes == 40
        assert ref_sig_full.n_cell_types == 4

    def test_repr(self, ref_sig_full):
        r = repr(ref_sig_full)
        assert "40g" in r
        assert "4k" in r

    def test_as_phi_direct(self, ref_sig_full):
        phi = ref_sig_full.as_phi()
        assert phi.shape == (40, 4)
        np.testing.assert_allclose(phi.sum(axis=0), np.ones(4), atol=1e-6)

    def test_as_phi_derived_from_R_cpm(self, ref_sig_cpm_only):
        phi = ref_sig_cpm_only.as_phi()
        assert phi.shape == (40, 4)
        np.testing.assert_allclose(phi.sum(axis=0), np.ones(4), atol=1e-5)

    def test_as_phi_raises_when_empty(self):
        ref = ReferenceSignature(gene_names=["G1"], cell_types=["A"])
        with pytest.raises(ValueError, match="neither phi nor R_cpm"):
            ref.as_phi()

    def test_as_R_cpm_direct(self, ref_sig_full):
        R = ref_sig_full.as_R_cpm()
        assert R.shape == (4, 40)

    def test_as_R_cpm_derived_from_phi(self, ref_sig_phi_only):
        R = ref_sig_phi_only.as_R_cpm()
        assert R.shape == (4, 40)
        # Derived from phi: rows should sum to ~1e6
        np.testing.assert_allclose(R.sum(axis=1), np.full(4, 1e6), rtol=1e-5)

    def test_as_R_cpm_raises_when_empty(self):
        ref = ReferenceSignature(gene_names=["G1"], cell_types=["A"])
        with pytest.raises(ValueError, match="neither R_cpm nor phi"):
            ref.as_R_cpm()

    def test_as_R_log(self, ref_sig_full):
        rl = ref_sig_full.as_R_log()
        assert rl.shape == (4, 40)
        assert (rl >= 0).all()

    def test_subset_genes(self, ref_sig_full):
        subset_genes = ref_sig_full.gene_names[:10]
        sub = ref_sig_full.subset_genes(subset_genes)
        assert sub.n_genes == 10
        assert sub.phi.shape == (10, 4)
        assert sub.R_cpm.shape == (4, 10)
        assert sub.phi_g.shape == (10,)
        assert sub.donor_cv.shape == (10, 4)

    def test_subset_genes_warns_on_missing(self, ref_sig_full):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sub = ref_sig_full.subset_genes(["GENE_000", "NONEXISTENT_GENE"])
        assert any("excluded" in str(warning.message).lower() for warning in w)
        assert sub.n_genes == 1

    def test_subset_genes_raises_on_all_missing(self, ref_sig_full):
        with pytest.raises(ValueError, match="no requested genes"):
            ref_sig_full.subset_genes(["FAKE_GENE_1", "FAKE_GENE_2"])

    def test_validate_passes_on_valid(self, ref_sig_full):
        ref_sig_full.validate()  # should not raise

    def test_validate_warns_on_empty(self):
        ref = ReferenceSignature(gene_names=["G1"], cell_types=["A"])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ref.validate()
        assert any("neither phi" in str(warning.message).lower() for warning in w)

    def test_validate_raises_on_bad_phi_shape(self, phi_matrix):
        ref = ReferenceSignature(
            gene_names=["G1", "G2"], cell_types=["A", "B"],
            phi=phi_matrix,  # wrong shape: (40, 4) not (2, 2)
        )
        with pytest.raises(ValueError, match="phi shape"):
            ref.validate()

    def test_save_load_roundtrip_full(self, ref_sig_full, tmp_path):
        save_dir = tmp_path / "ref_full"
        ref_sig_full.save(save_dir)
        loaded = ReferenceSignature.load(save_dir)

        assert loaded.gene_names == ref_sig_full.gene_names
        assert loaded.cell_types == ref_sig_full.cell_types
        assert loaded.genome == ref_sig_full.genome
        assert loaded.n_cells_per_type == ref_sig_full.n_cells_per_type
        np.testing.assert_array_equal(loaded.phi, ref_sig_full.phi)
        np.testing.assert_array_equal(loaded.R_cpm, ref_sig_full.R_cpm)
        np.testing.assert_array_equal(loaded.phi_g, ref_sig_full.phi_g)
        np.testing.assert_array_equal(loaded.donor_cv, ref_sig_full.donor_cv)

    def test_save_load_roundtrip_phi_only(self, ref_sig_phi_only, tmp_path):
        save_dir = tmp_path / "ref_phi"
        ref_sig_phi_only.save(save_dir)
        loaded = ReferenceSignature.load(save_dir)

        assert loaded.R_cpm is None
        assert loaded.phi_g is None
        np.testing.assert_array_equal(loaded.phi, ref_sig_phi_only.phi)

    def test_save_creates_metadata_json(self, ref_sig_full, tmp_path):
        save_dir = tmp_path / "ref_meta"
        ref_sig_full.save(save_dir)
        meta_path = save_dir / "metadata.json"
        assert meta_path.exists()
        with meta_path.open() as fh:
            meta = json.load(fh)
        assert meta["tissueresolve_object"] == "ReferenceSignature"
        assert meta["n_genes"] == 40
        assert meta["n_cell_types"] == 4
        assert meta["has_phi"] is True
        assert meta["has_R_cpm"] is True

    def test_load_raises_on_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            ReferenceSignature.load(tmp_path / "nonexistent")


# ---------------------------------------------------------------------------
# BulkDeconvResult
# ---------------------------------------------------------------------------

class TestBulkDeconvResult:

    def test_estimate_type_is_mrna_proportion(self, bulk_result):
        """DESIGN_SPEC: bulk outputs must be labelled as mRNA proportions."""
        assert bulk_result.ESTIMATE_TYPE == "mRNA_proportion"

    def test_estimate_type_cannot_be_overridden_at_init(self):
        """ESTIMATE_TYPE must not be settable to a misleading value at construction."""
        proportions = pd.DataFrame(
            np.ones((2, 2)) / 2,
            index=["s1", "s2"],
            columns=["A", "B"],
        )
        result = BulkDeconvResult(
            proportions=proportions,
            coverage_r2=pd.Series([0.9, 0.8], index=["s1", "s2"]),
            gene_panel=["G1", "G2"],
        )
        assert result.ESTIMATE_TYPE == "mRNA_proportion"

    def test_n_samples_and_n_cell_types(self, bulk_result):
        assert bulk_result.n_samples == 6
        assert bulk_result.n_cell_types == 4

    def test_summary(self, bulk_result):
        s = bulk_result.summary()
        assert set(s.columns) == {"mean", "sd", "median", "min", "max"}
        assert len(s) == 4

    def test_save_load_roundtrip_minimal(self, bulk_result, tmp_path):
        save_dir = tmp_path / "bulk_min"
        bulk_result.save(save_dir)
        loaded = BulkDeconvResult.load(save_dir)

        pd.testing.assert_frame_equal(loaded.proportions, bulk_result.proportions)
        pd.testing.assert_series_equal(loaded.coverage_r2, bulk_result.coverage_r2)
        assert loaded.gene_panel == bulk_result.gene_panel
        assert loaded.lower_ci is None
        assert loaded.upper_ci is None
        assert loaded.cell_fractions is None

    def test_save_load_roundtrip_with_ci(self, bulk_result_with_ci, tmp_path):
        save_dir = tmp_path / "bulk_ci"
        bulk_result_with_ci.save(save_dir)
        loaded = BulkDeconvResult.load(save_dir)

        assert loaded.lower_ci is not None
        assert loaded.upper_ci is not None
        pd.testing.assert_frame_equal(loaded.lower_ci, bulk_result_with_ci.lower_ci)

    def test_proportions_tsv_has_mrna_warning(self, bulk_result, tmp_path):
        """Saved TSV must contain the mRNA proportion warning comment."""
        save_dir = tmp_path / "bulk_warn"
        bulk_result.save(save_dir)
        content = (save_dir / "proportions.tsv").read_text()
        assert "mRNA_proportion" in content
        assert "NOT cell fractions" in content

    def test_metadata_json_has_estimate_type(self, bulk_result, tmp_path):
        save_dir = tmp_path / "bulk_meta"
        bulk_result.save(save_dir)
        with (save_dir / "metadata.json").open() as fh:
            meta = json.load(fh)
        assert meta["estimate_type"] == "mRNA_proportion"
        assert meta["tissueresolve_object"] == "BulkDeconvResult"

    def test_load_raises_on_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            BulkDeconvResult.load(tmp_path / "nonexistent")


# ---------------------------------------------------------------------------
# SpatialDeconvResult
# ---------------------------------------------------------------------------

class TestSpatialDeconvResult:

    def test_estimate_type_is_spot_rna_composition(self, spatial_result):
        """DESIGN_SPEC: spatial outputs must not be labelled as cell counts."""
        assert spatial_result.ESTIMATE_TYPE == "spot_rna_composition"

    def test_n_spots_and_n_cell_types(self, spatial_result):
        assert spatial_result.n_spots == 20
        assert spatial_result.n_cell_types == 4

    def test_save_load_roundtrip(self, spatial_result, tmp_path):
        save_dir = tmp_path / "spatial"
        spatial_result.save(save_dir)
        loaded = SpatialDeconvResult.load(save_dir)

        pd.testing.assert_frame_equal(loaded.proportions, spatial_result.proportions)
        assert loaded.marker_genes == spatial_result.marker_genes
        assert loaded.n_iter == spatial_result.n_iter
        assert loaded.converged == spatial_result.converged
        assert loaded.lambda_spatial == pytest.approx(spatial_result.lambda_spatial)
        np.testing.assert_array_equal(
            loaded.mismatch_factors, spatial_result.mismatch_factors
        )

    def test_proportions_tsv_records_lambda(self, spatial_result, tmp_path):
        """lambda_spatial must always be written to the proportions TSV header."""
        save_dir = tmp_path / "spatial_lambda"
        spatial_result.save(save_dir)
        content = (save_dir / "proportions.tsv").read_text()
        assert "lambda_spatial" in content

    def test_metadata_json_records_n_smooth(self, spatial_result, tmp_path):
        """n_smooth must be persisted even when None."""
        save_dir = tmp_path / "spatial_nsmooth"
        spatial_result.save(save_dir)
        with (save_dir / "metadata.json").open() as fh:
            meta = json.load(fh)
        # n_smooth is None for this fixture; it must still be recorded
        assert "n_smooth" in meta

    def test_metadata_json_has_estimate_type(self, spatial_result, tmp_path):
        save_dir = tmp_path / "spatial_et"
        spatial_result.save(save_dir)
        with (save_dir / "metadata.json").open() as fh:
            meta = json.load(fh)
        assert meta["estimate_type"] == "spot_rna_composition"

    def test_convergence_trace_roundtrip(self, spatial_result, tmp_path):
        save_dir = tmp_path / "spatial_trace"
        spatial_result.save(save_dir)
        loaded = SpatialDeconvResult.load(save_dir)
        assert len(loaded.convergence_trace) == len(spatial_result.convergence_trace)
        assert loaded.convergence_trace[0] == pytest.approx(
            spatial_result.convergence_trace[0], rel=1e-5
        )

    def test_load_raises_on_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            SpatialDeconvResult.load(tmp_path / "nonexistent")


# ---------------------------------------------------------------------------
# QCReport
# ---------------------------------------------------------------------------

class TestQCReport:

    def test_modality_validation(self):
        with pytest.raises(ValueError, match="modality must be"):
            QCReport(modality="invalid")

    def test_bulk_qc_save_load_roundtrip(self, bulk_qc_report, tmp_path):
        save_dir = tmp_path / "bulk_qc"
        bulk_qc_report.save(save_dir)
        loaded = QCReport.load(save_dir)

        assert loaded.modality == "bulk"
        assert loaded.condition_number == pytest.approx(bulk_qc_report.condition_number)
        assert loaded.protocol_risk_level == bulk_qc_report.protocol_risk_level
        pd.testing.assert_series_equal(
            loaded.recon_r2, bulk_qc_report.recon_r2, check_names=False
        )

    def test_spatial_qc_save_load_roundtrip(self, spatial_qc_report, tmp_path):
        save_dir = tmp_path / "spatial_qc"
        spatial_qc_report.save(save_dir)
        loaded = QCReport.load(save_dir)

        assert loaded.modality == "spatial"
        assert loaded.spot_qc is not None
        pd.testing.assert_series_equal(
            loaded.morans_i, spatial_qc_report.morans_i, check_names=False
        )

    def test_recommendations_roundtrip(self, bulk_qc_report, tmp_path):
        save_dir = tmp_path / "bulk_qc_recs"
        bulk_qc_report.save(save_dir)
        loaded = QCReport.load(save_dir)
        assert loaded.recommendations == bulk_qc_report.recommendations

    def test_heuristic_thresholds_in_per_sample_comment(self, bulk_qc_report, tmp_path):
        """Heuristic threshold notice must appear in the per-sample QC TSV."""
        save_dir = tmp_path / "bulk_qc_heuristic"
        bulk_qc_report.save(save_dir)
        content = (save_dir / "per_sample_qc.tsv").read_text()
        assert "heuristic" in content.lower()


# ---------------------------------------------------------------------------
# SeparabilityReport
# ---------------------------------------------------------------------------

class TestSeparabilityReport:

    def test_risk_level_counts(self, separability_report):
        assert separability_report.n_critical == 1  # BC=0.98
        assert separability_report.n_high == 1      # BC=0.91
        assert separability_report.has_problems is True

    def test_pair_risk_levels(self):
        p_ok = PairSeparability("A", "B", 0.70, 2.0, 0.60, 100)
        p_med = PairSeparability("A", "B", 0.85, 1.0, 0.80, 50)
        p_high = PairSeparability("A", "B", 0.93, 0.5, 0.90, 10)
        p_crit = PairSeparability("A", "B", 0.98, 0.1, 0.97, 2)
        assert p_ok.risk_level == "OK"
        assert p_med.risk_level == "MEDIUM"
        assert p_high.risk_level == "HIGH"
        assert p_crit.risk_level == "CRITICAL"

    def test_is_problematic(self):
        p = PairSeparability("A", "B", 0.95, 0.3, 0.93, 5)
        assert p.is_problematic is True
        p2 = PairSeparability("A", "B", 0.75, 1.5, 0.65, 80)
        assert p2.is_problematic is False

    def test_summary_str(self, separability_report):
        s = separability_report.summary_str()
        assert "CRITICAL" in s
        assert "HIGH" in s
        assert "SeparabilityReport" in s

    def test_save_load_roundtrip(self, separability_report, tmp_path):
        save_dir = tmp_path / "sep"
        separability_report.save(save_dir)
        loaded = SeparabilityReport.load(save_dir)

        assert len(loaded.pairs) == len(separability_report.pairs)
        assert loaded.n_critical == separability_report.n_critical
        assert loaded.n_high == separability_report.n_high
        assert loaded.pairs[0].bhattacharyya_coeff == pytest.approx(
            separability_report.pairs[0].bhattacharyya_coeff
        )

    def test_saved_tsv_contains_bc_warning(self, separability_report, tmp_path):
        save_dir = tmp_path / "sep_warn"
        separability_report.save(save_dir)
        content = (save_dir / "separability_pairs.tsv").read_text()
        assert "BC > 0.90" in content
        assert "unreliable" in content

    def test_load_raises_on_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            SeparabilityReport.load(tmp_path / "nonexistent")


# ---------------------------------------------------------------------------
# BenchmarkResult
# ---------------------------------------------------------------------------

class TestBenchmarkResult:

    def test_basic_construction(self, benchmark_result):
        assert benchmark_result.method == "SpatCAR"
        assert benchmark_result.scenario == "basic"
        assert "rmse" in benchmark_result.metrics

    def test_save_load_roundtrip(self, benchmark_result, tmp_path):
        save_dir = tmp_path / "bench"
        benchmark_result.save(save_dir)
        loaded = BenchmarkResult.load(save_dir)

        assert loaded.method == benchmark_result.method
        assert loaded.scenario == benchmark_result.scenario
        assert loaded.metrics["rmse"] == pytest.approx(benchmark_result.metrics["rmse"])
        assert loaded.n_samples_or_spots == benchmark_result.n_samples_or_spots
        assert loaded.run_time_s == pytest.approx(benchmark_result.run_time_s)

    def test_per_celltype_metrics_roundtrip(self, benchmark_result, tmp_path):
        save_dir = tmp_path / "bench_ct"
        benchmark_result.save(save_dir)
        loaded = BenchmarkResult.load(save_dir)
        assert loaded.per_celltype_metrics is not None
        pd.testing.assert_frame_equal(
            loaded.per_celltype_metrics, benchmark_result.per_celltype_metrics
        )

    def test_load_raises_on_missing_dir(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            BenchmarkResult.load(tmp_path / "nonexistent")
