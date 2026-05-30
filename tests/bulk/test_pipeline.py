"""
Tests for tissueresolve.bulk.pipeline (BulkPipeline).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.bulk.pipeline import BulkPipeline, BulkPipelineResult
from tissueresolve.bulk.solver import MRNAContentCorrector
from tissueresolve.config import TissueResolveConfig
from tissueresolve.results import BulkDeconvResult, QCReport, ReferenceSignature


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ref(n_genes: int = 100, n_types: int = 4, seed: int = 2) -> ReferenceSignature:
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    phi_raw = np.zeros((n_genes, n_types))
    for k in range(n_types):
        phi_raw[k * block:(k + 1) * block, k] = rng.uniform(2.0, 5.0, block)
    phi_raw += rng.uniform(0.0, 0.05, phi_raw.shape)
    col_sums = phi_raw.sum(axis=0, keepdims=True)
    phi = (phi_raw / col_sums).astype(np.float64)
    R_cpm = (phi.T * 1e6).astype(np.float32)
    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_types)]
    return ReferenceSignature(
        gene_names=gene_names, cell_types=cell_types,
        phi=phi, R_cpm=R_cpm,
        n_cells_per_type={ct: 50 for ct in cell_types},
    )


def _make_bulk(
    ref: ReferenceSignature,
    n_samples: int = 6,
    seed: int = 13,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    K = ref.n_cell_types
    props = rng.dirichlet(np.ones(K), size=n_samples)
    counts = (ref.as_phi() @ props.T) * 1e4
    return pd.DataFrame(counts, index=ref.gene_names,
                        columns=[f"S{i}" for i in range(n_samples)])


# ---------------------------------------------------------------------------
# End-to-end smoke test
# ---------------------------------------------------------------------------


class TestBulkPipelineEndToEnd:
    def test_runs_without_error(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        result = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert result is not None

    def test_returns_bulk_pipeline_result(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        result = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert isinstance(result, BulkPipelineResult)

    def test_deconv_is_bulk_deconv_result(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert isinstance(pr.deconv, BulkDeconvResult)

    def test_qc_is_qc_report(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert isinstance(pr.qc, QCReport)

    def test_proportions_sum_to_one(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        row_sums = pr.deconv.proportions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-9)

    def test_estimate_type_is_mrna_proportion(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.deconv.ESTIMATE_TYPE == "mRNA_proportion"

    def test_proportions_non_negative(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert (pr.deconv.proportions.values >= 0).all()

    def test_n_samples_matches_bulk(self):
        ref = _make_ref()
        bulk = _make_bulk(ref, n_samples=6)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.deconv.n_samples == 6

    def test_n_cell_types_matches_ref(self):
        ref = _make_ref(n_types=4)
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.deconv.n_cell_types == 4

    def test_recovers_proportions_within_tolerance(self):
        rng = np.random.default_rng(42)
        ref = _make_ref()
        K = ref.n_cell_types
        n_samples = 8
        true_props = rng.dirichlet(np.ones(K), size=n_samples)
        counts = (ref.as_phi() @ true_props.T) * 1e4
        bulk = pd.DataFrame(counts, index=ref.gene_names,
                            columns=[f"S{i}" for i in range(n_samples)])
        true_df = pd.DataFrame(true_props, index=bulk.columns, columns=ref.cell_types)

        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        err = (pr.deconv.proportions - true_df).abs().mean().mean()
        assert err < 0.05, f"Recovery error too high: {err:.4f}"


# ---------------------------------------------------------------------------
# Pre-selected gene panel
# ---------------------------------------------------------------------------


class TestBulkPipelinePreSelectedPanel:
    def test_pre_selected_panel_used(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        panel = ref.gene_names[:40]
        pr = BulkPipeline().run(bulk, ref, gene_panel=panel, n_bootstrap=0)
        assert set(pr.deconv.gene_panel).issubset(set(panel))

    def test_auto_marker_selection_when_no_panel(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.gene_selection is not None

    def test_no_marker_selection_when_panel_provided(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        panel = ref.gene_names[:40]
        pr = BulkPipeline().run(bulk, ref, gene_panel=panel, n_bootstrap=0)
        assert pr.gene_selection is None


# ---------------------------------------------------------------------------
# Bootstrap integration
# ---------------------------------------------------------------------------


class TestBulkPipelineBootstrap:
    def test_bootstrap_populates_ci(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=20)
        assert pr.deconv.lower_ci is not None
        assert pr.deconv.upper_ci is not None

    def test_no_bootstrap_when_zero(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.deconv.lower_ci is None
        assert pr.deconv.upper_ci is None

    def test_ci_shape_correct(self):
        ref = _make_ref()
        bulk = _make_bulk(ref, n_samples=6)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=20)
        assert pr.deconv.lower_ci.shape == (6, ref.n_cell_types)


# ---------------------------------------------------------------------------
# mRNA correction integration
# ---------------------------------------------------------------------------


class TestBulkPipelineMRNACorrection:
    def test_mrna_correction_populates_cell_fractions(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        content = pd.Series({ct: float(k + 1) for k, ct in enumerate(ref.cell_types)})
        corrector = MRNAContentCorrector(content)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0, mrna_corrector=corrector)
        assert pr.deconv.cell_fractions is not None

    def test_mrna_correction_does_not_modify_proportions(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)

        pr_no_corr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        props_no_corr = pr_no_corr.deconv.proportions.copy()

        content = pd.Series({ct: float(k + 1) for k, ct in enumerate(ref.cell_types)})
        corrector = MRNAContentCorrector(content)
        pr_corr = BulkPipeline().run(bulk, ref, n_bootstrap=0, mrna_corrector=corrector)

        # proportions unchanged
        assert np.allclose(pr_corr.deconv.proportions.values, props_no_corr.values)

    def test_cell_fractions_sum_to_one(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        content = pd.Series({ct: float(k + 1) for k, ct in enumerate(ref.cell_types)})
        corrector = MRNAContentCorrector(content)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0, mrna_corrector=corrector)
        row_sums = pr.deconv.cell_fractions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-9)

    def test_no_cell_fractions_without_corrector(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.deconv.cell_fractions is None


# ---------------------------------------------------------------------------
# QC integration
# ---------------------------------------------------------------------------


class TestBulkPipelineQC:
    def test_qc_run_by_default(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.qc.recon_r2 is not None

    def test_qc_skipped_when_disabled(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0, run_qc=False)
        assert pr.qc.recon_r2 is None


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


class TestBulkPipelineMetadata:
    def test_run_metadata_has_estimate_type(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert pr.deconv.run_metadata.get("estimate_type") == "mRNA_proportion"

    def test_run_metadata_has_n_genes_panel(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        pr = BulkPipeline().run(bulk, ref, n_bootstrap=0)
        assert "n_genes_panel" in pr.deconv.run_metadata or "n_genes_panel" in pr.run_metadata
