"""
Tests for tissueresolve.spatial.pipeline (SpatialPipeline) and
uncertainty.bootstrap.SpatialBootstrapCI.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature, SpatialDeconvResult
from tissueresolve.spatial.graph import build_hex_graph_from_arrays
from tissueresolve.spatial.pipeline import SpatialPipeline, SpatialPipelineResult


# ---------------------------------------------------------------------------
# Helpers — synthetic Visium dataset
# ---------------------------------------------------------------------------


def _hex_grid(n_rows=4, n_cols=5):
    rows, cols = [], []
    for r in range(n_rows):
        for c in range(n_cols):
            rows.append(r)
            cols.append(c * 2 + (r % 2))
    return np.array(rows, np.int32), np.array(cols, np.int32)


def _make_synthetic_dataset(n_spots=20, n_types=3, n_genes=30, seed=7):
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    R_cpm = np.zeros((n_types, n_genes), dtype=np.float32)
    for k in range(n_types):
        R_cpm[k, k * block:(k + 1) * block] = rng.uniform(300, 1500, block)
    R_cpm += rng.uniform(0, 10, R_cpm.shape)
    phi_g = np.full(n_genes, 5.0, dtype=np.float32)

    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_types)]
    ref = ReferenceSignature(
        gene_names=gene_names, cell_types=cell_types,
        R_cpm=R_cpm, phi_g=phi_g,
    )

    true_props = rng.dirichlet(np.ones(n_types), size=n_spots).astype(np.float32)
    R_lin = (R_cpm / 1e6)
    lib = rng.integers(500, 3000, size=n_spots).astype(np.float32)
    Mu = lib[:, None] * (true_props @ R_lin)
    phi = 5.0
    p = phi / (phi + Mu + 1e-8)
    p = np.clip(p, 1e-6, 1 - 1e-6)
    Y = rng.negative_binomial(phi * np.ones_like(Mu), p).astype(np.float32)

    r, c = _hex_grid(4, 5)
    r, c = r[:n_spots], c[:n_spots]

    return Y, ref, r, c, lib, gene_names, true_props


# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------


class TestSpatialPipelineEndToEnd:
    def _run(self, **kwargs):
        Y, ref, r, c, lib, gene_names, _ = _make_synthetic_dataset()
        pipeline = SpatialPipeline()
        return pipeline.run(
            Y, ref, r, c, lib, gene_names,
            marker_genes=gene_names,  # use all genes to speed up test
            **kwargs,
        )

    def test_returns_spatial_pipeline_result(self):
        result = self._run()
        assert isinstance(result, SpatialPipelineResult)

    def test_deconv_is_spatial_deconv_result(self):
        result = self._run()
        assert isinstance(result.deconv, SpatialDeconvResult)

    def test_estimate_type_is_spot_rna_composition(self):
        result = self._run()
        assert result.deconv.ESTIMATE_TYPE == "spot_rna_composition"

    def test_proportions_shape(self):
        Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
        result = SpatialPipeline().run(Y, ref, r, c, lib, genes, marker_genes=genes)
        assert result.deconv.proportions.shape == (len(r), ref.n_cell_types)

    def test_proportions_sum_to_one(self):
        result = self._run()
        row_sums = result.deconv.proportions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-4)

    def test_proportions_non_negative(self):
        result = self._run()
        assert (result.deconv.proportions.values >= 0).all()

    def test_cell_types_match_reference(self):
        Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
        result = SpatialPipeline().run(Y, ref, r, c, lib, genes, marker_genes=genes)
        assert list(result.deconv.cell_types) == list(ref.cell_types)

    def test_convergence_trace_recorded(self):
        result = self._run()
        assert len(result.deconv.convergence_trace) > 0

    def test_lambda_spatial_recorded(self):
        result = self._run()
        assert result.deconv.lambda_spatial is not None

    def test_qc_report_modality_is_spatial(self):
        result = self._run()
        assert result.qc.modality == "spatial"

    def test_spot_qc_dataframe_present(self):
        result = self._run()
        assert result.spot_qc is not None
        assert len(result.spot_qc) == 20  # n_spots

    def test_nb_loglik_non_nan_from_pipeline(self):
        """Pipeline must pass model arrays so nb_loglik is never NaN."""
        result = self._run()
        assert not result.spot_qc["nb_loglik"].isna().any(), (
            "Pipeline should pass model arrays to compute_spot_qc so nb_loglik is non-NaN"
        )

    def test_morans_i_series_present(self):
        result = self._run()
        assert result.morans_i is not None
        assert len(result.morans_i) == 3  # n_types

    def test_neighbourhood_none_by_default(self):
        result = self._run()
        assert result.neighbourhood is None

    def test_neighbourhood_present_when_requested(self):
        result = self._run(run_neighbourhood=True, n_neighbourhood_perm=9)
        assert result.neighbourhood is not None

    def test_run_metadata_has_estimate_type(self):
        result = self._run()
        assert result.run_metadata.get("estimate_type") == "spot_rna_composition"

    def test_run_metadata_has_lambda(self):
        result = self._run()
        assert "lambda_spatial" in result.run_metadata

    def test_mismatch_factors_attached(self):
        result = self._run()
        assert result.deconv.mismatch_factors is not None

    def test_mismatch_factors_shape(self):
        result = self._run()
        n_genes = 30
        assert len(result.deconv.mismatch_factors) == n_genes

    def test_marker_genes_stored(self):
        result = self._run()
        assert len(result.deconv.marker_genes) == 30


# ---------------------------------------------------------------------------
# SpatialBootstrapCI
# ---------------------------------------------------------------------------


class TestSpatialBootstrapCI:
    def _fit_model(self):
        from tissueresolve.spatial.model import SpatCARModel
        Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
        n = len(r)
        graph = build_hex_graph_from_arrays(r, c)
        model = SpatCARModel(max_iter=10, lambda_spatial=0.05, verbose=False)
        model.fit(Y, ref.subset_genes(genes), graph, lib)
        return model

    def test_returns_three_items(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        model = self._fit_model()
        result = SpatialBootstrapCI(n_bootstrap=3, n_iter_per_boot=2).compute(model)
        assert len(result) == 3

    def test_ci_lo_shape(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        model = self._fit_model()
        ci_lo, ci_hi, _ = SpatialBootstrapCI(n_bootstrap=3, n_iter_per_boot=2).compute(model)
        assert ci_lo.shape == model.proportions_.shape

    def test_ci_hi_shape(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        model = self._fit_model()
        ci_lo, ci_hi, _ = SpatialBootstrapCI(n_bootstrap=3, n_iter_per_boot=2).compute(model)
        assert ci_hi.shape == model.proportions_.shape

    def test_ci_lo_le_ci_hi(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        model = self._fit_model()
        ci_lo, ci_hi, _ = SpatialBootstrapCI(n_bootstrap=5, n_iter_per_boot=3).compute(model)
        assert (ci_hi >= ci_lo - 1e-6).all()

    def test_coverage_note_non_empty(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        model = self._fit_model()
        _, _, note = SpatialBootstrapCI(n_bootstrap=3, n_iter_per_boot=2).compute(model)
        assert isinstance(note, str)
        assert len(note) > 0

    def test_coverage_note_mentions_empirical(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        _, _, note = SpatialBootstrapCI(n_bootstrap=2, n_iter_per_boot=2).compute(
            self._fit_model()
        )
        assert "empirical" in note.lower() or "coverage" in note.lower()

    def test_attach_method(self):
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        model = self._fit_model()
        Y, ref, r, c, lib, genes, _ = _make_synthetic_dataset()
        result = SpatialPipeline().run(Y, ref, r, c, lib, genes, marker_genes=genes)
        deconv = result.deconv
        assert deconv.lower_ci is None
        SpatialBootstrapCI(n_bootstrap=3, n_iter_per_boot=2).attach(deconv, model)
        assert deconv.lower_ci is not None
        assert deconv.upper_ci is not None
        assert deconv.bootstrap_coverage_note is not None

    def test_ci_metadata_serialisable(self):
        import json
        from tissueresolve.uncertainty.bootstrap import SpatialBootstrapCI
        meta = SpatialBootstrapCI(n_bootstrap=5, n_iter_per_boot=3).ci_metadata
        json.dumps(meta)  # should not raise
