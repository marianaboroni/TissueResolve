"""
Tests for tissueresolve.bulk.solver (WNNLSSolver + MRNAContentCorrector).
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.bulk.solver import (
    MRNAContentCorrector,
    WNNLSSolver,
    _l1norm_cols,
    _solve_one,
)
from tissueresolve.results import ReferenceSignature


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------


def _make_ref(n_genes: int = 60, n_types: int = 4, seed: int = 1) -> ReferenceSignature:
    """Build a well-separated reference with disjoint marker blocks."""
    rng = np.random.default_rng(seed)
    phi_raw = np.zeros((n_genes, n_types))
    block = n_genes // n_types
    for k in range(n_types):
        phi_raw[k * block:(k + 1) * block, k] = rng.uniform(2.0, 4.0, block)
    phi_raw += rng.uniform(0.0, 0.05, phi_raw.shape)
    col_sums = phi_raw.sum(axis=0, keepdims=True)
    phi = phi_raw / col_sums

    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_types)]
    return ReferenceSignature(
        gene_names=gene_names,
        cell_types=cell_types,
        phi=phi.astype(np.float64),
    )


def _make_bulk(
    ref: ReferenceSignature,
    n_samples: int = 8,
    seed: int = 7,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create noiseless synthetic bulk = phi @ proportions and return true props."""
    rng = np.random.default_rng(seed)
    K = ref.n_cell_types
    G = ref.n_genes
    true_props = rng.dirichlet(np.ones(K), size=n_samples)   # (N, K)
    phi = ref.as_phi()
    counts = (phi @ true_props.T) * 1e4  # (G, N)

    gene_names = ref.gene_names
    sample_ids = [f"S{i}" for i in range(n_samples)]
    bulk = pd.DataFrame(counts, index=gene_names, columns=sample_ids)
    true_df = pd.DataFrame(true_props, index=sample_ids, columns=ref.cell_types)
    return bulk, true_df


# ---------------------------------------------------------------------------
# WNNLSSolver — basic correctness
# ---------------------------------------------------------------------------


class TestWNNLSSolverRecovery:
    def test_recovers_known_proportions(self):
        ref = _make_ref()
        bulk, true_props = _make_bulk(ref)
        solver = WNNLSSolver()
        result = solver.solve(bulk, ref, ref.gene_names)
        # For noiseless well-separated data, should be within 3 %
        err = (result.proportions - true_props).abs().mean().mean()
        assert err < 0.03, f"Mean absolute recovery error too high: {err:.4f}"

    def test_coefficients_are_non_negative(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        assert (result.proportions.values >= 0).all(), "Negative proportion detected"

    def test_proportions_sum_to_one(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        row_sums = result.proportions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-9), "Proportions do not sum to 1"

    def test_returns_bulk_deconv_result(self):
        from tissueresolve.results import BulkDeconvResult
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        assert isinstance(result, BulkDeconvResult)

    def test_estimate_type_is_mrna_proportion(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        assert result.ESTIMATE_TYPE == "mRNA_proportion"

    def test_gene_panel_stored(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        panel = ref.gene_names[:30]
        result = WNNLSSolver().solve(bulk, ref, panel)
        assert set(result.gene_panel).issubset(set(ref.gene_names))

    def test_coverage_r2_positive_for_good_fit(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        assert (result.coverage_r2 > 0).all(), "Expected positive R² for noiseless data"

    def test_r2_near_one_for_noiseless_data(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        assert result.coverage_r2.mean() > 0.90, "Expected high R² for noiseless data"


class TestWNNLSSolverWeights:
    def test_weighted_differs_from_unweighted(self):
        # Use a collinear reference so weights can change the solution.
        rng = np.random.default_rng(99)
        G, K = 40, 3
        # Make two columns of R nearly identical (collinear) to create ambiguity
        base = rng.uniform(0.5, 1.5, (G,))
        phi_raw = np.column_stack([
            base + rng.uniform(0, 0.3, G),
            base + rng.uniform(0, 0.3, G),    # near-duplicate of col 0
            rng.uniform(0.0, 0.1, G),          # sparse third type
        ]).astype(np.float64)
        col_sums = phi_raw.sum(axis=0, keepdims=True)
        phi = phi_raw / col_sums

        gene_names = [f"G{i:04d}" for i in range(G)]
        cell_types = ["CT0", "CT1", "CT2"]
        ref_col = ReferenceSignature(
            gene_names=gene_names, cell_types=cell_types, phi=phi
        )
        # Bulk is a mix skewed toward CT2
        true_p = np.array([0.3, 0.3, 0.4])
        counts = (phi @ true_p) * 1e4 + rng.uniform(0, 5, G)
        bulk = pd.DataFrame(counts[:, None], index=gene_names, columns=["S0"])

        result_unweighted = WNNLSSolver().solve(bulk, ref_col, gene_names)

        # Extreme weights: down-weight first 20 genes, up-weight last 20
        w = np.ones(G)
        w[:20] = 0.01
        w[20:] = 10.0
        weights = pd.Series(w, index=gene_names)
        result_weighted = WNNLSSolver().solve(bulk, ref_col, gene_names, gene_weights=weights)

        diff = (result_weighted.proportions - result_unweighted.proportions).abs().mean().mean()
        assert diff > 1e-6, (
            f"Weighted result identical to unweighted — weights have no effect (diff={diff})"
        )

    def test_uniform_weights_equivalent_to_no_weights(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        panel = ref.gene_names

        result_no_w = WNNLSSolver().solve(bulk, ref, panel)
        uniform_w = pd.Series(np.ones(len(panel)), index=panel)
        result_uniform_w = WNNLSSolver().solve(bulk, ref, panel, gene_weights=uniform_w)

        diff = (result_no_w.proportions - result_uniform_w.proportions).abs().max().max()
        assert diff < 1e-9, "Uniform weights should be equivalent to no weights"

    def test_gene_weights_stored_in_result(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        panel = ref.gene_names
        weights = pd.Series(np.ones(len(panel)), index=panel)
        result = WNNLSSolver().solve(bulk, ref, panel, gene_weights=weights)
        assert result.gene_weights is not None


class TestWNNLSSolverGeneOverlap:
    def test_empty_overlap_raises(self):
        ref = _make_ref()
        bulk, _ = _make_bulk(ref)
        with pytest.raises(ValueError, match="No genes overlap"):
            WNNLSSolver().solve(bulk, ref, ["NONEXISTENT_GENE"])

    def test_low_overlap_warns(self):
        ref = _make_ref(n_genes=60)
        bulk, _ = _make_bulk(ref)
        # Use only 5 genes — below MIN_GENES_WARN (20)
        tiny_panel = ref.gene_names[:5]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            WNNLSSolver().solve(bulk, ref, tiny_panel)
        texts = [str(w.message) for w in caught]
        assert any("only 5 genes" in t for t in texts), (
            "Expected warning about low gene count, got: " + str(texts)
        )

    def test_panel_intersection_with_bulk(self):
        ref = _make_ref(n_genes=60)
        bulk, _ = _make_bulk(ref)
        # Panel includes genes not in bulk; should still work on intersection
        extra = ref.gene_names[:30] + ["FAKE_GENE_1", "FAKE_GENE_2"]
        result = WNNLSSolver().solve(bulk, ref, extra)
        assert "FAKE_GENE_1" not in result.gene_panel


# ---------------------------------------------------------------------------
# MRNAContentCorrector
# ---------------------------------------------------------------------------


class TestMRNAContentCorrector:
    def _make_props(self, n_samples: int = 6) -> pd.DataFrame:
        rng = np.random.default_rng(42)
        raw = rng.dirichlet([2, 3, 1, 4], size=n_samples)
        return pd.DataFrame(raw, columns=["A", "B", "C", "D"])

    def test_correct_returns_dataframe(self):
        mrna = pd.Series({"A": 1.0, "B": 2.0, "C": 0.5, "D": 3.0})
        corrector = MRNAContentCorrector(mrna)
        props = self._make_props()
        result = corrector.correct(props)
        assert isinstance(result, pd.DataFrame)
        assert result.shape == props.shape

    def test_corrected_proportions_sum_to_one(self):
        mrna = pd.Series({"A": 1.0, "B": 2.0, "C": 0.5, "D": 3.0})
        corrector = MRNAContentCorrector(mrna)
        props = self._make_props()
        result = corrector.correct(props)
        assert np.allclose(result.sum(axis=1), 1.0, atol=1e-9)

    def test_corrected_values_non_negative(self):
        mrna = pd.Series({"A": 1.0, "B": 2.0, "C": 0.5, "D": 3.0})
        corrector = MRNAContentCorrector(mrna)
        props = self._make_props()
        result = corrector.correct(props)
        assert (result.values >= 0).all()

    def test_correction_differs_from_proportions_when_content_varies(self):
        mrna = pd.Series({"A": 1.0, "B": 8.0, "C": 0.5, "D": 1.0})
        corrector = MRNAContentCorrector(mrna)
        props = self._make_props()
        result = corrector.correct(props)
        diff = (result - props).abs().mean().mean()
        assert diff > 1e-4, "Correction should change estimates when content varies"

    def test_non_positive_content_raises(self):
        with pytest.raises(ValueError, match="strictly positive"):
            MRNAContentCorrector(pd.Series({"A": 1.0, "B": 0.0, "C": 0.5}))

    def test_uniform_content_warns(self):
        uniform = pd.Series({"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            MRNAContentCorrector(uniform)
        texts = [str(w.message) for w in caught]
        assert any("nearly uniform" in t or "negligible effect" in t for t in texts), (
            "Expected uniformity warning, got: " + str(texts)
        )

    def test_allow_uniform_suppresses_warning(self):
        uniform = pd.Series({"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            MRNAContentCorrector(uniform, allow_uniform=True)
        texts = [str(w.message) for w in caught]
        uniform_warnings = [t for t in texts if "nearly uniform" in t or "negligible effect" in t]
        assert len(uniform_warnings) == 0

    def test_missing_cell_type_warns_uses_default(self):
        mrna = pd.Series({"A": 1.0, "B": 2.0})
        corrector = MRNAContentCorrector(mrna)
        props = self._make_props()  # has A, B, C, D
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = corrector.correct(props)
        texts = [str(w.message) for w in caught]
        assert any("no mRNA content" in t.lower() or "content=1.0" in t for t in texts)
        # Result should still be valid
        assert np.allclose(result.sum(axis=1), 1.0, atol=1e-9)

    def test_correction_summary_shape(self):
        mrna = pd.Series({"A": 1.0, "B": 2.0, "C": 0.5, "D": 3.0})
        corrector = MRNAContentCorrector(mrna)
        props = self._make_props()
        summary = corrector.correction_summary(props)
        assert "cell_type" in summary.columns
        assert "mean_delta" in summary.columns
        assert len(summary) == 4

    def test_from_reference_counts(self):
        rng = np.random.default_rng(42)
        n_cells, n_genes = 100, 20
        counts = pd.DataFrame(
            rng.negative_binomial(5, 0.5, (n_genes, n_cells)).astype(float),
            index=[f"G{i}" for i in range(n_genes)],
            columns=[f"cell_{i}" for i in range(n_cells)],
        )
        labels = pd.Series(
            ["TypeA"] * 50 + ["TypeB"] * 50,
            index=counts.columns,
        )
        corrector = MRNAContentCorrector.from_reference_counts(counts, labels)
        assert isinstance(corrector, MRNAContentCorrector)
        assert "TypeA" in corrector.mrna_content.index
        assert "TypeB" in corrector.mrna_content.index


class TestSolveOneHelper:
    def test_non_negative_output(self):
        rng = np.random.default_rng(0)
        G, K = 20, 3
        R = rng.uniform(0, 1, (G, K))
        b = R @ rng.dirichlet(np.ones(K))
        w = np.ones(G)
        theta = _solve_one(b, R, w)
        assert (theta >= 0).all()

    def test_sums_to_one(self):
        rng = np.random.default_rng(1)
        G, K = 20, 3
        R = rng.uniform(0, 1, (G, K))
        b = R @ rng.dirichlet(np.ones(K))
        w = np.ones(G)
        theta = _solve_one(b, R, w)
        assert abs(theta.sum() - 1.0) < 1e-9
