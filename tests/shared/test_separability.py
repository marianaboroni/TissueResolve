"""
Tests: separability diagnostics and merge_nonseparable_types (including
regression test for the P0-3 dispersion bug).
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from tissueresolve.reference.separability import (
    SeparabilityWarning,
    compute_separability,
    merge_nonseparable_types,
    separability_heatmap_data,
)
from tissueresolve.results import PairSeparability, ReferenceSignature, SeparabilityReport


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ref(
    phi: np.ndarray,
    *,
    phi_g: np.ndarray | None = None,
    n_cells: dict | None = None,
) -> ReferenceSignature:
    G, K = phi.shape
    R_cpm = (phi.T * 1e6 + 1.0).astype(np.float32)
    gene_names = [f"G{i}" for i in range(G)]
    cell_types = [f"CT{k}" for k in range(K)]
    return ReferenceSignature(
        gene_names=gene_names,
        cell_types=cell_types,
        phi=phi.astype(np.float64),
        R_cpm=R_cpm,
        R_log=np.log1p(R_cpm).astype(np.float32),
        phi_g=phi_g,
        n_cells_per_type=n_cells or {f"CT{k}": 50 for k in range(K)},
    )


def _make_identical_ref() -> ReferenceSignature:
    """Reference where CT0 and CT1 are identical (BC = 1.0)."""
    rng = np.random.default_rng(0)
    G = 30
    profile = rng.exponential(1.0, G)
    profile /= profile.sum()
    # CT0 = CT1 = same profile; CT2 = different
    other = rng.exponential(1.0, G)
    other /= other.sum()
    phi = np.column_stack([profile, profile, other])
    return _make_ref(phi)


def _make_well_separated_ref() -> ReferenceSignature:
    """Reference where all pairs are well separated."""
    G = 30
    K = 3
    rng = np.random.default_rng(42)
    phi = np.zeros((G, K))
    genes_per = G // K
    for k in range(K):
        phi[k * genes_per:(k + 1) * genes_per, k] = 1.0
    phi += 1e-6
    col_sums = phi.sum(axis=0, keepdims=True)
    phi /= col_sums
    return _make_ref(phi)


# ---------------------------------------------------------------------------
# compute_separability
# ---------------------------------------------------------------------------

class TestComputeSeparability:

    def test_returns_separability_report(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        assert isinstance(report, SeparabilityReport)

    def test_n_pairs_is_K_choose_2(self):
        ref = _make_well_separated_ref()  # K=3
        report = compute_separability(ref)
        assert len(report.pairs) == 3  # 3*(3-1)/2

    def test_bc_in_zero_one_range(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        for p in report.pairs:
            assert 0.0 <= p.bhattacharyya_coeff <= 1.0

    def test_identical_types_have_bc_near_one(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        max_bc = max(p.bhattacharyya_coeff for p in report.pairs)
        assert max_bc > 0.99

    def test_well_separated_has_low_bc(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        max_bc = max(p.bhattacharyya_coeff for p in report.pairs)
        assert max_bc < 0.80

    def test_pairs_sorted_worst_first(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        bcs = [p.bhattacharyya_coeff for p in report.pairs]
        assert bcs == sorted(bcs, reverse=True)

    def test_risk_level_counts(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        total = report.n_critical + report.n_high + report.n_medium + report.n_ok
        assert total == len(report.pairs)

    def test_separability_warning_emitted(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compute_separability(ref, warn_threshold=0.90)
        sep_warnings = [x for x in w if issubclass(x.category, SeparabilityWarning)]
        assert len(sep_warnings) > 0

    def test_no_warning_for_well_separated(self):
        ref = _make_well_separated_ref()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compute_separability(ref, warn_threshold=0.90)
        sep_warnings = [x for x in w if issubclass(x.category, SeparabilityWarning)]
        assert len(sep_warnings) == 0

    def test_raise_on_critical(self):
        ref = _make_identical_ref()
        with pytest.raises(ValueError, match="critically"):
            compute_separability(ref, warn_threshold=0.90, raise_on_critical=True)

    def test_discriminability_score_is_complement_of_bc(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        for p in report.pairs:
            expected = max(0.0, 1.0 - p.bhattacharyya_coeff)
            assert abs(p.discriminability_score - expected) < 1e-9

    def test_jeffreys_divergence_non_negative(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        for p in report.pairs:
            assert p.jeffreys_divergence >= 0.0

    def test_n_discriminating_genes_non_negative(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        for p in report.pairs:
            assert p.n_discriminating_genes >= 0

    def test_pearson_r_range(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        for p in report.pairs:
            assert -1.0 <= p.pearson_r <= 1.0


# ---------------------------------------------------------------------------
# merge_nonseparable_types — P0-3 regression test
# ---------------------------------------------------------------------------

class TestMergeNonseparableTypes:

    def test_merge_reduces_n_cell_types(self):
        ref = _make_identical_ref()  # CT0 and CT1 should merge
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        merged, _ = merge_nonseparable_types(ref, report, merge_threshold=0.90)
        assert merged.n_cell_types < ref.n_cell_types

    def test_merge_map_covers_all_original_types(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        _, merge_map = merge_nonseparable_types(ref, report)
        assert set(merge_map.keys()) == set(ref.cell_types)

    def test_merged_ref_has_valid_shapes(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        merged, _ = merge_nonseparable_types(ref, report)
        merged.validate()

    def test_merged_phi_columns_sum_to_one(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        merged, _ = merge_nonseparable_types(ref, report)
        np.testing.assert_allclose(
            merged.phi.sum(axis=0), np.ones(merged.n_cell_types), atol=1e-5
        )

    def test_p03_dispersion_bug_regression(self):
        """Regression test for P0-3: merge must NOT apply global maximum dispersion.

        Original SpatCAR bug: the loop ran np.maximum(phi_new, ref.phi_g) on
        EVERY original type, effectively setting phi_new = phi_g regardless
        of group membership.  The correct fix: phi_g.copy() — a single direct
        assignment.

        This test verifies that phi_g is preserved exactly (not inflated by a
        spurious global max loop).
        """
        rng = np.random.default_rng(99)
        G = 20
        # Build two types that are nearly identical → will be merged
        profile = rng.exponential(1.0, G)
        profile /= profile.sum()
        phi = np.column_stack([profile, profile + 1e-4, rng.exponential(1.0, G)])
        phi = phi / phi.sum(axis=0, keepdims=True)
        # Deliberately set phi_g with specific values we can verify
        phi_g_orig = rng.uniform(1.0, 50.0, G).astype(np.float32)

        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=["CT0", "CT1", "CT2"],
            phi=phi.astype(np.float64),
            R_cpm=(phi.T * 1e6 + 1.0).astype(np.float32),
            phi_g=phi_g_orig.copy(),
            n_cells_per_type={"CT0": 50, "CT1": 50, "CT2": 50},
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        merged, _ = merge_nonseparable_types(ref, report, merge_threshold=0.80)

        # phi_g must be the original phi_g, not a global-max-corrupted version
        assert merged.phi_g is not None
        np.testing.assert_array_equal(
            merged.phi_g, phi_g_orig,
            err_msg=(
                "P0-3 regression: merged phi_g must equal original phi_g. "
                "The loop-based global-max bug would have produced a different result."
            ),
        )

    def test_no_merge_when_all_well_separated(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        merged, merge_map = merge_nonseparable_types(ref, report, merge_threshold=0.90)
        # Nothing should be merged
        assert merged.n_cell_types == ref.n_cell_types
        # All types map to themselves
        for orig, merged_label in merge_map.items():
            assert merged_label == orig

    def test_gene_names_preserved_after_merge(self):
        ref = _make_identical_ref()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        merged, _ = merge_nonseparable_types(ref, report)
        assert merged.gene_names == ref.gene_names

    def test_donor_cv_none_after_merge(self):
        """donor_cv cannot be meaningfully merged; must be set to None."""
        rng = np.random.default_rng(100)
        G = 15
        profile = rng.exponential(1.0, G)
        profile /= profile.sum()
        phi = np.column_stack([profile, profile, rng.exponential(1.0, G)])
        phi = phi / phi.sum(axis=0, keepdims=True)
        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=["CT0", "CT1", "CT2"],
            phi=phi.astype(np.float64),
            R_cpm=(phi.T * 1e6 + 1.0).astype(np.float32),
            donor_cv=np.zeros((G, 3)),
            n_cells_per_type={"CT0": 40, "CT1": 40, "CT2": 40},
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeparabilityWarning)
            report = compute_separability(ref)
        merged, _ = merge_nonseparable_types(ref, report)
        # donor_cv is not preserved through merge
        assert merged.donor_cv is None


# ---------------------------------------------------------------------------
# separability_heatmap_data
# ---------------------------------------------------------------------------

class TestSeparabilityHeatmapData:

    def test_shape_is_K_by_K(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        M = separability_heatmap_data(report, ref.cell_types)
        assert M.shape == (ref.n_cell_types, ref.n_cell_types)

    def test_diagonal_is_one(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        M = separability_heatmap_data(report, ref.cell_types)
        np.testing.assert_allclose(np.diag(M), np.ones(ref.n_cell_types))

    def test_symmetric(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        M = separability_heatmap_data(report, ref.cell_types)
        np.testing.assert_array_equal(M, M.T)

    def test_values_in_zero_one(self):
        ref = _make_well_separated_ref()
        report = compute_separability(ref)
        M = separability_heatmap_data(report, ref.cell_types)
        assert (M >= 0.0).all() and (M <= 1.0).all()
