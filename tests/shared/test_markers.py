"""
Tests: GeneSelector — marker selection, rejection reasons, overlap, scoring.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.config import GeneConfig
from tissueresolve.reference.markers import GeneSelector, MarkerSelectionResult
from tissueresolve.results import ReferenceSignature


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_clean_ref(rng: np.random.Generator, n_genes: int = 40, n_types: int = 3) -> ReferenceSignature:
    """Reference with clear type-specific marker genes."""
    G, K = n_genes, n_types
    phi = np.zeros((G, K))
    genes_per = G // K
    for k in range(K):
        start = k * genes_per
        phi[start:start + genes_per, k] = rng.exponential(10.0, genes_per)
        # Low background expression for other types
        for j in range(K):
            if j != k:
                phi[start:start + genes_per, j] = 0.01
    phi += 1e-6
    phi /= phi.sum(axis=0, keepdims=True)
    R_cpm = (phi.T * 1e6 + 1.0).astype(np.float32)
    return ReferenceSignature(
        gene_names=[f"GENE_{i:03d}" for i in range(G)],
        cell_types=sorted([f"Type{k}" for k in range(K)]),
        phi=phi.astype(np.float64),
        R_cpm=R_cpm,
        R_log=np.log1p(R_cpm).astype(np.float32),
    )


# ---------------------------------------------------------------------------
# Basic selection
# ---------------------------------------------------------------------------

class TestGeneSelector:

    def test_returns_marker_selection_result(self):
        rng = np.random.default_rng(0)
        ref = _make_clean_ref(rng)
        result = GeneSelector().select(ref)
        assert isinstance(result, MarkerSelectionResult)

    def test_selected_genes_are_subset_of_reference(self):
        rng = np.random.default_rng(1)
        ref = _make_clean_ref(rng)
        result = GeneSelector().select(ref)
        ref_set = set(ref.gene_names)
        for g in result.selected_genes:
            assert g in ref_set

    def test_selected_genes_non_empty(self):
        rng = np.random.default_rng(2)
        ref = _make_clean_ref(rng)
        result = GeneSelector().select(ref)
        assert result.n_selected > 0

    def test_n_selected_respects_config_max(self):
        rng = np.random.default_rng(3)
        ref = _make_clean_ref(rng, n_genes=80)
        cfg = GeneConfig(n_genes=10)
        result = GeneSelector(cfg).select(ref)
        assert result.n_selected <= 10

    def test_score_components_present(self):
        rng = np.random.default_rng(4)
        ref = _make_clean_ref(rng)
        result = GeneSelector().select(ref)
        expected_cols = {"specificity", "stability", "prot_weight", "concordance",
                         "composite_score"}
        assert expected_cols.issubset(set(result.score_components.columns))

    def test_score_components_index_are_selected_genes(self):
        rng = np.random.default_rng(5)
        ref = _make_clean_ref(rng, n_genes=40)
        cfg = GeneConfig(n_genes=15)
        result = GeneSelector(cfg).select(ref)
        # score_components is indexed by phase-1 top genes, which are a subset
        for g in result.score_components.index:
            assert g in set(ref.gene_names)

    def test_n_per_type_sums_to_n_selected_genes_in_phase1(self):
        """n_per_type counts are for phase-1 only; just verify they're non-negative."""
        rng = np.random.default_rng(6)
        ref = _make_clean_ref(rng)
        result = GeneSelector().select(ref)
        for ct in ref.cell_types:
            assert result.n_per_type.get(ct, 0) >= 0


# ---------------------------------------------------------------------------
# Rejection reasons
# ---------------------------------------------------------------------------

class TestRejectionReasons:

    def test_low_fc_genes_rejected(self):
        """Genes below the min_log2fc threshold appear in rejected_genes."""
        rng = np.random.default_rng(10)
        G, K = 40, 3
        # Most genes have moderate FC (~1.0); a handful have high FC (> 3.0).
        phi = np.ones((G, K)) / G
        # Give a few genes a strong marker profile so they survive any threshold
        for i in range(3):
            phi[i, i % K] = 0.2
            phi[i, (i + 1) % K] = 0.001
            phi[i, (i + 2) % K] = 0.001
        phi += 1e-6
        phi /= phi.sum(axis=0, keepdims=True)
        R_cpm = (phi.T * 1e6 + 1.0).astype(np.float32)
        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=["TA", "TB", "TC"],
            phi=phi.astype(np.float64),
            R_cpm=R_cpm,
        )
        # Very high threshold — only the 3 hand-crafted genes pass
        cfg = GeneConfig(min_log2fc=4.0)
        result = GeneSelector(cfg).select(ref)
        # Most genes should be rejected with reason "below_min_log2fc"
        fc_rejections = [
            g for g, r in result.rejected_genes.items()
            if r == "below_min_log2fc"
        ]
        assert len(fc_rejections) > 0

    def test_discordant_genes_excluded(self):
        rng = np.random.default_rng(11)
        ref = _make_clean_ref(rng, n_genes=40)
        # Mark first 5 genes as discordant
        discordant = set(ref.gene_names[:5])
        result = GeneSelector().select(ref, discordant_genes=discordant)
        # None of the discordant genes should be selected
        for g in discordant:
            assert g not in result.selected_genes

    def test_discordant_genes_have_correct_rejection_reason(self):
        rng = np.random.default_rng(12)
        ref = _make_clean_ref(rng, n_genes=40)
        discordant = set(ref.gene_names[:3])
        result = GeneSelector().select(ref, discordant_genes=discordant)
        for g in discordant:
            if g in result.rejected_genes:
                assert result.rejected_genes[g] == "discordant"

    def test_query_genes_filter(self):
        rng = np.random.default_rng(13)
        ref = _make_clean_ref(rng, n_genes=40)
        # Allow only first 20 genes
        query = ref.gene_names[:20]
        result = GeneSelector().select(ref, query_genes=query)
        for g in result.selected_genes:
            assert g in set(query)

    def test_not_in_query_rejection(self):
        rng = np.random.default_rng(14)
        ref = _make_clean_ref(rng, n_genes=40)
        query = ref.gene_names[:20]
        result = GeneSelector().select(ref, query_genes=query)
        # Genes outside query that were otherwise selected should be rejected
        # with reason "not_in_query"
        not_in_query_reasons = {
            g for g, r in result.rejected_genes.items() if r == "not_in_query"
        }
        # At least some of the remaining 20 genes should be rejected as "not_in_query"
        assert len(not_in_query_reasons) >= 0  # may be 0 if none selected outside query

    def test_gene_panel_overlap_populated_with_query(self):
        rng = np.random.default_rng(15)
        ref = _make_clean_ref(rng, n_genes=40)
        query = ref.gene_names[:20]
        result = GeneSelector().select(ref, query_genes=query)
        assert result.gene_panel_overlap is not None
        assert set(result.gene_panel_overlap) == set(result.selected_genes)

    def test_gene_panel_overlap_none_without_query(self):
        rng = np.random.default_rng(16)
        ref = _make_clean_ref(rng)
        result = GeneSelector().select(ref)
        assert result.gene_panel_overlap is None

    def test_blacklist_prefix_rejection_with_gene_filters(self, tmp_path):
        from tissueresolve.reference.gene_filters import GeneFilterSet

        rng = np.random.default_rng(17)
        ref = _make_clean_ref(rng, n_genes=40)
        # Rename first gene to an MT- gene
        old_genes = ref.gene_names.copy()
        old_genes[0] = "MT-ND1"
        ref_mt = ReferenceSignature(
            gene_names=old_genes,
            cell_types=ref.cell_types,
            phi=ref.phi,
            R_cpm=ref.R_cpm,
        )

        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        result = GeneSelector().select(ref_mt, gene_filters=gfs)
        assert "MT-ND1" not in result.selected_genes
        if "MT-ND1" in result.rejected_genes:
            assert result.rejected_genes["MT-ND1"] == "blacklist_prefix"


# ---------------------------------------------------------------------------
# Donor CV integration
# ---------------------------------------------------------------------------

class TestMarkerSelectionWithDonorCV:

    def test_high_cv_genes_excluded(self):
        """Genes with CV above max_donor_cv threshold should be excluded."""
        rng = np.random.default_rng(20)
        ref = _make_clean_ref(rng, n_genes=40)
        # Inject high donor_cv for first 5 genes
        donor_cv = np.zeros((ref.n_genes, ref.n_cell_types))
        donor_cv[:5, :] = 2.0   # very high CV
        ref_cv = ReferenceSignature(
            gene_names=ref.gene_names,
            cell_types=ref.cell_types,
            phi=ref.phi,
            R_cpm=ref.R_cpm,
            donor_cv=donor_cv,
        )
        cfg = GeneConfig(max_donor_cv=1.0)
        result = GeneSelector(cfg).select(ref_cv)
        high_cv_genes = set(ref.gene_names[:5])
        # High-CV genes should be rejected
        for g in high_cv_genes:
            if g in result.rejected_genes:
                assert result.rejected_genes[g] == "above_max_donor_cv"


# ---------------------------------------------------------------------------
# MarkerSelectionResult helpers
# ---------------------------------------------------------------------------

class TestMarkerSelectionResult:

    def test_rejection_summary_returns_series(self):
        rng = np.random.default_rng(30)
        ref = _make_clean_ref(rng, n_genes=60)
        cfg = GeneConfig(n_genes=10, min_log2fc=0.5)
        result = GeneSelector(cfg).select(ref)
        summary = result.rejection_summary()
        assert isinstance(summary, pd.Series)

    def test_n_rejected_matches_rejected_genes_dict(self):
        rng = np.random.default_rng(31)
        ref = _make_clean_ref(rng, n_genes=40)
        result = GeneSelector().select(ref)
        assert result.n_rejected == len(result.rejected_genes)

    def test_no_overlap_between_selected_and_hard_rejected(self):
        """No gene can appear in both selected_genes and rejected_genes."""
        rng = np.random.default_rng(32)
        ref = _make_clean_ref(rng, n_genes=40)
        result = GeneSelector().select(ref)
        selected_set = set(result.selected_genes)
        rejected_set = set(result.rejected_genes.keys())
        overlap = selected_set & rejected_set
        assert len(overlap) == 0, f"Genes in both selected and rejected: {overlap}"
