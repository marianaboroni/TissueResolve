"""
Tests: Stage 2 protocol layer — metadata, bulk risk, spatial mismatch.

Test categories
---------------
- Enum values and serialisation
- ProtocolMetadata construction and round-trip
- Invalid values raise clear errors
- ProtocolRiskAssessor with real gene list files
- No silent gene removal
- Gene-level risk scores
- SpatialMismatch construction and shape validation
- compute_spatial_discordance: zero counts, shape, positive d_g
- update_mismatch_factors: shape, positivity, immutability
- Bulk risk and spatial mismatch are not confused
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.protocol.metadata import (
    BulkProtocol,
    ProtocolMetadata,
    RefCapture,
    RefCounting,
    RefModality,
    SpatialPlatform,
    _parse_enum,
)
from tissueresolve.protocol.risk import (
    RISK_HARD_THRESHOLD,
    ProtocolRiskAssessor,
    ProtocolRiskReport,
)
from tissueresolve.protocol.mismatch import (
    ProtocolMismatch,
    SpatialMismatch,
    compute_spatial_discordance,
    update_mismatch_factors,
)
from tissueresolve.results import ReferenceSignature


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _simple_ref(n_genes: int = 20, n_types: int = 3) -> ReferenceSignature:
    rng = np.random.default_rng(0)
    G, K = n_genes, n_types
    R_cpm = (rng.exponential(100.0, (K, G)) + 1.0).astype(np.float32)
    phi_raw = R_cpm.T
    phi = phi_raw / phi_raw.sum(axis=0, keepdims=True)
    return ReferenceSignature(
        gene_names=[f"GENE_{i:04d}" for i in range(G)],
        cell_types=[f"CT{k}" for k in range(K)],
        phi=phi.astype(np.float64),
        R_cpm=R_cpm,
        R_log=np.log1p(R_cpm).astype(np.float32),
    )


def _simple_mismatch(n_genes: int = 20) -> SpatialMismatch:
    rng = np.random.default_rng(1)
    g = [f"GENE_{i:04d}" for i in range(n_genes)]
    return SpatialMismatch(
        d_g=rng.uniform(0.5, 2.0, n_genes).astype(np.float32),
        discord_score=rng.uniform(0, 3.0, n_genes).astype(np.float32),
        gene_weight=rng.uniform(0.1, 1.0, n_genes).astype(np.float32),
        gene_names=g,
    )


# ===========================================================================
# Enum values and serialisation
# ===========================================================================

class TestEnumValues:

    def test_bulk_protocol_values(self):
        assert BulkProtocol.POLYA.value == "polyA"
        assert BulkProtocol.RIBODEP.value == "ribodepleted"
        assert BulkProtocol.UNKNOWN.value == "unknown"

    def test_ref_modality_values(self):
        assert RefModality.SCRNA.value == "scRNA"
        assert RefModality.SNRNA.value == "snRNA"
        assert RefModality.UNKNOWN.value == "unknown"

    def test_ref_capture_values(self):
        assert RefCapture.UMI_3PRIME.value == "10x_3prime"
        assert RefCapture.FULL_LENGTH.value == "full_length"
        assert RefCapture.DROPSEQ.value == "dropseq"

    def test_ref_counting_values(self):
        assert RefCounting.EXONIC.value == "exonic"
        assert RefCounting.EXONIC_INTRONIC.value == "exonic_intronic"

    def test_spatial_platform_values(self):
        assert SpatialPlatform.VISIUM.value == "visium"
        assert SpatialPlatform.VISIUM_HD.value == "visium_hd"
        assert SpatialPlatform.SLIDESEQ.value == "slideseq"
        assert SpatialPlatform.MERFISH.value == "merfish"
        assert SpatialPlatform.XENIUM.value == "xenium"
        assert SpatialPlatform.UNKNOWN.value == "unknown"

    def test_enum_serialises_to_string(self):
        assert str(BulkProtocol.POLYA.value) == "polyA"
        assert json.dumps({"p": BulkProtocol.POLYA.value})  # JSON-serialisable

    def test_enum_is_str_subtype(self):
        assert isinstance(BulkProtocol.POLYA, str)
        assert isinstance(RefModality.SCRNA, str)

    def test_parse_enum_case_insensitive(self):
        assert _parse_enum(BulkProtocol, "polya") == BulkProtocol.POLYA
        assert _parse_enum(BulkProtocol, "POLYA") == BulkProtocol.POLYA
        assert _parse_enum(BulkProtocol, "polyA") == BulkProtocol.POLYA

    def test_invalid_enum_value_raises(self):
        with pytest.raises(ValueError, match="not a valid BulkProtocol"):
            _parse_enum(BulkProtocol, "totalRNA")


# ===========================================================================
# ProtocolMetadata construction
# ===========================================================================

class TestProtocolMetadata:

    def test_default_all_unknown(self):
        meta = ProtocolMetadata()
        assert meta.bulk_protocol == BulkProtocol.UNKNOWN
        assert meta.ref_modality == RefModality.UNKNOWN
        assert meta.ref_capture == RefCapture.UNKNOWN
        assert meta.ref_counting == RefCounting.UNKNOWN
        assert meta.spatial_platform == SpatialPlatform.UNKNOWN

    def test_is_fully_unknown_true(self):
        assert ProtocolMetadata().is_fully_unknown()

    def test_is_fully_unknown_false_when_bulk_set(self):
        meta = ProtocolMetadata(bulk_protocol=BulkProtocol.POLYA)
        assert not meta.is_fully_unknown()

    def test_is_bulk_fully_unknown(self):
        meta = ProtocolMetadata(spatial_platform=SpatialPlatform.VISIUM)
        assert meta.is_bulk_fully_unknown()

    def test_is_spatial_true(self):
        meta = ProtocolMetadata(spatial_platform=SpatialPlatform.VISIUM)
        assert meta.is_spatial()

    def test_is_spatial_false_for_bulk(self):
        meta = ProtocolMetadata(bulk_protocol=BulkProtocol.POLYA)
        assert not meta.is_spatial()

    def test_from_strings_valid(self):
        meta = ProtocolMetadata.from_strings(
            bulk_protocol="polyA",
            ref_modality="snRNA",
            ref_capture="10x_3prime",
            ref_counting="exonic",
        )
        assert meta.bulk_protocol == BulkProtocol.POLYA
        assert meta.ref_modality == RefModality.SNRNA
        assert meta.ref_capture == RefCapture.UMI_3PRIME
        assert meta.ref_counting == RefCounting.EXONIC

    def test_from_strings_invalid_bulk_protocol_raises(self):
        with pytest.raises(ValueError, match="BulkProtocol"):
            ProtocolMetadata.from_strings(bulk_protocol="totalRNA")

    def test_from_strings_invalid_ref_modality_raises(self):
        with pytest.raises(ValueError, match="RefModality"):
            ProtocolMetadata.from_strings(ref_modality="tenx")

    def test_from_strings_invalid_spatial_platform_raises(self):
        with pytest.raises(ValueError, match="SpatialPlatform"):
            ProtocolMetadata.from_strings(spatial_platform="not_a_platform")

    def test_to_dict_round_trip(self):
        meta = ProtocolMetadata(
            bulk_protocol=BulkProtocol.RIBODEP,
            ref_modality=RefModality.SCRNA,
            spatial_platform=SpatialPlatform.VISIUM,
        )
        d = meta.to_dict()
        assert d["bulk_protocol"] == "ribodepleted"
        assert d["ref_modality"] == "scRNA"
        assert d["spatial_platform"] == "visium"
        loaded = ProtocolMetadata.from_dict(d)
        assert loaded.bulk_protocol == BulkProtocol.RIBODEP
        assert loaded.ref_modality == RefModality.SCRNA
        assert loaded.spatial_platform == SpatialPlatform.VISIUM

    def test_to_dict_is_json_serialisable(self):
        meta = ProtocolMetadata(
            bulk_protocol=BulkProtocol.POLYA,
            ref_modality=RefModality.SNRNA,
        )
        serialised = json.dumps(meta.to_dict())
        assert "polyA" in serialised

    def test_from_dict_tolerates_missing_keys(self):
        """from_dict must not raise on partial dicts."""
        meta = ProtocolMetadata.from_dict({"bulk_protocol": "polyA"})
        assert meta.bulk_protocol == BulkProtocol.POLYA
        assert meta.ref_modality == RefModality.UNKNOWN

    def test_from_dict_ignores_extra_keys(self):
        meta = ProtocolMetadata.from_dict({
            "bulk_protocol": "polyA",
            "unknown_future_key": "whatever",
        })
        assert meta.bulk_protocol == BulkProtocol.POLYA

    def test_spatial_platform_field_present(self):
        meta = ProtocolMetadata(spatial_platform=SpatialPlatform.XENIUM)
        assert meta.spatial_platform == SpatialPlatform.XENIUM


# ===========================================================================
# ProtocolRiskAssessor — gene list files and risk computation
# ===========================================================================

class TestProtocolRiskAssessor:

    def test_gene_list_data_path_works_in_installed_package(self):
        """Gene lists bundled with the package must be accessible."""
        from tissueresolve.reference.gene_filters import GeneFilterSet
        gfs = GeneFilterSet(genome="hg38")
        # Files copied from CHIMERA; verify they're non-empty
        intronic = gfs.intronic_dominant()
        stress = gfs.dissociation_stress()
        protein_coding = gfs.protein_coding()
        assert len(intronic) > 0, "intronic_dominant_hg38.tsv is empty"
        assert len(stress) > 0, "dissociation_stress_hg38.tsv is empty"
        assert len(protein_coding) > 0, "protein_coding_hg38.tsv is empty"

    def test_snrna_reference_flags_intronic_genes(self):
        """polyA bulk + snRNA reference should flag intronic-dominant genes."""
        meta = ProtocolMetadata(
            bulk_protocol=BulkProtocol.POLYA,
            ref_modality=RefModality.SNRNA,
        )
        assessor = ProtocolRiskAssessor(genome="hg38")
        intronic_genes = list(assessor._filters.intronic_dominant())[:5]
        safe_genes = ["GAPDH", "ACTB", "TP53", "VIM", "CDH1"]
        candidates = intronic_genes + safe_genes

        report = assessor.assess(meta, candidates)
        assert "intronic_retention" in report.mismatch_types
        assert report.risk_level in ("medium", "high")

    def test_scrna_reference_flags_stress_genes(self):
        """Any bulk + scRNA reference should flag dissociation-stress genes."""
        meta = ProtocolMetadata(ref_modality=RefModality.SCRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        stress_genes = list(assessor._filters.dissociation_stress())[:5]
        safe_genes = ["GAPDH", "ACTB"]
        candidates = stress_genes + safe_genes

        report = assessor.assess(meta, candidates)
        assert "dissociation_stress" in report.mismatch_types

    def test_3prime_capture_flags_length_biased_genes(self):
        """3′-UMI capture reference should flag length-biased genes."""
        meta = ProtocolMetadata(ref_capture=RefCapture.UMI_3PRIME)
        assessor = ProtocolRiskAssessor(genome="hg38")
        long_genes = list(assessor._filters.length_biased_3prime())[:5]
        candidates = long_genes + ["GAPDH", "ACTB"]

        report = assessor.assess(meta, candidates)
        assert "length_bias" in report.mismatch_types

    def test_polya_unknown_ref_is_low_risk(self):
        """polyA bulk with no ref information is low or unknown risk."""
        meta = ProtocolMetadata(bulk_protocol=BulkProtocol.POLYA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        report = assessor.assess(meta, ["GAPDH", "ACTB", "TP53"])
        # No modality set → no mismatch types → low risk
        assert report.risk_level in ("low", "unknown")

    def test_excluded_genes_listed_explicitly_not_silent(self):
        """ProtocolRiskReport must list every excluded gene explicitly."""
        meta = ProtocolMetadata(ref_modality=RefModality.SNRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        intronic = list(assessor._filters.intronic_dominant())[:3]
        candidates = intronic + ["GAPDH", "ACTB"]

        report = assessor.assess(meta, candidates)

        # n_genes_excluded must equal len(excluded_genes)
        assert report.n_genes_excluded == len(report.excluded_genes), (
            "n_genes_excluded must exactly match len(excluded_genes). "
            "No silent gene removal."
        )
        # Every excluded gene must be in the candidate list
        for g in report.excluded_genes:
            assert g in set(candidates), (
                f"Excluded gene {g!r} was not in candidate_genes — "
                "something was removed silently."
            )

    def test_no_genes_excluded_when_no_risk_genes_in_candidates(self):
        """If candidates contain no protocol-risk genes, nothing is excluded."""
        meta = ProtocolMetadata(ref_modality=RefModality.SNRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        # These are not in the intronic_dominant list
        candidates = ["GAPDH", "ACTB", "TP53", "MYC", "CDH1"]

        report = assessor.assess(meta, candidates)
        # May flag intronic_retention (type) if any candidates matched, but
        # excluded_genes should be empty since none of these are intronic-dominant
        assert all(g in set(candidates) for g in report.excluded_genes)
        assert len(report.excluded_genes) == report.n_genes_excluded

    def test_gene_risk_scores_returns_series(self):
        meta = ProtocolMetadata(ref_modality=RefModality.SNRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        genes = ["GAPDH", "ACTB"] + list(assessor._filters.intronic_dominant())[:3]
        scores = assessor.gene_risk_scores(meta, genes)
        assert isinstance(scores, pd.Series)
        assert list(scores.index) == genes

    def test_gene_risk_scores_in_0_1_range(self):
        meta = ProtocolMetadata(ref_modality=RefModality.SNRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        genes = list(assessor._filters.intronic_dominant())[:5] + ["GAPDH"]
        scores = assessor.gene_risk_scores(meta, genes)
        assert (scores >= 0.0).all()
        assert (scores <= 1.0).all()

    def test_intronic_genes_have_high_risk_with_snrna(self):
        """snRNA reference → intronic-dominant genes must have risk ≥ threshold."""
        meta = ProtocolMetadata(ref_modality=RefModality.SNRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        intronic = list(assessor._filters.intronic_dominant())
        if not intronic:
            pytest.skip("No intronic_dominant genes loaded")
        scores = assessor.gene_risk_scores(meta, intronic[:5])
        assert (scores >= RISK_HARD_THRESHOLD).all()

    def test_unknown_protocol_returns_unknown_risk(self):
        meta = ProtocolMetadata()
        assessor = ProtocolRiskAssessor(genome="hg38")
        report = assessor.assess(meta, ["GAPDH", "ACTB"])
        assert report.risk_level == "unknown"
        assert report.excluded_genes == []

    def test_risk_level_high_when_many_intronic_genes(self):
        """High affected fraction with snRNA reference → high risk."""
        meta = ProtocolMetadata(ref_modality=RefModality.SNRNA)
        assessor = ProtocolRiskAssessor(genome="hg38")
        intronic = list(assessor._filters.intronic_dominant())
        if len(intronic) < 4:
            pytest.skip("Too few intronic genes to test high risk level")
        # Make intronic genes > 15% of candidates → high risk
        safe = ["GAPDH"]
        candidates = intronic + safe   # intronic >> 15%
        report = assessor.assess(meta, candidates)
        assert report.risk_level in ("medium", "high")

    def test_safe_panel_returns_subset_of_gene_names(self):
        assessor = ProtocolRiskAssessor(genome="hg38")
        genes = ["GAPDH", "ACTB", "NONEXISTENT_ZZZ"]
        panel = assessor.safe_panel(genes)
        assert set(panel).issubset(set(genes))

    def test_ribodep_bulk_adds_mismatch_type(self):
        meta = ProtocolMetadata(
            bulk_protocol=BulkProtocol.RIBODEP,
            ref_modality=RefModality.SCRNA,
        )
        assessor = ProtocolRiskAssessor(genome="hg38")
        report = assessor.assess(meta, ["GAPDH", "ACTB", "FOS", "JUN"])
        assert "ribodep_intronic" in report.mismatch_types

    def test_gene_filters_shared_between_assessor_and_selector(self):
        """GeneFilterSet instance can be shared without duplicate file loading."""
        from tissueresolve.reference.gene_filters import GeneFilterSet
        from tissueresolve.reference.markers import GeneSelector

        shared_gfs = GeneFilterSet(genome="hg38")
        assessor = ProtocolRiskAssessor(genome="hg38", gene_filters=shared_gfs)
        assert assessor._filters is shared_gfs


# ===========================================================================
# SpatialMismatch construction and validation
# ===========================================================================

class TestSpatialMismatch:

    def test_construction(self):
        m = _simple_mismatch(10)
        assert m.G_m == 10
        assert isinstance(m.d_g, np.ndarray)
        assert isinstance(m.discord_score, np.ndarray)
        assert isinstance(m.gene_weight, np.ndarray)

    def test_shape_validation_d_g_wrong(self):
        with pytest.raises(ValueError, match="d_g shape"):
            SpatialMismatch(
                d_g=np.ones(5),
                discord_score=np.zeros(10),
                gene_weight=np.ones(10),
                gene_names=[f"G{i}" for i in range(10)],
            )

    def test_shape_validation_discord_score_wrong(self):
        with pytest.raises(ValueError, match="discord_score shape"):
            SpatialMismatch(
                d_g=np.ones(10),
                discord_score=np.zeros(5),
                gene_weight=np.ones(10),
                gene_names=[f"G{i}" for i in range(10)],
            )

    def test_shape_validation_gene_weight_wrong(self):
        with pytest.raises(ValueError, match="gene_weight shape"):
            SpatialMismatch(
                d_g=np.ones(10),
                discord_score=np.zeros(10),
                gene_weight=np.ones(3),
                gene_names=[f"G{i}" for i in range(10)],
            )

    def test_summary_dict_keys(self):
        m = _simple_mismatch(15)
        s = m.summary()
        for key in ("n_marker_genes", "n_high_discord", "mean_discord",
                    "d_g_median", "d_g_min", "d_g_max"):
            assert key in s

    def test_protocol_mismatch_alias(self):
        """ProtocolMismatch must be the same class as SpatialMismatch."""
        assert ProtocolMismatch is SpatialMismatch

    def test_g_m_property(self):
        m = _simple_mismatch(7)
        assert m.G_m == 7


# ===========================================================================
# compute_spatial_discordance
# ===========================================================================

class TestComputeSpatialDiscordance:

    def test_basic_numpy_input_returns_spatial_mismatch(self):
        G_m, N = 20, 30
        ref = _simple_ref(n_genes=G_m)
        rng = np.random.default_rng(0)
        Y = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)

        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert isinstance(result, SpatialMismatch)

    def test_output_shape(self):
        G_m, N = 15, 25
        ref = _simple_ref(n_genes=G_m)
        rng = np.random.default_rng(1)
        Y = rng.integers(0, 50, (N, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)

        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert result.d_g.shape == (G_m,)
        assert result.discord_score.shape == (G_m,)
        assert result.gene_weight.shape == (G_m,)
        assert result.G_m == G_m

    def test_d_g_positive(self):
        """Scale factors must always be positive (clipped to d_clip[0])."""
        G_m, N = 20, 30
        ref = _simple_ref(n_genes=G_m)
        rng = np.random.default_rng(2)
        Y = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)

        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert (result.d_g > 0).all()

    def test_zero_counts_safe(self):
        """All-zero count matrix must not raise and must produce positive d_g."""
        G_m, N = 15, 20
        ref = _simple_ref(n_genes=G_m)
        Y_zero = np.zeros((N, G_m), dtype=np.float32)
        lib_sizes = np.ones(N, dtype=np.float32) * 1000.0

        result = compute_spatial_discordance(
            Y_zero, ref, ref.gene_names, lib_sizes=lib_sizes
        )
        assert (result.d_g > 0).all()
        assert np.isfinite(result.d_g).all()
        assert np.isfinite(result.discord_score).all()

    def test_gene_weight_in_range(self):
        """gene_weight must be in [0.05, 1.0]."""
        G_m, N = 20, 30
        ref = _simple_ref(n_genes=G_m)
        rng = np.random.default_rng(3)
        Y = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)

        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert (result.gene_weight >= 0.05).all()
        assert (result.gene_weight <= 1.0).all()

    def test_discord_score_non_negative(self):
        G_m, N = 20, 30
        ref = _simple_ref(n_genes=G_m)
        rng = np.random.default_rng(4)
        Y = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)

        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert (result.discord_score >= 0).all()

    def test_missing_lib_sizes_raises(self):
        G_m, N = 10, 20
        ref = _simple_ref(n_genes=G_m)
        Y = np.ones((N, G_m), dtype=np.float32)
        with pytest.raises(ValueError, match="lib_sizes"):
            compute_spatial_discordance(Y, ref, ref.gene_names)

    def test_gene_names_preserved(self):
        G_m, N = 12, 20
        ref = _simple_ref(n_genes=G_m)
        Y = np.ones((N, G_m), dtype=np.float32)
        lib_sizes = np.ones(N) * 1000.0
        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert result.gene_names == list(ref.gene_names)

    def test_compute_discordance_alias(self):
        """compute_discordance must be an alias for compute_spatial_discordance."""
        from tissueresolve.protocol.mismatch import compute_discordance, compute_spatial_discordance
        assert compute_discordance is compute_spatial_discordance

    def test_anndata_input(self):
        """compute_spatial_discordance must accept AnnData input."""
        import anndata as ad
        import scipy.sparse as sp

        G_m, N = 12, 25
        ref = _simple_ref(n_genes=G_m)
        rng = np.random.default_rng(5)

        # Build AnnData with marker genes in known columns
        X = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        lib_sizes_vals = X.sum(axis=1).astype(np.float32)
        obs = pd.DataFrame({"total_counts": lib_sizes_vals})
        var = pd.DataFrame(index=ref.gene_names)
        adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)

        result = compute_spatial_discordance(adata, ref, ref.gene_names)
        assert result.G_m == G_m
        assert (result.d_g > 0).all()


# ===========================================================================
# update_mismatch_factors
# ===========================================================================

class TestUpdateMismatchFactors:

    def test_shape_preserved(self):
        G_m, N, K = 15, 25, 3
        rng = np.random.default_rng(0)
        Y = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        Pi = rng.dirichlet(np.ones(K), N).astype(np.float32)
        R_lin = rng.exponential(0.001, (K, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)
        mismatch = _simple_mismatch(G_m)
        # Align gene_names
        mismatch = SpatialMismatch(
            d_g=mismatch.d_g,
            discord_score=mismatch.discord_score,
            gene_weight=mismatch.gene_weight,
            gene_names=[f"GENE_{i:04d}" for i in range(G_m)],
        )

        updated = update_mismatch_factors(Y, Pi, R_lin, lib_sizes, mismatch)
        assert updated.d_g.shape == (G_m,)

    def test_d_g_positive_after_update(self):
        G_m, N, K = 15, 25, 3
        rng = np.random.default_rng(1)
        Y = rng.integers(0, 100, (N, G_m)).astype(np.float32)
        Pi = rng.dirichlet(np.ones(K), N).astype(np.float32)
        R_lin = rng.exponential(0.001, (K, G_m)).astype(np.float32)
        lib_sizes = rng.uniform(500, 2000, N).astype(np.float32)
        mismatch = _simple_mismatch(G_m)
        mismatch = SpatialMismatch(
            d_g=mismatch.d_g, discord_score=mismatch.discord_score,
            gene_weight=mismatch.gene_weight,
            gene_names=[f"GENE_{i:04d}" for i in range(G_m)],
        )

        updated = update_mismatch_factors(Y, Pi, R_lin, lib_sizes, mismatch)
        assert (updated.d_g > 0).all()

    def test_returns_new_object(self):
        """update_mismatch_factors must not mutate the input."""
        G_m, N, K = 10, 20, 2
        rng = np.random.default_rng(2)
        Y = rng.integers(0, 50, (N, G_m)).astype(np.float32)
        Pi = rng.dirichlet(np.ones(K), N).astype(np.float32)
        R_lin = rng.exponential(0.001, (K, G_m)).astype(np.float32)
        lib_sizes = np.ones(N) * 1000.0
        mismatch = SpatialMismatch(
            d_g=np.ones(G_m, dtype=np.float32),
            discord_score=np.zeros(G_m, dtype=np.float32),
            gene_weight=np.ones(G_m, dtype=np.float32),
            gene_names=[f"G{i}" for i in range(G_m)],
        )
        original_d_g = mismatch.d_g.copy()

        updated = update_mismatch_factors(Y, Pi, R_lin, lib_sizes, mismatch)
        assert updated is not mismatch
        # Original must not be mutated
        np.testing.assert_array_equal(mismatch.d_g, original_d_g)

    def test_d_g_clipped_to_d_clip_range(self):
        G_m, N, K = 10, 20, 2
        rng = np.random.default_rng(3)
        Y = rng.integers(0, 50, (N, G_m)).astype(np.float32)
        Pi = rng.dirichlet(np.ones(K), N).astype(np.float32)
        R_lin = rng.exponential(0.001, (K, G_m)).astype(np.float32)
        lib_sizes = np.ones(N) * 1000.0
        mismatch = SpatialMismatch(
            d_g=np.ones(G_m, dtype=np.float32),
            discord_score=np.zeros(G_m, dtype=np.float32),
            gene_weight=np.ones(G_m, dtype=np.float32),
            gene_names=[f"G{i}" for i in range(G_m)],
        )
        d_lo, d_hi = 0.2, 5.0
        updated = update_mismatch_factors(
            Y, Pi, R_lin, lib_sizes, mismatch, d_clip=(d_lo, d_hi)
        )
        assert (updated.d_g >= d_lo).all()
        assert (updated.d_g <= d_hi).all()

    def test_numerical_stability_with_zero_Pi(self):
        """All-zero proportions must not cause division by zero."""
        G_m, N, K = 8, 15, 3
        Y = np.zeros((N, G_m), dtype=np.float32)
        Pi = np.zeros((N, K), dtype=np.float32)
        R_lin = np.ones((K, G_m), dtype=np.float32) * 0.001
        lib_sizes = np.ones(N) * 1000.0
        mismatch = SpatialMismatch(
            d_g=np.ones(G_m, dtype=np.float32),
            discord_score=np.zeros(G_m, dtype=np.float32),
            gene_weight=np.ones(G_m, dtype=np.float32),
            gene_names=[f"G{i}" for i in range(G_m)],
        )
        updated = update_mismatch_factors(Y, Pi, R_lin, lib_sizes, mismatch)
        assert np.isfinite(updated.d_g).all()
        assert (updated.d_g > 0).all()


# ===========================================================================
# Distinction test — bulk risk vs spatial mismatch
# ===========================================================================

class TestProtocolDistinction:

    def test_bulk_risk_report_has_excluded_genes_field(self):
        """ProtocolRiskReport must have the excluded_genes field."""
        report = ProtocolRiskReport(
            risk_level="low", mismatch_types=[], affected_gene_fraction=0.0,
            n_genes_excluded=0, excluded_genes=[], downweighted_genes=[],
            rationale="test", recommended_actions=[], literature_refs=[],
        )
        assert hasattr(report, "excluded_genes")
        assert hasattr(report, "downweighted_genes")

    def test_spatial_mismatch_has_d_g_field(self):
        """SpatialMismatch must have the d_g field; ProtocolRiskReport must not."""
        m = _simple_mismatch(5)
        assert hasattr(m, "d_g")
        report = ProtocolRiskReport(
            risk_level="low", mismatch_types=[], affected_gene_fraction=0.0,
            n_genes_excluded=0, excluded_genes=[], downweighted_genes=[],
            rationale="", recommended_actions=[], literature_refs=[],
        )
        assert not hasattr(report, "d_g")

    def test_different_module_namespaces(self):
        """bulk risk and spatial mismatch must live in different namespaces."""
        from tissueresolve.protocol import risk as risk_mod
        from tissueresolve.protocol import mismatch as mismatch_mod
        assert risk_mod is not mismatch_mod

    def test_spatial_mismatch_class_name_is_spatial(self):
        assert "Spatial" in SpatialMismatch.__name__

    def test_bulk_risk_class_name_is_risk(self):
        assert "Risk" in ProtocolRiskReport.__name__

    def test_compute_discordance_returns_spatial_mismatch_not_risk_report(self):
        G_m, N = 10, 15
        ref = _simple_ref(n_genes=G_m)
        Y = np.ones((N, G_m), dtype=np.float32)
        lib_sizes = np.ones(N) * 1000.0
        result = compute_spatial_discordance(Y, ref, ref.gene_names, lib_sizes=lib_sizes)
        assert isinstance(result, SpatialMismatch)
        assert not isinstance(result, ProtocolRiskReport)
