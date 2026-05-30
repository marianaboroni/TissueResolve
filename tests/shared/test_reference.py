"""
Tests: ReferenceSignature (orientation + summary) and ReferenceBuilder.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.reference.build import ReferenceBuilder
from tissueresolve.results import ReferenceSignature
from tissueresolve.config import ReferenceConfig


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_GENES = 30
N_CELLS = 60
N_TYPES = 3
GENE_NAMES = [f"GENE_{i:03d}" for i in range(N_GENES)]
CELL_TYPES_LIST = ["TypeA", "TypeB", "TypeC"]
CELL_TYPES_SORTED = sorted(CELL_TYPES_LIST)


def _make_counts(rng: np.random.Generator) -> pd.DataFrame:
    """Synthetic genes × cells count DataFrame."""
    X = rng.integers(0, 50, size=(N_GENES, N_CELLS)).astype(float)
    cells = [f"cell_{i}" for i in range(N_CELLS)]
    return pd.DataFrame(X, index=GENE_NAMES, columns=cells)


def _make_metadata(
    counts: pd.DataFrame,
    *,
    n_donors: int = 0,
) -> pd.DataFrame:
    """Cell metadata with optional donor column.

    Donor assignment interleaves donors across cell types so that every cell
    type has ``n_donors`` unique donors (required for donor_cv computation).
    """
    n = counts.shape[1]
    types = np.tile(CELL_TYPES_LIST, n // len(CELL_TYPES_LIST) + 1)[:n]
    meta = pd.DataFrame({"cell_type": types}, index=counts.columns)
    if n_donors > 0:
        # Assign donor by block-of-n_types so each type gets all donors.
        # E.g. cells 0,1,2 = (TypeA/B/C, donor_0),
        #      cells 3,4,5 = (TypeA/B/C, donor_1), ...
        K = len(CELL_TYPES_LIST)
        donors = [f"donor_{(i // K) % n_donors}" for i in range(n)]
        meta["donor"] = donors
    return meta


# ---------------------------------------------------------------------------
# ReferenceSignature orientation and shape
# ---------------------------------------------------------------------------

class TestReferenceSignatureOrientation:

    def test_phi_is_G_by_K(self):
        rng = np.random.default_rng(0)
        G, K = 20, 4
        raw = rng.exponential(1.0, (G, K))
        raw /= raw.sum(axis=0, keepdims=True)
        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=[f"CT{k}" for k in range(K)],
            phi=raw,
        )
        assert ref.phi.shape == (G, K), f"phi must be (G={G}, K={K})"

    def test_R_cpm_is_K_by_G(self):
        rng = np.random.default_rng(1)
        G, K = 20, 4
        R = rng.exponential(100.0, (K, G)).astype(np.float32)
        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=[f"CT{k}" for k in range(K)],
            R_cpm=R,
        )
        assert ref.R_cpm.shape == (K, G), f"R_cpm must be (K={K}, G={G})"

    def test_phi_columns_sum_to_one(self):
        rng = np.random.default_rng(2)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        builder = ReferenceBuilder()
        ref = builder.build_from_df(counts, meta)
        col_sums = ref.phi.sum(axis=0)
        np.testing.assert_allclose(col_sums, np.ones(ref.n_cell_types), atol=1e-6)

    def test_R_cpm_rows_are_K_first(self):
        rng = np.random.default_rng(3)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        builder = ReferenceBuilder()
        ref = builder.build_from_df(counts, meta)
        assert ref.R_cpm.shape == (ref.n_cell_types, ref.n_genes)

    def test_as_phi_from_R_cpm_is_G_by_K(self):
        rng = np.random.default_rng(4)
        G, K = 15, 3
        R = rng.exponential(100.0, (K, G)).astype(np.float32)
        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=[f"CT{k}" for k in range(K)],
            R_cpm=R,
        )
        phi = ref.as_phi()
        assert phi.shape == (G, K)
        np.testing.assert_allclose(phi.sum(axis=0), np.ones(K), atol=1e-5)

    def test_as_R_cpm_from_phi_is_K_by_G(self):
        rng = np.random.default_rng(5)
        G, K = 15, 3
        raw = rng.exponential(1.0, (G, K))
        raw /= raw.sum(axis=0, keepdims=True)
        ref = ReferenceSignature(
            gene_names=[f"G{i}" for i in range(G)],
            cell_types=[f"CT{k}" for k in range(K)],
            phi=raw,
        )
        R = ref.as_R_cpm()
        assert R.shape == (K, G)

    def test_as_R_log_shape(self):
        rng = np.random.default_rng(6)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        rl = ref.as_R_log()
        assert rl.shape == (ref.n_cell_types, ref.n_genes)

    def test_phi_g_shape(self):
        rng = np.random.default_rng(7)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta, estimate_overdispersion=True)
        assert ref.phi_g is not None
        assert ref.phi_g.shape == (ref.n_genes,)


# ---------------------------------------------------------------------------
# ReferenceSignature.summary()
# ---------------------------------------------------------------------------

class TestReferenceSignatureSummary:

    def test_summary_shape(self):
        rng = np.random.default_rng(10)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        s = ref.summary()
        assert s.shape == (ref.n_cell_types, 4)
        assert list(s.columns) == ["n_cells", "mean_cpm", "max_cpm", "profile_entropy"]

    def test_summary_index_is_cell_types(self):
        rng = np.random.default_rng(11)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        s = ref.summary()
        assert list(s.index) == ref.cell_types

    def test_profile_entropy_positive(self):
        rng = np.random.default_rng(12)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        s = ref.summary()
        assert (s["profile_entropy"] > 0).all()

    def test_n_cells_matches_metadata(self):
        rng = np.random.default_rng(13)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        s = ref.summary()
        for ct in ref.cell_types:
            expected = (meta["cell_type"] == ct).sum()
            assert s.loc[ct, "n_cells"] == expected


# ---------------------------------------------------------------------------
# ReferenceBuilder: basic construction
# ---------------------------------------------------------------------------

class TestReferenceBuilderBasic:

    def test_build_from_df_returns_reference_signature(self):
        rng = np.random.default_rng(20)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert isinstance(ref, ReferenceSignature)

    def test_cell_types_sorted_alphabetically(self):
        rng = np.random.default_rng(21)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.cell_types == sorted(ref.cell_types)

    def test_gene_names_preserved(self):
        rng = np.random.default_rng(22)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.gene_names == list(counts.index)

    def test_n_cell_types(self):
        rng = np.random.default_rng(23)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.n_cell_types == N_TYPES

    def test_n_genes(self):
        rng = np.random.default_rng(24)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.n_genes == N_GENES

    def test_metadata_as_series(self):
        """build_from_df should accept a Series as metadata."""
        rng = np.random.default_rng(25)
        counts = _make_counts(rng)
        n = counts.shape[1]
        types = np.tile(CELL_TYPES_LIST, n // len(CELL_TYPES_LIST) + 1)[:n]
        meta_series = pd.Series(types, index=counts.columns, name="cell_type")
        ref = ReferenceBuilder().build_from_df(counts, meta_series)
        assert ref.n_cell_types == N_TYPES

    def test_n_cells_per_type_populated(self):
        rng = np.random.default_rng(26)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert len(ref.n_cells_per_type) == N_TYPES
        total = sum(ref.n_cells_per_type.values())
        assert total == N_CELLS

    def test_genome_passed_through(self):
        rng = np.random.default_rng(27)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        cfg = ReferenceConfig(genome="mm10")
        ref = ReferenceBuilder(cfg).build_from_df(counts, meta)
        assert ref.genome == "mm10"

    def test_validate_passes_after_build(self):
        rng = np.random.default_rng(28)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        ref.validate()  # should not raise


# ---------------------------------------------------------------------------
# ReferenceBuilder: input validation errors
# ---------------------------------------------------------------------------

class TestReferenceBuilderValidation:

    def test_duplicate_genes_raises(self):
        rng = np.random.default_rng(30)
        counts = _make_counts(rng)
        # Introduce duplicate gene names
        counts.index = list(counts.index[:-1]) + [counts.index[0]]
        meta = _make_metadata(counts)
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            ReferenceBuilder().build_from_df(counts, meta)

    def test_missing_celltype_col_raises(self):
        rng = np.random.default_rng(31)
        counts = _make_counts(rng)
        meta = pd.DataFrame({"other_col": np.zeros(N_CELLS)}, index=counts.columns)
        with pytest.raises(KeyError, match="cell_type"):
            ReferenceBuilder().build_from_df(counts, meta)

    def test_no_shared_cells_raises(self):
        rng = np.random.default_rng(32)
        counts = _make_counts(rng)
        meta = pd.DataFrame(
            {"cell_type": ["TypeA"] * N_CELLS},
            index=[f"other_{i}" for i in range(N_CELLS)],
        )
        with pytest.raises(ValueError, match="[Nn]o shared cell"):
            ReferenceBuilder().build_from_df(counts, meta)

    def test_all_types_below_min_cells_raises(self):
        rng = np.random.default_rng(33)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        # Require more cells than any type has
        cfg = ReferenceConfig(min_cells=N_CELLS + 1)
        with pytest.raises(ValueError, match="No cell types survived"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ReferenceBuilder(cfg).build_from_df(counts, meta)

    def test_rare_type_dropped_with_warning(self):
        rng = np.random.default_rng(34)
        counts = _make_counts(rng)
        # Make TypeC very rare (1 cell)
        meta = _make_metadata(counts)
        # Override: set most TypeC cells to TypeA
        n = counts.shape[1]
        types = list(meta["cell_type"])
        type_c_indices = [i for i, t in enumerate(types) if t == "TypeC"]
        for idx in type_c_indices[1:]:
            types[idx] = "TypeA"
        meta["cell_type"] = types

        cfg = ReferenceConfig(min_cells=2)
        with pytest.warns(UserWarning, match="TypeC"):
            ref = ReferenceBuilder(cfg).build_from_df(counts, meta)
        # TypeC was dropped — only 2 types remain
        assert "TypeC" not in ref.cell_types

    def test_negative_counts_raises(self):
        rng = np.random.default_rng(35)
        counts = _make_counts(rng)
        counts.iloc[0, 0] = -1.0
        meta = _make_metadata(counts)
        with pytest.raises(ValueError, match="negative"):
            ReferenceBuilder().build_from_df(counts, meta)


# ---------------------------------------------------------------------------
# ReferenceBuilder: donor-aware aggregation
# ---------------------------------------------------------------------------

class TestDonorAwareAggregation:

    def test_donor_cv_present_with_two_donors(self):
        rng = np.random.default_rng(40)
        counts = _make_counts(rng)
        meta = _make_metadata(counts, n_donors=2)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.donor_cv is not None
        assert ref.donor_cv.shape == (N_GENES, N_TYPES)

    def test_donor_cv_none_with_one_donor(self):
        rng = np.random.default_rng(41)
        counts = _make_counts(rng)
        meta = _make_metadata(counts, n_donors=1)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.donor_cv is None

    def test_donor_cv_none_without_donor_col(self):
        rng = np.random.default_rng(42)
        counts = _make_counts(rng)
        meta = _make_metadata(counts, n_donors=0)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.donor_cv is None

    def test_donor_aware_differs_from_simple_mean(self):
        """With imbalanced donors, donor-aware aggregation differs from simple mean."""
        rng = np.random.default_rng(43)
        G, C = 10, 40
        X_raw = rng.integers(0, 100, (G, C)).astype(float)
        # Create 2 donors with very different expression for TypeA
        # Donor 0: 30 TypeA cells;  Donor 1: 2 TypeA cells (much smaller group)
        cell_types = ["TypeA"] * 32 + ["TypeB"] * 8
        donors = ["D0"] * 30 + ["D1"] * 2 + ["D1"] * 8
        # Boost D0's TypeA cells on gene 0
        X_raw[0, :30] = 1000.0
        X_raw[0, 30:32] = 1.0

        genes = [f"G{i}" for i in range(G)]
        cells = [f"c{i}" for i in range(C)]
        counts = pd.DataFrame(X_raw, index=genes, columns=cells)
        meta = pd.DataFrame(
            {"cell_type": cell_types, "donor": donors}, index=cells
        )

        cfg = ReferenceConfig(min_cells=2)
        ref_donor = ReferenceBuilder(cfg).build_from_df(counts, meta)
        ref_simple = ReferenceBuilder(
            ReferenceConfig(donor_col=None, min_cells=2)
        ).build_from_df(counts, meta)

        # TypeA, gene 0: donor-aware averages D0 and D1 equally → lower than
        # the simple mean which is dominated by the 30 D0 cells.
        type_a_idx = ref_donor.cell_types.index("TypeA")
        phi_donor = ref_donor.phi[0, type_a_idx]
        phi_simple = ref_simple.phi[0, type_a_idx]
        # They should differ (not necessarily by direction, just differ)
        assert phi_donor != pytest.approx(phi_simple, rel=0.01)

    def test_donor_cv_non_negative(self):
        rng = np.random.default_rng(44)
        counts = _make_counts(rng)
        meta = _make_metadata(counts, n_donors=3)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        assert ref.donor_cv is not None
        assert (ref.donor_cv >= 0).all()


# ---------------------------------------------------------------------------
# ReferenceBuilder: NB overdispersion
# ---------------------------------------------------------------------------

class TestOverdispersionEstimation:

    def test_phi_g_shape_after_estimation(self):
        rng = np.random.default_rng(50)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta, estimate_overdispersion=True)
        assert ref.phi_g is not None
        assert ref.phi_g.shape == (N_GENES,)

    def test_phi_g_in_valid_range(self):
        rng = np.random.default_rng(51)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta, estimate_overdispersion=True)
        assert (ref.phi_g >= 0.5).all()
        assert (ref.phi_g <= 100.0).all()

    def test_phi_g_none_when_not_requested(self):
        rng = np.random.default_rng(52)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta, estimate_overdispersion=False)
        assert ref.phi_g is None

    def test_phi_g_dtype_float32(self):
        rng = np.random.default_rng(53)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta, estimate_overdispersion=True)
        assert ref.phi_g.dtype == np.float32


# ---------------------------------------------------------------------------
# ReferenceBuilder: build_from_adata
# ---------------------------------------------------------------------------

class TestBuildFromAdata:

    def test_build_from_adata_basic(self):
        """build_from_adata should produce same shapes as build_from_df."""
        import anndata as ad
        import scipy.sparse as sp

        rng = np.random.default_rng(60)
        X = rng.integers(0, 50, (N_CELLS, N_GENES)).astype(np.float32)
        n = N_CELLS
        types = np.tile(CELL_TYPES_LIST, n // len(CELL_TYPES_LIST) + 1)[:n]
        obs = pd.DataFrame({"cell_type": types})
        var = pd.DataFrame(index=GENE_NAMES)
        adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)

        ref = ReferenceBuilder().build_from_adata(adata)
        assert ref.n_genes == N_GENES
        assert ref.n_cell_types == N_TYPES
        assert ref.cell_types == sorted(CELL_TYPES_LIST)

    def test_build_from_adata_phi_L1_normalised(self):
        import anndata as ad
        import scipy.sparse as sp

        rng = np.random.default_rng(61)
        X = rng.integers(0, 50, (N_CELLS, N_GENES)).astype(np.float32)
        n = N_CELLS
        types = np.tile(CELL_TYPES_LIST, n // len(CELL_TYPES_LIST) + 1)[:n]
        obs = pd.DataFrame({"cell_type": types})
        var = pd.DataFrame(index=GENE_NAMES)
        adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)

        ref = ReferenceBuilder().build_from_adata(adata)
        np.testing.assert_allclose(
            ref.phi.sum(axis=0), np.ones(ref.n_cell_types), atol=1e-6
        )

    def test_build_from_adata_missing_celltype_col_raises(self):
        import anndata as ad

        rng = np.random.default_rng(62)
        X = rng.integers(0, 50, (N_CELLS, N_GENES)).astype(np.float32)
        obs = pd.DataFrame({"wrong_col": np.zeros(N_CELLS)})
        var = pd.DataFrame(index=GENE_NAMES)
        adata = ad.AnnData(X=X, obs=obs, var=var)
        with pytest.raises(KeyError, match="cell_type"):
            ReferenceBuilder().build_from_adata(adata)


# ---------------------------------------------------------------------------
# ReferenceBuilder: save/load round-trip after build
# ---------------------------------------------------------------------------

class TestBuiltReferenceSaveLoad:

    def test_save_load_preserves_phi(self, tmp_path):
        rng = np.random.default_rng(70)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        save_dir = tmp_path / "built_ref"
        ref.save(save_dir)
        loaded = ReferenceSignature.load(save_dir)
        np.testing.assert_array_almost_equal(loaded.phi, ref.phi)

    def test_save_load_preserves_R_cpm(self, tmp_path):
        rng = np.random.default_rng(71)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        save_dir = tmp_path / "built_ref2"
        ref.save(save_dir)
        loaded = ReferenceSignature.load(save_dir)
        np.testing.assert_array_almost_equal(loaded.R_cpm, ref.R_cpm, decimal=4)

    def test_save_load_with_overdispersion(self, tmp_path):
        rng = np.random.default_rng(72)
        counts = _make_counts(rng)
        meta = _make_metadata(counts)
        ref = ReferenceBuilder().build_from_df(counts, meta, estimate_overdispersion=True)
        save_dir = tmp_path / "ref_od"
        ref.save(save_dir)
        loaded = ReferenceSignature.load(save_dir)
        assert loaded.phi_g is not None
        np.testing.assert_array_equal(loaded.phi_g, ref.phi_g)

    def test_save_load_with_donor_cv(self, tmp_path):
        rng = np.random.default_rng(73)
        counts = _make_counts(rng)
        meta = _make_metadata(counts, n_donors=3)
        ref = ReferenceBuilder().build_from_df(counts, meta)
        save_dir = tmp_path / "ref_cv"
        ref.save(save_dir)
        loaded = ReferenceSignature.load(save_dir)
        assert loaded.donor_cv is not None
        np.testing.assert_array_almost_equal(loaded.donor_cv, ref.donor_cv)
