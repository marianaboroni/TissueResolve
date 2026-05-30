"""
Tests for tissueresolve.uncertainty.bootstrap (BulkBootstrapCI).
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.config import BootstrapConfig
from tissueresolve.results import BulkDeconvResult, ReferenceSignature
from tissueresolve.uncertainty.bootstrap import BulkBootstrapCI, MIN_PANEL_GENES


# ---------------------------------------------------------------------------
# Helpers (local — do not rely on conftest for novel bulk fixtures)
# ---------------------------------------------------------------------------


def _make_ref(n_genes: int = 80, n_types: int = 3, seed: int = 5) -> ReferenceSignature:
    rng = np.random.default_rng(seed)
    phi_raw = np.zeros((n_genes, n_types))
    block = n_genes // n_types
    for k in range(n_types):
        phi_raw[k * block:(k + 1) * block, k] = rng.uniform(1.5, 3.0, block)
    phi_raw += rng.uniform(0.0, 0.05, phi_raw.shape)
    col_sums = phi_raw.sum(axis=0, keepdims=True)
    phi = (phi_raw / col_sums).astype(np.float64)
    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_types)]
    return ReferenceSignature(gene_names=gene_names, cell_types=cell_types, phi=phi)


def _make_bulk(
    ref: ReferenceSignature,
    n_samples: int = 6,
    seed: int = 11,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    K = ref.n_cell_types
    props = rng.dirichlet(np.ones(K), size=n_samples)
    counts = (ref.as_phi() @ props.T) * 1e4
    return pd.DataFrame(counts, index=ref.gene_names,
                        columns=[f"S{i}" for i in range(n_samples)])


# ---------------------------------------------------------------------------
# Shape and type checks
# ---------------------------------------------------------------------------


class TestBulkBootstrapCIShapes:
    def test_returns_two_dataframes(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=10, bootstrap_frac=0.8, ci_level=0.95, seed=42)
        bs = BulkBootstrapCI(config=cfg)
        lo, hi = bs.compute(bulk, ref, ref.gene_names)
        assert isinstance(lo, pd.DataFrame)
        assert isinstance(hi, pd.DataFrame)

    def test_lower_shape_matches_bulk(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=10, seed=42)
        lo, hi = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert lo.shape == (bulk.shape[1], ref.n_cell_types)

    def test_upper_shape_matches_lower(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=10, seed=42)
        lo, hi = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert lo.shape == hi.shape

    def test_sample_index_matches_bulk_columns(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=10, seed=42)
        lo, hi = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert list(lo.index) == list(bulk.columns)
        assert list(hi.index) == list(bulk.columns)

    def test_columns_match_cell_types(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=10, seed=42)
        lo, hi = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert list(lo.columns) == list(ref.cell_types)

    def test_ci_lower_le_upper(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=20, seed=42)
        lo, hi = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert (hi.values >= lo.values - 1e-12).all()

    def test_ci_values_in_zero_one(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=20, seed=42)
        lo, hi = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert (lo.values >= -1e-9).all()
        assert (hi.values <= 1.0 + 1e-9).all()


# ---------------------------------------------------------------------------
# Warning for too few genes
# ---------------------------------------------------------------------------


class TestBulkBootstrapCIWarnings:
    def test_warns_on_too_few_genes(self):
        ref = _make_ref(n_genes=80)
        bulk = _make_bulk(ref)
        tiny_panel = ref.gene_names[:MIN_PANEL_GENES - 1]  # below threshold
        cfg = BootstrapConfig(n_bootstrap=5, seed=42)
        bs = BulkBootstrapCI(cfg)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            bs.compute(bulk, ref, tiny_panel)
        texts = [str(w.message) for w in caught]
        assert any("only" in t and "genes" in t for t in texts), (
            "Expected warning about too few genes, got: " + str(texts)
        )

    def test_no_warning_for_sufficient_genes(self):
        ref = _make_ref(n_genes=80)
        bulk = _make_bulk(ref)
        panel = ref.gene_names  # all 80 genes
        cfg = BootstrapConfig(n_bootstrap=5, seed=42)
        bs = BulkBootstrapCI(cfg)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            bs.compute(bulk, ref, panel)
        texts = [str(w.message) for w in caught if "genes" in str(w.message).lower()]
        # No gene-count warnings expected
        panel_size_warnings = [t for t in texts if "too few" in t or "only" in t]
        assert len(panel_size_warnings) == 0

    def test_no_overlap_raises(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=5, seed=42)
        with pytest.raises(ValueError, match="No gene overlap"):
            BulkBootstrapCI(cfg).compute(bulk, ref, ["NONEXISTENT"])


# ---------------------------------------------------------------------------
# attach() method
# ---------------------------------------------------------------------------


class TestBulkBootstrapCIAttach:
    def test_attach_populates_result_ci(self):
        from tissueresolve.bulk.solver import WNNLSSolver
        ref = _make_ref()
        bulk = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        assert result.lower_ci is None
        assert result.upper_ci is None

        cfg = BootstrapConfig(n_bootstrap=10, seed=42)
        BulkBootstrapCI(cfg).attach(result, bulk, ref)

        assert result.lower_ci is not None
        assert result.upper_ci is not None

    def test_attach_ci_shape_correct(self):
        from tissueresolve.bulk.solver import WNNLSSolver
        ref = _make_ref()
        bulk = _make_bulk(ref)
        result = WNNLSSolver().solve(bulk, ref, ref.gene_names)
        cfg = BootstrapConfig(n_bootstrap=10, seed=42)
        BulkBootstrapCI(cfg).attach(result, bulk, ref)

        assert result.lower_ci.shape == result.proportions.shape
        assert result.upper_ci.shape == result.proportions.shape


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


class TestBulkBootstrapCIMetadata:
    def test_ci_metadata_contains_n_bootstrap(self):
        cfg = BootstrapConfig(n_bootstrap=50, seed=0)
        bs = BulkBootstrapCI(cfg)
        meta = bs.ci_metadata
        assert meta["n_bootstrap"] == 50

    def test_ci_metadata_contains_ci_level(self):
        cfg = BootstrapConfig(ci_level=0.90, seed=0)
        bs = BulkBootstrapCI(cfg)
        assert bs.ci_metadata["ci_level"] == 0.90

    def test_ci_metadata_serialisable(self):
        import json
        cfg = BootstrapConfig(n_bootstrap=20, seed=1)
        meta = BulkBootstrapCI(cfg).ci_metadata
        json.dumps(meta)  # should not raise


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------


class TestBulkBootstrapCIReproducibility:
    def test_same_seed_gives_same_ci(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        cfg = BootstrapConfig(n_bootstrap=30, seed=99)
        lo1, hi1 = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        lo2, hi2 = BulkBootstrapCI(cfg).compute(bulk, ref, ref.gene_names)
        assert np.allclose(lo1.values, lo2.values)
        assert np.allclose(hi1.values, hi2.values)

    def test_different_seeds_differ(self):
        ref = _make_ref()
        bulk = _make_bulk(ref)
        lo1, _ = BulkBootstrapCI(BootstrapConfig(n_bootstrap=30, seed=1)).compute(bulk, ref, ref.gene_names)
        lo2, _ = BulkBootstrapCI(BootstrapConfig(n_bootstrap=30, seed=2)).compute(bulk, ref, ref.gene_names)
        assert not np.allclose(lo1.values, lo2.values)
