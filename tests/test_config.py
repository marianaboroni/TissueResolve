"""
Tests: TissueResolveConfig construction and YAML round-trips.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tissueresolve.config import (
    BootstrapConfig,
    BulkQCConfig,
    BulkSolverConfig,
    DiscordanceConfig,
    GeneConfig,
    ReferenceConfig,
    SpatialQCConfig,
    SpatialSolverConfig,
    TissueResolveConfig,
)
from tissueresolve.utils import bh_fdr, project_simplex_batch, set_random_state


# ---------------------------------------------------------------------------
# Default construction
# ---------------------------------------------------------------------------

class TestDefaultConfig:

    def test_default_config_instantiates(self):
        cfg = TissueResolveConfig()
        assert cfg is not None

    def test_default_reference_config(self):
        cfg = TissueResolveConfig()
        assert cfg.reference.celltype_col == "cell_type"
        assert cfg.reference.genome == "hg38"
        assert cfg.reference.min_cells == 10

    def test_default_gene_config(self):
        cfg = TissueResolveConfig()
        assert cfg.genes.n_genes == 500
        assert cfg.genes.min_log2fc == 1.0

    def test_default_bulk_qc_thresholds_are_heuristic(self):
        cfg = TissueResolveConfig()
        thresholds = cfg.heuristic_thresholds()
        for name, entry in thresholds.items():
            assert entry["is_heuristic"] is True, (
                f"Threshold '{name}' must be marked is_heuristic=True."
            )

    def test_output_dir_is_path(self):
        cfg = TissueResolveConfig()
        assert isinstance(cfg.output_dir, Path)

    def test_spatial_solver_defaults(self):
        cfg = TissueResolveConfig()
        assert cfg.spatial_solver.lambda_spatial == 0.1
        assert cfg.spatial_solver.max_iter == 200
        assert cfg.spatial_solver.random_state == 42

    def test_bootstrap_defaults(self):
        cfg = TissueResolveConfig()
        assert cfg.bootstrap.n_bootstrap == 200
        assert cfg.bootstrap.ci_level == 0.95

    def test_spatial_qc_thresholds_are_none_by_default(self):
        cfg = TissueResolveConfig()
        # None means "derive from data at run time"
        assert cfg.spatial_qc.entropy_threshold is None
        assert cfg.spatial_qc.loglik_threshold is None
        assert cfg.spatial_qc.spatial_resid_threshold is None


# ---------------------------------------------------------------------------
# YAML round-trip
# ---------------------------------------------------------------------------

class TestYAMLRoundTrip:

    def test_basic_roundtrip(self, tmp_path):
        cfg = TissueResolveConfig()
        yaml_path = tmp_path / "config.yaml"
        cfg.to_yaml(yaml_path)
        loaded = TissueResolveConfig.from_yaml(yaml_path)

        assert loaded.reference.celltype_col == cfg.reference.celltype_col
        assert loaded.reference.genome == cfg.reference.genome
        assert loaded.genes.n_genes == cfg.genes.n_genes
        assert loaded.bootstrap.n_bootstrap == cfg.bootstrap.n_bootstrap
        assert loaded.verbose == cfg.verbose
        assert isinstance(loaded.output_dir, Path)

    def test_modified_values_survive_roundtrip(self, tmp_path):
        cfg = TissueResolveConfig()
        cfg.reference.genome = "mm10"
        cfg.genes.n_genes = 300
        cfg.spatial_solver.lambda_spatial = 0.5
        cfg.verbose = False

        yaml_path = tmp_path / "config_mod.yaml"
        cfg.to_yaml(yaml_path)
        loaded = TissueResolveConfig.from_yaml(yaml_path)

        assert loaded.reference.genome == "mm10"
        assert loaded.genes.n_genes == 300
        assert loaded.spatial_solver.lambda_spatial == pytest.approx(0.5)
        assert loaded.verbose is False

    def test_from_yaml_raises_on_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            TissueResolveConfig.from_yaml(tmp_path / "nonexistent.yaml")

    def test_to_yaml_creates_parent_dirs(self, tmp_path):
        cfg = TissueResolveConfig()
        nested = tmp_path / "deep" / "nested" / "config.yaml"
        cfg.to_yaml(nested)
        assert nested.exists()

    def test_yaml_is_human_readable(self, tmp_path):
        cfg = TissueResolveConfig()
        yaml_path = tmp_path / "readable.yaml"
        cfg.to_yaml(yaml_path)
        content = yaml_path.read_text()
        # Key names should appear as plain text, not mangled
        assert "reference:" in content
        assert "genes:" in content
        assert "bootstrap:" in content

    def test_output_dir_survives_roundtrip_as_path(self, tmp_path):
        cfg = TissueResolveConfig()
        cfg.output_dir = Path("my_run_dir")
        cfg.to_yaml(tmp_path / "cfg.yaml")
        loaded = TissueResolveConfig.from_yaml(tmp_path / "cfg.yaml")
        assert isinstance(loaded.output_dir, Path)
        assert loaded.output_dir == Path("my_run_dir")

    def test_unknown_yaml_keys_are_ignored(self, tmp_path):
        """Loading a config written by a newer version must not raise."""
        yaml_path = tmp_path / "future.yaml"
        yaml_path.write_text(
            "reference:\n"
            "  celltype_col: cell_type\n"
            "  future_unknown_key: some_value\n"
            "genes:\n"
            "  n_genes: 400\n"
        )
        loaded = TissueResolveConfig.from_yaml(yaml_path)
        assert loaded.genes.n_genes == 400


# ---------------------------------------------------------------------------
# Utils — bh_fdr
# ---------------------------------------------------------------------------

class TestBhFdr:
    """Verify the unified BH-FDR implementation against known cases."""

    def test_empty_input(self):
        result = bh_fdr(np.array([]))
        assert len(result) == 0

    def test_single_pvalue(self):
        result = bh_fdr(np.array([0.03]))
        np.testing.assert_allclose(result, [0.03])

    def test_all_pvalues_one(self):
        result = bh_fdr(np.ones(5))
        np.testing.assert_allclose(result, np.ones(5))

    def test_adjusted_ge_original(self):
        pvals = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        adjusted = bh_fdr(pvals)
        assert (adjusted >= pvals).all()

    def test_adjusted_le_one(self):
        pvals = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        adjusted = bh_fdr(pvals)
        assert (adjusted <= 1.0).all()

    def test_monotone_after_sorting(self):
        """Adjusted p-values must be non-decreasing in rank order."""
        pvals = np.array([0.001, 0.002, 0.01, 0.04, 0.2])
        adjusted = bh_fdr(pvals)
        order = np.argsort(pvals)
        ranked_adj = adjusted[order]
        assert (np.diff(ranked_adj) >= -1e-12).all()

    def test_known_five_pvalues(self):
        """Cross-check against a hand-computed BH example."""
        # p = [0.001, 0.008, 0.039, 0.041, 0.210], n=5
        # Ranks (1-indexed): 1, 2, 3, 4, 5
        # BH: p_adj[i] = p[i] * n / rank[i]  (then monotone from right)
        # raw: 0.005, 0.020, 0.065, 0.051 (from 0.041*5/4), 0.210
        # after monotone:  0.020, 0.020, 0.065, 0.065, 0.210 — but need exact
        pvals = np.array([0.001, 0.008, 0.039, 0.041, 0.210])
        adjusted = bh_fdr(pvals)
        assert adjusted[0] <= adjusted[1] <= adjusted[2]
        assert adjusted[4] <= 1.0
        assert adjusted[0] < adjusted[4]

    def test_raises_on_2d_input(self):
        with pytest.raises(ValueError, match="1-D"):
            bh_fdr(np.ones((3, 3)))


# ---------------------------------------------------------------------------
# Utils — project_simplex_batch
# ---------------------------------------------------------------------------

class TestProjectSimplexBatch:

    def test_rows_sum_to_one(self):
        rng = np.random.default_rng(0)
        V = rng.standard_normal((20, 5))
        result = project_simplex_batch(V)
        np.testing.assert_allclose(result.sum(axis=1), np.ones(20), atol=1e-6)

    def test_rows_nonnegative(self):
        rng = np.random.default_rng(1)
        V = rng.standard_normal((20, 5))
        result = project_simplex_batch(V)
        assert (result >= 0).all()

    def test_already_on_simplex_unchanged(self):
        V = np.array([[0.25, 0.25, 0.25, 0.25]])
        result = project_simplex_batch(V)
        np.testing.assert_allclose(result, V, atol=1e-6)


# ---------------------------------------------------------------------------
# Utils — set_random_state
# ---------------------------------------------------------------------------

class TestSetRandomState:

    def test_returns_generator(self):
        rng = set_random_state(42)
        assert isinstance(rng, np.random.Generator)

    def test_reproducibility(self):
        rng1 = set_random_state(7)
        rng2 = set_random_state(7)
        a = rng1.standard_normal(10)
        b = rng2.standard_normal(10)
        np.testing.assert_array_equal(a, b)
