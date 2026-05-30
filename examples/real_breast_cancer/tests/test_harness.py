"""
Offline unit tests for the harness logic in ``_harness.py``.

No network, no downloads, tiny synthetic data only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Folders & flags
# ---------------------------------------------------------------------------


def test_ensure_dirs_creates_all(harness):
    harness.ensure_dirs()
    for d in harness.ALL_DIRS:
        assert d.exists() and d.is_dir()


def test_real_data_disabled_by_default(harness, monkeypatch):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    assert harness.real_data_enabled(False) is False


def test_real_data_enabled_by_flag(harness, monkeypatch):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    assert harness.real_data_enabled(True) is True


def test_real_data_enabled_by_env(harness, monkeypatch):
    monkeypatch.setenv(harness.REAL_DATA_ENV, "1")
    assert harness.real_data_enabled(False) is True


# ---------------------------------------------------------------------------
# Manifest schema
# ---------------------------------------------------------------------------


def test_default_manifest_is_valid(harness):
    m = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    harness.validate_manifest(m)  # no raise
    assert set(m["datasets"]) == {"reference", "spatial"}


def test_manifest_missing_key_raises(harness):
    m = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    del m["datasets"]["reference"]["access_date"]
    with pytest.raises(ValueError, match="missing keys"):
        harness.validate_manifest(m)


def test_manifest_roundtrip(harness, tmp_path):
    m = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    p = harness.write_manifest(m, tmp_path / "manifest.json")
    assert p.exists()
    harness.validate_manifest(harness.read_manifest(p))


# ---------------------------------------------------------------------------
# Cell-type column detection
# ---------------------------------------------------------------------------


def test_detect_cell_type_col_default(harness):
    obs = pd.DataFrame({"cell_type": ["A", "B"]})
    assert harness.detect_cell_type_col(obs) == "cell_type"


def test_detect_cell_type_col_priority(harness):
    obs = pd.DataFrame({"cell_type": ["A"], "cell_type_major": ["X"]})
    assert harness.detect_cell_type_col(obs) == "cell_type_major"


def test_detect_cell_type_col_override(harness):
    obs = pd.DataFrame({"cell_type": ["A"], "my_col": ["Z"]})
    assert harness.detect_cell_type_col(obs, "my_col") == "my_col"


def test_detect_cell_type_col_missing_raises(harness):
    obs = pd.DataFrame({"nope": ["A"]})
    with pytest.raises(KeyError):
        harness.detect_cell_type_col(obs)


# ---------------------------------------------------------------------------
# Reference preparation
# ---------------------------------------------------------------------------


def test_prepare_reference_builds_signature(harness, tiny_sc_adata):
    prep = harness.prepare_reference(tiny_sc_adata, min_cells=50)
    assert prep.cell_type_col == "cell_type"
    # Rare type (5 cells) dropped by min_cells=50; 3 majors survive.
    assert prep.reference.n_cell_types == 3
    assert "Rare" not in prep.reference.cell_types
    # NB overdispersion present so the reference also serves spatial.
    assert prep.reference.phi_g is not None
    assert prep.reference.n_genes == tiny_sc_adata.n_vars


def test_prepare_reference_summary_columns(harness, tiny_sc_adata):
    prep = harness.prepare_reference(tiny_sc_adata, min_cells=50)
    assert set(prep.summary["metric"]) >= {"n_cells", "n_genes", "n_cell_types"}
    assert prep.cell_type_counts.sum() == tiny_sc_adata.n_obs


# ---------------------------------------------------------------------------
# Pseudobulk generation
# ---------------------------------------------------------------------------


def test_generate_pseudobulk_shapes_and_simplex(harness, tiny_sc_adata):
    counts, props, meta = harness.generate_pseudobulk(
        tiny_sc_adata, "cell_type", n_per_regime=4, n_cells=80, seed=0
    )
    assert counts.shape[1] >= 12           # >= 12 mixtures
    assert counts.shape[0] == tiny_sc_adata.n_vars
    np.testing.assert_allclose(props.to_numpy().sum(axis=1), 1.0, atol=1e-6)
    assert (counts.to_numpy() >= 0).all()
    assert np.array_equal(counts.to_numpy(), np.round(counts.to_numpy()))
    assert set(meta["regime"]) == {"easy", "medium", "hard"}


def test_generate_pseudobulk_deterministic(harness, tiny_sc_adata):
    a = harness.generate_pseudobulk(tiny_sc_adata, "cell_type", seed=7)[0]
    b = harness.generate_pseudobulk(tiny_sc_adata, "cell_type", seed=7)[0]
    pd.testing.assert_frame_equal(a, b)


# ---------------------------------------------------------------------------
# Metrics (hand-checked)
# ---------------------------------------------------------------------------


def test_metrics_perfect_recovery(harness):
    true = pd.DataFrame({"A": [0.6, 0.2], "B": [0.4, 0.8]}, index=["s0", "s1"])
    m = harness.compute_bulk_metrics(true, true.copy())
    assert m["rmse"] == pytest.approx(0.0, abs=1e-12)
    assert m["mae"] == pytest.approx(0.0, abs=1e-12)
    assert m["signed_bias"] == pytest.approx(0.0, abs=1e-12)
    assert m["pearson"] == pytest.approx(1.0, abs=1e-9)


def test_metrics_known_error(harness):
    true = pd.DataFrame({"A": [1.0, 0.0], "B": [0.0, 1.0]}, index=["s0", "s1"])
    est = pd.DataFrame({"A": [0.8, 0.1], "B": [0.2, 0.9]}, index=["s0", "s1"])
    m = harness.compute_bulk_metrics(true, est)
    # errors are all magnitude 0.1 / 0.2 / 0.1 / 0.2 ... mae = 0.15, rmse = sqrt(0.025)
    assert m["mae"] == pytest.approx(0.15, abs=1e-9)
    assert m["rmse"] == pytest.approx(np.sqrt(0.025), abs=1e-9)
    assert m["signed_bias"] == pytest.approx(0.0, abs=1e-9)


def test_per_celltype_metrics(harness):
    true = pd.DataFrame({"A": [1.0, 0.0], "B": [0.0, 1.0]}, index=["s0", "s1"])
    est = pd.DataFrame({"A": [0.9, 0.1], "B": [0.1, 0.9]}, index=["s0", "s1"])
    df = harness.per_celltype_metrics(true, est)
    assert set(df.index) == {"A", "B"}
    assert df.loc["A", "rmse"] == pytest.approx(0.1, abs=1e-9)


def test_gene_overlap(harness):
    o = harness.gene_overlap(["g1", "g2", "g3"], ["g2", "g3", "g4"])
    assert o == {"n_query": 3, "n_reference": 3, "n_shared": 2,
                 "n_query_only": 1, "n_reference_only": 1}


def test_align_proportions_no_overlap_raises(harness):
    true = pd.DataFrame({"A": [1.0]}, index=["s0"])
    est = pd.DataFrame({"B": [1.0]}, index=["x9"])
    with pytest.raises(ValueError):
        harness.align_proportions(true, est)
