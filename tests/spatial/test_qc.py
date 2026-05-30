"""
Tests for tissueresolve.spatial.qc.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.spatial.graph import build_hex_graph_from_arrays
from tissueresolve.spatial.model import SpatCARModel
from tissueresolve.spatial.qc import (
    boundary_sharpness,
    compute_model_qc,
    compute_morans_i,
    compute_spot_qc,
    flag_low_quality_spots,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_hex_coords(n_rows=4, n_cols=5):
    rows, cols = [], []
    for r in range(n_rows):
        for c in range(n_cols):
            rows.append(r)
            cols.append(c * 2 + (r % 2))
    return np.array(rows, np.int32), np.array(cols, np.int32)


def _make_props(n_spots=20, n_types=3, seed=42):
    rng = np.random.default_rng(seed)
    return rng.dirichlet(np.ones(n_types), size=n_spots).astype(np.float32)


def _make_fitted_model(n_spots=20, n_types=3, n_genes=24, seed=1):
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    R_cpm = np.zeros((n_types, n_genes), dtype=np.float32)
    for k in range(n_types):
        R_cpm[k, k * block:(k + 1) * block] = rng.uniform(200, 1000, block)
    R_cpm += rng.uniform(0, 5, R_cpm.shape)
    phi_g = np.full(n_genes, 5.0, dtype=np.float32)
    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_types)]
    ref = ReferenceSignature(
        gene_names=gene_names, cell_types=cell_types,
        R_cpm=R_cpm, phi_g=phi_g,
    )
    true_props = rng.dirichlet(np.ones(n_types), size=n_spots).astype(np.float32)
    R_lin = (R_cpm / 1e6).astype(np.float32)
    lib = rng.integers(300, 1500, size=n_spots).astype(np.float32)
    Mu = lib[:, None] * (true_props @ R_lin)
    phi = 5.0
    p_nb = phi / (phi + Mu + 1e-8)
    p_nb = np.clip(p_nb, 1e-6, 1 - 1e-6)
    Y = rng.negative_binomial(phi * np.ones_like(Mu), p_nb).astype(np.float32)

    r, c = _make_hex_coords(4, 5)
    r, c = r[:n_spots], c[:n_spots]
    graph = build_hex_graph_from_arrays(r, c)
    model = SpatCARModel(max_iter=10, lambda_spatial=0.05, verbose=False)
    model.fit(Y, ref, graph, lib)
    return model, graph, Y, lib


# ---------------------------------------------------------------------------
# compute_spot_qc
# ---------------------------------------------------------------------------


class TestComputeSpotQC:
    def test_returns_dataframe(self):
        Pi = _make_props()
        g = build_hex_graph_from_arrays(*_make_hex_coords(4, 5)[:])
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        df = compute_spot_qc(Pi, g)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        df = compute_spot_qc(Pi, g)
        for col in ("prop_entropy", "nb_loglik", "dominant_frac",
                    "n_types_gt05", "spatial_resid"):
            assert col in df.columns

    def test_nb_loglik_nan_without_arrays(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        df = compute_spot_qc(Pi, g)
        assert df["nb_loglik"].isna().all(), "Expected NaN nb_loglik without model arrays"

    def test_nb_loglik_non_nan_with_arrays(self):
        model, graph, Y, lib = _make_fitted_model()
        Pi = model.proportions_
        d_g = model._mismatch.d_g
        R_d = (model._R_lin * d_g[np.newaxis, :]).astype(np.float32)
        df = compute_spot_qc(
            Pi, graph,
            cell_types=model.cell_types_,
            Y_marker=Y,
            R_d=R_d,
            phi_g=model._phi_g,
            lib_sizes=lib,
        )
        assert not df["nb_loglik"].isna().any(), "Expected non-NaN nb_loglik with arrays"

    def test_nb_loglik_negative(self):
        model, graph, Y, lib = _make_fitted_model()
        Pi = model.proportions_
        d_g = model._mismatch.d_g
        R_d = (model._R_lin * d_g[np.newaxis, :]).astype(np.float32)
        df = compute_spot_qc(Pi, graph, Y_marker=Y, R_d=R_d,
                             phi_g=model._phi_g, lib_sizes=lib)
        # NB log-likelihood should be negative (it's a log-prob)
        assert (df["nb_loglik"] < 0).all()

    def test_spatial_resid_nan_without_graph(self):
        Pi = _make_props()
        df = compute_spot_qc(Pi, graph=None)
        assert df["spatial_resid"].isna().all()

    def test_spatial_resid_non_negative(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        df = compute_spot_qc(Pi, g)
        valid = df["spatial_resid"].dropna()
        assert (valid >= 0).all()

    def test_entropy_positive(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        df = compute_spot_qc(Pi, g)
        assert (df["prop_entropy"] >= 0).all()

    def test_dominant_frac_in_zero_one(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        df = compute_spot_qc(Pi, g)
        assert (df["dominant_frac"] >= 0).all()
        assert (df["dominant_frac"] <= 1.0 + 1e-5).all()

    def test_row_count_matches_n_spots(self):
        n = 15
        Pi = _make_props(n_spots=n)
        r, c = _make_hex_coords(3, 5)
        g = build_hex_graph_from_arrays(r[:n], c[:n])
        df = compute_spot_qc(Pi, g)
        assert len(df) == n


# ---------------------------------------------------------------------------
# compute_morans_i
# ---------------------------------------------------------------------------


class TestComputeMonarsI:
    def test_returns_series(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        result = compute_morans_i(Pi, g, cell_types=["A", "B", "C"])
        assert isinstance(result, pd.Series)

    def test_length_equals_n_types(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        result = compute_morans_i(Pi, g)
        assert len(result) == Pi.shape[1]

    def test_high_for_spatially_structured_proportions(self):
        """Left-side spots get CT0; right-side get CT1 → high Moran's I."""
        rng = np.random.default_rng(0)
        r, c = _make_hex_coords(6, 6)
        N = len(r)
        K = 2
        Pi = np.zeros((N, K), dtype=np.float32)
        left = c < c.mean()
        Pi[left, 0] = 0.85
        Pi[left, 1] = 0.15
        Pi[~left, 0] = 0.15
        Pi[~left, 1] = 0.85
        g = build_hex_graph_from_arrays(r, c)
        mi = compute_morans_i(Pi, g)
        # Both types should have high positive Moran's I
        assert mi.min() > 0.0, f"Expected positive Moran's I, got {mi.values}"

    def test_low_for_random_proportions(self):
        """Random proportions should give Moran's I near zero."""
        rng = np.random.default_rng(42)
        r, c = _make_hex_coords(6, 6)
        Pi = rng.dirichlet(np.ones(3), size=len(r)).astype(np.float32)
        g = build_hex_graph_from_arrays(r, c)
        mi = compute_morans_i(Pi, g)
        # Random should have low absolute Moran's I
        assert abs(float(mi.mean())) < 0.5


# ---------------------------------------------------------------------------
# flag_low_quality_spots
# ---------------------------------------------------------------------------


class TestFlagLowQualitySpots:
    def test_returns_boolean_series(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        qc_df = compute_spot_qc(Pi, g)
        flag = flag_low_quality_spots(qc_df)
        assert isinstance(flag, pd.Series)
        assert flag.dtype == bool

    def test_flags_some_spots(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        qc_df = compute_spot_qc(Pi, g)
        flag = flag_low_quality_spots(qc_df)
        # Should flag some but not all spots by default (95th-percentile threshold)
        assert flag.any()
        assert not flag.all()

    def test_explicit_threshold(self):
        Pi = _make_props()
        r, c = _make_hex_coords(4, 5)
        g = build_hex_graph_from_arrays(r[:20], c[:20])
        qc_df = compute_spot_qc(Pi, g)
        # Very high threshold → nothing flagged
        flag_none = flag_low_quality_spots(qc_df, entropy_threshold=100.0)
        assert not flag_none.any()
        # Very low threshold → everything flagged
        flag_all = flag_low_quality_spots(qc_df, entropy_threshold=0.0)
        assert flag_all.all()


# ---------------------------------------------------------------------------
# compute_model_qc
# ---------------------------------------------------------------------------


class TestComputeModelQC:
    def test_returns_dict(self):
        result = compute_model_qc(50, True, [0.1, 0.05, 0.01])
        assert isinstance(result, dict)

    def test_n_iter_present(self):
        result = compute_model_qc(42, True, [0.1])
        assert result["n_iter"] == 42

    def test_converged_flag(self):
        result = compute_model_qc(10, False, [0.5])
        assert result["converged"] is False

    def test_d_g_stats_present(self):
        d_g = np.array([0.9, 1.1, 1.2, 0.8], dtype=np.float32)
        result = compute_model_qc(5, True, [0.1], d_g=d_g)
        assert "d_g_median" in result
        assert "d_g_min" in result
        assert "d_g_max" in result

    def test_n_marker_genes_present(self):
        result = compute_model_qc(5, True, [0.1], marker_genes=["A", "B", "C"])
        assert result["n_marker_genes"] == 3


# ---------------------------------------------------------------------------
# boundary_sharpness
# ---------------------------------------------------------------------------


class TestBoundarySharpness:
    def test_sharp_boundary(self):
        N = 20
        Pi = np.zeros((N, 2), dtype=np.float32)
        array_col = np.arange(N, dtype=np.int32)
        Pi[:10, 0] = 0.9; Pi[:10, 1] = 0.1
        Pi[10:, 0] = 0.1; Pi[10:, 1] = 0.9
        array_row = np.zeros(N, dtype=np.int32)
        sharpness = boundary_sharpness(Pi, array_row, array_col, boundary_col=10)
        assert sharpness > 0.5, f"Expected high sharpness, got {sharpness}"

    def test_gradual_boundary_low_sharpness(self):
        N = 20
        Pi = np.ones((N, 2), dtype=np.float32) * 0.5  # uniform
        array_row = np.zeros(N, dtype=np.int32)
        array_col = np.arange(N, dtype=np.int32)
        sharpness = boundary_sharpness(Pi, array_row, array_col, boundary_col=10)
        assert sharpness < 0.1

    def test_confident_type_flip_scores_higher_than_same_type(self):
        # Regression: a confident dominant-TYPE flip across the boundary is
        # the signature of a sharp boundary and must score higher than a
        # boundary where the same type dominates both sides at equal
        # confidence.  The buggy implementation collapsed each side to the
        # dominant-fraction magnitude and scored the type flip as 0.0.
        N = 20
        array_row = np.zeros(N, dtype=np.int32)
        array_col = np.arange(N, dtype=np.int32)

        flip = np.zeros((N, 2), dtype=np.float32)
        flip[:10, 0] = 0.9; flip[:10, 1] = 0.1
        flip[10:, 0] = 0.1; flip[10:, 1] = 0.9

        same = np.zeros((N, 2), dtype=np.float32)
        same[:, 0] = 0.9; same[:, 1] = 0.1  # type 0 dominant at 0.9 everywhere

        s_flip = boundary_sharpness(flip, array_row, array_col, boundary_col=10)
        s_same = boundary_sharpness(same, array_row, array_col, boundary_col=10)
        assert s_flip > s_same
        assert s_same == 0.0  # no compositional change -> no boundary

    def test_sharpness_bounded_in_unit_interval(self):
        N = 20
        array_row = np.zeros(N, dtype=np.int32)
        array_col = np.arange(N, dtype=np.int32)
        rng = np.random.default_rng(0)
        Pi = rng.dirichlet(np.ones(3), size=N).astype(np.float32)
        s = boundary_sharpness(Pi, array_row, array_col, boundary_col=10)
        assert 0.0 <= s <= 1.0
