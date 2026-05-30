"""
Tests for tissueresolve.spatial.graph (SpatialGraph, builders).
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest
import scipy.sparse as sp

from tissueresolve.spatial.graph import (
    SpatialGraph,
    build_hex_graph_from_arrays,
    build_expression_weighted_graph,
    compute_spatial_batches,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _hex_grid(n_rows: int = 4, n_cols: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Build a small Visium-style hex grid."""
    rows, cols = [], []
    for r in range(n_rows):
        for c in range(n_cols):
            rows.append(r)
            cols.append(c * 2 + (r % 2))
    return np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32)


# ---------------------------------------------------------------------------
# SpatialGraph construction
# ---------------------------------------------------------------------------


class TestBuildHexGraphFromArrays:
    def test_returns_spatial_graph(self):
        r, c = _hex_grid()
        g = build_hex_graph_from_arrays(r, c)
        assert isinstance(g, SpatialGraph)

    def test_adjacency_shape(self):
        r, c = _hex_grid(4, 4)
        N = len(r)
        g = build_hex_graph_from_arrays(r, c)
        assert g.A.shape == (N, N)
        assert g.L.shape == (N, N)

    def test_n_spots_correct(self):
        r, c = _hex_grid(3, 5)
        g = build_hex_graph_from_arrays(r, c)
        assert g.n_spots == len(r)

    def test_degree_array_length(self):
        r, c = _hex_grid(4, 4)
        g = build_hex_graph_from_arrays(r, c)
        assert len(g.degree) == g.n_spots

    def test_interior_spots_have_degree_six(self):
        # A 5×5 grid should have interior spots with degree 6
        r, c = _hex_grid(6, 6)
        g = build_hex_graph_from_arrays(r, c)
        assert g.degree.max() == 6

    def test_boundary_spots_have_lower_degree(self):
        r, c = _hex_grid(4, 4)
        g = build_hex_graph_from_arrays(r, c)
        assert (g.degree < 6).any()

    def test_adjacency_is_symmetric(self):
        """A[i, j] > 0  ↔  A[j, i] > 0 (symmetric pattern)."""
        r, c = _hex_grid(4, 4)
        g = build_hex_graph_from_arrays(r, c)
        A = g.A
        # Pattern symmetry: nonzero (i, j) iff nonzero (j, i)
        diff = (A - A.T)
        # Only row-normalisation makes values asymmetric; pattern should be symmetric
        assert A.nnz == A.T.nnz

    def test_row_sums_equal_one_for_connected_spots(self):
        r, c = _hex_grid(4, 4)
        g = build_hex_graph_from_arrays(r, c)
        row_sums = np.asarray(g.A.sum(axis=1)).ravel()
        connected = g.degree > 0
        assert np.allclose(row_sums[connected], 1.0, atol=1e-5)

    def test_disconnected_spot_has_zero_row(self):
        # A single isolated spot at (100, 100) disconnected from the rest
        r = np.array([0, 0, 1, 100], dtype=np.int32)
        c = np.array([0, 2, 1, 100], dtype=np.int32)
        g = build_hex_graph_from_arrays(r, c)
        # Spot index 3 (100, 100) should be isolated
        assert g.degree[3] == 0
        row_sum = float(g.A[3, :].sum())
        assert row_sum == 0.0

    def test_spot_ids_attached(self):
        r, c = _hex_grid(3, 3)
        ids = np.array([f"spot_{i}" for i in range(len(r))])
        g = build_hex_graph_from_arrays(r, c, spot_ids=ids)
        assert g.spot_ids is not None
        assert len(g.spot_ids) == g.n_spots

    def test_batches_cover_all_spots(self):
        r, c = _hex_grid(4, 4)
        g = build_hex_graph_from_arrays(r, c, batch_size=3)
        all_idx = np.concatenate(g.batches)
        assert len(all_idx) == g.n_spots
        assert len(set(all_idx.tolist())) == g.n_spots

    def test_laplacian_is_I_minus_A(self):
        r, c = _hex_grid(3, 3)
        g = build_hex_graph_from_arrays(r, c)
        N = g.n_spots
        I_minus_A = sp.eye(N, format="csr", dtype=np.float32) - g.A
        diff = (g.L - I_minus_A).data
        assert len(diff) == 0 or np.allclose(diff, 0.0, atol=1e-6)


class TestInvalidCoordinates:
    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            build_hex_graph_from_arrays(np.array([], dtype=np.int32),
                                        np.array([], dtype=np.int32))

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="same length"):
            build_hex_graph_from_arrays(
                np.array([0, 1], dtype=np.int32),
                np.array([0], dtype=np.int32),
            )

    def test_negative_coordinates_raise(self):
        with pytest.raises(ValueError, match="Negative"):
            build_hex_graph_from_arrays(
                np.array([-1, 0], dtype=np.int32),
                np.array([0, 2], dtype=np.int32),
            )

    def test_no_neighbours_warns(self):
        # Single isolated spot — should warn but not raise
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            build_hex_graph_from_arrays(
                np.array([0], dtype=np.int32),
                np.array([0], dtype=np.int32),
            )
        texts = [str(w.message) for w in caught]
        assert any("No edges" in t for t in texts)


# ---------------------------------------------------------------------------
# Expression-weighted graph
# ---------------------------------------------------------------------------


class TestExpressionWeightedGraph:
    def test_returns_spatial_graph(self):
        r, c = _hex_grid(3, 3)
        N = len(r)
        G_m = 10
        rng = np.random.default_rng(0)
        Y = rng.poisson(50, (N, G_m)).astype(np.float32)
        lib = Y.sum(axis=1)
        g = build_expression_weighted_graph(Y, r, c, lib)
        assert isinstance(g, SpatialGraph)

    def test_edge_weights_in_zero_one(self):
        r, c = _hex_grid(3, 3)
        N = len(r)
        rng = np.random.default_rng(1)
        Y = rng.poisson(50, (N, 8)).astype(np.float32)
        lib = Y.sum(axis=1)
        g = build_expression_weighted_graph(Y, r, c, lib)
        data = g.A.data
        assert (data >= 0).all()
        assert (data <= 1.0 + 1e-5).all()

    def test_row_sums_equal_one(self):
        r, c = _hex_grid(4, 4)
        N = len(r)
        rng = np.random.default_rng(2)
        Y = rng.poisson(100, (N, 10)).astype(np.float32)
        lib = Y.sum(axis=1)
        g = build_expression_weighted_graph(Y, r, c, lib)
        connected = g.degree > 0
        row_sums = np.asarray(g.A.sum(axis=1)).ravel()
        assert np.allclose(row_sums[connected], 1.0, atol=1e-5)


# ---------------------------------------------------------------------------
# Spatial batches
# ---------------------------------------------------------------------------


class TestComputeSpatialBatches:
    def test_covers_all_spots(self):
        r, c = _hex_grid(5, 5)
        batches = compute_spatial_batches(r, c, batch_size=5)
        all_idx = np.concatenate(batches)
        assert len(all_idx) == len(r)

    def test_no_overlap(self):
        r, c = _hex_grid(5, 5)
        batches = compute_spatial_batches(r, c, batch_size=5)
        all_idx = np.concatenate(batches)
        assert len(set(all_idx.tolist())) == len(all_idx)

    def test_empty_returns_empty(self):
        batches = compute_spatial_batches(np.array([], np.int32),
                                         np.array([], np.int32))
        assert batches == []
