"""
Tests for tissueresolve.spatial.model (SpatCARModel).
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.spatial.graph import build_hex_graph_from_arrays
from tissueresolve.spatial.model import SpatCARModel, _nb_multiplicative_update


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_hex_coords(n_rows: int = 4, n_cols: int = 4):
    rows, cols = [], []
    for r in range(n_rows):
        for c in range(n_cols):
            rows.append(r)
            cols.append(c * 2 + (r % 2))
    return np.array(rows, np.int32), np.array(cols, np.int32)


def _make_synthetic_data(
    n_spots: int = 20,
    n_cell_types: int = 3,
    n_genes: int = 30,
    seed: int = 42,
) -> tuple[np.ndarray, ReferenceSignature, np.ndarray, np.ndarray, np.ndarray]:
    """Returns Y_marker, ref_marker, graph, lib_sizes, true_props."""
    rng = np.random.default_rng(seed)

    # Build reference with block-diagonal structure for separability
    block = n_genes // n_cell_types
    R_cpm = np.zeros((n_cell_types, n_genes), dtype=np.float32)
    for k in range(n_cell_types):
        R_cpm[k, k * block:(k + 1) * block] = rng.uniform(500, 2000, block)
    R_cpm += rng.uniform(0, 10, R_cpm.shape)
    phi_g = np.full(n_genes, 5.0, dtype=np.float32)

    gene_names = [f"G{i:04d}" for i in range(n_genes)]
    cell_types = [f"CT{k}" for k in range(n_cell_types)]
    ref = ReferenceSignature(
        gene_names=gene_names,
        cell_types=cell_types,
        R_cpm=R_cpm,
        phi_g=phi_g,
    )

    # True proportions
    true_props = rng.dirichlet(np.ones(n_cell_types), size=n_spots).astype(np.float32)

    # Synthetic NB counts
    R_lin = (R_cpm / 1e6).astype(np.float32)
    lib_sizes = rng.integers(500, 2000, size=n_spots).astype(np.float32)
    Mu = lib_sizes[:, None] * (true_props @ R_lin)
    phi = 5.0
    p_nb = phi / (phi + Mu + 1e-8)
    p_nb = np.clip(p_nb, 1e-6, 1 - 1e-6)
    Y_marker = rng.negative_binomial(
        phi * np.ones_like(Mu), p_nb
    ).astype(np.float32)

    return Y_marker, ref, lib_sizes, true_props


# ---------------------------------------------------------------------------
# NNLS warm start
# ---------------------------------------------------------------------------


class TestNNLSWarmStart:
    def test_nnls_returns_non_negative_proportions(self):
        Y, ref, lib, _ = _make_synthetic_data()
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)

        model = SpatCARModel(lambda_spatial=0.0, max_iter=1, verbose=False)
        model.fit(Y, ref, graph, lib)
        assert (model.proportions_ >= 0).all()

    def test_nnls_proportions_sum_to_one(self):
        Y, ref, lib, _ = _make_synthetic_data()
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)

        model = SpatCARModel(lambda_spatial=0.0, max_iter=1, verbose=False)
        model.fit(Y, ref, graph, lib)
        row_sums = model.proportions_.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-5)


# ---------------------------------------------------------------------------
# Model fit
# ---------------------------------------------------------------------------


class TestSpatCARModelFit:
    def _fit(self, lambda_spatial=0.0, max_iter=30, seed=42):
        Y, ref, lib, _ = _make_synthetic_data(seed=seed)
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)
        model = SpatCARModel(
            lambda_spatial=lambda_spatial,
            max_iter=max_iter,
            random_state=seed,
            verbose=False,
        )
        model.fit(Y, ref, graph, lib)
        return model

    def test_fit_returns_self(self):
        Y, ref, lib, _ = _make_synthetic_data()
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)
        model = SpatCARModel(max_iter=5, verbose=False)
        ret = model.fit(Y, ref, graph, lib)
        assert ret is model

    def test_proportions_shape(self):
        model = self._fit()
        n_spots = 20
        n_cell_types = 3
        assert model.proportions_.shape == (n_spots, n_cell_types)

    def test_proportions_non_negative(self):
        model = self._fit()
        assert (model.proportions_ >= 0).all()

    def test_proportions_sum_to_one(self):
        model = self._fit()
        row_sums = model.proportions_.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=1e-4)

    def test_convergence_trace_recorded(self):
        model = self._fit(max_iter=20)
        assert len(model.convergence_trace_) > 0

    def test_n_iter_recorded(self):
        model = self._fit(max_iter=10)
        assert model.n_iter_ > 0
        assert model.n_iter_ <= 10

    def test_lambda_spatial_recorded(self):
        model = self._fit(lambda_spatial=0.2)
        assert model.lambda_spatial == 0.2

    def test_alpha_recorded_when_lambda_zero(self):
        model = self._fit(lambda_spatial=0.0)
        assert model.alpha == 0.0

    def test_alpha_positive_when_lambda_positive(self):
        model = self._fit(lambda_spatial=0.3)
        assert model.alpha > 0.0

    def test_alpha_max_half(self):
        # alpha = min(0.5, lambda) so never > 0.5
        model = self._fit(lambda_spatial=10.0)
        assert model.alpha <= 0.5

    def test_cell_types_stored(self):
        model = self._fit()
        assert len(model.cell_types_) == 3

    def test_marker_genes_stored(self):
        model = self._fit()
        assert len(model.marker_genes_) == 30

    def test_convergence_trace_is_decreasing_approximately(self):
        model = self._fit(max_iter=50)
        trace = model.convergence_trace_
        if len(trace) >= 3:
            # First few steps should generally decrease
            assert trace[-1] < trace[0] * 2, (
                "Trace should roughly decrease over 50 iterations"
            )

    def test_no_convergence_warning_when_converged(self):
        # With enough iterations and very loose tol, should converge
        Y, ref, lib, _ = _make_synthetic_data()
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            SpatCARModel(max_iter=200, tol=1e-2, verbose=False).fit(Y, ref, graph, lib)
        conv_warns = [
            w for w in caught
            if "converge" in str(w.message).lower()
            and "did not converge" in str(w.message).lower()
        ]
        # Should converge with loose tol — no "did not converge" warning
        assert len(conv_warns) == 0

    def test_non_convergence_emits_warning(self):
        Y, ref, lib, _ = _make_synthetic_data()
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            SpatCARModel(max_iter=1, tol=1e-12, verbose=False).fit(Y, ref, graph, lib)
        texts = [str(w.message) for w in caught]
        assert any("converge" in t.lower() for t in texts), (
            "Expected non-convergence warning, got: " + str(texts)
        )

    def test_smoothing_lambda_changes_result(self):
        """Model with lambda > 0 should differ from lambda=0."""
        Y, ref, lib, _ = _make_synthetic_data(seed=7)
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)

        m0 = SpatCARModel(lambda_spatial=0.0, max_iter=20, verbose=False)
        m0.fit(Y, ref, graph, lib)
        m1 = SpatCARModel(lambda_spatial=0.5, max_iter=20, verbose=False)
        m1.fit(Y, ref, graph, lib)

        diff = np.abs(m0.proportions_ - m1.proportions_).mean()
        assert diff > 1e-6, "Spatial smoothing should change proportions"

    def test_no_phi_g_warns(self):
        """When ReferenceSignature has no phi_g, a warning is emitted."""
        Y, ref_full, lib, _ = _make_synthetic_data()
        n_spots = Y.shape[0]
        r, c = _make_hex_coords(4, 5)
        r, c = r[:n_spots], c[:n_spots]
        graph = build_hex_graph_from_arrays(r, c)

        # Build ref without phi_g
        ref_no_phi = ReferenceSignature(
            gene_names=ref_full.gene_names,
            cell_types=ref_full.cell_types,
            R_cpm=ref_full.R_cpm,
            phi_g=None,
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            SpatCARModel(max_iter=2, verbose=False).fit(Y, ref_no_phi, graph, lib)
        texts = [str(w.message) for w in caught]
        assert any("phi_g" in t.lower() or "dispersion" in t.lower() for t in texts)


# ---------------------------------------------------------------------------
# NB multiplicative update helper
# ---------------------------------------------------------------------------


class TestNBMultiplicativeUpdate:
    def test_output_non_negative(self):
        rng = np.random.default_rng(0)
        B, K, G_m = 5, 3, 10
        Y_b = rng.poisson(50, (B, G_m)).astype(np.float32)
        Pi_b = rng.dirichlet(np.ones(K), size=B).astype(np.float32)
        R_d = rng.uniform(1e-4, 1e-3, (K, G_m)).astype(np.float32)
        phi_g = np.ones(G_m, dtype=np.float32) * 5.0
        lib_b = np.full(B, 1000.0, dtype=np.float32)

        result = _nb_multiplicative_update(Y_b, Pi_b, R_d, phi_g, lib_b)
        assert (result >= 0).all()

    def test_output_sums_to_one(self):
        rng = np.random.default_rng(1)
        B, K, G_m = 4, 3, 10
        Y_b = rng.poisson(30, (B, G_m)).astype(np.float32)
        Pi_b = rng.dirichlet(np.ones(K), size=B).astype(np.float32)
        R_d = rng.uniform(1e-4, 1e-3, (K, G_m)).astype(np.float32)
        phi_g = np.ones(G_m, dtype=np.float32) * 3.0
        lib_b = np.full(B, 800.0, dtype=np.float32)

        result = _nb_multiplicative_update(Y_b, Pi_b, R_d, phi_g, lib_b)
        assert np.allclose(result.sum(axis=1), 1.0, atol=1e-5)
