"""
Tests for benchmark.spillover and reference.pairwise_markers.

Synthetic, offline.  Covers spillover-matrix shape/normalisation, main leaking
partner identification, the expression-only proxy, and pairwise marker
selection/scoring/augmentation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.benchmark.spillover import (
    compute_spillover_risk,
    estimate_spillover_matrix,
    expression_spillover_proxy,
    simulate_pairwise_mixtures,
    simulate_pure_profiles,
    simulate_realistic_mixtures,
    summarize_spillover_partners,
)
from tissueresolve.reference.pairwise_markers import (
    augment_marker_panel_for_confusable_pairs,
    score_pairwise_markers,
    select_pairwise_discriminative_genes,
)
from tissueresolve.results import ReferenceSignature


def _block_ref(n_types=3, block=12):
    G = n_types * block
    R = np.full((n_types, G), 5.0)
    for k in range(n_types):
        R[k, k * block:(k + 1) * block] = 1000.0
    R = R.astype(np.float32)
    cts = [f"T{k}" for k in range(n_types)]
    return ReferenceSignature(
        gene_names=[f"G{i:03d}" for i in range(G)],
        cell_types=cts, R_cpm=R, R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={ct: 100 for ct in cts},
    )


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def test_simulate_pure_profiles_shape():
    ref = _block_ref()
    pure = simulate_pure_profiles(ref)
    assert pure.shape == (ref.n_genes, ref.n_cell_types)
    assert (pure.to_numpy() >= 0).all()
    assert list(pure.columns) == list(ref.cell_types)


def test_simulate_pairwise_mixtures_truth_sums_to_one():
    ref = _block_ref()
    counts, truth = simulate_pairwise_mixtures(ref)
    assert counts.shape[0] == ref.n_genes
    np.testing.assert_allclose(truth.to_numpy().sum(axis=1), 1.0, atol=1e-9)


def test_simulate_realistic_mixtures():
    ref = _block_ref()
    counts, truth = simulate_realistic_mixtures(ref, n=10, seed=0)
    assert counts.shape == (ref.n_genes, 10)
    np.testing.assert_allclose(truth.to_numpy().sum(axis=1), 1.0, atol=1e-9)


# ---------------------------------------------------------------------------
# Spillover matrix
# ---------------------------------------------------------------------------


def test_spillover_matrix_shape_and_rows_sum_to_one():
    ref = _block_ref()
    M = estimate_spillover_matrix(ref)
    assert M.shape == (ref.n_cell_types, ref.n_cell_types)
    np.testing.assert_allclose(M.to_numpy().sum(axis=1), 1.0, atol=1e-4)


def test_spillover_matrix_separable_is_near_identity():
    ref = _block_ref()  # disjoint blocks → easily identifiable
    M = estimate_spillover_matrix(ref)
    # diagonal dominates for a clearly separable reference
    assert (np.diag(M.to_numpy()) > 0.8).all()


def test_compute_spillover_risk_identifies_partner():
    # Two near-identical types (T0,T1) + a distinct one (T2): T0 should leak to T1.
    G = 30
    R = np.full((3, G), 5.0)
    R[0, :10] = 1000.0
    R[1, :10] = 990.0           # almost identical to T0
    R[2, 20:] = 1000.0          # distinct
    R = R.astype(np.float32)
    cts = ["T0", "T1", "T2"]
    ref = ReferenceSignature(
        gene_names=[f"G{i:03d}" for i in range(G)], cell_types=cts,
        R_cpm=R, R_log=np.log1p(R).astype(np.float32),
        phi_g=np.full(G, 5.0, dtype=np.float32),
        n_cells_per_type={ct: 100 for ct in cts})
    M = estimate_spillover_matrix(ref)
    risk = compute_spillover_risk(M)
    assert set(risk.columns) == {"self_retention", "spillover_risk",
                                 "main_partner", "main_partner_fraction"}
    # T0's main leaking partner should be the near-identical T1 (not T2).
    assert risk.loc["T0", "main_partner"] == "T1"


def test_summarize_spillover_partners_long_form():
    ref = _block_ref()
    M = estimate_spillover_matrix(ref)
    long = summarize_spillover_partners(M, min_fraction=0.0)
    assert set(long.columns) == {"true_type", "leaks_into", "fraction"}


def test_expression_proxy_rows_sum_to_one():
    ref = _block_ref()
    proxy = expression_spillover_proxy(ref)
    assert proxy.shape == (ref.n_cell_types, ref.n_cell_types)
    np.testing.assert_allclose(proxy.to_numpy().sum(axis=1), 1.0, atol=1e-6)


# ---------------------------------------------------------------------------
# Pairwise markers
# ---------------------------------------------------------------------------


def test_select_pairwise_discriminative_genes():
    ref = _block_ref()  # T0 high in G000-G011, T1 in G012-G023
    genes = select_pairwise_discriminative_genes(ref, "T0", "T1", top_n=8)
    assert len(genes) == 8
    # discriminative genes should come from the two blocks
    block_genes = {f"G{i:03d}" for i in range(24)}
    assert all(g in block_genes for g in genes)


def test_score_pairwise_markers_columns():
    ref = _block_ref()
    scored = score_pairwise_markers(ref, "T0", "T1")
    for col in ("log2fc", "abs_log2fc", "donor_stability",
                "detectability", "leakage", "score"):
        assert col in scored.columns
    assert scored["score"].is_monotonic_decreasing


def test_augment_marker_panel_adds_pair_genes():
    ref = _block_ref()
    base = ["G000", "G001"]
    panel, added = augment_marker_panel_for_confusable_pairs(
        ref, base, [("T0", "T1")], top_n=5)
    assert base[0] in panel  # base preserved
    assert len(panel) > len(base)
    assert added[("T0", "T1")]  # some genes were added
