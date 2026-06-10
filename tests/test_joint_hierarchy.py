"""Invariant + synthetic-recovery tests for the soft joint hierarchy (Phase 2B).

Offline, deterministic.
"""
import numpy as np
import pandas as pd
import pytest

from tissueresolve.experimental.soft_hierarchy import (
    fit_alternating_soft_hierarchy, fit_joint_soft_hierarchy,
    reconcile_broad_fine, compute_hierarchy_consistency, compute_error_decomposition)

FAMILY_MAP = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B"}
SUBTYPES = ["A1", "A2", "A3", "B1", "B2"]


def _toy_sig(G=60, seed=0):
    rng = np.random.default_rng(seed)
    # distinct-ish subtype signatures (G genes × 5 subtypes), non-negative
    S = rng.gamma(2.0, 1.0, size=(G, len(SUBTYPES)))
    return S


def _synth_samples(S, n=8, seed=1):
    rng = np.random.default_rng(seed)
    theta = rng.dirichlet(np.ones(len(SUBTYPES)), size=n)   # known fine truth
    Y = theta @ S.T                                         # y = S θ (noiseless)
    return Y, theta


def test_invariants_nonneg_simplex_and_aggregation():
    S = _toy_sig()
    Y, _ = _synth_samples(S)
    r = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=1.0)
    pi, th = r.broad_proportions.to_numpy(), r.reconciled_fine_proportions.to_numpy()
    assert (pi >= -1e-9).all() and (th >= -1e-9).all()
    assert np.allclose(pi.sum(1), 1.0, atol=1e-6)
    assert np.allclose(th.sum(1), 1.0, atol=1e-6)
    # reconciled fine aggregates to broad (hierarchy consistency, mass conserved)
    cons = compute_hierarchy_consistency(r.broad_proportions, r.reconciled_fine_proportions, FAMILY_MAP)
    assert float(cons.max()) < 1e-6


def test_synthetic_recovery_noiseless():
    """y = S θ_true with distinct signatures → recover θ_true reasonably."""
    S = _toy_sig(G=120, seed=3)
    Y, theta_true = _synth_samples(S, n=10, seed=4)
    r = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=0.5, max_iter=200)
    est = r.reconciled_fine_proportions.to_numpy()
    err = np.abs(est - theta_true).mean()
    assert err < 0.05, f"mean abs recovery error {err:.3f} too high"


def test_lambda_h_zero_decouples_to_flat_like():
    """λ_h=0 → fine estimate ≈ unconstrained flat NNLS (decoupled)."""
    from scipy.optimize import nnls
    S = _toy_sig(seed=5); Y, _ = _synth_samples(S, n=4, seed=6)
    r = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=0.0)
    for s in range(Y.shape[0]):
        flat, _ = nnls(S, Y[s]); flat = flat / flat.sum()
        assert np.allclose(r.raw_fine_proportions.to_numpy()[s], flat, atol=1e-4)


def test_strong_lambda_h_increases_consistency():
    S = _toy_sig(seed=7); Y, _ = _synth_samples(S, n=6, seed=8)
    weak = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=0.0)
    strong = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=50.0)
    # raw consistency ‖π−Aθ‖ should be smaller (or equal) with strong coupling
    assert strong.hierarchy_consistency_error.mean() <= weak.hierarchy_consistency_error.mean() + 1e-9


def test_no_sparsity_collapse():
    """Joint solver must not collapse to ~1 population (no strong-L1 behaviour)."""
    S = _toy_sig(seed=9); Y, _ = _synth_samples(S, n=6, seed=10)
    r = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=1.0, lambda_2=0.01)
    th = r.reconciled_fine_proportions.to_numpy()
    eff_n = np.exp(-(np.where(th > 0, th * np.log(np.clip(th, 1e-12, None)), 0)).sum(1)).mean()
    assert eff_n > 2.0, f"effective-N collapsed to {eff_n:.2f}"


def test_flat_initialization_and_convergence():
    S = _toy_sig(seed=11); Y, _ = _synth_samples(S, n=4, seed=12)
    r = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=1.0,
                                       init="flat_nnls", max_iter=200, tol=1e-6)
    assert all(n >= 1 for n in r.n_iterations)
    assert any(r.convergence_status)  # at least some samples converge


def test_reconcile_mass_conservation():
    broad = pd.DataFrame({"A": [0.6, 0.3], "B": [0.4, 0.7]}, index=["s0", "s1"])
    fine = pd.DataFrame({"A1": [0.3, 0.1], "A2": [0.3, 0.2], "A3": [0.0, 0.0],
                         "B1": [0.2, 0.4], "B2": [0.2, 0.3]}, index=["s0", "s1"])
    rec = reconcile_broad_fine(broad, fine, FAMILY_MAP)
    for fam, mem in {"A": ["A1", "A2", "A3"], "B": ["B1", "B2"]}.items():
        assert np.allclose(rec[mem].sum(1), broad[fam], atol=1e-9)


def test_joint_slsqp_runs_and_conserves():
    S = _toy_sig(G=40, seed=13); Y, _ = _synth_samples(S, n=3, seed=14)
    r = fit_joint_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=1.0)
    assert np.allclose(r.broad_proportions.sum(1), 1.0, atol=1e-4)
    assert np.allclose(r.raw_fine_proportions.sum(1), 1.0, atol=1e-4)


def test_error_decomposition_keys():
    S = _toy_sig(seed=15); Y, theta = _synth_samples(S, n=5, seed=16)
    r = fit_alternating_soft_hierarchy(Y, S, FAMILY_MAP, SUBTYPES, lambda_h=1.0)
    tdf = pd.DataFrame(theta, index=r.reconciled_fine_proportions.index, columns=SUBTYPES)
    dec = compute_error_decomposition(tdf, r.broad_proportions, r.reconciled_fine_proportions, FAMILY_MAP)
    for k in ["broad_error", "conditional_fine_error", "consistency_error", "total_fine_error"]:
        assert k in dec and np.isfinite(dec[k])
