"""Cohort covariance estimator — synthetic proofs (Stage 2A, Experiment B).

Central claim: with a KNOWN reference R, cohort cross-sample covariance cannot recover θ-directions
in null(R) — not with, not without, cross-sample variation. It only reduces variance on the
already-identifiable (row-space) directions. These tests encode exactly that.
"""
from __future__ import annotations

import numpy as np
import pytest

from tissueresolve.experimental.cohort_covariance import cohort_deconvolve, null_direction


def _rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def _confounded_R(G=120):
    """A and B share an identical profile (confounded); C is distinct."""
    colA = np.zeros(G); colA[:40] = 1.0
    colC = np.zeros(G); colC[40:80] = 1.0
    return np.vstack([colA, colA.copy(), colC]).T          # (G, 3), cols A==B


def test_null_direction_identifies_confounded_axis():
    nd = null_direction(_confounded_R())
    assert nd is not None
    # the near-null θ-direction is A − B (up to sign); C is uninvolved
    assert abs(abs(nd[0]) - abs(nd[1])) < 1e-6 and abs(nd[2]) < 1e-6
    assert np.sign(nd[0]) != np.sign(nd[1])


def test_cohort_covariance_cannot_create_information_without_variation():
    """Proportions constant across the cohort ⇒ cross-covariance adds nothing on the null axis."""
    rng = np.random.default_rng(0)
    R = _confounded_R(); K, N = 3, 50
    nd = null_direction(R)
    Th = np.tile(np.array([0.5, 0.3, 0.2])[:, None], (1, N))   # identical every sample
    Y = R @ Th + rng.normal(0, 0.05, (R.shape[0], N))
    s = cohort_deconvolve(Y, R, "single").theta
    cc = cohort_deconvolve(Y, R, "cross_covariance").theta
    assert abs(_rms((cc - Th).T @ nd) - _rms((s - Th).T @ nd)) < 1e-8


def test_cohort_covariance_no_gain_on_nullspace_even_with_variation():
    """Even when A and B vary independently across samples, a KNOWN R fixes the null space:
    cross-covariance error on the confounded axis equals single-sample error."""
    rng = np.random.default_rng(1)
    R = _confounded_R(); K, N = 3, 60
    nd = null_direction(R)
    Th = np.abs(rng.normal(0.4, 0.15, (K, N)))                 # A, B vary independently
    Y = R @ Th + rng.normal(0, 0.05, (R.shape[0], N))
    s = cohort_deconvolve(Y, R, "single").theta
    cc = cohort_deconvolve(Y, R, "cross_covariance").theta
    assert abs(_rms((cc - Th).T @ nd) - _rms((s - Th).T @ nd)) < 1e-8


def test_cohort_covariance_gain_only_where_mathematically_supported():
    """Positive control: on identifiable directions with high noise and clustered proportions,
    shrinkage (first_moment) and empirical Bayes (cross_covariance) reduce error — but NOT on a
    confounded axis. Gain appears only where the linear algebra supports it."""
    rng = np.random.default_rng(2)
    G, K, N = 100, 3, 60
    R = np.zeros((G, K)); R[:33, 0] = 1; R[33:66, 1] = 1; R[66:, 2] = 1   # orthogonal, well-conditioned
    Th = np.array([0.5, 0.3, 0.2])[:, None] + rng.normal(0, 0.03, (K, N)) # clustered
    Y = R @ Th + rng.normal(0, 0.5, (G, N))                               # high noise
    single = cohort_deconvolve(Y, R, "single", sigma2=0.25).theta
    fm = cohort_deconvolve(Y, R, "first_moment", sigma2=0.25).theta
    cc = cohort_deconvolve(Y, R, "cross_covariance", sigma2=0.25).theta
    # variance reduction on identifiable directions
    assert _rms(fm - Th) < _rms(single - Th)
    assert _rms(cc - Th) < _rms(single - Th)


def test_oracle_returns_truth():
    R = _confounded_R()
    Th = np.abs(np.random.default_rng(3).normal(0.3, 0.1, (3, 10)))
    Y = R @ Th
    out = cohort_deconvolve(Y, R, "oracle", theta_true=Th).theta
    np.testing.assert_array_equal(out, Th)


def test_oracle_requires_truth_and_bad_method_rejected():
    R = _confounded_R(); Y = R @ np.ones((3, 5))
    with pytest.raises(ValueError):
        cohort_deconvolve(Y, R, "oracle")
    with pytest.raises(ValueError):
        cohort_deconvolve(Y, R, "not_a_method")


def test_deterministic():
    rng = np.random.default_rng(4)
    R = _confounded_R(); Y = R @ np.abs(rng.normal(0.3, 0.1, (3, 20))) + rng.normal(0, 0.05, (120, 20))
    a = cohort_deconvolve(Y, R, "cross_covariance").theta
    b = cohort_deconvolve(Y, R, "cross_covariance").theta
    np.testing.assert_array_equal(a, b)
