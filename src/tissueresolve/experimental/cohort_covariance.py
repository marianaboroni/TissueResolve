"""Cohort covariance deconvolution — EXPERIMENT ONLY (Stage 2A, Experiment B).

Central question: given a cohort of bulk samples y_1..y_N that share a KNOWN reference R,
does modelling the cross-sample covariance of proportions recover information that per-sample
deconvolution cannot?

This module implements four estimators purely to *answer* that question empirically. It is NOT a
production solver: it is not registered in the solver registry, changes no default, and is
feature-flagged experimental. The estimators operate in the linear-Gaussian regime y ≈ R θ + ε
(the regime in which cohort covariance is most favourable — so a negative result here is strong).

Estimators (all take known R):
  single           per-sample min-norm least squares  θ̂ = R⁺ y                (baseline)
  first_moment     per-sample LS shrunk toward the cohort mean (James–Stein-like variance cut)
  cross_covariance empirical-Bayes with prior θ_n ~ N(μ, C), C estimated across the cohort
  oracle           returns the true θ (achievability ceiling)

Key theoretical fact these estimators expose: with R FIXED and KNOWN, pooling samples cannot
change the null space of R. Any θ-direction in null(R) is unidentifiable per sample and remains
unidentifiable for the whole cohort — cohort covariance cannot recover it. Cohort methods only
reduce *variance* on the already-identifiable (row-space) directions. See tests.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = ["cohort_deconvolve", "CohortResult", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "cohort_covariance-0.1.0"
METHODS = ("single", "first_moment", "cross_covariance", "oracle")


@dataclass
class CohortResult:
    theta: np.ndarray                 # (K, N) estimated proportions
    method: str
    sigma2: float
    cohort_cov: Optional[np.ndarray] = None   # (K, K) estimated proportion covariance (cross_covariance)
    metadata: dict = field(default_factory=dict)


def _nonneg_normalize(theta: np.ndarray) -> np.ndarray:
    theta = np.clip(theta, 0.0, None)
    s = theta.sum(axis=0, keepdims=True)
    s[s == 0] = 1.0
    return theta / s


def cohort_deconvolve(
    Y: np.ndarray,                    # (G, N) bulk (CPM or counts) — one column per sample
    R: np.ndarray,                    # (G, K) reference gene × type profiles (same scale as Y)
    method: str = "single",
    *,
    theta_true: Optional[np.ndarray] = None,   # (K, N) required only for method="oracle"
    sigma2: Optional[float] = None,            # noise variance; estimated from residuals if None
    shrink: float = 0.5,                       # first_moment shrinkage weight toward cohort mean
    debias_cov: bool = True,                    # subtract LS noise floor from estimated cohort cov
    nonneg: bool = False,                       # clip to ≥0 and renormalise to proportions
) -> CohortResult:
    """Deconvolve a cohort with one of the four estimators. R is assumed known and fixed."""
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    Y = np.asarray(Y, dtype=float)
    R = np.asarray(R, dtype=float)
    G, N = Y.shape
    K = R.shape[1]
    if R.shape[0] != G:
        raise ValueError(f"R has {R.shape[0]} genes but Y has {G}")

    if method == "oracle":
        if theta_true is None:
            raise ValueError("method='oracle' requires theta_true")
        return CohortResult(theta=np.asarray(theta_true, float), method=method, sigma2=0.0,
                            metadata={"algorithm_version": ALGORITHM_VERSION})

    Rpinv = np.linalg.pinv(R)                      # (K, G) min-norm inverse
    single = Rpinv @ Y                             # (K, N) — null(R) component is exactly 0
    resid = Y - R @ single
    if sigma2 is None:
        dof = max(G - K, 1)
        sigma2 = float((resid ** 2).sum() / (dof * N))
    sigma2 = max(sigma2, 1e-12)

    if method == "single":
        theta = single
    elif method == "first_moment":
        mu = single.mean(axis=1, keepdims=True)
        theta = (1.0 - shrink) * single + shrink * mu
    else:  # cross_covariance — empirical Bayes with prior N(mu, C)
        mu = single.mean(axis=1)
        C = np.cov(single, rowvar=True) if N > 1 else np.zeros((K, K))
        C = np.atleast_2d(C)
        if debias_cov:
            # single-sample estimates inflate cov by the LS noise floor σ²(RᵀR)⁺
            C = C - sigma2 * np.linalg.pinv(R.T @ R)
            C = _psd_clip(C)
        RC = R @ C                                 # (G, K)
        S = RC @ R.T + sigma2 * np.eye(G)          # (G, G) marginal cov of y given prior
        M = C @ R.T @ np.linalg.pinv(S)            # (K, G) posterior gain
        theta = mu[:, None] + M @ (Y - (R @ mu)[:, None])
        cohort_cov = C

    if nonneg:
        theta = _nonneg_normalize(theta)
    meta = {"algorithm_version": ALGORITHM_VERSION, "n_samples": N, "n_types": K, "n_genes": G}
    return CohortResult(theta=theta, method=method, sigma2=sigma2,
                        cohort_cov=(cohort_cov if method == "cross_covariance" else None),
                        metadata=meta)


def _psd_clip(C: np.ndarray) -> np.ndarray:
    """Project a symmetric matrix onto the PSD cone (clip negative eigenvalues to 0)."""
    C = 0.5 * (C + C.T)
    w, Vt = np.linalg.eigh(C)
    w = np.clip(w, 0.0, None)
    return (Vt * w) @ Vt.T


def null_direction(R: np.ndarray, tol: float = 1e-6) -> Optional[np.ndarray]:
    """Return a unit vector spanning R's (near-)null space in θ-space, or None if full rank."""
    R = np.asarray(R, float)
    _, S, Vt = np.linalg.svd(R, full_matrices=False)
    smax = S[0] if len(S) else 0.0
    for i in range(len(S) - 1, -1, -1):
        if S[i] <= tol * smax:
            return Vt[i]
    return None
