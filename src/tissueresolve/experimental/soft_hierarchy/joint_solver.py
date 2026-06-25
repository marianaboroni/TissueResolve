"""Soft joint broad–fine hierarchy solver (Phase 2B, experimental).

Tests the hypothesis that the residual hierarchical-vs-flat gap is caused by the
*sequential* factorisation (broad → freeze → fine).  Estimates broad ``π`` and
fine ``θ`` with a SOFT consistency coupling ``π ≈ Aθ`` instead of freezing broad
mass.  Per-sample constrained objective (no strong L1; rule 6):

    min_{π,θ≥0}  ‖W_b(y − Bπ)‖² + α‖W_f(y − Sθ)‖² + λ_h‖π − Aθ‖² + λ₂‖θ‖²
    s.t.  Σπ = 1,  Σθ = 1

B (family signatures) = per-family mean of member fine signatures; A (F×K) is the
hierarchy aggregation (A[f,k]=1 iff k∈f).  Simplex constraints are enforced by
renormalisation after a non-negative solve (matching the package's flat solver).

Two solvers: **alternating** (NNLS sub-problems; primary) and **joint** (SLSQP;
secondary).  Pure scipy/numpy — no PyTorch/JAX/Pyro (rules 7,8).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Optional

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ALGORITHM_VERSION = "joint_soft_hierarchy-0.1.0"
FEATURE_STATUS = "experimental"


@dataclass
class JointHierarchyResult:
    broad_proportions: pd.DataFrame          # samples × families
    raw_fine_proportions: pd.DataFrame       # samples × subtypes (joint estimate)
    reconciled_fine_proportions: pd.DataFrame  # samples × subtypes (after reconcile)
    hierarchy_consistency_error: pd.Series   # per-sample ‖π − Aθ‖₂
    broad_reconstruction_error: pd.Series
    fine_reconstruction_error: pd.Series
    convergence_status: list
    n_iterations: list
    runtime: float
    warnings: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)


def _aggregation(subtypes, families, family_map):
    A = np.zeros((len(families), len(subtypes)))
    fidx = {f: i for i, f in enumerate(families)}
    for j, k in enumerate(subtypes):
        A[fidx[family_map.get(k, k)], j] = 1.0
    return A


def _family_signatures(S, A):
    """B[:,f] = mean of member fine signatures (rows of A normalised)."""
    counts = A.sum(axis=1, keepdims=True)
    counts[counts == 0] = 1.0
    return S @ (A.T / counts.T)              # G×F


def _simplex_nnls(M, b):
    """min‖Mθ − b‖² s.t. θ≥0, then renormalise to Σθ=1 (package convention)."""
    theta, _ = nnls(M, b)
    s = theta.sum()
    return theta / s if s > 0 else np.full(M.shape[1], 1.0 / M.shape[1])


def _alt_one(y, S, B, A, alpha, lam_h, lam_2, max_iter, tol, init):
    G, K = S.shape
    F = B.shape[1]
    sa, sh, s2 = np.sqrt(alpha), np.sqrt(lam_h), np.sqrt(lam_2)
    # init θ
    if init == "uniform":
        theta = np.full(K, 1.0 / K)
    else:                                    # flat NNLS init (default) / sequential / broad-first
        theta = _simplex_nnls(S, y)
    pi = A @ theta
    pi = pi / pi.sum() if pi.sum() > 0 else np.full(F, 1.0 / F)
    n_it, conv = 0, False
    for n_it in range(1, max_iter + 1):
        theta_prev = theta
        # π update: min ‖B π − y‖² + λ_h‖π − Aθ‖²  (π≥0, renorm)
        if lam_h > 0:
            Mp = np.vstack([B, sh * np.eye(F)])
            bp = np.concatenate([y, sh * (A @ theta)])
        else:
            Mp, bp = B, y
        pi = _simplex_nnls(Mp, bp)
        # θ update: min α‖Sθ − y‖² + λ_h‖Aθ − π‖² + λ₂‖θ‖²  (θ≥0, renorm)
        blocks_M = [sa * S]
        blocks_b = [sa * y]
        if lam_h > 0:
            blocks_M.append(sh * A); blocks_b.append(sh * pi)
        if lam_2 > 0:
            blocks_M.append(s2 * np.eye(K)); blocks_b.append(np.zeros(K))
        theta = _simplex_nnls(np.vstack(blocks_M), np.concatenate(blocks_b))
        if np.linalg.norm(theta - theta_prev) < tol:
            conv = True
            break
    return pi, theta, n_it, conv


def fit_alternating_soft_hierarchy(
    y: np.ndarray, S: np.ndarray, family_map: Mapping[str, str],
    subtypes, *, alpha: float = 1.0, lambda_h: float = 1.0, lambda_2: float = 0.0,
    max_iter: int = 100, tol: float = 1e-5, init: str = "flat_nnls",
    W_f: Optional[np.ndarray] = None, W_b: Optional[np.ndarray] = None,
    sample_ids=None,
) -> JointHierarchyResult:
    """Alternating constrained solver. ``y`` = samples×G, ``S`` = G×K."""
    import time
    t0 = time.perf_counter()
    y = np.atleast_2d(np.asarray(y, float))
    S = np.asarray(S, float)
    subtypes = [str(c) for c in subtypes]
    families = sorted({str(family_map.get(k, k)) for k in subtypes})
    A = _aggregation(subtypes, families, family_map)
    Wf = np.ones(S.shape[0]) if W_f is None else np.asarray(W_f, float)
    Wb = np.ones(S.shape[0]) if W_b is None else np.asarray(W_b, float)
    Sw, Bw = (np.sqrt(Wf)[:, None] * S), None
    B = _family_signatures(S, A)
    Bw = np.sqrt(Wb)[:, None] * B
    PIs, THs, hcons, brec, frec, conv, nit = [], [], [], [], [], [], []
    for s in range(y.shape[0]):
        ys_f = np.sqrt(Wf) * y[s]
        ys_b = np.sqrt(Wb) * y[s]
        pi, th, n, c = _alt_one(ys_f, Sw, Bw, A, alpha, lambda_h, lambda_2,
                                max_iter, tol, init)
        # reconstruction errors use the (weighted) signatures consistently
        PIs.append(pi); THs.append(th)
        hcons.append(float(np.linalg.norm(pi - A @ th)))
        brec.append(float(np.linalg.norm(ys_b - Bw @ pi)))
        frec.append(float(np.linalg.norm(ys_f - Sw @ th)))
        conv.append(bool(c)); nit.append(int(n))
    idx = list(sample_ids) if sample_ids is not None else [f"s{i}" for i in range(y.shape[0])]
    PI = pd.DataFrame(PIs, index=idx, columns=families)
    TH = pd.DataFrame(THs, index=idx, columns=subtypes)
    # reconcile: scale fine within family so Σ_{k∈f} θ_k == π_f (soft→hard at output)
    TH_rec = reconcile_broad_fine(PI, TH, family_map)
    return JointHierarchyResult(
        broad_proportions=PI, raw_fine_proportions=TH,
        reconciled_fine_proportions=TH_rec,
        hierarchy_consistency_error=pd.Series(hcons, index=idx),
        broad_reconstruction_error=pd.Series(brec, index=idx),
        fine_reconstruction_error=pd.Series(frec, index=idx),
        convergence_status=conv, n_iterations=nit,
        runtime=time.perf_counter() - t0,
        metadata={"algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
                  "solver": "alternating", "alpha": alpha, "lambda_h": lambda_h,
                  "lambda_2": lambda_2, "init": init, "estimate_type": "rna_proportion"})


def fit_joint_soft_hierarchy(y, S, family_map, subtypes, *, alpha=1.0, lambda_h=1.0,
                             lambda_2=0.0, sample_ids=None, **kw) -> JointHierarchyResult:
    """Joint SLSQP solver over [π,θ] with simplex constraints (secondary; slower).
    Falls back to the alternating result's structure."""
    import time
    from scipy.optimize import minimize
    t0 = time.perf_counter()
    y = np.atleast_2d(np.asarray(y, float)); S = np.asarray(S, float)
    subtypes = [str(c) for c in subtypes]
    families = sorted({str(family_map.get(k, k)) for k in subtypes})
    A = _aggregation(subtypes, families, family_map); B = _family_signatures(S, A)
    F, K = len(families), len(subtypes)
    PIs, THs, hcons, brec, frec, conv, nit = [], [], [], [], [], [], []
    for s in range(y.shape[0]):
        ys = y[s]
        def obj(z):
            pi, th = z[:F], z[F:]
            return (np.sum((ys - B @ pi) ** 2) + alpha * np.sum((ys - S @ th) ** 2)
                    + lambda_h * np.sum((pi - A @ th) ** 2) + lambda_2 * np.sum(th ** 2))
        th0 = _simplex_nnls(S, ys); z0 = np.concatenate([A @ th0, th0])
        cons = [{"type": "eq", "fun": lambda z: z[:F].sum() - 1},
                {"type": "eq", "fun": lambda z: z[F:].sum() - 1}]
        bnds = [(0, None)] * (F + K)
        r = minimize(obj, z0, method="SLSQP", bounds=bnds, constraints=cons,
                     options={"maxiter": 200, "ftol": 1e-8})
        pi, th = np.clip(r.x[:F], 0, None), np.clip(r.x[F:], 0, None)
        pi = pi / pi.sum() if pi.sum() > 0 else np.full(F, 1 / F)
        th = th / th.sum() if th.sum() > 0 else np.full(K, 1 / K)
        PIs.append(pi); THs.append(th); hcons.append(float(np.linalg.norm(pi - A @ th)))
        brec.append(float(np.linalg.norm(ys - B @ pi))); frec.append(float(np.linalg.norm(ys - S @ th)))
        conv.append(bool(r.success)); nit.append(int(r.nit))
    idx = list(sample_ids) if sample_ids is not None else [f"s{i}" for i in range(y.shape[0])]
    PI = pd.DataFrame(PIs, index=idx, columns=families); TH = pd.DataFrame(THs, index=idx, columns=subtypes)
    return JointHierarchyResult(
        broad_proportions=PI, raw_fine_proportions=TH,
        reconciled_fine_proportions=reconcile_broad_fine(PI, TH, family_map),
        hierarchy_consistency_error=pd.Series(hcons, index=idx),
        broad_reconstruction_error=pd.Series(brec, index=idx),
        fine_reconstruction_error=pd.Series(frec, index=idx),
        convergence_status=conv, n_iterations=nit, runtime=time.perf_counter() - t0,
        metadata={"algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
                  "solver": "joint_slsqp", "alpha": alpha, "lambda_h": lambda_h, "lambda_2": lambda_2})


def reconcile_broad_fine(broad: pd.DataFrame, fine: pd.DataFrame,
                         family_map: Mapping[str, str]) -> pd.DataFrame:
    """Scale fine within each family so Σ_{k∈f} θ_k == π_f (mass-consistent output)."""
    out = fine.copy().astype(float)
    fam_of = {c: family_map.get(c, c) for c in fine.columns}
    for fam in broad.columns:
        mem = [c for c in fine.columns if fam_of[c] == fam]
        if not mem:
            continue
        s = fine[mem].sum(axis=1)
        scale = np.where(s > 0, broad[fam] / s.replace(0, np.nan), 0.0)
        for c in mem:
            out[c] = fine[c] * pd.Series(scale, index=fine.index).fillna(0.0)
        # if family signal is zero, distribute uniformly
        zero = s <= 0
        if zero.any():
            for c in mem:
                out.loc[zero, c] = broad.loc[zero, fam] / len(mem)
    return out


def compute_hierarchy_consistency(broad, fine, family_map) -> pd.Series:
    A_fine = fine.copy()
    fam_of = {c: family_map.get(c, c) for c in fine.columns}
    agg = pd.DataFrame({f: A_fine[[c for c in fine.columns if fam_of[c] == f]].sum(axis=1)
                        for f in broad.columns})
    return (broad - agg).pow(2).sum(axis=1).pow(0.5)


def compute_error_decomposition(truth_fine, broad_pred, fine_pred, family_map) -> dict:
    """Per-sample error components vs truth (RMSE units)."""
    cols = [c for c in truth_fine.columns]
    fam_of = {c: family_map.get(c, c) for c in cols}
    fams = sorted(set(fam_of.values()))
    tfam = pd.DataFrame({f: truth_fine[[c for c in cols if fam_of[c] == f]].sum(1) for f in fams})
    pfam = broad_pred.reindex(columns=fams).fillna(0.0)
    broad_err = float(np.sqrt(((tfam - pfam) ** 2).to_numpy().mean()))
    # conditional within-family error
    ct, cp = truth_fine * 0.0, fine_pred.reindex(columns=cols).fillna(0.0) * 0.0
    fp = fine_pred.reindex(columns=cols).fillna(0.0)
    for f in fams:
        mem = [c for c in cols if fam_of[c] == f]
        st_, sp_ = truth_fine[mem].sum(1), fp[mem].sum(1)
        for c in mem:
            ct[c] = np.where(st_ > 0, truth_fine[c] / st_.replace(0, np.nan), 0)
            cp[c] = np.where(sp_ > 0, fp[c] / sp_.replace(0, np.nan), 0)
    cond_err = float(np.sqrt(((ct.fillna(0) - cp.fillna(0)) ** 2).to_numpy().mean()))
    total_err = float(np.sqrt(((truth_fine[cols] - fp[cols]) ** 2).to_numpy().mean()))
    consist = float(compute_hierarchy_consistency(pfam, fp, family_map).mean())
    return {"broad_error": broad_err, "conditional_fine_error": cond_err,
            "consistency_error": consist, "total_fine_error": total_err}
