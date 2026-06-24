"""Experimental in-solver state-similarity regularization (Option A; opt-in).

A **separate** experimental optimiser (the production NB-CAR ``SpatCARModel`` is NOT
modified). It minimises an explicit joint objective by projected gradient descent on
a Gaussian (least-squares) surrogate of the reconstruction term:

    min_Theta  ||Y_norm - Theta @ R||^2
             + lambda_spatial * tr(Theta^T L_sym Theta)        # spot smoothing
             + lambda_state   * P_state(Theta)                  # state coupling
             + lambda_sparse  * sum_ik Theta_ik log Theta_ik    # neg-entropy sparsity
        s.t. Theta >= 0, rows on the simplex (or per-family mass)

State-penalty variants:
  * ``"competition"`` (default): tr(Theta A_state Theta^T) — discourages assigning
    mass to several *similar* states in the same spot (targets diffuse effective-N /
    false-positive subtype mass);
  * ``"laplacian"``: tr(Theta L_state Theta^T) — smooths abundances of similar
    states across spots (can worsen rare niches).

Inspired by the *concept* of state-aware regularization; NOT Redeconve, not
Redeconve-equivalent, no external code copied. See
``docs/IN_SOLVER_STATE_REGULARIZATION_REPORT.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.sparse as sp

from tissueresolve.experimental.state_similarity_regularization import StateSimilarityGraph

__all__ = ["StateRegularizedSolverResult", "fit_state_regularized_spatial", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"


@dataclass
class StateRegularizedSolverResult:
    """Output of :func:`fit_state_regularized_spatial`."""

    theta: np.ndarray                 # (N, K) non-negative, valid compositions
    cell_types: list
    n_iter: int
    converged: bool
    loss_trace: list
    recon_trace: list
    final_loss: float
    recon_loss: float
    spatial_penalty: float
    state_penalty: float
    sparsity_penalty: float
    loss_monotonic: bool
    metadata: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Projections
# ---------------------------------------------------------------------------


def _project_simplex_rows(V: np.ndarray, z) -> np.ndarray:
    """Euclidean projection of each row of ``V`` onto {x >= 0, sum(x) = z}.

    ``z`` is a scalar or a per-row ``(N,)`` / ``(N, 1)`` array of positive targets.
    Duchi et al. (2008), vectorised over rows.
    """
    V = np.asarray(V, dtype=np.float64)
    n, k = V.shape
    z_arr = np.full((n, 1), float(z)) if np.isscalar(z) else np.asarray(z, dtype=np.float64).reshape(n, 1)
    z_arr = np.maximum(z_arr, 1e-12)
    U = np.sort(V, axis=1)[:, ::-1]
    css = np.cumsum(U, axis=1) - z_arr
    ind = np.arange(1, k + 1)
    cond = U - css / ind > 0
    rho = np.where(cond.any(axis=1), cond.sum(axis=1), 1)
    theta = css[np.arange(n), rho - 1] / rho
    return np.maximum(V - theta[:, None], 0.0)


def _apply_rare_protection(theta, cell_types, rare_protection, z):
    """Floor protected state columns at their per-state minimum, re-balance the rest.

    ``rare_protection`` maps ``state -> min fraction`` (of the row mass ``z``). Only
    spots where the floor exceeds the current value are adjusted; the remaining
    (non-protected) mass is rescaled so each row still sums to ``z``.
    """
    if not rare_protection:
        return theta
    idx = {c: i for i, c in enumerate(cell_types)}
    z_arr = np.full(theta.shape[0], float(z)) if np.isscalar(z) else np.asarray(z).reshape(-1)
    prot_cols = {idx[s]: float(f) for s, f in rare_protection.items() if s in idx}
    if not prot_cols:
        return theta
    out = theta.copy()
    cols = list(prot_cols)
    floor = np.zeros((theta.shape[0], len(cols)))
    for j, c in enumerate(cols):
        floor[:, j] = prot_cols[c] * z_arr
    cur = out[:, cols]
    new_prot = np.maximum(cur, floor)
    prot_total = new_prot.sum(axis=1)
    out[:, cols] = new_prot
    non_mask = np.ones(theta.shape[1], dtype=bool)
    non_mask[cols] = False
    non = out[:, non_mask]
    non_sum = non.sum(axis=1)
    target_non = np.maximum(z_arr - prot_total, 0.0)
    scale = np.divide(target_non, np.maximum(non_sum, 1e-12))
    out[:, non_mask] = non * scale[:, None]
    # if a row's protected floors exceeded z (degenerate), renormalise whole row
    bad = prot_total > z_arr + 1e-9
    if bad.any():
        rs = out[bad].sum(axis=1, keepdims=True)
        out[bad] = out[bad] / np.maximum(rs, 1e-12) * z_arr[bad, None]
    return out


# ---------------------------------------------------------------------------
# Penalty matrices
# ---------------------------------------------------------------------------


def _state_matrix(state_graph: StateSimilarityGraph, cell_types: list, kind: str) -> np.ndarray:
    """Build the (K, K) state-penalty matrix aligned to ``cell_types`` order.

    ``kind="competition"`` → symmetric adjacency ``A`` (off-diagonal weights);
    ``kind="laplacian"``   → ``L = diag(A.sum) - A``.
    """
    K = len(cell_types)
    A = np.zeros((K, K), dtype=np.float64)
    if state_graph is not None and state_graph.n_edges:
        pos = {c: i for i, c in enumerate(cell_types)}
        for e in range(state_graph.n_edges):
            a = state_graph.cell_types[int(state_graph.state_i[e])]
            b = state_graph.cell_types[int(state_graph.state_j[e])]
            if a in pos and b in pos:
                w = float(state_graph.weights[e])
                A[pos[a], pos[b]] = w
                A[pos[b], pos[a]] = w
    if kind == "competition":
        return A
    if kind == "laplacian":
        return np.diag(A.sum(axis=1)) - A
    raise ValueError(f"Unknown state_penalty {kind!r}; use 'competition' or 'laplacian'.")


def _spatial_laplacian_sym(spot_graph) -> Optional[sp.csr_matrix]:
    """Symmetric spot Laplacian ``L_sym = I - (A + A^T)/2`` from the graph adjacency."""
    if spot_graph is None:
        return None
    A = spot_graph.A
    n = A.shape[0]
    As = (A + A.T) * 0.5
    return (sp.identity(n, format="csr") - As).tocsr()


# ---------------------------------------------------------------------------
# Objective + gradient
# ---------------------------------------------------------------------------


def _objective(theta, Yn, R, Lsym, M, lam_spatial, lam_state, lam_sparse, kind):
    resid = theta @ R - Yn
    recon = float(np.sum(resid * resid))
    spatial = 0.0
    if Lsym is not None and lam_spatial > 0:
        spatial = float(np.sum(theta * (Lsym @ theta)))
    state = 0.0
    if M is not None and lam_state > 0:
        state = float(np.sum(theta * (theta @ M)))  # tr(Θ M Θᵀ) for both kinds
    sparse = 0.0
    if lam_sparse > 0:
        # Shannon entropy H = -Σ Θ logΘ (≥0). Adding +λ·H to the minimised objective
        # drives H down → concentration → lower effective-N (a smooth sparsity prior,
        # not hard thresholding, so rare states are not destroyed).
        sparse = -float(np.sum(theta * np.log(theta + 1e-12)))
    total = recon + lam_spatial * spatial + lam_state * state + lam_sparse * sparse
    return total, recon, lam_spatial * spatial, lam_state * state, lam_sparse * sparse


def _gradient(theta, Yn, R, Rt, Lsym, M, lam_spatial, lam_state, lam_sparse):
    g = 2.0 * ((theta @ R - Yn) @ Rt)
    if Lsym is not None and lam_spatial > 0:
        g = g + 2.0 * lam_spatial * (Lsym @ theta)
    if M is not None and lam_state > 0:
        g = g + 2.0 * lam_state * (theta @ M)
    if lam_sparse > 0:
        # ∇(λ·H) = λ·∇(-Σ Θ logΘ) = -λ(logΘ + 1)
        g = g - lam_sparse * (np.log(theta + 1e-12) + 1.0)
    return g


# ---------------------------------------------------------------------------
# Public solver
# ---------------------------------------------------------------------------


def fit_state_regularized_spatial(
    Y: np.ndarray,
    reference_profiles: np.ndarray,
    spot_graph=None,
    state_graph: Optional[StateSimilarityGraph] = None,
    cell_types: Optional[list] = None,
    init_theta: Optional[np.ndarray] = None,
    lib_sizes: Optional[np.ndarray] = None,
    lambda_spatial: float = 0.02,
    lambda_state: float = 0.01,
    lambda_sparse: float = 0.001,
    state_penalty: str = "competition",
    within_family_only: bool = True,
    preserve_broad_mass: bool = False,
    family_map: Optional[dict] = None,
    rare_protection: Optional[dict] = None,
    max_iter: int = 500,
    tol: float = 1e-5,
    optimizer: str = "projected_gradient",
) -> StateRegularizedSolverResult:
    """Projected-gradient optimiser for the joint state-regularized objective.

    ``Y`` is ``(N, G)`` counts (normalised internally by ``lib_sizes`` or row sums);
    ``reference_profiles`` is ``(K, G)`` in proportion scale (rows ~ sum to 1, e.g.
    CPM/1e6). ``init_theta`` is an optional ``(N, K)`` warm start (e.g. the
    production NB-CAR fit). Returns a :class:`StateRegularizedSolverResult` whose
    ``theta`` satisfies ``theta >= 0`` and valid compositions (or conserved
    per-family mass when ``preserve_broad_mass``).
    """
    if optimizer != "projected_gradient":
        raise ValueError(f"Unsupported optimizer {optimizer!r} (only 'projected_gradient').")

    Y = np.asarray(Y, dtype=np.float64)
    R = np.asarray(reference_profiles, dtype=np.float64)
    N, G = Y.shape
    K = R.shape[0]
    if R.shape[1] != G:
        raise ValueError(f"reference_profiles has {R.shape[1]} genes but Y has {G}.")
    cell_types = list(cell_types) if cell_types is not None else [f"state_{k}" for k in range(K)]

    lib = (np.asarray(lib_sizes, dtype=np.float64) if lib_sizes is not None
           else np.maximum(Y.sum(axis=1), 1.0))
    Yn = Y / np.maximum(lib[:, None], 1e-12)

    # init
    if init_theta is not None:
        theta = np.asarray(init_theta, dtype=np.float64).copy()
    else:
        theta = np.full((N, K), 1.0 / K, dtype=np.float64)
    theta = np.maximum(theta, 1e-10)

    # per-family targets for preserve_broad_mass
    fam_blocks = None
    z_target = 1.0
    if preserve_broad_mass and family_map is not None:
        fam_of = [str(family_map.get(c, c)) for c in cell_types]
        fams = {}
        for j, f in enumerate(fam_of):
            fams.setdefault(f, []).append(j)
        fam_blocks = list(fams.values())
        # baseline per-family mass per spot (from init)
        z_target = {tuple(cols): theta[:, cols].sum(axis=1, keepdims=True) for cols in fam_blocks}

    Rt = R.T.copy()
    Lsym = _spatial_laplacian_sym(spot_graph)
    M = _state_matrix(state_graph, cell_types, state_penalty) if (
        state_graph is not None and lambda_state > 0) else None

    def project(th):
        if fam_blocks is not None:
            out = th.copy()
            for cols in fam_blocks:
                out[:, cols] = _project_simplex_rows(th[:, cols], z_target[tuple(cols)])
            z_full = None
        else:
            out = _project_simplex_rows(th, 1.0)
            z_full = 1.0
        if rare_protection:
            zz = 1.0 if fam_blocks is None else out.sum(axis=1)
            out = _apply_rare_protection(out, cell_types, rare_protection, zz)
        return np.maximum(out, 0.0)

    theta = project(theta)
    f, recon, sp_pen, st_pen, sps_pen = _objective(
        theta, Yn, R, Lsym, M, lambda_spatial, lambda_state, lambda_sparse, state_penalty)
    loss_trace = [f]
    recon_trace = [recon]

    eta = 1.0
    converged = False
    monotonic = True
    n_iter = 0
    for it in range(1, max_iter + 1):
        n_iter = it
        g = _gradient(theta, Yn, R, Rt, Lsym, M, lambda_spatial, lambda_state, lambda_sparse)
        # backtracking projected-gradient: ensure non-increasing objective
        stepped = False
        for _ in range(30):
            cand = project(theta - eta * g)
            fc, rc, spc, stc, spsc = _objective(
                cand, Yn, R, Lsym, M, lambda_spatial, lambda_state, lambda_sparse, state_penalty)
            if fc <= f + 1e-12:
                theta, f_new, recon = cand, fc, rc
                sp_pen, st_pen, sps_pen = spc, stc, spsc
                eta *= 1.2
                stepped = True
                break
            eta *= 0.5
        if not stepped:
            break  # cannot decrease further → converged at a stationary point
        rel = (f - f_new) / max(abs(f), 1e-12)
        f = f_new
        loss_trace.append(f)
        recon_trace.append(recon)
        if len(loss_trace) >= 2 and loss_trace[-1] > loss_trace[-2] + 1e-9:
            monotonic = False
        if rel < tol or eta < 1e-12:
            converged = True
            break

    metadata = {
        "feature_status": FEATURE_STATUS,
        "optimizer": optimizer,
        "state_penalty": state_penalty,
        "within_family_only": bool(within_family_only),
        "preserve_broad_mass": bool(preserve_broad_mass),
        "lambda_spatial": float(lambda_spatial),
        "lambda_state": float(lambda_state),
        "lambda_sparse": float(lambda_sparse),
        "n_iter": n_iter,
        "converged": bool(converged),
        "loss_monotonic": bool(monotonic),
        "rare_protection_used": bool(rare_protection),
        "state_graph_edges": int(state_graph.n_edges) if state_graph is not None else 0,
        "n_states": K,
        "n_spots": N,
    }
    return StateRegularizedSolverResult(
        theta=theta.astype(np.float32), cell_types=cell_types, n_iter=n_iter,
        converged=converged, loss_trace=loss_trace, recon_trace=recon_trace,
        final_loss=float(f), recon_loss=float(recon), spatial_penalty=float(sp_pen),
        state_penalty=float(st_pen), sparsity_penalty=float(sps_pen),
        loss_monotonic=monotonic, metadata=metadata,
    )
