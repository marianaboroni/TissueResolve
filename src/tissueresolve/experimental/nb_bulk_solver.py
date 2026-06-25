"""Experimental count-likelihood (Poisson / negative-binomial) bulk solver (P1).

Opt-in, non-default. Tests whether a noise model matched to RNA-seq counts beats
the default weighted-NNLS (Gaussian/L2 on L1-normalised profiles), especially for
low-abundance and within-family subtypes. Outputs remain **RNA-derived (mRNA)
proportions**, never cell fractions; the default wNNLS solver is unchanged.

Model (per sample n), design ``X_{gk} = ℓ_n · Φ_{gk}`` with Φ in proportion scale
(``CPM/1e6``), expected counts ``μ_{gn} = Σ_k X_{gk} θ_{kn}``:

    loss="poisson":  min_θ≥0  Σ_g w_g [ μ_g − y_g log μ_g ]
    loss="nb":       min_θ≥0  Σ_g w_g [ (y_g+φ_g) log(μ_g+φ_g) − y_g log μ_g ]

Solved by a multiplicative majorise–minimise (MM) update — monotone non-increasing
NLL, keeps ``θ ≥ 0`` — then renormalised to the simplex (mRNA proportions):

    NB:       θ_k ← θ_k · [Σ_g w_g X_{gk} y_g/μ_g] / [Σ_g w_g X_{gk}(y_g+φ_g)/(μ_g+φ_g)]
    Poisson:  θ_k ← θ_k · [Σ_g w_g X_{gk} y_g/μ_g] / [Σ_g w_g X_{gk}]

The MM form mirrors the spatial NB update mathematically but is re-derived here in a
separate module — the production spatial and bulk solvers are not modified.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.results import BulkDeconvResult, ReferenceSignature

__all__ = ["fit_bulk_nb_glm", "NBGLMBulkSolver", "NBGLMResult", "FEATURE_STATUS"]

logger = logging.getLogger("tissueresolve.experimental.nb_bulk_solver")
FEATURE_STATUS = "experimental"
_MIN_GENES_WARN = 20
_EPS = 1e-10


@dataclass
class NBGLMResult:
    """Array-level result of :func:`fit_bulk_nb_glm`."""

    theta: np.ndarray          # (K, N) mRNA proportions (columns sum to 1)
    theta_raw: np.ndarray      # (K, N) unnormalised non-negative abundances
    n_iter: int
    converged: bool
    loss_trace: list
    final_loss: float
    loss_monotonic: bool
    metadata: dict = field(default_factory=dict)


def _nll(Y, Phi, Theta, lib, w, phi_g, loss):
    """θ-dependent part of the (weighted) Poisson / NB negative log-likelihood."""
    M_nolib = Phi @ Theta                      # (G, N)
    mu = np.maximum(lib[None, :] * M_nolib, _EPS)
    if loss == "poisson":
        per_g = mu - Y * np.log(mu)
    else:
        per_g = (Y + phi_g[:, None]) * np.log(mu + phi_g[:, None]) - Y * np.log(mu)
    return float(np.sum(w[:, None] * per_g))


def fit_bulk_nb_glm(
    bulk_counts: np.ndarray,
    reference_profiles: np.ndarray,
    library_sizes: Optional[np.ndarray] = None,
    gene_dispersion: Optional[np.ndarray] = None,
    init_theta: Optional[np.ndarray] = None,
    gene_weights: Optional[np.ndarray] = None,
    max_iter: int = 500,
    tol: float = 1e-6,
    loss: str = "nb",
    optimizer: str = "multiplicative_update",
) -> NBGLMResult:
    """Fit Poisson / NB bulk deconvolution by multiplicative MM updates.

    Parameters
    ----------
    bulk_counts:
        ``(G, N)`` non-negative counts (genes × samples), panel-aligned to
        ``reference_profiles``.
    reference_profiles:
        ``(G, K)`` reference in proportion scale (``CPM/1e6``; columns ~ per-type
        gene distributions).
    library_sizes:
        ``(N,)`` per-sample library sizes. ``None`` → column sums of ``bulk_counts``.
    gene_dispersion:
        ``(G,)`` NB dispersion ``φ_g``. ``None`` with ``loss="nb"`` → Poisson
        fallback (recorded).
    gene_weights:
        ``(G,)`` non-negative per-gene weights. ``None`` → uniform.
    loss:
        ``"poisson"`` or ``"nb"``.
    optimizer:
        Only ``"multiplicative_update"`` is implemented.
    """
    if optimizer != "multiplicative_update":
        raise ValueError(f"Unsupported optimizer {optimizer!r} (only 'multiplicative_update').")
    if loss not in ("poisson", "nb"):
        raise ValueError(f"Unsupported loss {loss!r} (use 'poisson' or 'nb').")

    Y = np.asarray(bulk_counts, dtype=np.float64)
    Phi = np.asarray(reference_profiles, dtype=np.float64)
    G, N = Y.shape
    if Phi.shape[0] != G:
        raise ValueError(f"reference_profiles has {Phi.shape[0]} genes but bulk has {G}.")
    K = Phi.shape[1]

    lib = (np.asarray(library_sizes, dtype=np.float64) if library_sizes is not None
           else np.maximum(Y.sum(axis=0), 1.0))
    lib = np.maximum(lib, 1.0)
    w = (np.clip(np.asarray(gene_weights, dtype=np.float64), 0, None) if gene_weights is not None
         else np.ones(G))
    if w.sum() > 0:
        w = w / w.mean()

    dispersion_source = "none"
    effective_loss = loss
    if loss == "nb":
        if gene_dispersion is not None:
            phi_g = np.clip(np.asarray(gene_dispersion, dtype=np.float64), 1e-6, None)
            dispersion_source = "provided"
        else:
            effective_loss = "poisson"
            dispersion_source = "fallback_poisson"
            warnings.warn(
                "fit_bulk_nb_glm: loss='nb' requested but no gene_dispersion was "
                "provided; falling back to Poisson (recorded as "
                "dispersion_source='fallback_poisson'). Build the reference with "
                "estimate_overdispersion=True for a true NB fit.", stacklevel=2)
            phi_g = np.full(G, np.inf)
    else:
        phi_g = np.full(G, np.inf)

    if init_theta is not None:
        Theta = np.asarray(init_theta, dtype=np.float64).copy()
    elif effective_loss == "nb":
        # NB warm start from a Poisson fit: with small dispersion the NB
        # multiplicative ratio flattens to ≈1 (uninformative from a uniform init),
        # so we start NB from the well-behaved Poisson optimum and let it refine.
        pois = fit_bulk_nb_glm(Y, Phi, library_sizes=lib, gene_weights=gene_weights,
                               max_iter=max_iter, tol=tol, loss="poisson")
        Theta = pois.theta_raw.copy()
    else:
        Theta = np.full((K, N), 1.0 / K)
    Theta = np.maximum(Theta, _EPS)

    # Poisson denominator is constant: Σ_g w_g X_gk = lib_n · Σ_g w_g Φ_gk
    wPhi = (w[:, None] * Phi).sum(axis=0)       # (K,)
    loss_trace = [_nll(Y, Phi, Theta, lib, w, phi_g, effective_loss)]
    converged = False
    monotonic = True
    n_iter = 0
    for it in range(1, max_iter + 1):
        n_iter = it
        M_nolib = np.maximum(Phi @ Theta, _EPS)             # (G, N)
        mu = np.maximum(lib[None, :] * M_nolib, _EPS)        # (G, N)
        # numerator (K, N): lib cancels → Φᵀ (w·Y/M_nolib)
        num = Phi.T @ (w[:, None] * Y / M_nolib)
        if effective_loss == "poisson":
            den = wPhi[:, None] * lib[None, :]               # (K, N)
        else:
            ratio = (Y + phi_g[:, None]) / (mu + phi_g[:, None])
            den = lib[None, :] * (Phi.T @ (w[:, None] * ratio))
        Theta = np.maximum(Theta * num / np.maximum(den, _EPS), _EPS)
        f = _nll(Y, Phi, Theta, lib, w, phi_g, effective_loss)
        loss_trace.append(f)
        if f > loss_trace[-2] + 1e-6 * max(abs(loss_trace[-2]), 1.0):
            monotonic = False
        rel = (loss_trace[-2] - f) / max(abs(loss_trace[-2]), 1e-12)
        if abs(rel) < tol:
            converged = True
            break

    theta_raw = Theta.copy()
    col = Theta.sum(axis=0, keepdims=True)
    theta = Theta / np.where(col > 0, col, 1.0)
    # pathological all-zero column → uniform
    bad = (col.ravel() <= _EPS)
    if bad.any():
        theta[:, bad] = 1.0 / K

    metadata = {
        "feature_status": FEATURE_STATUS,
        "solver": f"{loss}_glm_experimental",
        "likelihood": effective_loss,
        "optimizer": optimizer,
        "n_iter": n_iter,
        "converged": bool(converged),
        "loss_monotonic": bool(monotonic),
        "final_loss": float(loss_trace[-1]),
        "n_genes": int(G),
        "n_samples": int(N),
        "n_cell_types": int(K),
        "dispersion_source": dispersion_source,
        "fallback_to_poisson": effective_loss != loss,
        "library_size_normalization": True,
        "estimate_type": "mRNA_proportion",
    }
    return NBGLMResult(theta=theta, theta_raw=theta_raw, n_iter=n_iter,
                       converged=converged, loss_trace=loss_trace,
                       final_loss=float(loss_trace[-1]), loss_monotonic=monotonic,
                       metadata=metadata)


class NBGLMBulkSolver:
    """Pipeline-compatible wrapper returning a :class:`BulkDeconvResult`.

    Mirrors ``WNNLSSolver.solve`` (same inputs/outputs) so reports and the
    pipeline consume it unchanged. ``loss`` selects Poisson vs NB.
    """

    def __init__(self, *, loss: str = "nb", seed: int = 42, max_iter: int = 500,
                 tol: float = 1e-6) -> None:
        self.loss = loss
        self.seed = seed
        self.max_iter = max_iter
        self.tol = tol

    def solve(self, bulk: pd.DataFrame, ref: ReferenceSignature, gene_panel: list,
              gene_weights: Optional[pd.Series] = None) -> BulkDeconvResult:
        common = sorted(set(gene_panel) & set(bulk.index) & set(ref.gene_names))
        if not common:
            raise ValueError(
                "No genes overlap between bulk, reference, and gene_panel. "
                "Check gene ID conventions (symbol vs Ensembl, capitalisation).")
        if len(common) < _MIN_GENES_WARN:
            warnings.warn(
                f"NBGLMBulkSolver: only {len(common)} genes in the final panel "
                f"(reliable threshold: {_MIN_GENES_WARN}). Results may be unreliable.",
                stacklevel=2)

        # panel-aligned arrays (no silent gene drop: 'common' is the explicit overlap)
        ct = list(ref.cell_types)
        r_idx = {g: i for i, g in enumerate(ref.gene_names)}
        R_cpm = ref.as_R_cpm()                               # (K, G_ref)
        panel_cols = [r_idx[g] for g in common]
        Phi = (R_cpm[:, panel_cols].T / 1e6)                 # (G, K) proportion scale
        Y = bulk.loc[common].to_numpy(dtype=np.float64)      # (G, N)

        phi_g = None
        if self.loss == "nb" and ref.phi_g is not None:
            phi_g = np.asarray(ref.phi_g, dtype=np.float64)[panel_cols]

        w = (gene_weights.reindex(common).fillna(1.0).to_numpy(float)
             if gene_weights is not None else None)

        res = fit_bulk_nb_glm(Y, Phi, gene_dispersion=phi_g, gene_weights=w,
                              max_iter=self.max_iter, tol=self.tol, loss=self.loss)

        sample_ids = list(bulk.columns)
        props = pd.DataFrame(res.theta.T, index=sample_ids, columns=ct)
        props.index.name = "sample"

        # reconstruction R² on L1-normalised profiles (schema parity with wNNLS)
        Phi_norm = Phi / np.where(Phi.sum(0, keepdims=True) > 0, Phi.sum(0, keepdims=True), 1.0)
        Bn = Y / np.where(Y.sum(0, keepdims=True) > 0, Y.sum(0, keepdims=True), 1.0)
        recon = Phi_norm @ res.theta                          # (G, N)
        r2 = np.empty(len(sample_ids))
        for n in range(len(sample_ids)):
            y, yh = Bn[:, n], recon[:, n]
            ss_res = float(np.dot(y - yh, y - yh))
            ss_tot = float(np.dot(y - y.mean(), y - y.mean()))
            r2[n] = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        run_metadata = dict(res.metadata)
        run_metadata["n_genes_panel"] = len(common)
        logger.info("NB/Poisson GLM bulk: %d samples, %d types, %d genes, likelihood=%s, "
                    "mean R²=%.4f.", len(sample_ids), len(ct), len(common),
                    res.metadata["likelihood"], float(r2.mean()))
        return BulkDeconvResult(
            proportions=props,
            coverage_r2=pd.Series(r2, index=sample_ids, name="coverage_r2"),
            gene_panel=common,
            gene_weights=gene_weights.reindex(common) if gene_weights is not None else None,
            run_metadata=run_metadata,
        )
