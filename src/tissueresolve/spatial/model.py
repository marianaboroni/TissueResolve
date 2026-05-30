"""
SpatCAR spatial deconvolution model for TissueResolve.

Algorithm
---------
Alternating block coordinate descent over spot-level proportion matrix Π (N × K):

1. **Multiplicative NB update** per batch::

       π_{sk}^{NB} ∝ π_{sk} · (Σ_g R^d_{kg} y_{sg} / Mu_no_lib_{sg})
                              / (Σ_g R^d_{kg} (y_{sg}+φ_g) / (Mu_{sg}+φ_g))

   Guaranteed non-negative and non-increasing in NB-NLL per step.

2. **Spatial proximal mixing** (blending with row-normalised neighbour mean)::

       π_{s}^{mix} = (1 − α) π_s^{NB} + α (A π)_s

   where α = min(0.5, λ) and A is the row-normalised adjacency.
   **Smoothing is always recorded** via ``lambda_spatial`` and ``alpha``.

3. **Mismatch factor update** (every ``update_mismatch_every`` iterations).

Non-convergence warning
-----------------------
If the model does not converge within ``max_iter`` iterations, a
non-suppressible ``UserWarning`` is emitted in addition to the INFO log.

Estimate type
-------------
All fitted proportions are **spot-level RNA-derived cellular composition
estimates**, not direct single-cell counts.  This is enforced in
:class:`~tissueresolve.results.SpatialDeconvResult`.
"""
from __future__ import annotations

import logging
import time
import warnings
from typing import Optional

import numpy as np
import scipy.optimize as opt
import scipy.sparse as sp

from tissueresolve.protocol.mismatch import SpatialMismatch, compute_spatial_discordance, update_mismatch_factors
from tissueresolve.results import ReferenceSignature
from tissueresolve.spatial.graph import SpatialGraph
from tissueresolve.utils import project_simplex_batch, set_random_state

__all__ = ["SpatCARModel"]

logger = logging.getLogger("tissueresolve.spatial.model")


class SpatCARModel:
    """Spatial CAR deconvolution model (array-based interface).

    Parameters
    ----------
    lambda_spatial:
        Spatial regularisation strength λ ≥ 0.  ``0`` = pure NB-MAP (no spatial
        smoothing).  **Always recorded.** Never smoothed silently.
    max_iter:
        Maximum outer optimisation iterations.
    tol:
        Convergence tolerance on ``‖ΔΠ‖_F / N``.
    update_mismatch_every:
        Update mismatch scale factors every this many iterations.
    n_jobs:
        Parallel workers for NNLS warm start.  ``1`` = sequential.
    random_state:
        Integer seed for reproducibility.
    verbose:
        Show tqdm progress bar.

    Fitted attributes (available after :meth:`fit`)
    ------------------------------------------------
    proportions_   : (N, K) float32
    cell_types_    : list[str]
    marker_genes_  : list[str]
    n_iter_        : int
    converged_     : bool
    convergence_trace_ : list[float]
    lambda_spatial, alpha  : smoothing parameters recorded always
    """

    def __init__(
        self,
        *,
        lambda_spatial: float = 0.1,
        max_iter: int = 200,
        tol: float = 5e-5,
        update_mismatch_every: int = 5,
        n_jobs: int = 1,
        random_state: int = 42,
        verbose: bool = False,
    ) -> None:
        self.lambda_spatial = lambda_spatial
        self.max_iter = max_iter
        self.tol = tol
        self.update_mismatch_every = update_mismatch_every
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose

        # Fitted state
        self.proportions_: Optional[np.ndarray] = None
        self.cell_types_: list[str] = []
        self.marker_genes_: list[str] = []
        self.n_iter_: int = 0
        self.converged_: bool = False
        self.convergence_trace_: list[float] = []
        self.alpha: float = 0.0        # smoothing mix weight — always recorded

        # Stored for bootstrap
        self._R_lin: Optional[np.ndarray] = None
        self._phi_g: Optional[np.ndarray] = None
        self._lib_sizes: Optional[np.ndarray] = None
        self._Y_marker: Optional[np.ndarray] = None
        self._mismatch: Optional[SpatialMismatch] = None
        self._graph_ref: Optional[SpatialGraph] = None
        self._rng: Optional[np.random.Generator] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fit(
        self,
        Y_marker: np.ndarray,
        ref_marker: ReferenceSignature,
        graph: SpatialGraph,
        lib_sizes: np.ndarray,
        *,
        mismatch: Optional[SpatialMismatch] = None,
    ) -> "SpatCARModel":
        """Fit the SpatCAR model.

        Parameters
        ----------
        Y_marker:
            Dense marker-gene count matrix, shape ``(N, G_m)`` float32.
        ref_marker:
            :class:`~tissueresolve.results.ReferenceSignature` **already
            subsetted** to the G_m marker genes.
        graph:
            :class:`~tissueresolve.spatial.graph.SpatialGraph`.
        lib_sizes:
            Per-spot library sizes, shape ``(N,)`` float32.
        mismatch:
            Pre-computed :class:`~tissueresolve.protocol.mismatch.SpatialMismatch`.
            When ``None``, discordance is computed from ``Y_marker`` and
            ``ref_marker`` automatically.

        Returns
        -------
        SpatCARModel
            ``self`` (fitted in-place).
        """
        t0 = time.perf_counter()
        self._rng = set_random_state(self.random_state)
        self._graph_ref = graph

        N, G_m = Y_marker.shape
        K = ref_marker.n_cell_types
        self.cell_types_ = list(ref_marker.cell_types)
        self.marker_genes_ = list(ref_marker.gene_names)

        logger.info(
            "Fitting SpatCAR: N=%d spots, K=%d cell types, G_m=%d marker genes, λ=%.3f.",
            N, K, G_m, self.lambda_spatial,
        )

        # Reference in proportion scale: R_lin[k, g] = CPM/1e6
        R_cpm = ref_marker.as_R_cpm()          # (K, G_m)
        R_lin = (R_cpm / 1e6).astype(np.float32)

        phi_g = ref_marker.phi_g
        if phi_g is None:
            # Fallback: uniform dispersion when not estimated
            phi_g = np.ones(G_m, dtype=np.float32)
            warnings.warn(
                "SpatCARModel.fit: phi_g (NB dispersion) not present in "
                "ReferenceSignature.  Using uniform dispersion (phi_g=1).  "
                "Build the reference with estimate_overdispersion=True for "
                "better-calibrated NB likelihoods.",
                stacklevel=2,
            )

        Y_marker = np.asarray(Y_marker, dtype=np.float32)
        lib_sizes = np.asarray(lib_sizes, dtype=np.float32)

        # Mismatch initialisation
        if mismatch is None:
            mismatch = compute_spatial_discordance(
                Y_marker, ref_marker, list(ref_marker.gene_names),
                lib_sizes=lib_sizes,
            )
        self._mismatch = mismatch

        # NNLS warm start
        logger.info("NNLS warm start …")
        Pi = self._nnls_init(Y_marker, R_lin, lib_sizes)
        logger.info("NNLS init done.  Max proportion: %.3f.", float(Pi.max()))

        # Compute alpha (smoothing weight) — recorded always
        if self.lambda_spatial > 0:
            self.alpha = min(0.5, self.lambda_spatial)
        else:
            self.alpha = 0.0
        logger.info(
            "Spatial smoothing: lambda=%.3f, alpha=%.3f.", self.lambda_spatial, self.alpha
        )

        # Main optimisation
        Pi, trace, n_iter, converged = self._coordinate_descent(
            Y_marker, Pi, R_lin, phi_g, lib_sizes, graph
        )

        # Store fitted state
        self.proportions_ = Pi.astype(np.float32)
        self.convergence_trace_ = trace
        self.n_iter_ = n_iter
        self.converged_ = converged
        self._R_lin = R_lin.copy()
        self._phi_g = phi_g.copy()
        self._lib_sizes = lib_sizes.copy()
        self._Y_marker = Y_marker

        elapsed = time.perf_counter() - t0
        logger.info(
            "SpatCAR fit complete: %d iters, converged=%s, %.1f s.",
            n_iter, converged, elapsed,
        )

        if not converged:
            warnings.warn(
                f"SpatCARModel did not converge within {self.max_iter} iterations "
                f"(final δΠ = {trace[-1]:.2e}, tol = {self.tol:.2e}).  "
                "Consider increasing max_iter or loosening tol.  "
                "Estimates may be suboptimal.",
                stacklevel=2,
            )
        return self

    # ------------------------------------------------------------------
    # NNLS warm start
    # ------------------------------------------------------------------

    def _nnls_init(
        self,
        Y_marker: np.ndarray,   # (N, G_m)
        R_lin: np.ndarray,      # (K, G_m)
        lib_sizes: np.ndarray,  # (N,)
    ) -> np.ndarray:
        lib_safe = np.maximum(lib_sizes, 1.0)
        Y_norm = (Y_marker / lib_safe[:, None]).astype(np.float64)
        R_nnls = R_lin.T.astype(np.float64)   # (G_m, K)

        N = Y_norm.shape[0]
        K = R_lin.shape[0]

        Pi = np.zeros((N, K), dtype=np.float32)
        for s in range(N):
            pi_s, _ = opt.nnls(R_nnls, Y_norm[s])
            Pi[s] = pi_s.astype(np.float32)

        eps = 1e-4
        Pi = Pi + eps
        Pi = project_simplex_batch(Pi)
        return Pi

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------

    def _coordinate_descent(
        self,
        Y_marker: np.ndarray,
        Pi: np.ndarray,
        R_lin: np.ndarray,
        phi_g: np.ndarray,
        lib_sizes: np.ndarray,
        graph: SpatialGraph,
    ) -> tuple[np.ndarray, list[float], int, bool]:
        batches = graph.batches
        A = graph.A
        alpha = self.alpha

        d_g = self._mismatch.d_g
        R_d = (R_lin * d_g[np.newaxis, :]).astype(np.float32)

        trace: list[float] = []
        converged = False

        for iteration in range(1, self.max_iter + 1):
            Pi_prev = Pi.copy()
            Pi = self._one_pass(Y_marker, Pi, R_d, phi_g, lib_sizes, A, batches, alpha)

            if iteration % self.update_mismatch_every == 0:
                self._mismatch = update_mismatch_factors(
                    Y_marker, Pi, R_lin, lib_sizes, self._mismatch
                )
                d_g = self._mismatch.d_g
                R_d = (R_lin * d_g[np.newaxis, :]).astype(np.float32)

            delta = float(np.linalg.norm(Pi - Pi_prev, "fro")) / Pi.shape[0]
            trace.append(delta)

            max_elem = float(np.abs(Pi - Pi_prev).max())
            if delta < self.tol and max_elem < self.tol * 50:
                logger.info(
                    "Converged at iteration %d (δΠ=%.2e, max_elem=%.2e).",
                    iteration, delta, max_elem,
                )
                converged = True
                break

        if not converged:
            logger.warning(
                "Did not converge within %d iterations (final δΠ=%.2e).",
                self.max_iter,
                trace[-1] if trace else float("nan"),
            )

        return Pi, trace, len(trace), converged

    def _one_pass(
        self,
        Y_marker: np.ndarray,
        Pi: np.ndarray,
        R_d: np.ndarray,
        phi_g: np.ndarray,
        lib_sizes: np.ndarray,
        A: sp.csr_matrix,
        batches: list[np.ndarray],
        alpha: float,
    ) -> np.ndarray:
        """One full pass over all spatial batches (Jacobi-style update)."""
        Pi_next = Pi.copy()
        for batch_idx in batches:
            Y_b = Y_marker[batch_idx]
            Pi_b = Pi[batch_idx]
            lib_b = lib_sizes[batch_idx]

            Pi_b_nb = _nb_multiplicative_update(Y_b, Pi_b, R_d, phi_g, lib_b)

            if alpha > 0.0:
                A_b = A[batch_idx, :]
                Pi_neigh = A_b @ Pi
                Pi_b_mix = (1.0 - alpha) * Pi_b_nb + alpha * Pi_neigh
                row_sums = Pi_b_mix.sum(axis=1, keepdims=True)
                Pi_b_mix = Pi_b_mix / np.maximum(row_sums, 1e-10)
                Pi_b_mix = np.maximum(Pi_b_mix, 1e-10)
                Pi_next[batch_idx] = Pi_b_mix
            else:
                Pi_next[batch_idx] = Pi_b_nb

        return Pi_next


# ---------------------------------------------------------------------------
# NB multiplicative update
# ---------------------------------------------------------------------------


def _nb_multiplicative_update(
    Y_b: np.ndarray,    # (B, G_m)
    Pi_b: np.ndarray,   # (B, K)
    R_d: np.ndarray,    # (K, G_m)
    phi_g: np.ndarray,  # (G_m,)
    lib_b: np.ndarray,  # (B,)
) -> np.ndarray:
    """Multiplicative NB update for a batch of spots.

    Guaranteed non-negative per-element.  Each row is renormalised to the
    simplex after the update.
    """
    eps = 1e-8
    Pi_b = np.maximum(Pi_b, eps).astype(np.float32)

    Mu_no_lib = np.maximum(Pi_b @ R_d, eps)   # (B, G_m)
    Mu = np.maximum(lib_b[:, None] * Mu_no_lib, eps)  # (B, G_m)

    # Numerator: Σ_g R_d[k,g] · y[b,g] / Mu_no_lib[b,g]
    num = (Y_b / Mu_no_lib) @ R_d.T   # (B, K)

    # Denominator: Σ_g R_d[k,g] · (y[b,g]+φ_g) / (Mu[b,g]+φ_g)
    den_ratio = (Y_b + phi_g[None, :]) / (Mu + phi_g[None, :])  # (B, G_m)
    den = lib_b[:, None] * (den_ratio @ R_d.T)  # (B, K)
    den = np.maximum(den, eps)

    Pi_new = np.maximum(Pi_b * (num / den), eps)
    row_sums = Pi_new.sum(axis=1, keepdims=True)
    return (Pi_new / np.maximum(row_sums, eps)).astype(np.float32)
