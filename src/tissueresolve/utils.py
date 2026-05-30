"""
Shared utilities for TissueResolve.

Provides:
- Logging configuration (setup_logging)
- Benjamini-Hochberg FDR correction (bh_fdr)
- Random state management (set_random_state)
- Numeric helpers: safe_log, project_simplex_batch, nb_loglik_batch
- Sparse matrix helpers: dense_subset
- Memory and timing helpers: MemoryTracker, timer
"""
from __future__ import annotations

import logging
import time
import tracemalloc
from contextlib import contextmanager
from typing import Generator

import numpy as np

__all__ = [
    "setup_logging",
    "bh_fdr",
    "set_random_state",
    "safe_log",
    "project_simplex_batch",
    "nb_loglik_batch",
    "dense_subset",
    "MemoryTracker",
    "timer",
]


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------


def setup_logging(level: str = "INFO", *, logger_name: str = "tissueresolve") -> None:
    """Configure the TissueResolve package logger.

    Parameters
    ----------
    level:
        Logging level string, e.g. ``"DEBUG"``, ``"INFO"``, ``"WARNING"``.
    logger_name:
        Name of the root logger to configure.  Defaults to ``"tissueresolve"``.
    """
    numeric = getattr(logging, level.upper(), logging.INFO)
    log = logging.getLogger(logger_name)
    if not log.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                "[%(asctime)s] %(levelname)-8s %(name)s — %(message)s",
                datefmt="%H:%M:%S",
            )
        )
        log.addHandler(handler)
    log.setLevel(numeric)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Unified implementation replacing independent copies in both legacy
    packages (CHIMERA ``discordance._bh_correction`` and SpatCAR
    ``neighborhood_stats._bh_fdr``).

    Parameters
    ----------
    pvalues:
        1-D array of raw p-values.

    Returns
    -------
    np.ndarray
        BH-adjusted p-values, same shape.  Values are clipped to [0, 1].

    Notes
    -----
    For a vector of length *n*, the BH procedure sets the adjusted p-value
    for rank *i* (1-indexed, ascending) to ``min(p_i * n / i, 1)``, then
    enforces monotonicity from right to left so that lower-ranked (more
    significant) tests never have larger adjusted p-values than higher-ranked
    tests.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    if pvalues.ndim != 1:
        raise ValueError(f"pvalues must be 1-D, got shape {pvalues.shape}.")
    n = len(pvalues)
    if n == 0:
        return pvalues.copy()

    order = np.argsort(pvalues)
    ranked = pvalues[order]
    # scale by n / rank (rank is 1-indexed)
    adjusted = ranked * (n / np.arange(1, n + 1))
    # enforce monotone non-decrease from right to left
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    result = np.empty(n)
    result[order] = adjusted
    return result


# ---------------------------------------------------------------------------
# Random state
# ---------------------------------------------------------------------------


def set_random_state(seed: int) -> np.random.Generator:
    """Seed NumPy's global state and return a ``np.random.Generator``.

    Parameters
    ----------
    seed:
        Integer random seed.

    Returns
    -------
    np.random.Generator
    """
    np.random.seed(seed)
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# Numeric helpers (algorithms implemented in later stages)
# ---------------------------------------------------------------------------


def safe_log(x: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """Element-wise natural log with a floor at *eps* to prevent -inf."""
    return np.log(np.maximum(x, eps))


def project_simplex_batch(V: np.ndarray) -> np.ndarray:
    """Project rows of *V* onto the probability simplex.

    Uses the O(K log K) algorithm of Duchi et al. (2008).

    Parameters
    ----------
    V:
        Array of shape ``(B, K)``.  Each row is projected independently.

    Returns
    -------
    np.ndarray
        Shape ``(B, K)``, rows non-negative and summing to 1.
    """
    V = np.asarray(V, dtype=float)
    if V.ndim == 1:
        V = V[np.newaxis, :]
        return _project_rows(V)[0]
    return _project_rows(V)


def _project_rows(V: np.ndarray) -> np.ndarray:
    B, K = V.shape
    U = np.sort(V, axis=1)[:, ::-1]
    cssv = np.cumsum(U, axis=1)
    ind = np.arange(1, K + 1, dtype=np.float64)
    cond = U * ind > (cssv - 1.0)
    rho = K - 1 - np.argmax(cond[:, ::-1], axis=1)
    theta = (cssv[np.arange(B), rho] - 1.0) / (rho + 1.0)
    return np.maximum(V - theta[:, np.newaxis], 0.0)


def nb_loglik_batch(
    Y: np.ndarray,
    Mu: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    """Per-row NB log-likelihood (sum over columns), returning shape ``(B,)``.

    Parameterisation: NB(μ, φ) where Var = μ + μ²/φ.
    Only terms depending on μ are included (gamma terms omitted).

    Parameters
    ----------
    Y:
        Observed counts, shape ``(B, G)``.
    Mu:
        Expected means, shape ``(B, G)``.
    phi:
        Per-gene dispersion, shape ``(G,)``.
    """
    eps = 1e-8
    Mu_ = np.maximum(Mu, eps)
    phi_ = phi[np.newaxis, :]
    ll = Y * safe_log(Mu_) - (Y + phi_) * safe_log(Mu_ + phi_)
    return ll.sum(axis=1)


def dense_subset(X: "sp.spmatrix | np.ndarray", col_indices: np.ndarray) -> np.ndarray:  # type: ignore[name-defined]
    """Extract *col_indices* columns from *X* as a dense float32 array.

    Parameters
    ----------
    X:
        Sparse or dense matrix ``(N, G)``.
    col_indices:
        Integer column indices to extract, shape ``(G_m,)``.

    Returns
    -------
    np.ndarray
        Dense array of shape ``(N, G_m)``, dtype float32.
    """
    import scipy.sparse as sp  # lazy import — only needed in spatial workflow

    col_indices = np.asarray(col_indices, dtype=np.intp)
    if sp.issparse(X):
        csc = X.tocsc()
        return np.asarray(csc[:, col_indices].todense(), dtype=np.float32)
    return np.asarray(X[:, col_indices], dtype=np.float32)


# ---------------------------------------------------------------------------
# Memory and timing
# ---------------------------------------------------------------------------


class MemoryTracker:
    """Context manager that reports peak tracemalloc allocation.

    Example
    -------
    >>> with MemoryTracker("NNLS init") as mt:
    ...     run_expensive_function()
    >>> print(f"{mt.peak_mb:.1f} MB peak")
    """

    def __init__(self, label: str = "", logger: logging.Logger | None = None) -> None:
        self.label = label
        self._log = logger or logging.getLogger("tissueresolve.memory")
        self.peak_mb: float = 0.0

    def __enter__(self) -> "MemoryTracker":
        tracemalloc.start()
        self._start = time.perf_counter()
        return self

    def __exit__(self, *_: object) -> None:
        _current, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        elapsed = time.perf_counter() - self._start
        self.peak_mb = peak_bytes / 1024**2
        if self.label:
            self._log.debug(
                "%s: %.1f MB peak, %.2f s elapsed",
                self.label, self.peak_mb, elapsed,
            )


@contextmanager
def timer(label: str, log: logging.Logger) -> Generator[None, None, None]:
    """Wall-clock timer context manager."""
    t0 = time.perf_counter()
    yield
    log.info("%s finished in %.2f s", label, time.perf_counter() - t0)
