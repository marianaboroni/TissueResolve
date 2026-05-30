"""
Spatial deconvolution QC metrics for TissueResolve.

Per-spot metrics
----------------
- ``prop_entropy``  : Shannon entropy of the proportion vector.
- ``nb_loglik``     : Per-spot NB log-likelihood at the fitted proportions.
                      **Non-NaN only when Y_marker, R_d, phi_g, lib_sizes are
                      provided.** The pipeline always passes these arrays.
- ``dominant_type`` : Cell type with the largest proportion.
- ``dominant_frac`` : Proportion of the dominant type.
- ``n_types_gt05``  : Number of types with proportion > 0.05.
- ``spatial_resid`` : ‖π_s − mean(π_t, t∈N(s))‖₂.

Model-level
-----------
- ``compute_morans_i``   : Moran's I spatial autocorrelation per cell type.
- ``compute_model_qc``   : Convergence, mismatch factor summary.
- ``flag_low_quality_spots`` : Boolean mask of flagged spots.
- ``boundary_sharpness`` : Mean proportion gradient across a known boundary.
"""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.utils import nb_loglik_batch

__all__ = [
    "compute_spot_qc",
    "compute_morans_i",
    "flag_low_quality_spots",
    "compute_model_qc",
    "boundary_sharpness",
]

logger = logging.getLogger("tissueresolve.spatial.qc")


# ---------------------------------------------------------------------------
# Per-spot QC
# ---------------------------------------------------------------------------


def compute_spot_qc(
    Pi: np.ndarray,
    graph: Optional["SpatialGraph"],  # type: ignore[name-defined]
    cell_types: Optional[list[str]] = None,
    *,
    Y_marker: Optional[np.ndarray] = None,
    R_d: Optional[np.ndarray] = None,
    phi_g: Optional[np.ndarray] = None,
    lib_sizes: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Compute per-spot QC metrics from fitted proportions.

    Parameters
    ----------
    Pi:
        Fitted proportion matrix, shape ``(N, K)``.
    graph:
        :class:`~tissueresolve.spatial.graph.SpatialGraph`.
        Required for spatial residuals; set to ``None`` to skip.
    cell_types:
        Cell-type names for ``dominant_type`` column.
    Y_marker:
        Dense marker-gene count matrix ``(N, G_m)``.  Required for
        non-NaN ``nb_loglik``.
    R_d:
        Scaled reference ``(K, G_m)`` (R_lin × d_g).  Required for
        non-NaN ``nb_loglik``.
    phi_g:
        Per-gene NB dispersion ``(G_m,)``.  Required for non-NaN
        ``nb_loglik``.
    lib_sizes:
        Per-spot library sizes ``(N,)``.  Required for non-NaN
        ``nb_loglik``.

    Returns
    -------
    pd.DataFrame
        Rows indexed 0 … N−1.  Columns: ``prop_entropy``, ``nb_loglik``,
        ``dominant_type``, ``dominant_frac``, ``n_types_gt05``,
        ``spatial_resid``.
    """
    N, K = Pi.shape
    metrics: dict[str, np.ndarray] = {}

    # Shannon entropy
    Pi_safe = np.maximum(Pi, 1e-10)
    Pi_norm = Pi_safe / Pi_safe.sum(axis=1, keepdims=True)
    entropy = -np.sum(Pi_norm * np.log(Pi_norm), axis=1).astype(np.float32)
    metrics["prop_entropy"] = entropy

    # Dominant type
    dom_idx = Pi.argmax(axis=1).astype(np.int32)
    dom_frac = Pi[np.arange(N), dom_idx].astype(np.float32)
    if cell_types is not None and len(cell_types) == K:
        ct_arr = np.array(cell_types)
        metrics["dominant_type"] = ct_arr[dom_idx]
    else:
        metrics["dominant_type"] = dom_idx
    metrics["dominant_frac"] = dom_frac

    # Number of types > 5 %
    metrics["n_types_gt05"] = (Pi > 0.05).sum(axis=1).astype(np.int32)

    # NB log-likelihood — non-NaN only when model arrays are provided
    if (
        Y_marker is not None
        and R_d is not None
        and phi_g is not None
        and lib_sizes is not None
    ):
        Mu = lib_sizes[:, None] * (Pi @ R_d)
        Mu = np.maximum(Mu, 1e-8)
        ll = nb_loglik_batch(
            Y_marker.astype(np.float64),
            Mu.astype(np.float64),
            phi_g.astype(np.float64),
        )
        metrics["nb_loglik"] = ll.astype(np.float32)
    else:
        metrics["nb_loglik"] = np.full(N, np.nan, dtype=np.float32)

    # Spatial residual ‖π_s − π̄_neighbours‖₂
    if graph is not None:
        Pi_neigh = (graph.A @ Pi).astype(np.float32)
        spatial_resid = np.linalg.norm(Pi - Pi_neigh, axis=1).astype(np.float32)
        metrics["spatial_resid"] = spatial_resid
    else:
        metrics["spatial_resid"] = np.full(N, np.nan, dtype=np.float32)

    logger.info(
        "Spot QC: entropy mean=%.2f, n_uncertain(>0.9*log(K))=%d, "
        "nb_loglik_nan=%d.",
        float(entropy.mean()),
        int((entropy > 0.9 * np.log(max(K, 2))).sum()),
        int(np.isnan(metrics["nb_loglik"]).sum()),
    )
    return pd.DataFrame(metrics)


# ---------------------------------------------------------------------------
# Moran's I
# ---------------------------------------------------------------------------


def compute_morans_i(
    Pi: np.ndarray,
    graph: "SpatialGraph",  # type: ignore[name-defined]
    cell_types: Optional[list[str]] = None,
) -> pd.Series:
    """Compute Moran's I spatial autocorrelation per cell-type proportion map.

    Uses the formula::

        I_k = (z_k^T A z_k) / (z_k^T z_k)

    where z_k = π_{·k} − π̄_k (mean-centred) and A is the row-normalised
    adjacency.  Because A is row-normalised, the sum of weights equals N,
    simplifying the standard formula.

    Returns
    -------
    pd.Series
        Moran's I values indexed by cell-type name.
    """
    A = graph.A
    N, K = Pi.shape
    if cell_types is None:
        cell_types = [f"type_{k}" for k in range(K)]

    morans = np.empty(K, dtype=np.float32)
    for k in range(K):
        z = Pi[:, k] - Pi[:, k].mean()
        denom = float(z @ z)
        if denom < 1e-12:
            morans[k] = 0.0
            continue
        morans[k] = float(z @ (A @ z)) / denom

    result = pd.Series(morans, index=cell_types, name="morans_i", dtype=np.float32)
    logger.info(
        "Moran's I: median=%.3f, min=%.3f, max=%.3f.",
        float(result.median()), float(result.min()), float(result.max()),
    )
    return result


# ---------------------------------------------------------------------------
# Flagging
# ---------------------------------------------------------------------------


def flag_low_quality_spots(
    qc_df: pd.DataFrame,
    *,
    entropy_threshold: Optional[float] = None,
    loglik_threshold: Optional[float] = None,
    spatial_resid_threshold: Optional[float] = None,
) -> pd.Series:
    """Return a boolean mask of spots that fail QC.

    Thresholds are auto-derived from the data distribution when not provided
    (95th-percentile for entropy; 5th-percentile for loglik; 3× median for
    spatial residuals).  All thresholds are heuristic.

    Returns
    -------
    pd.Series
        ``True`` = flagged as low quality.
    """
    flag = pd.Series(False, index=qc_df.index)

    if "prop_entropy" in qc_df.columns:
        thresh_e = (
            entropy_threshold
            if entropy_threshold is not None
            else float(np.percentile(qc_df["prop_entropy"].dropna(), 95))
        )
        flag |= qc_df["prop_entropy"] > thresh_e

    if "nb_loglik" in qc_df.columns and qc_df["nb_loglik"].notna().any():
        thresh_l = (
            loglik_threshold
            if loglik_threshold is not None
            else float(np.percentile(qc_df["nb_loglik"].dropna(), 5))
        )
        flag |= qc_df["nb_loglik"] < thresh_l

    if "spatial_resid" in qc_df.columns and qc_df["spatial_resid"].notna().any():
        thresh_s = (
            spatial_resid_threshold
            if spatial_resid_threshold is not None
            else float(np.median(qc_df["spatial_resid"].dropna()) * 3.0)
        )
        flag |= qc_df["spatial_resid"] > thresh_s

    n_flagged = int(flag.sum())
    logger.info(
        "Flagged %d / %d spots as low quality (%.1f%%).",
        n_flagged, len(flag), 100.0 * n_flagged / max(len(flag), 1),
    )
    return flag


# ---------------------------------------------------------------------------
# Model-level QC
# ---------------------------------------------------------------------------


def compute_model_qc(
    n_iter: int,
    converged: bool,
    convergence_trace: list[float],
    d_g: Optional[np.ndarray] = None,
    marker_genes: Optional[list[str]] = None,
) -> dict:
    """Assemble a model-level QC summary dict.

    Returns
    -------
    dict
        Flat dict suitable for ``SpatialDeconvResult.run_metadata`` or AnnData.
    """
    summary: dict = {
        "n_iter": n_iter,
        "converged": converged,
        "final_delta": float(convergence_trace[-1]) if convergence_trace else None,
    }
    if d_g is not None:
        summary["d_g_median"] = round(float(np.median(d_g)), 4)
        summary["d_g_min"] = round(float(d_g.min()), 4)
        summary["d_g_max"] = round(float(d_g.max()), 4)
        summary["n_extreme_dg"] = int((d_g < 0.2).sum() + (d_g > 5.0).sum())
    if marker_genes is not None:
        summary["n_marker_genes"] = len(marker_genes)
    return summary


# ---------------------------------------------------------------------------
# Boundary sharpness
# ---------------------------------------------------------------------------


def boundary_sharpness(
    Pi: np.ndarray,
    array_row: np.ndarray,
    array_col: np.ndarray,
    boundary_col: float,
) -> float:
    """Mean dominant-type proportion gradient across a known tissue boundary.

    Parameters
    ----------
    Pi:
        Proportion matrix ``(N, K)``.
    array_row / array_col:
        Spot coordinates ``(N,)``.
    boundary_col:
        Column coordinate of the tissue boundary.

    Returns
    -------
    float
        Total-variation distance between the mean composition vectors on
        either side of the boundary, in ``[0, 1]``.  Higher = sharper
        boundary recovery.  A value of 1.0 means the two sides have
        completely disjoint dominant compositions; 0.0 means identical
        mean composition on both sides.

    Notes
    -----
    The metric compares the *full composition vector* across the boundary,
    not the dominant-type fraction magnitude.  A confident type flip
    (e.g. ``[0.9, 0.1]`` → ``[0.1, 0.9]``) is the signature of a sharp
    boundary and must score highly even though the dominant fraction
    (0.9) is unchanged.  Collapsing each side to ``argmax`` proportion
    would discard which type dominates and report 0.0 for a perfect
    boundary, so we use the L1/2 (total-variation) distance between the
    mean proportion vectors instead.
    """
    left = array_col < boundary_col
    right = array_col >= boundary_col
    if not left.any() or not right.any():
        return 0.0

    mean_l = Pi[left].mean(axis=0)   # (K,) mean composition left of boundary
    mean_r = Pi[right].mean(axis=0)  # (K,) mean composition right of boundary
    return float(0.5 * np.abs(mean_l - mean_r).sum())
