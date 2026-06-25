"""Experimental edge-aware + family-specific spatial smoothing (opt-in).

Post-fit, interpretable re-smoothing of spot composition estimates. It does NOT
modify the NB-CAR solver, the default ``lambda_spatial`` (0.1), soft gating, or
unresolved mass; it operates on an already-computed proportion matrix plus the
spot coordinates (and optionally the expression matrix). This is the "wrapper
that reweights smoothing after an initial fit" route — chosen to avoid any
solver rewrite.

Two strategies, both opt-in:

* **edge-aware** — smoothing is attenuated across likely tissue boundaries /
  rare niches (where neighbouring spots disagree in expression or composition)
  and strengthened within homogeneous neighbourhoods.
* **family-specific lambda** — structural families (stroma/endothelium/smooth
  muscle/broad epithelium) get more smoothing; immune/rare/fine/high-collinearity
  families get less.

No Redeconve-like state-graph regularization; no new model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "spatial_adaptive_smoothing-0.1.0"

# Default family classes for the conservative family-lambda rules.
STRUCTURAL_FAMILIES = {
    "Stromal/Fibroblast", "Stromal", "Fibroblast", "Fibroblast lineage",
    "Endothelial", "Blood vessels", "Lymphatic EC", "Mural", "Smooth muscle",
    "Epithelial", "Airway epithelium", "Alveolar epithelium", "Mesothelium",
}
IMMUNE_RARE_FAMILIES = {
    "T/NK", "Lymphoid", "Myeloid", "Immune", "B", "Plasma",
}


# ---------------------------------------------------------------------------
# Neighbour graph helpers (numpy-only; no hard sklearn/scipy dependency)
# ---------------------------------------------------------------------------


def _knn(coords: np.ndarray, k: int) -> np.ndarray:
    """Return (N, k) integer neighbour indices by Euclidean coordinate distance."""
    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    k = int(min(k, max(n - 1, 1)))
    d2 = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    return np.argsort(d2, axis=1)[:, :k]


def _row_norm(df_or_arr):
    arr = np.asarray(df_or_arr, dtype=float)
    s = arr.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return arr / s


# ---------------------------------------------------------------------------
# Phase 1 — edge-aware weights
# ---------------------------------------------------------------------------


@dataclass
class EdgeAwareWeights:
    neighbors: np.ndarray              # (N, k) neighbour indices
    edge_weights: np.ndarray           # (N, k) per-edge weight in [min, max]
    spot_smoothness: np.ndarray        # (N,) mean edge weight per spot
    metadata: dict[str, Any] = field(default_factory=dict)


def compute_edge_aware_spatial_weights(
    expression,
    coordinates,
    initial_proportions=None,
    k_neighbors: int = 6,
    expression_weight: float = 0.5,
    composition_weight: float = 0.5,
    min_edge_weight: float = 0.1,
    max_edge_weight: float = 1.0,
) -> EdgeAwareWeights:
    """Edge weights for the spot-neighbour graph: LOW across likely boundaries.

    A high expression/composition distance between neighbouring spots => low
    edge weight => weaker smoothing across that edge. Weights are bounded to
    ``[min_edge_weight, max_edge_weight]`` and finite.
    """
    coords = np.asarray(coordinates, dtype=float)
    n = coords.shape[0]
    nbr = _knn(coords, k_neighbors)
    X = np.asarray(expression.toarray() if hasattr(expression, "toarray") else expression,
                   dtype=float)
    # robust per-feature scaling so distances are comparable
    Xz = (X - X.mean(0)) / (X.std(0) + 1e-8)
    P = _row_norm(initial_proportions.values if isinstance(initial_proportions, pd.DataFrame)
                  else initial_proportions) if initial_proportions is not None else None

    ew = np.full((n, nbr.shape[1]), np.nan)
    for i in range(n):
        for jj, j in enumerate(nbr[i]):
            d_expr = np.linalg.norm(Xz[i] - Xz[j]) / np.sqrt(Xz.shape[1])
            d_comp = (np.abs(P[i] - P[j]).sum() * 0.5) if P is not None else 0.0
            wsum = expression_weight + (composition_weight if P is not None else 0.0)
            wsum = wsum if wsum > 0 else 1.0
            dist = (expression_weight * d_expr +
                    (composition_weight * d_comp if P is not None else 0.0)) / wsum
            ew[i, jj] = np.exp(-3.0 * dist)  # decay; tuneable but fixed/conservative
    # normalise to [0,1] across all edges, then rescale to [min,max]
    lo, hi = np.nanmin(ew), np.nanmax(ew)
    ew01 = (ew - lo) / (hi - lo) if hi > lo else np.ones_like(ew)
    ew_scaled = min_edge_weight + ew01 * (max_edge_weight - min_edge_weight)
    ew_scaled = np.clip(ew_scaled, min_edge_weight, max_edge_weight)
    spot_smoothness = ew_scaled.mean(axis=1)
    meta = {
        "edge_aware": True,
        "k_neighbors": int(nbr.shape[1]),
        "expression_weight": float(expression_weight),
        "composition_weight": float(composition_weight if initial_proportions is not None else 0.0),
        "used_initial_proportions": initial_proportions is not None,
        "min_edge_weight": float(min_edge_weight),
        "max_edge_weight": float(max_edge_weight),
        "edge_weight_summary": {
            "mean": float(np.nanmean(ew_scaled)), "min": float(np.nanmin(ew_scaled)),
            "max": float(np.nanmax(ew_scaled)), "p10": float(np.nanpercentile(ew_scaled, 10)),
        },
    }
    return EdgeAwareWeights(neighbors=nbr, edge_weights=ew_scaled,
                            spot_smoothness=spot_smoothness, metadata=meta)


def apply_edge_aware_smoothing(
    initial_proportions: pd.DataFrame,
    coordinates,
    expression=None,
    edge_weights: Optional[EdgeAwareWeights] = None,
    base_lambda: float = 0.1,
    k_neighbors: int = 6,
    **kw,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Post-fit edge-aware smoothing of a (spots × types) proportion matrix.

    For each spot, a convex blend toward its edge-weighted neighbour mean:
    ``θ_i ← (1 − s_i)·θ_i + s_i·(Σ_j w_ij θ_j / Σ_j w_ij)`` where the smoothing
    fraction ``s_i = base_lambda · spot_smoothness_i`` is attenuated at edges
    (low ``spot_smoothness``). Rows remain valid compositions. Default smoothing
    is untouched (this is a separate, opt-in operator).
    """
    if edge_weights is None:
        if expression is None:
            raise ValueError("apply_edge_aware_smoothing needs expression or precomputed edge_weights")
        edge_weights = compute_edge_aware_spatial_weights(
            expression, coordinates, initial_proportions=initial_proportions,
            k_neighbors=k_neighbors)
    P = initial_proportions.to_numpy(dtype=float)
    nbr, ew, smooth = edge_weights.neighbors, edge_weights.edge_weights, edge_weights.spot_smoothness
    # normalise base_lambda into a [0,1) blend; cap to keep it conservative
    blend = float(min(max(base_lambda, 0.0), 0.5))
    out = P.copy()
    for i in range(P.shape[0]):
        w = ew[i]
        wsum = w.sum()
        if wsum <= 0:
            continue
        nbr_mean = (w[:, None] * P[nbr[i]]).sum(0) / wsum
        s_i = blend * float(smooth[i])           # less smoothing where edges are sharp
        out[i] = (1.0 - s_i) * P[i] + s_i * nbr_mean
    out = _row_norm(out)
    smoothed = pd.DataFrame(out, index=initial_proportions.index, columns=initial_proportions.columns)
    meta = {
        "adaptive_smoothing_used": True, "strategy": "edge_aware",
        "base_lambda": float(base_lambda), "default_lambda_unchanged": True,
        **edge_weights.metadata,
    }
    return smoothed, meta


# ---------------------------------------------------------------------------
# Phase 2 — family-specific lambda
# ---------------------------------------------------------------------------


def _family_class(family: str) -> str:
    if family in STRUCTURAL_FAMILIES:
        return "structural"
    if family in IMMUNE_RARE_FAMILIES:
        return "immune_rare"
    return "other"


def compute_family_lambda(
    family: str,
    family_metrics: Optional[dict] = None,
    base_lambda: float = 0.1,
    min_lambda: float = 0.01,
    max_lambda: float = 0.1,
) -> tuple[float, str]:
    """Conservative per-family smoothing strength. Returns ``(lambda, rule)``.

    Structural families → toward ``max_lambda``; immune/rare-niche → toward
    ``min_lambda``; high-collinearity / rare-subtype-bearing families are NOT
    given more smoothing (prefer unresolved mass instead). Bounded to
    ``[min_lambda, max_lambda]``.
    """
    fm = family_metrics or {}
    cls = _family_class(str(family))
    if cls == "structural":
        lam, rule = max_lambda, "structural→high"
    elif cls == "immune_rare":
        lam, rule = min_lambda, "immune_rare→low"
    else:
        lam, rule = base_lambda, "other→base"
    # evidence overrides (conservative): rare subtype present or high collinearity
    if fm.get("has_rare_subtype") or fm.get("signature_collinearity", 0.0) >= 0.9:
        lam = min(lam, (min_lambda + base_lambda) / 2.0)
        rule += "+rare/collinear→reduce"
    return float(np.clip(lam, min_lambda, max_lambda)), rule


def family_adaptive_smoothing(
    initial_proportions: pd.DataFrame,
    coordinates,
    mapping: dict[str, str],
    family_metrics: Optional[dict[str, dict]] = None,
    base_lambda: float = 0.1,
    min_lambda: float = 0.01,
    max_lambda: float = 0.1,
    k_neighbors: int = 6,
    fine_le_broad: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Smooth each fine column with its family's lambda (coordinate kNN mean blend).

    ``fine_le_broad`` enforces that a fine subtype's lambda does not exceed its
    family (broad) lambda. Rows remain valid compositions.
    """
    coords = np.asarray(coordinates, dtype=float)
    nbr = _knn(coords, k_neighbors)
    P = initial_proportions.to_numpy(dtype=float)
    cols = list(initial_proportions.columns)
    fam_metrics = family_metrics or {}

    fam_lambda: dict[str, float] = {}
    fam_rule: dict[str, str] = {}
    for fam in {str(mapping.get(c, c)) for c in cols}:
        lam, rule = compute_family_lambda(fam, fam_metrics.get(fam), base_lambda,
                                          min_lambda, max_lambda)
        fam_lambda[fam], fam_rule[fam] = lam, rule

    out = P.copy()
    col_lambda = {}
    for ci, c in enumerate(cols):
        fam = str(mapping.get(c, c))
        lam = fam_lambda.get(fam, base_lambda)
        if fine_le_broad:
            lam = min(lam, fam_lambda.get(fam, base_lambda))  # fine ≤ broad
        col_lambda[c] = lam
        blend = float(min(max(lam, 0.0), 0.5))
        col = P[:, ci]
        nbr_mean = col[nbr].mean(axis=1)
        out[:, ci] = (1.0 - blend) * col + blend * nbr_mean
    out = _row_norm(out)
    smoothed = pd.DataFrame(out, index=initial_proportions.index, columns=cols)
    meta = {
        "adaptive_smoothing_used": True, "strategy": "family_adaptive",
        "base_lambda": float(base_lambda), "min_lambda": float(min_lambda),
        "max_lambda": float(max_lambda), "default_lambda_unchanged": True,
        "family_lambdas": fam_lambda, "family_rules": fam_rule,
        "families_low_lambda": [f for f, l in fam_lambda.items() if l <= (min_lambda + base_lambda) / 2.0],
        "level": "fine(≤broad)" if fine_le_broad else "fine",
    }
    return smoothed, meta
