"""Discriminative-marker-reweighted within-family conditional split (experimental).

Motivation (see docs/dev/DISTRIBUTION_AWARE_TNK_PILOT.md): collinear fine states
(e.g. T/NK subtypes, mean BC > 0.98) are actually ~95% separable at the cell level —
the few discriminating genes are real but get **diluted** in the full mean profile, so
the standard within-family split fails. This module reweights the within-family
conditional deconvolution by each gene's **between-sibling discriminative power**, so
the genes that separate the states drive the split instead of being drowned by shared
pan-family genes.

Opt-in, experimental: a within-family refinement that conserves the family mass and
leaves the rest of the pipeline (and defaults) unchanged.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from scipy.optimize import nnls

__all__ = ["discriminative_gene_weights", "solve_conditional_split",
           "refine_family_split", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
_EPS = 1e-12


def discriminative_gene_weights(
    state_profiles: np.ndarray,
    power: float = 1.0,
    min_weight: float = 0.05,
    log_transform: bool = True,
) -> np.ndarray:
    """Per-gene weight ∝ between-state (between-sibling) variance.

    Parameters
    ----------
    state_profiles:
        ``(S, G)`` profiles of the family's states (CPM / proportion scale, S ≥ 2).
    power:
        Sharpening exponent on the variance (1 = linear; >1 emphasises the top
        discriminating genes more). ``power=0`` → uniform weights (baseline).
    min_weight:
        Floor (relative to the mean) so shared genes are down-weighted but not zeroed
        (keeps the system identifiable / well-posed).
    log_transform:
        Use ``log1p`` of the profiles (scale-stabilised) before the variance.

    Returns
    -------
    ``(G,)`` non-negative weights, mean-normalised to 1.
    """
    P = np.asarray(state_profiles, dtype=np.float64)
    if P.ndim != 2 or P.shape[0] < 2:
        return np.ones(P.shape[-1])
    X = np.log1p(P) if log_transform else P
    v = X.var(axis=0)                                  # between-state variance per gene
    if power == 0:
        return np.ones(P.shape[1])
    w = np.power(np.maximum(v, 0.0), power)
    m = w.mean()
    w = w / m if m > 0 else np.ones_like(w)
    # floor relative to mean, then renormalise to mean 1
    w = np.maximum(w, min_weight)
    w = w / w.mean()
    return w


def solve_conditional_split(
    b: np.ndarray,
    Phi_states: np.ndarray,
    weights: Optional[np.ndarray] = None,
    total_mass: float = 1.0,
) -> np.ndarray:
    """Weighted-NNLS within-family split conditional on a fixed family mass.

    Solves ``min_{θ≥0} ‖√w ⊙ (b − Φ θ)‖²`` then renormalises θ to ``total_mass``.

    Parameters
    ----------
    b:
        ``(G,)`` observed (family) signal (e.g. normalised bulk on the gene panel).
    Phi_states:
        ``(G, S)`` state profiles (genes × states).
    weights:
        ``(G,)`` per-gene weights (e.g. from :func:`discriminative_gene_weights`).
        ``None`` → uniform (ordinary NNLS).
    total_mass:
        Sum the returned proportions should equal (the family mass).
    """
    b = np.asarray(b, dtype=np.float64)
    Phi = np.asarray(Phi_states, dtype=np.float64)
    G, S = Phi.shape
    w = np.ones(G) if weights is None else np.clip(np.asarray(weights, float), 0, None)
    sw = np.sqrt(w)
    theta, _ = nnls(Phi * sw[:, None], b * sw)
    tot = theta.sum()
    if tot <= _EPS:
        theta = np.full(S, 1.0 / S)
        tot = 1.0
    return theta / tot * total_mass


def refine_family_split(
    bulk_norm: np.ndarray,
    family_state_profiles: np.ndarray,
    family_mass: float = 1.0,
    *,
    power: float = 1.0,
    min_weight: float = 0.05,
) -> tuple[np.ndarray, dict]:
    """Convenience: compute discriminative weights and solve one family's split.

    ``bulk_norm`` is ``(G,)``; ``family_state_profiles`` is ``(S, G)``. Returns
    ``(theta_states summing to family_mass, metadata)``.
    """
    w = discriminative_gene_weights(family_state_profiles, power=power, min_weight=min_weight)
    theta = solve_conditional_split(bulk_norm, family_state_profiles.T, w, total_mass=family_mass)
    meta = {"feature_status": FEATURE_STATUS, "power": float(power),
            "min_weight": float(min_weight), "n_states": int(family_state_profiles.shape[0]),
            "weight_max": float(w.max()), "weight_top_frac":
            float((w > 1.0).mean())}
    return theta, meta
