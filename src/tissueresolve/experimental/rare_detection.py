"""Rare-subtype detection / calibration layer (P2; opt-in, post-hoc decision layer).

Count-likelihood deconvolution (Poisson/NB GLM) improves rare-subtype *recall* but
can raise the false-positive rate / lower precision for rare types (collinear or
weakly-supported states get assigned mass). This module addresses that **outside the
solver**, as a decision/detection layer on top of any solver's output:

  1. marker-support evidence  — does the query actually express the state's markers?
  2. rare detection probability — a calibrated P(present) from (predicted mass,
     marker support), optionally fit on labelled mixtures (isotonic / logistic);
  3. detection gate            — call a rare state present only if probability ≥
     threshold AND marker support ≥ a floor; otherwise move its mass to
     ``unresolved_<family>`` (mass-conserving) — never silently deleted;
  4. precision–recall curve    — diagnostic over thresholds.

It does NOT modify the solver, defaults, or the produced proportions unless the gate
is explicitly applied; the gated output is clearly labelled and conserves mass.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["RareDetectionResult", "compute_marker_support", "rare_detection_probability",
           "calibrate_detection", "precision_recall_curve", "apply_rare_detection_gate",
           "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
_EPS = 1e-12


@dataclass
class RareDetectionResult:
    probabilities: pd.DataFrame      # samples × states, calibrated P(present) in [0,1]
    marker_support: pd.DataFrame     # samples × states, fraction of markers detected
    calls: pd.DataFrame              # samples × states, bool (present)
    pr_curve: pd.DataFrame           # threshold, precision, recall, fpr (if labels given)
    metadata: dict = field(default_factory=dict)


def _top_markers(ref, states: Sequence[str], n_markers: int) -> dict:
    """Top specificity markers per state (genes where the state dominates CPM)."""
    R = ref.as_R_cpm().astype(np.float64)               # (K, G)
    R = R / (R.sum(axis=1, keepdims=True) + _EPS)       # row-normalise to per-gene share
    genes = list(ref.gene_names)
    others_max = np.zeros_like(R)
    for k in range(R.shape[0]):
        mask = np.ones(R.shape[0], bool); mask[k] = False
        others_max[k] = R[mask].max(axis=0)
    spec = R / (others_max + _EPS)                       # specificity ratio
    ct_idx = {str(c): i for i, c in enumerate(ref.cell_types)}
    out = {}
    for st in states:
        if str(st) not in ct_idx:
            out[str(st)] = []
            continue
        k = ct_idx[str(st)]
        order = np.argsort(spec[k])[::-1][:n_markers]
        out[str(st)] = [genes[i] for i in order]
    return out


def compute_marker_support(
    bulk: pd.DataFrame,
    ref,
    states: Sequence[str],
    *,
    n_markers: int = 25,
    detect_quantile: float = 0.5,
) -> pd.DataFrame:
    """Per (sample, state) fraction of the state's top markers expressed in the sample.

    ``bulk`` is genes × samples (counts or normalised). A marker is "detected" in a
    sample if its expression exceeds the sample's ``detect_quantile`` over the marker
    panel (robust, scale-free). Returns samples × states in [0, 1].
    """
    markers = _top_markers(ref, states, n_markers)
    samples = list(bulk.columns)
    out = pd.DataFrame(0.0, index=samples, columns=[str(s) for s in states])
    for st in out.columns:
        gset = [g for g in markers[st] if g in bulk.index]
        if not gset:
            out[st] = np.nan
            continue
        sub = bulk.loc[gset].to_numpy(float)            # (m, N)
        # per-sample detection threshold from the marker panel distribution
        thr = np.quantile(sub, detect_quantile, axis=0, keepdims=True)
        detected = (sub > np.maximum(thr, _EPS)).mean(axis=0)
        out[st] = detected
    return out


def rare_detection_probability(
    proportions: pd.DataFrame,
    marker_support: pd.DataFrame,
    states: Sequence[str],
    *,
    prop_scale: float = 0.02,
) -> pd.DataFrame:
    """Uncalibrated detection score in [0,1] from predicted mass × marker support.

    ``prop_scale`` sets the proportion at which mass alone gives ~0.5 evidence
    (``tanh(p/prop_scale)``); multiplied by marker support so unsupported mass is
    down-weighted. Use :func:`calibrate_detection` to map to a true probability.
    """
    states = [str(s) for s in states]
    P = proportions.reindex(columns=states).fillna(0.0)
    S = marker_support.reindex(index=P.index, columns=states).fillna(0.0)
    mass_ev = np.tanh(P.to_numpy(float) / max(prop_scale, _EPS))
    score = mass_ev * S.to_numpy(float)
    return pd.DataFrame(np.clip(score, 0.0, 1.0), index=P.index, columns=states)


def calibrate_detection(scores: np.ndarray, labels: np.ndarray, *, kind: str = "isotonic"):
    """Fit score→P(present) on labelled mixtures. Returns a callable score→prob.

    Caller guarantees calibration/evaluation split separation (no tuning on test).
    """
    x = np.asarray(scores, float).ravel()
    y = np.asarray(labels, float).ravel()
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size < 3 or len(np.unique(y)) < 2:
        return lambda s: np.clip(np.asarray(s, float), 0.0, 1.0)  # identity fallback
    if kind == "logistic":
        from sklearn.linear_model import LogisticRegression
        m = LogisticRegression(max_iter=1000).fit(x.reshape(-1, 1), y)
        return lambda s: m.predict_proba(np.asarray(s, float).reshape(-1, 1))[:, 1]
    from sklearn.isotonic import IsotonicRegression
    m = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(x, y)
    return lambda s: np.clip(m.predict(np.asarray(s, float).ravel()), 0.0, 1.0)


def precision_recall_curve(scores: np.ndarray, labels: np.ndarray,
                           thresholds: Optional[Sequence[float]] = None) -> pd.DataFrame:
    """Precision / recall / FPR over thresholds for a binary detection score."""
    x = np.asarray(scores, float).ravel()
    y = np.asarray(labels, float).ravel().astype(bool)
    mask = np.isfinite(x)
    x, y = x[mask], y[mask]
    if thresholds is None:
        thresholds = np.linspace(0.0, 1.0, 21)
    rows = []
    for t in thresholds:
        pred = x >= t
        tp = int((pred & y).sum()); fp = int((pred & ~y).sum())
        fn = int((~pred & y).sum()); tn = int((~pred & ~y).sum())
        prec = tp / (tp + fp) if (tp + fp) else float("nan")
        rec = tp / (tp + fn) if (tp + fn) else float("nan")
        fpr = fp / (fp + tn) if (fp + tn) else float("nan")
        rows.append({"threshold": float(t), "precision": prec, "recall": rec,
                     "fpr": fpr, "tp": tp, "fp": fp, "fn": fn, "tn": tn})
    return pd.DataFrame(rows)


def apply_rare_detection_gate(
    proportions: pd.DataFrame,
    probabilities: pd.DataFrame,
    marker_support: pd.DataFrame,
    states: Sequence[str],
    *,
    family_map: Optional[Mapping[str, str]] = None,
    threshold: float = 0.5,
    min_marker_support: float = 0.1,
    protected: Optional[Iterable[str]] = None,
) -> tuple[pd.DataFrame, dict]:
    """Gate rare-state calls: drop mass of un-supported / low-probability calls.

    A state is *kept present* when ``probability ≥ threshold`` AND
    ``marker_support ≥ min_marker_support`` (or it is ``protected``). Dropped mass is
    moved to ``unresolved_<family>`` (mass-conserving) — never silently deleted. Only
    the listed rare ``states`` are gated; all other columns are untouched.
    """
    states = [str(s) for s in states]
    protected = {str(s) for s in (protected or [])}
    df = proportions.copy()
    prob = probabilities.reindex(index=df.index, columns=states)
    sup = marker_support.reindex(index=df.index, columns=states)
    n_dropped = 0
    moved_mass = 0.0
    for st in states:
        if st not in df.columns:
            continue
        keep = (prob[st].fillna(0.0) >= threshold) & (sup[st].fillna(0.0) >= min_marker_support)
        if st in protected:
            keep[:] = True
        drop = ~keep
        if not drop.any():
            continue
        fam = str(family_map.get(st, st)) if family_map else st
        ucol = f"unresolved_{fam}"
        if ucol not in df.columns:
            df[ucol] = 0.0
        moved = df.loc[drop, st].to_numpy(float)
        moved_mass += float(moved.sum())
        df.loc[drop, ucol] = df.loc[drop, ucol].to_numpy(float) + moved
        df.loc[drop, st] = 0.0
        n_dropped += int(drop.sum())
    meta = {"feature_status": FEATURE_STATUS, "threshold": float(threshold),
            "min_marker_support": float(min_marker_support),
            "n_calls_dropped": n_dropped, "mass_moved_to_unresolved": moved_mass,
            "protected": sorted(protected), "gated_states": states,
            "mass_conserving": True}
    return df, meta
