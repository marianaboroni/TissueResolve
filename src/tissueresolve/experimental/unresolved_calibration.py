"""Experimental family-specific unresolved-mass calibration (opt-in).

Rebalances, **within each broad family**, how much mass is reported at the fine
level vs routed to ``unresolved_<family>`` — using per-family resolution
decisions, confidence, and support evidence. It is mass-conserving: each
family's total mass (its fine members + its unresolved bucket) is preserved, so
broad-level mass and overall mass are unchanged.

The goal is NOT to reduce unresolved mass everywhere, but to improve the balance
between false-resolution and false-abstention without zeroing rare,
marker-supported subtypes. Default behaviour is unchanged (``none``).

No solver changes; no soft-gating default change; unresolved mass is preserved
as a concept (never removed).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "unresolved_calibration-0.1.0"

UNRESOLVED_PREFIX = "unresolved_"


@dataclass
class CalibrationParams:
    confidence_threshold: float = 0.5      # min confidence to keep a subtype at fine level
    rare_abundance: float = 0.01           # below this mean abundance => "rare"
    collinearity_unresolved: float = 0.9   # >= => treat unsupported family conservatively
    keep_marker_supported_rare: bool = True


@dataclass
class CalibrationResult:
    predictions: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)


def _family_members(columns, mapping):
    fams: dict[str, list[str]] = {}
    for c in columns:
        if str(c).startswith(UNRESOLVED_PREFIX):
            continue
        fams.setdefault(str(mapping.get(c, c)), []).append(c)
    return fams


def calibrate_unresolved_by_family(
    raw_fine_predictions: pd.DataFrame,
    confidence_scores: Optional[dict] = None,
    resolution_decisions: Optional[dict] = None,
    family_metrics: Optional[dict] = None,
    calibration_params: Optional[CalibrationParams] = None,
    mapping: Optional[dict] = None,
    mode: str = "family_calibrated",
) -> CalibrationResult:
    """Return mass-conserving, family-calibrated predictions.

    Parameters mirror the spec. ``mode='none'`` returns the input unchanged.
    ``confidence_scores`` and ``family_metrics`` are keyed by fine label / family
    respectively; ``resolution_decisions`` maps family -> one of
    ``broad_only`` | ``selected_fine`` | ``full_fine``.
    """
    if mode == "none":
        return CalibrationResult(raw_fine_predictions.copy(),
                                 {"unresolved_calibration": "none"})
    params = calibration_params or CalibrationParams()
    conf = confidence_scores or {}
    decisions = resolution_decisions or {}
    fmetrics = family_metrics or {}
    mapping = mapping or {}

    df = raw_fine_predictions.copy().astype(float)
    fams = _family_members(df.columns, mapping)
    actions: dict[str, dict] = {}

    for fam, members in fams.items():
        members = [m for m in members if m in df.columns]
        if not members:
            continue
        ucol = f"{UNRESOLVED_PREFIX}{fam}"
        if ucol not in df.columns:
            df[ucol] = 0.0
        # total family mass to conserve (fine members + existing unresolved)
        fam_total = df[members].sum(axis=1) + df[ucol]
        decision = decisions.get(fam, "selected_fine" if len(members) > 1 else "full_fine")
        collinear = float(fmetrics.get(fam, {}).get("signature_collinearity", 0.0))
        mean_ab = df[members].mean(axis=0)

        keep_cols: list[str] = []
        if decision == "broad_only":
            keep_cols = []                                   # fine -> unresolved (diagnostic only)
        elif decision == "full_fine":
            keep_cols = list(members)
        else:  # selected_fine
            for m in members:
                c = float(conf.get(m, 1.0))
                rare = float(mean_ab[m]) < params.rare_abundance
                marker = bool(fmetrics.get(fam, {}).get("marker_supported", {}).get(m, False))
                if c >= params.confidence_threshold:
                    keep_cols.append(m)
                elif rare and marker and params.keep_marker_supported_rare:
                    keep_cols.append(m)                      # don't zero rare marker-supported
            if collinear >= params.collinearity_unresolved and not keep_cols:
                keep_cols = []                               # conservative: stay unresolved

        drop_cols = [m for m in members if m not in keep_cols]
        # route dropped mass to unresolved; renormalise kept members to (family - dropped) share
        dropped_mass = df[drop_cols].sum(axis=1) if drop_cols else pd.Series(0.0, index=df.index)
        for m in drop_cols:
            df[m] = 0.0
        df[ucol] = df[ucol] + dropped_mass
        # mass conservation guard: re-pin family total exactly
        new_total = df[members].sum(axis=1) + df[ucol]
        gap = (fam_total - new_total)
        df[ucol] = df[ucol] + gap                            # absorb rounding into unresolved
        df[ucol] = df[ucol].clip(lower=0.0)
        actions[fam] = {
            "decision": decision, "kept": keep_cols, "dropped": drop_cols,
            "collinearity": collinear,
            "mean_unresolved": float(df[ucol].mean()),
        }

    meta = {
        "unresolved_calibration": mode,
        "n_families": len(fams),
        "actions": actions,
        "confidence_threshold": params.confidence_threshold,
        "mass_conserving": True,
    }
    return CalibrationResult(df, meta)
