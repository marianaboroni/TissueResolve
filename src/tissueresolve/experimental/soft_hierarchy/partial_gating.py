"""Partial confidence-weighted unresolved mass (Phase 2A core).

EXPERIMENTAL. Replaces the binary all-or-nothing subtype gate
(``estimate_partial_subtype_resolution``, which keeps a subtype's *full* mass iff
confidence ≥ threshold else zeroes it) with a **continuous** weight
``c_{k,s} ∈ [0,1]``:

    resolved_θ[k,s] = raw_θ[k,s] · c_{k,s}
    unresolved[f,s] = broad_mass[f,s] − Σ_{k∈f} resolved_θ[k,s]

This is a strict generalisation: ``c ∈ {0,1}`` reproduces the current hard gate;
calibrated ``c`` recovers the mass the threshold over-discards.  Mass is conserved
by construction (``c ≤ 1`` and the raw fine estimates sum to the family mass), so
``resolved_f + unresolved_f = broad_mass_f`` and Σ_f broad_mass_f is preserved.

Pure / deterministic / no core-pipeline imports beyond dataclasses.
"""
from __future__ import annotations

from typing import Mapping, Optional, Union

import numpy as np
import pandas as pd

from .config import SoftGatingConfig, SoftGatingResult


def _as_sample_subtype(conf, index, columns) -> pd.DataFrame:
    """Coerce confidence (scalar / Series[subtype] / DataFrame) to samples×subtypes."""
    if isinstance(conf, pd.DataFrame):
        return conf.reindex(index=index, columns=columns)
    if isinstance(conf, pd.Series):
        row = conf.reindex(columns)
        return pd.DataFrame(np.tile(row.to_numpy(float), (len(index), 1)),
                            index=index, columns=columns)
    # scalar
    return pd.DataFrame(float(conf), index=index, columns=columns)


def apply_partial_confidence_gating(
    broad_mass: pd.DataFrame,
    fine_estimates: pd.DataFrame,
    family_map: Mapping[str, str],
    *,
    confidence: Union[pd.DataFrame, pd.Series, float, None] = None,
    confidence_features: Optional[pd.DataFrame] = None,
    confidence_model=None,
    config: Optional[SoftGatingConfig] = None,
) -> SoftGatingResult:
    """Apply continuous confidence-weighted resolution with unresolved mass.

    Parameters
    ----------
    broad_mass : samples × families  (broad-family composition; rows sum to ~1)
    fine_estimates : samples × subtypes  (absolute fine = family × conditional,
        ALL subtypes — the *pre-gating* estimate; within each family the columns
        should sum to that family's broad mass)
    family_map : ``{subtype: family}``
    confidence : explicit ``c_{k,s}`` in [0,1] (DataFrame/Series/scalar). Mutually
        exclusive with (confidence_features + confidence_model).
    confidence_features, confidence_model : if given, ``confidence =
        confidence_model.predict(confidence_features)`` (samples × subtypes in [0,1]).
    """
    cfg = config or SoftGatingConfig()
    warnings: list = []
    subtypes = [str(c) for c in fine_estimates.columns]
    families = [str(c) for c in broad_mass.columns]
    idx = fine_estimates.index
    group_of = {st: str(family_map.get(st, st)) for st in subtypes}

    # ---- resolve confidence ----
    if confidence_model is not None and confidence_features is not None:
        conf = confidence_model.predict(confidence_features)
        conf = _as_sample_subtype(conf, idx, subtypes)
    elif confidence is not None:
        conf = _as_sample_subtype(confidence, idx, subtypes)
    else:
        raise ValueError("provide either `confidence` or "
                         "(`confidence_features` + `confidence_model`)")
    # missing-feature / NaN handling: missing confidence ⇒ 0 (fully unresolved),
    # which is the conservative (no-false-precision) default. Recorded as a warning.
    n_missing = int(conf.isna().to_numpy().sum())
    if n_missing:
        warnings.append(f"{n_missing} missing confidence value(s) set to 0 "
                        "(conservative: mass routed to unresolved).")
    conf = conf.fillna(0.0).clip(cfg.min_confidence, cfg.max_confidence)

    raw = fine_estimates.copy().astype(float)
    if cfg.clip_negative:
        raw = raw.clip(lower=0.0)
    resolved = raw * conf

    row_total_before = broad_mass.sum(axis=1)
    unresolved_cols = {}
    fam_mass_err = 0.0
    for fam in families:
        members = [st for st in subtypes if group_of[st] == fam]
        if not members:
            continue
        res_sum = resolved[members].sum(axis=1)
        fam_mass = broad_mass[fam].astype(float)
        # numerical guard: resolved within a family can't exceed family mass
        over = (res_sum - fam_mass).clip(lower=0.0)
        if float(over.max()) > cfg.mass_tol:
            # rescale resolved within the family down to family mass
            scale = np.where(res_sum > 0, fam_mass / res_sum.replace(0, np.nan), 1.0)
            scale = np.clip(scale, 0.0, 1.0)
            for m in members:
                resolved[m] = resolved[m] * pd.Series(scale, index=idx).fillna(1.0)
            res_sum = resolved[members].sum(axis=1)
            warnings.append(f"family '{fam}': resolved mass exceeded family mass; "
                            "rescaled to conserve mass.")
        u = (fam_mass - resolved[members].sum(axis=1)).clip(lower=0.0)
        unresolved_cols[f"unresolved_{fam}"] = u
        # track family-mass conservation error
        fam_mass_err = max(fam_mass_err,
                           float((resolved[members].sum(axis=1) + u - fam_mass).abs().max()))

    unresolved = (pd.DataFrame(unresolved_cols, index=idx)
                  if unresolved_cols else pd.DataFrame(index=idx))
    row_total_after = resolved.sum(axis=1) + (unresolved.sum(axis=1) if not unresolved.empty else 0.0)
    mass_err = float((row_total_after - row_total_before).abs().max())

    meta = {
        "algorithm_version": cfg.algorithm_version,
        "feature_status": cfg.feature_status,
        "n_subtypes": len(subtypes), "n_families": len(families),
        "n_samples": int(len(idx)),
        "mean_confidence": float(conf.to_numpy().mean()),
        "mean_unresolved_fraction": float(unresolved.sum(axis=1).mean())
        if not unresolved.empty else 0.0,
        "family_mass_error": fam_mass_err,
        "estimate_type": "rna_proportion",
    }
    return SoftGatingResult(
        raw_fine_estimates=raw, confidence_by_subtype=conf,
        resolved_fine_estimates=resolved, unresolved_by_family=unresolved,
        mass_conservation_error=mass_err, confidence_features=confidence_features,
        warnings=warnings, metadata=meta)
