"""Experimental diagnostic-driven auto-mode recommendation (opt-in, advisory).

``recommend_mode`` inspects reference/query QC, per-family resolution decisions,
and (for spatial) spatial diagnostics, and returns a conservative
:class:`ModeRecommendation`. It is **advisory only**: it never overrides a
user-specified mode/preset (it only warns), and with incomplete evidence it
recommends the safest option. Nothing here changes default behaviour.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "diagnostic_auto_mode-0.1.0"


@dataclass
class ModeRecommendation:
    modality: str
    recommended_mode: str
    recommended_preset: Optional[str]
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)
    is_experimental: bool = True


def _g(d: Optional[dict], key: str, default=None):
    return (d or {}).get(key, default)


def recommend_mode(
    modality: str,
    reference_qc: Optional[dict] = None,
    query_qc: Optional[dict] = None,
    resolution_decisions: Optional[dict] = None,
    spatial_diagnostics: Optional[dict] = None,
    benchmark_priors: Optional[dict] = None,
    config: Optional[dict] = None,
) -> ModeRecommendation:
    """Recommend a conservative mode/preset from diagnostics. Advisory only."""
    reasons: list[str] = []
    warnings: list[str] = []
    metrics: dict[str, float] = {}
    cfg = config or {}
    user_mode = cfg.get("user_mode")
    user_preset = cfg.get("user_preset")

    # ---- evidence completeness ----
    have_ref = bool(reference_qc)
    have_query = bool(query_qc)
    decisions = resolution_decisions or {}
    incomplete = not (have_ref and have_query) or (modality == "bulk" and not decisions)

    if modality == "bulk":
        rec_mode, rec_preset = "hierarchical_soft", None
        if incomplete:
            rec_mode = "broad_only"
            reasons.append("incomplete evidence → conservative broad_only")
        else:
            statuses = list(decisions.values())
            n_full = sum(1 for s in statuses if s == "full_fine")
            n_broad = sum(1 for s in statuses if s == "broad_only")
            frac_unsupported = (n_broad / len(statuses)) if statuses else 1.0
            metrics["frac_broad_only_families"] = float(frac_unsupported)
            if frac_unsupported >= 0.5:
                rec_mode = "broad_only"
                reasons.append("most families lack fine support → broad_only / unresolved")
            elif n_full == len(statuses) and len(statuses) > 0:
                rec_mode = "flat"
                reasons.append("all families fine-supported → flat acceptable")
            else:
                rec_mode = "hierarchical_soft"
                reasons.append("mixed support → hierarchical_soft with unresolved mass")
            if _g(reference_qc, "spillover_risk", 0.0) >= 0.5:
                warnings.append("high spillover risk; fine subtypes may be unreliable")
        return ModeRecommendation("bulk", rec_mode, rec_preset, reasons, warnings, metrics, True)

    # ---- spatial ----
    sd = spatial_diagnostics or {}
    boundary = float(_g(sd, "boundary_score", 0.0))
    rare = float(_g(sd, "rare_niche_signal", 0.0))
    oversmooth_risk = float(_g(sd, "oversmoothing_risk", 0.0))
    effn_inflation = float(_g(sd, "effective_n_inflation", 0.0))
    graph_quality = float(_g(sd, "graph_quality", 1.0))
    metrics.update({"boundary_score": boundary, "rare_niche_signal": rare,
                    "oversmoothing_risk": oversmooth_risk,
                    "effective_n_inflation": effn_inflation, "graph_quality": graph_quality})

    rec_mode, rec_preset = "default", "default"
    if not sd:
        reasons.append("no spatial diagnostics → default smoothing (conservative)")
        return ModeRecommendation("spatial", rec_mode, rec_preset, reasons, warnings, metrics, True)

    if graph_quality < 0.5:
        rec_preset = "default"
        warnings.append("poor spatial graph quality → do not rely heavily on smoothing")
        reasons.append("low graph quality → keep default, avoid aggressive presets")
    elif boundary >= 0.5 or oversmooth_risk >= 0.5:
        rec_preset = "edge_aware_smoothing"
        reasons.append("high boundary/oversmoothing risk → edge-aware smoothing (experimental)")
    elif rare >= 0.5:
        rec_preset = "weak_smoothing"
        reasons.append("strong rare-niche signal → weaker smoothing (experimental)")
    elif effn_inflation >= 0.5:
        rec_preset = "weak_smoothing"
        reasons.append("effective-N inflation → weaker smoothing (experimental)")
    else:
        reasons.append("no boundary/rare/effN risk → default smoothing")

    rec = ModeRecommendation("spatial", rec_mode, rec_preset, reasons, warnings, metrics, True)

    # ---- never override an explicit user choice ----
    if user_mode is not None and user_mode != rec.recommended_mode:
        rec.warnings.append(
            f"user mode '{user_mode}' kept; recommendation was '{rec.recommended_mode}' (not applied)")
        rec.recommended_mode = user_mode
    if user_preset is not None and user_preset != rec.recommended_preset:
        rec.warnings.append(
            f"user preset '{user_preset}' kept; recommendation was '{rec.recommended_preset}' (not applied)")
        rec.recommended_preset = user_preset
    return rec
