"""Resolution Decision Layer — decide *trusted* fine resolution per broad family
BEFORE fine-subpopulation predictions are interpreted.

This is an organizational/decision layer, not a new algorithm and not a new solver.
It formalises the pipeline principle that three distinct identifiability levels must
not be conflated:

    1. cell-level classifiability        (may be reported; NOT a decision criterion)
    2. pairwise mixture recoverability   (supports interpretation; cannot override #3)
    3. full-panel deconvolution reliability  (DOMINATES the trusted-resolution status)

For each family the trusted resolution is one of:

    * ``broad_only``    — subtypes are not reliably deconvolvable in the full panel;
                          show broad mass as trusted, fine only as diagnostic.
    * ``selected_fine`` — only some subtypes are supported; show those, route the rest
                          to unresolved / mark not-trusted.
    * ``full_fine``     — full fine composition is supportable (still soft-gated).

The decision is conservative under missing evidence and is recorded in run metadata so
it affects the *outputs* (status + report), not only visualisation. Soft gating remains
the final uncertainty layer and is applied exactly once, after fine prediction — this
layer never gates and never deconvolves.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional

try:                                   # Literal is 3.8+, the project targets 3.9+
    from typing import Literal
    _STATUS = Literal["broad_only", "selected_fine", "full_fine"]
except Exception:                      # pragma: no cover
    _STATUS = str

import numpy as np
import pandas as pd

ALGORITHM_VERSION = "resolution_decision-1.0"
BROAD_ONLY = "broad_only"
SELECTED_FINE = "selected_fine"
FULL_FINE = "full_fine"


@dataclass
class ResolutionDecisionConfig:
    """Thresholds for the trusted-resolution decision (full-panel-reliability first)."""
    # signature/identifiability (within-family) evidence
    sep_min: float = 0.10            # below → broad_only
    sep_full: float = 0.30           # at/above (with the rest) → full_fine eligible
    spillover_max: float = 0.30      # above → broad_only
    spillover_full: float = 0.20     # at/below → full_fine eligible
    disc_min: int = 10               # discriminating genes for any fine resolution
    disc_full: int = 20              # discriminating genes for full_fine
    # per-subtype confidence (within-family) for "supported" subtypes
    confidence_min: float = 0.10
    # query compatibility
    gene_overlap_min: float = 0.30   # below → broad_only (query too incompatible)
    query_marker_min: float = 0.10   # per-subtype query-detectable marker support
    # full-panel benchmark reliability (DOMINANT when available)
    benchmark_cond_rmse_broadonly: float = 0.35   # above → force broad_only
    benchmark_cond_rmse_full: float = 0.20        # at/below → full_fine eligible
    benchmark_spillover_max: float = 0.15         # above → cap at selected_fine
    feature_status: str = "stable"
    version: str = ALGORITHM_VERSION


@dataclass
class ResolutionDecision:
    family: str
    status: _STATUS
    reasons: list = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    thresholds: dict = field(default_factory=dict)
    supported_subtypes: list = field(default_factory=list)
    unsupported_subtypes: list = field(default_factory=list)
    warnings: list = field(default_factory=list)


def _members_by_family(family_map: Mapping[str, str]) -> dict[str, list]:
    fams: dict[str, list] = {}
    for sub, fam in family_map.items():
        fams.setdefault(str(fam), []).append(str(sub))
    return fams


def decide_trusted_resolution(
    reference_qc: Optional[Mapping[str, Any]] = None,
    query_qc: Optional[Mapping[str, Any]] = None,
    identifiability_metrics: Optional[pd.DataFrame] = None,
    benchmark_metrics: Optional[Mapping[str, Mapping[str, float]]] = None,
    family_map: Optional[Mapping[str, str]] = None,
    config: Optional[ResolutionDecisionConfig] = None,
    *,
    subtype_confidence: Optional[Mapping[str, Mapping[str, float]]] = None,
) -> dict[str, ResolutionDecision]:
    """Classify the trusted fine resolution per broad family.

    Parameters
    ----------
    reference_qc : optional ``{family: {"n_cells", "n_donors", ...}}`` reference support.
    query_qc : optional dict with ``"gene_overlap_frac"`` and optionally
        ``{"query_marker_support": {subtype: frac}}`` and ``"modality_compatible": bool``.
    identifiability_metrics : the within-family resolvability table (columns
        ``n_subtypes, mean_separability, min_discriminating_genes, mean_spillover,
        resolvable``) — the **full-panel-relevant signature evidence**.
    benchmark_metrics : optional ``{family: {"cond_rmse", "pairwise_spillover", ...}}``
        from prior validation — **full-panel reliability dominates** when present.
    family_map : ``{fine: broad}``.
    subtype_confidence : optional ``{subtype: {"confidence": float}}`` for selecting
        supported subtypes.

    Notes
    -----
    Cell-level AUROC is intentionally **not** an input — by rule it must not drive the
    decision.  Missing evidence is treated conservatively (never upgrades to full_fine).
    """
    cfg = config or ResolutionDecisionConfig()
    fams = _members_by_family(family_map or {})
    if identifiability_metrics is not None:
        for f in identifiability_metrics.index:
            fams.setdefault(str(f), [])
    overlap = float((query_qc or {}).get("gene_overlap_frac", 1.0))
    qmark = (query_qc or {}).get("query_marker_support", {})
    out: dict[str, ResolutionDecision] = {}

    for fam, members in fams.items():
        reasons, warns, metrics = [], [], {}
        row = (identifiability_metrics.loc[fam].to_dict()
               if identifiability_metrics is not None and fam in identifiability_metrics.index else {})
        n_sub = int(row.get("n_subtypes", len(members)) or len(members))
        sep = row.get("mean_separability")
        spill = row.get("mean_spillover")
        disc = row.get("min_discriminating_genes")
        resolvable = bool(row.get("resolvable", True))
        metrics.update({"n_subtypes": n_sub, "mean_separability": _f(sep),
                        "mean_spillover": _f(spill), "min_discriminating_genes": _f(disc),
                        "gene_overlap_frac": overlap})

        # supported subtypes by within-family confidence + query-detectable markers
        supported, unsupported = [], []
        for m in members:
            c = float((subtype_confidence or {}).get(m, {}).get("confidence", np.nan))
            qm = float(qmark.get(m, np.nan))
            ok_c = (not np.isfinite(c)) or (c >= cfg.confidence_min)
            ok_q = (not np.isfinite(qm)) or (qm >= cfg.query_marker_min)
            (supported if (ok_c and ok_q) else unsupported).append(m)
        if not supported and members:               # never strand a whole family
            supported = list(members); unsupported = []

        # ----- single-subtype family: trivially the family itself -----
        if n_sub <= 1:
            out[fam] = ResolutionDecision(
                fam, FULL_FINE, ["single subtype (fine == family)"], metrics,
                _thresholds(cfg), list(members), [], warns)
            continue

        # ----- benchmark full-panel reliability DOMINATES when available -----
        bench = (benchmark_metrics or {}).get(fam, {})
        b_rmse = bench.get("cond_rmse")
        b_spill = bench.get("pairwise_spillover")
        force_broad = force_selected = False
        if b_rmse is not None:
            metrics["benchmark_cond_rmse"] = float(b_rmse)
            if float(b_rmse) > cfg.benchmark_cond_rmse_broadonly:
                force_broad = True
                reasons.append(f"full-panel cond RMSE {b_rmse:.3f} > "
                               f"{cfg.benchmark_cond_rmse_broadonly} → broad_only")
        if b_spill is not None:
            metrics["benchmark_pairwise_spillover"] = float(b_spill)
            if float(b_spill) > cfg.benchmark_spillover_max:
                force_selected = True
                reasons.append(f"full-panel spillover {b_spill:.3f} > "
                               f"{cfg.benchmark_spillover_max} → cap at selected_fine")

        # ----- signature/query gates -----
        broad_reasons = []
        if not resolvable:
            broad_reasons.append("family flagged non-resolvable (signature)")
        if sep is not None and np.isfinite(sep) and sep < cfg.sep_min:
            broad_reasons.append(f"separability {sep:.3f} < {cfg.sep_min}")
        if spill is not None and np.isfinite(spill) and spill > cfg.spillover_max:
            broad_reasons.append(f"spillover {spill:.3f} > {cfg.spillover_max}")
        if disc is not None and np.isfinite(disc) and 0 <= disc < cfg.disc_min:
            broad_reasons.append(f"discriminating genes {int(disc)} < {cfg.disc_min}")
        if overlap < cfg.gene_overlap_min:
            broad_reasons.append(f"query gene overlap {overlap:.2f} < {cfg.gene_overlap_min}")
        if not supported:
            broad_reasons.append("no subtype passes confidence/query support")

        # missing signature evidence → conservative (cannot reach full_fine)
        evidence_complete = all(v is not None and np.isfinite(_f(v))
                                for v in (sep, spill, disc))
        if not evidence_complete:
            warns.append("incomplete signature evidence → conservative resolution")

        # ----- classify -----
        if force_broad or broad_reasons:
            status = BROAD_ONLY
            reasons += broad_reasons or ["benchmark full-panel reliability too low"]
        else:
            full_eligible = (
                resolvable
                and (sep is not None and np.isfinite(sep) and sep >= cfg.sep_full)
                and (spill is not None and np.isfinite(spill) and spill <= cfg.spillover_full)
                and (disc is not None and np.isfinite(disc) and disc >= cfg.disc_full)
                and overlap >= cfg.gene_overlap_min
                and len(unsupported) == 0
                and evidence_complete
                and not force_selected
                and (b_rmse is None or float(b_rmse) <= cfg.benchmark_cond_rmse_full))
            if full_eligible:
                status = FULL_FINE
                reasons.append("all signature/query/benchmark gates pass → full_fine")
            else:
                status = SELECTED_FINE
                reasons.append("partial support → selected_fine (supported subtypes only)")

        out[fam] = ResolutionDecision(
            fam, status, reasons, metrics, _thresholds(cfg),
            supported if status != BROAD_ONLY else [],
            (unsupported if status == SELECTED_FINE else
             (list(members) if status == BROAD_ONLY else [])),
            warns)
    return out


def decisions_from_estimates(estimates, *, query_qc=None, benchmark_metrics=None,
                             config=None) -> dict[str, ResolutionDecision]:
    """Derive trusted-resolution decisions from an already-computed
    :class:`HierarchicalEstimates` (no extra deconvolution).

    Uses the estimates' ``resolvability`` table, ``mapping`` and any per-subtype
    confidence recorded in metadata.  Cheap to call inside the pipeline.
    """
    conf = (estimates.metadata or {}).get("subtype_confidence")
    return decide_trusted_resolution(
        reference_qc=None, query_qc=query_qc,
        identifiability_metrics=estimates.resolvability,
        benchmark_metrics=benchmark_metrics, family_map=estimates.mapping,
        config=config, subtype_confidence=conf)


def decisions_to_frame(decisions: Mapping[str, ResolutionDecision]) -> pd.DataFrame:
    """Tidy per-family table (for the QC-first report's trusted-resolution section)."""
    rows = []
    for fam, d in sorted(decisions.items()):
        rows.append({
            "broad_family": fam, "trusted_resolution": d.status,
            "n_supported_subtypes": len(d.supported_subtypes),
            "n_unsupported_subtypes": len(d.unsupported_subtypes),
            "mean_separability": d.metrics.get("mean_separability"),
            "mean_spillover": d.metrics.get("mean_spillover"),
            "min_discriminating_genes": d.metrics.get("min_discriminating_genes"),
            "benchmark_cond_rmse": d.metrics.get("benchmark_cond_rmse"),
            "reason": "; ".join(d.reasons[:3]),
        })
    return pd.DataFrame(rows)


def decisions_metadata(decisions: Mapping[str, ResolutionDecision]) -> dict:
    """Compact JSON-able summary of the decisions for run metadata."""
    return {
        "version": ALGORITHM_VERSION,
        "trusted_resolution": {f: d.status for f, d in decisions.items()},
        "supported_subtypes": {f: list(d.supported_subtypes) for f, d in decisions.items()},
        "unsupported_subtypes": {f: list(d.unsupported_subtypes) for f, d in decisions.items()},
        "n_broad_only": sum(1 for d in decisions.values() if d.status == BROAD_ONLY),
        "n_selected_fine": sum(1 for d in decisions.values() if d.status == SELECTED_FINE),
        "n_full_fine": sum(1 for d in decisions.values() if d.status == FULL_FINE),
        "note": ("trusted resolution decided from full-panel reliability + signature + "
                 "query compatibility; cell-level AUROC is NOT a decision criterion"),
    }


# ===========================================================================
# Modality-aware gene weighting (shared compatibility layer + modality criteria)
# ===========================================================================
def compute_modality_aware_gene_weights(
    reference, query, modality: str, *, family_map=None, config=None,
    protocol_meta=None, gene_panel=None,
):
    """Modality-aware gene weights as a shared compatibility layer.

    Thin wrapper over the existing :class:`GeneSelector`: it makes the modality
    explicit and records the weighting mode, rather than introducing a new
    weighting method.  Bulk uses protocol-aware weighted NNLS weights; spatial
    uses unweighted marker selection (the NB-CAR model handles spot noise), so
    spatial returns ``None`` weights with the recorded rationale.

    Returns ``(weights_or_None, info)`` where *info* records the modality, mode,
    and the shared vs modality-specific criteria considered.
    """
    if modality not in ("bulk", "spatial"):
        raise ValueError(f"modality must be 'bulk' or 'spatial' (got {modality!r})")
    shared = ["reference_detection", "query_detection", "donor_stability",
              "subtype/broad_specificity", "separability", "protocol_compatibility",
              "technical_risk", "redundancy"]
    if modality == "bulk":
        info = {"modality": "bulk", "gene_weighting_mode": "protocol_aware_weighted_nnls",
                "shared_criteria": shared,
                "bulk_criteria": ["bulk_gene_detection", "platform_mismatch",
                                  "library_normalization_compatibility", "batch_sensitivity",
                                  "rna_content_sensitivity", "sc_vs_bulk_unstable_genes"],
                "weights_applied": True}
        weights = None
        try:
            from tissueresolve.reference.markers import GeneSelector
            sel = GeneSelector()
            panel = gene_panel
            if panel is None and hasattr(sel, "select"):
                panel = getattr(sel.select(reference, query, protocol_meta=protocol_meta),
                                "genes", None)
            if panel is not None and hasattr(sel, "compute_gene_weights"):
                weights = sel.compute_gene_weights(reference, panel, protocol_meta=protocol_meta)
                info["n_weighted_genes"] = int(len(weights)) if weights is not None else 0
        except Exception as exc:                       # never fail the wrapper
            info["weights_note"] = f"delegated weighting unavailable: {exc}"
        return weights, info
    # spatial
    info = {"modality": "spatial", "gene_weighting_mode": "marker_unweighted_nb_car",
            "shared_criteria": shared,
            "spatial_criteria": ["spot_level_detection", "spatial_dropout", "ambient_rna_risk",
                                 "spot_mixing_risk", "spatially_artefactual_genes",
                                 "oversmoothing_sensitivity", "coordinate_he_compatibility"],
            "weights_applied": False,
            "rationale": "spatial uses marker selection + NB-CAR spot noise model (no per-gene WNNLS weights)"}
    return None, info


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _thresholds(cfg: ResolutionDecisionConfig) -> dict:
    return {"sep_min": cfg.sep_min, "sep_full": cfg.sep_full,
            "spillover_max": cfg.spillover_max, "spillover_full": cfg.spillover_full,
            "disc_min": cfg.disc_min, "disc_full": cfg.disc_full,
            "gene_overlap_min": cfg.gene_overlap_min,
            "benchmark_cond_rmse_broadonly": cfg.benchmark_cond_rmse_broadonly,
            "benchmark_cond_rmse_full": cfg.benchmark_cond_rmse_full}
