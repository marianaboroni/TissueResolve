"""
Reference suitability scoring.

Before deconvolution, evaluate whether a reference is suitable for a given query
and annotation scheme.  Each component returns a score in [0, 1] and a status;
the overall score is classified PASS / CAUTION / WARNING / FAIL.  Components that
cannot be evaluated (e.g. no batch metadata in a saved ReferenceSignature) are
reported as ``unknown`` (CAUTION) — never silently treated as PASS.

This module is read-only with respect to the reference and core algorithms.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = [
    "ComponentScore", "SuitabilityResult", "classify",
    "score_gene_overlap", "score_protocol_compatibility",
    "score_library_compatibility", "score_batch_confounding",
    "score_marker_stability", "score_celltype_balance",
    "score_fine_label_separability", "score_hierarchy_quality",
    "compute_reference_suitability_score", "summarize_reference_suitability",
    "save_reference_suitability",
]

# component weights (sum need not be 1; normalised internally)
_WEIGHTS = {
    "gene_overlap": 2.0, "protocol_compatibility": 1.0,
    "library_compatibility": 1.0, "batch_confounding": 1.0,
    "marker_stability": 1.0, "celltype_balance": 1.0,
    "fine_label_separability": 1.5, "hierarchy_quality": 1.0,
}


@dataclass
class ComponentScore:
    name: str
    score: float            # [0,1]; NaN when unknown
    status: str             # PASS | CAUTION | WARNING | FAIL | UNKNOWN
    detail: str = ""


@dataclass
class SuitabilityResult:
    overall_score: float
    classification: str
    components: list[ComponentScore] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def components_frame(self) -> pd.DataFrame:
        return pd.DataFrame([{"component": c.name, "score": c.score,
                              "status": c.status, "detail": c.detail}
                             for c in self.components]).set_index("component")


def classify(score: float) -> str:
    if np.isnan(score):
        return "CAUTION"
    if score >= 0.75:
        return "PASS"
    if score >= 0.55:
        return "CAUTION"
    if score >= 0.35:
        return "WARNING"
    return "FAIL"


def _status_from_score(score: float, *, pass_at=0.75, caution_at=0.55,
                       warn_at=0.35) -> str:
    if np.isnan(score):
        return "UNKNOWN"
    if score >= pass_at:
        return "PASS"
    if score >= caution_at:
        return "CAUTION"
    if score >= warn_at:
        return "WARNING"
    return "FAIL"


# --- components --------------------------------------------------------------


def score_gene_overlap(ref: ReferenceSignature,
                       query_genes: Optional[list[str]] = None) -> ComponentScore:
    if not query_genes:
        return ComponentScore("gene_overlap", float("nan"), "UNKNOWN",
                              "no query genes provided")
    rg = set(map(str, ref.gene_names))
    qg = set(map(str, query_genes))
    shared = rg & qg
    frac = len(shared) / max(len(rg), 1)
    s = float(min(1.0, frac / 0.6))  # ≥60% overlap → full marks
    return ComponentScore("gene_overlap", round(s, 3), _status_from_score(s),
                          f"{len(shared)} shared / {len(rg)} reference genes "
                          f"({100*frac:.0f}%)")


def score_protocol_compatibility(protocol_info: Optional[dict] = None) -> ComponentScore:
    if not protocol_info:
        return ComponentScore("protocol_compatibility", float("nan"), "UNKNOWN",
                              "protocol metadata not provided")
    risk = str(protocol_info.get("risk", "unknown")).lower()
    s = {"low": 1.0, "medium": 0.6, "high": 0.3}.get(risk, float("nan"))
    return ComponentScore("protocol_compatibility", s, _status_from_score(s)
                          if not np.isnan(s) else "UNKNOWN",
                          f"protocol risk={risk}")


def score_library_compatibility(library_info: Optional[dict] = None) -> ComponentScore:
    if not library_info:
        return ComponentScore("library_compatibility", float("nan"), "UNKNOWN",
                              "library-type metadata not provided")
    overall = str(library_info.get("overall", "unknown")).lower()
    if overall == "mixed":
        return ComponentScore("library_compatibility", 0.6, "CAUTION",
                              "mixed scRNA/snRNA reference — check confounding")
    if overall == "unknown":
        return ComponentScore("library_compatibility", float("nan"), "UNKNOWN",
                              "library type unknown")
    return ComponentScore("library_compatibility", 0.9, "PASS",
                          f"single library type: {overall}")


def score_batch_confounding(obs: Optional[pd.DataFrame] = None,
                            celltype_col: Optional[str] = None,
                            batch_col: Optional[str] = None) -> ComponentScore:
    if obs is None or not celltype_col or not batch_col:
        return ComponentScore("batch_confounding", float("nan"), "UNKNOWN",
                              "no batch/donor metadata available")
    try:
        from benchmarks.shared.batch_effects import compute_celltype_batch_confounding
        conf = compute_celltype_batch_confounding(obs, celltype_col, batch_col)
        frac = float(conf["confounded"].mean())
        s = float(1.0 - frac)
        return ComponentScore("batch_confounding", round(s, 3), _status_from_score(s),
                              f"{int(conf['confounded'].sum())}/{len(conf)} cell types "
                              "confounded with batch")
    except Exception as exc:  # noqa: BLE001
        return ComponentScore("batch_confounding", float("nan"), "UNKNOWN", str(exc))


def score_marker_stability(ref: ReferenceSignature) -> ComponentScore:
    cv = getattr(ref, "donor_cv", None)
    if cv is None:
        return ComponentScore("marker_stability", float("nan"), "UNKNOWN",
                              "no cross-donor CV (single donor or not recorded)")
    mean_cv = float(np.nanmean(cv))
    s = float(np.clip(1.0 - mean_cv, 0.0, 1.0))
    return ComponentScore("marker_stability", round(s, 3), _status_from_score(s),
                          f"mean cross-donor CV={mean_cv:.3f}")


def score_celltype_balance(ref: ReferenceSignature) -> ComponentScore:
    counts = np.array([ref.n_cells_per_type.get(c, 0) for c in ref.cell_types],
                      dtype=float)
    if counts.sum() == 0:
        return ComponentScore("celltype_balance", float("nan"), "UNKNOWN",
                              "no per-cell-type counts")
    min_cells = int(counts.min())
    # balance: 1 - Gini; plus a penalty if any type has very few cells
    s = float(np.clip(counts.min() / (counts.mean() + 1e-9), 0.0, 1.0))
    note = f"min cells/type={min_cells}, min/mean ratio={s:.2f}"
    if min_cells < 25:
        s = min(s, 0.5)
        note += " (some subtypes have <25 cells)"
    return ComponentScore("celltype_balance", round(s, 3), _status_from_score(s), note)


def score_fine_label_separability(ref: ReferenceSignature) -> ComponentScore:
    if ref.n_cell_types < 2:
        return ComponentScore("fine_label_separability", 1.0, "PASS",
                              "single cell type")
    try:
        from tissueresolve.reference.separability import compute_separability
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sep = compute_separability(ref, warn_threshold=2.0, raise_on_critical=False)
        seps = [1.0 - p.bhattacharyya_coeff for p in sep.pairs]
        frac_sep = float(np.mean([s >= 0.10 for s in seps])) if seps else 1.0
        s = frac_sep
        n_bad = int(sum(1 for p in sep.pairs if p.bhattacharyya_coeff > 0.90))
        return ComponentScore("fine_label_separability", round(s, 3),
                              _status_from_score(s),
                              f"{n_bad}/{len(sep.pairs)} pairs poorly separable "
                              f"(BC>0.9); {100*frac_sep:.0f}% separable")
    except Exception as exc:  # noqa: BLE001
        return ComponentScore("fine_label_separability", float("nan"), "UNKNOWN", str(exc))


def score_hierarchy_quality(ref: ReferenceSignature,
                            mapping: Optional[dict] = None) -> ComponentScore:
    if not mapping:
        return ComponentScore("hierarchy_quality", float("nan"), "UNKNOWN",
                              "no broad/fine hierarchy provided")
    cts = [str(c) for c in ref.cell_types]
    mapped = [c for c in cts if c in mapping]
    coverage = len(mapped) / max(len(cts), 1)
    fams: dict[str, int] = {}
    for c in mapped:
        fams[mapping[c]] = fams.get(mapping[c], 0) + 1
    n_other = sum(1 for c in mapped if str(mapping[c]).lower() == "other")
    other_frac = n_other / max(len(mapped), 1)
    s = float(np.clip(coverage * (1.0 - 0.5 * other_frac), 0.0, 1.0))
    return ComponentScore("hierarchy_quality", round(s, 3), _status_from_score(s),
                          f"{len(fams)} families, coverage={100*coverage:.0f}%, "
                          f"'Other'={100*other_frac:.0f}%")


# --- aggregate ---------------------------------------------------------------


def compute_reference_suitability_score(
    ref: ReferenceSignature, *, query_genes: Optional[list[str]] = None,
    mapping: Optional[dict] = None, protocol_info: Optional[dict] = None,
    library_info: Optional[dict] = None, obs: Optional[pd.DataFrame] = None,
    celltype_col: Optional[str] = None, batch_col: Optional[str] = None,
) -> SuitabilityResult:
    """Compute the overall reference suitability score and classification."""
    comps = [
        score_gene_overlap(ref, query_genes),
        score_protocol_compatibility(protocol_info),
        score_library_compatibility(library_info),
        score_batch_confounding(obs, celltype_col, batch_col),
        score_marker_stability(ref),
        score_celltype_balance(ref),
        score_fine_label_separability(ref),
        score_hierarchy_quality(ref, mapping),
    ]
    # weighted mean over the components that could be evaluated
    num = den = 0.0
    for c in comps:
        if not np.isnan(c.score):
            w = _WEIGHTS.get(c.name, 1.0)
            num += w * c.score
            den += w
    overall = float(num / den) if den > 0 else float("nan")
    cls = classify(overall)
    # worst-component override: a FAILing component must not be hidden behind a
    # high mean.  gene_overlap is critical (a FAIL there fails the reference);
    # any other FAIL caps the headline at WARNING.
    by_name = {c.name: c for c in comps}
    if by_name.get("gene_overlap") and by_name["gene_overlap"].status == "FAIL":
        cls = "FAIL"
    elif any(c.status == "FAIL" for c in comps):
        cls = "WARNING" if cls in ("PASS", "CAUTION") else cls
    n_unknown = sum(1 for c in comps if c.status == "UNKNOWN")
    if n_unknown >= 3 and cls == "PASS":
        cls = "CAUTION"   # too many unevaluated components to claim PASS
    warns = [f"{c.name}: {c.detail}" for c in comps
             if c.status in ("WARNING", "FAIL")]
    warns += [f"{c.name} could not be evaluated ({c.detail})"
              for c in comps if c.status == "UNKNOWN"]
    return SuitabilityResult(round(overall, 3) if not np.isnan(overall) else overall,
                             cls, comps, warns)


def summarize_reference_suitability(result: SuitabilityResult) -> str:
    lines = [f"# Reference suitability — **{result.classification}** "
             f"(score={result.overall_score})", ""]
    lines.append("| component | score | status | detail |")
    lines.append("|---|---|---|---|")
    for c in result.components:
        sc = "—" if np.isnan(c.score) else f"{c.score:.3f}"
        lines.append(f"| {c.name} | {sc} | {c.status} | {c.detail} |")
    if result.warnings:
        lines += ["", "## Warnings / not-evaluated", ""]
        lines += [f"- {w}" for w in result.warnings]
    lines += ["", "_PASS ≥ 0.75, CAUTION ≥ 0.55, WARNING ≥ 0.35, else FAIL. "
              "Unknown components are reported, not assumed to pass._"]
    return "\n".join(lines)


def save_reference_suitability(result: SuitabilityResult, out_dir: Path | str) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = {}
    score_df = pd.DataFrame([{"overall_score": result.overall_score,
                              "classification": result.classification}])
    p = out / "reference_suitability_score.tsv"
    score_df.to_csv(p, sep="\t", index=False)
    written["score"] = p
    p = out / "reference_suitability_components.tsv"
    result.components_frame().to_csv(p, sep="\t")
    written["components"] = p
    p = out / "reference_suitability_report.md"
    p.write_text(summarize_reference_suitability(result), encoding="utf-8")
    written["report"] = p
    import json
    p = out / "reference_suitability_warnings.json"
    p.write_text(json.dumps(result.warnings, indent=2), encoding="utf-8")
    written["warnings"] = p
    return written
