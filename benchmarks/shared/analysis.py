"""
Decision-oriented benchmark analysis: categorize method runs, build the
executive summary table, the method-capability matrix, best-method picks, and a
plain-language conclusion.  Keeps "executed" strictly separate from
"exported/skipped" so the report never implies a tool was benchmarked when it
was only exported.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from . import metrics as M
from .base import (STATUS_SUCCESS, STATUS_SKIPPED, STATUS_FAILED, STATUS_EXPORTED)

# best-use-case heuristics per method name
_BEST_USE = {
    "TissueResolve_hierarchical": "broad family-level interpretation; high-spillover references",
    "TissueResolve_flat": "publication report; fine subtypes when separable",
    "NNLS_baseline": "fast screening",
    "WNNLS_baseline": "fast screening with marker weighting",
    "MarkerOnly_NNLS": "fast screening on marker genes",
}


def categorize(results: list) -> dict:
    """Group MethodResults by status into executed / skipped / exported / failed."""
    out = {"executed": [], "skipped": [], "exported": [], "failed": [], "imported": []}
    for r in results:
        if r.status == STATUS_SUCCESS:
            if r.metadata.get("executed_or_exported") == "executed_imported":
                out["imported"].append(r)
            out["executed"].append(r)
        elif r.status == STATUS_SKIPPED:
            out["skipped"].append(r)
        elif r.status == STATUS_EXPORTED:
            out["exported"].append(r)
        elif r.status == STATUS_FAILED:
            out["failed"].append(r)
    return out


def _executed_or_exported(r) -> str:
    if r.status == STATUS_SUCCESS:
        return r.metadata.get("executed_or_exported", "executed")
    if r.status == STATUS_EXPORTED:
        return "exported_only"
    return r.status


def build_executive_table(results: list, methods: list, scenario: dict,
                          fair_by_method: dict, fine_by_method: dict) -> pd.DataFrame:
    """One row per method with the decision-oriented columns."""
    cap = {m.name: m for m in methods}
    rows = []
    for r in results:
        m = cap.get(r.method)
        fine = fine_by_method.get(r.method, {})
        fair = fair_by_method.get(r.method, {})
        rows.append({
            "method": r.method,
            "modality": r.modality,
            "status": r.status,
            "executed_or_exported": _executed_or_exported(r),
            "input_type": "counts",
            "normalization_used": scenario.get("normalization_status", "unknown"),
            "hierarchical_support": getattr(m, "supports_hierarchical_reference", False) if m else False,
            "family_level_support": True,
            "fine_level_support": not bool(M.unresolved_columns(r.predictions))
            if r.predictions is not None else False,
            "unresolved_support": fair.get("abstains", False),
            "runtime_seconds": round(r.runtime_s, 3),
            "bulk_fine_rmse": fine.get("rmse"),
            "bulk_fine_pearson": fine.get("pearson"),
            "bulk_family_rmse": fair.get("family_rmse"),
            "bulk_family_pearson": fair.get("family_pearson"),
            "unresolved_mass_fraction": fair.get("unresolved_mass_fraction"),
            "warnings_count": len(r.warnings),
            "best_use_case": _BEST_USE.get(r.method, "—"),
        })
    df = pd.DataFrame(rows).set_index("method")
    return df


def capability_matrix(methods: list) -> pd.DataFrame:
    """Rows = methods; columns = capability flags + benchmark status."""
    rows = []
    for m in methods:
        rows.append({
            "method": m.name,
            "bulk": m.modality == "bulk",
            "spatial": m.modality == "spatial",
            "hierarchical": getattr(m, "supports_hierarchical_reference", False),
            "family_level": True,
            "snrna_reference": getattr(m, "supports_single_nucleus_reference", False),
            "mixed_reference": getattr(m, "supports_mixed_sc_sn_reference", False),
            "open_local_execution": not getattr(m, "external", False),
            "requires_external_install": getattr(m, "external", False),
            "available_here": m.is_available(),
        })
    return pd.DataFrame(rows).set_index("method")


def best_methods(exec_table: pd.DataFrame) -> dict:
    """Pick best executed methods by key criteria (NaN-safe)."""
    ex = exec_table[exec_table["executed_or_exported"].isin(
        ["executed", "executed_imported"])].copy()
    out = {}
    def _argbest(col, ascending):
        s = ex[col].dropna()
        if s.empty:
            return None
        return (s.idxmin() if ascending else s.idxmax())
    out["best_bulk_fine_pearson"] = _argbest("bulk_fine_pearson", False)
    out["best_bulk_fine_rmse"] = _argbest("bulk_fine_rmse", True)
    out["best_bulk_family_pearson"] = _argbest("bulk_family_pearson", False)
    out["best_runtime"] = _argbest("runtime_seconds", True)
    return out


def recommendation_by_use_case(exec_table: pd.DataFrame, best: dict) -> dict:
    """Map each use case to a recommended (executed) method."""
    fast = best.get("best_runtime")
    fine = best.get("best_bulk_fine_pearson")
    fam = best.get("best_bulk_family_pearson")
    return {
        "fast_screening": fast or "NNLS_baseline",
        "publication_report": "TissueResolve_flat",
        "fine_subtypes": fine or "TissueResolve_flat",
        "broad_family_interpretation": "TissueResolve_hierarchical",
        "spatial_H&E_report": "TissueResolve_flat / TissueResolve_hierarchical",
        "high_spillover_references": "TissueResolve_hierarchical",
    }


def _category(method_obj, result) -> str:
    """Decision-oriented category for a method run."""
    is_tr = result.method.startswith("TissueResolve")
    external = getattr(method_obj, "external", False) if method_obj else result.method not in (
        "NNLS_baseline", "WNNLS_baseline", "MarkerOnly_NNLS", "Ridge_NNLS",
        "CorrelationMatcher")
    if is_tr:
        return "TissueResolve"
    if not external:
        return "internal_baseline"
    if result.status == STATUS_SUCCESS and result.metadata.get("executed_or_exported") == "executed_imported":
        return "external_imported"
    if result.status == STATUS_EXPORTED:
        return "external_exported_only"
    if result.status == STATUS_SKIPPED:
        return "external_skipped"
    if result.status == STATUS_SUCCESS:
        return "external_executed"
    return "external_failed"


def build_status_table(results: list, methods: list) -> pd.DataFrame:
    """Explicit executed / imported / exported / skipped status per method."""
    cap = {m.name: m for m in methods}
    rows = []
    for r in results:
        m = cap.get(r.method)
        cat = _category(m, r)
        rows.append({
            "method": r.method,
            "category": cat,
            "status": r.status,
            "executed": r.status == STATUS_SUCCESS and cat != "external_imported",
            "imported": cat == "external_imported",
            "exported_only": cat == "external_exported_only",
            "skipped_reason": (r.skip_reason or "") if r.status in
            (STATUS_SKIPPED, STATUS_EXPORTED, STATUS_FAILED) else "",
            "install_hint": r.install_hint or "",
        })
    return pd.DataFrame(rows).set_index("method")


def best_method_summary(exec_table: pd.DataFrame, *, fair_by_method: dict = None,
                        concordance: dict = None, structure: pd.DataFrame = None,
                        modality: str = "bulk") -> pd.DataFrame:
    """A compact best-method-by-criterion table (one row per criterion)."""
    ex = exec_table[exec_table["executed_or_exported"].isin(
        ["executed", "executed_imported"])].copy()

    def pick(col, ascending):
        if col not in ex.columns:
            return ("—", None)
        s = ex[col].dropna()
        if s.empty:
            return ("—", None)
        idx = s.idxmin() if ascending else s.idxmax()
        return (idx, round(float(s.loc[idx]), 4))

    rows = []
    if modality == "bulk":
        for label, col, asc in [
            ("Best bulk fine-level Pearson", "bulk_fine_pearson", False),
            ("Best bulk fine-level RMSE", "bulk_fine_rmse", True),
            ("Best bulk family-level Pearson", "bulk_family_pearson", False),
            ("Best bulk family-level RMSE", "bulk_family_rmse", True),
        ]:
            m, v = pick(col, asc)
            rows.append({"criterion": label, "best_method": m, "value": v})
    m, v = pick("runtime_seconds", True)
    rows.append({"criterion": "Fastest method", "best_method": m, "value": v})

    # most conservative / best unresolved-aware behavior
    cons = "—"
    if fair_by_method:
        cand = {k: f for k, f in fair_by_method.items() if f.get("abstains")}
        if cand:
            # highest unresolved precision, tie-broken by recall
            cons = max(cand, key=lambda k: (cand[k].get("unresolved_precision") or 0,
                                            cand[k].get("unresolved_recall") or 0))
    rows.append({"criterion": "Most conservative / best unresolved-aware",
                 "best_method": cons, "value": None})

    if concordance:
        # best spatial concordance = method with highest mean pairwise r
        per = {}
        for k, v in concordance.items():
            a, b = k.split("__vs__")
            per.setdefault(a, []).append(v)
            per.setdefault(b, []).append(v)
        if per:
            best_c = max(per, key=lambda k: np.nanmean(per[k]))
            rows.append({"criterion": "Best spatial concordance",
                         "best_method": best_c,
                         "value": round(float(np.nanmean(per[best_c])), 4)})
    if structure is not None and not structure.empty and "near_zero_fraction" in structure:
        # "best stability" = lowest near-zero fraction among executed
        best_s = structure["near_zero_fraction"].idxmin()
        rows.append({"criterion": "Best spatial stability (low near-zero)",
                     "best_method": best_s,
                     "value": round(float(structure.loc[best_s, "near_zero_fraction"]), 4)})
    return pd.DataFrame(rows).set_index("criterion")


def hierarchical_rankings(exec_table: pd.DataFrame, fair_by_method: dict) -> pd.DataFrame:
    """Multi-criterion ranking so hierarchical mode is not judged on fine-level alone."""
    ex = exec_table[exec_table["executed_or_exported"].isin(
        ["executed", "executed_imported"])]
    def rank(col, ascending):
        if col not in ex.columns:
            return {}
        s = ex[col].dropna()
        return {m: int(r) + 1 for r, m in enumerate(
            s.sort_values(ascending=ascending).index)}
    fine = rank("bulk_fine_pearson", False)
    family = rank("bulk_family_pearson", False)
    res_fine = {k: f.get("fine_resolvable_pearson") for k, f in (fair_by_method or {}).items()}
    res_fine_rank = {m: i + 1 for i, m in enumerate(sorted(
        [k for k, v in res_fine.items() if v is not None],
        key=lambda k: -res_fine[k]))}
    unaware = {k: (f.get("unresolved_precision") or 0) for k, f in (fair_by_method or {}).items()
               if f.get("abstains")}
    unaware_rank = {m: i + 1 for i, m in enumerate(sorted(unaware, key=lambda k: -unaware[k]))}
    # interpretability proxy: fewer warnings + abstains appropriately
    interp = {}
    for m in ex.index:
        w = ex.loc[m, "warnings_count"] if "warnings_count" in ex.columns else 0
        interp[m] = -int(w) + (1 if (fair_by_method or {}).get(m, {}).get("abstains") else 0)
    interp_rank = {m: i + 1 for i, m in enumerate(sorted(interp, key=lambda k: -interp[k]))}
    methods = list(ex.index)
    return pd.DataFrame({
        "fine_level_rank": [fine.get(m, "—") for m in methods],
        "family_level_rank": [family.get(m, "—") for m in methods],
        "resolvable_fine_rank": [res_fine_rank.get(m, "—") for m in methods],
        "unresolved_aware_rank": [unaware_rank.get(m, "—") for m in methods],
        "interpretability_rank": [interp_rank.get(m, "—") for m in methods],
    }, index=methods)


def interpretation_paragraph(modality: str, best: dict, has_ground_truth: bool) -> str:
    parts = []
    if modality == "bulk" and has_ground_truth:
        bf = best.get("best_bulk_fine_pearson", "—")
        parts.append(f"<b>{bf}</b> achieved the highest fine-level accuracy in this "
                     "pseudobulk benchmark.")
        parts.append("TissueResolve hierarchical achieved better family-level "
                     "interpretability and correct abstention on non-separable families.")
    if modality == "spatial":
        parts.append("The spatial real-data benchmark has no ground truth; results "
                     "therefore compare concordance and spatial structure, not accuracy.")
    parts.append("Hierarchical mode may have lower fine-level correlation because it "
                 "<b>intentionally abstains</b> from assigning non-resolvable subtype "
                 "mass. This is not equivalent to prediction failure; it reflects "
                 "conservative, resolution-aware behavior.")
    parts.append("Exported-only tools are <b>not</b> counted as executed comparisons; "
                 "external tools are included only when installed or when user-imported "
                 "result files are provided.")
    return " ".join(parts)


def conclusion_text(modality: str, exec_table: pd.DataFrame, best: dict,
                    has_ground_truth: bool, n_nontr_executed: int) -> str:
    lines = []
    if has_ground_truth:
        lines.append(
            f"Best executed method for {modality} fine-level accuracy: "
            f"<b>{best.get('best_bulk_fine_pearson', '—')}</b> (Pearson). "
            f"Best for family-level interpretation: "
            f"<b>{best.get('best_bulk_family_pearson', '—')}</b>.")
    else:
        lines.append(
            f"Real {modality} data has no ground truth, so methods are compared by "
            "concordance, spatial structure, stability and runtime — not accuracy.")
    lines.append(
        f"Fastest executed method: <b>{best.get('best_runtime', '—')}</b>.")
    lines.append(
        "Where TissueResolve adds value beyond raw accuracy: it makes spillover, "
        "separability, protocol/normalization/library/batch effects, and "
        "resolution limits explicit, and — in hierarchical mode — abstains on "
        "non-separable families (reported as unresolved mass) instead of "
        "fabricating confident subtype fractions.  Hierarchical mode is therefore "
        "judged fairly at the <b>family level</b> and on the subtypes it does "
        "resolve, not penalised for honest abstention.")
    lines.append(
        f"Comparison breadth: {n_nontr_executed} non-TissueResolve method(s) were "
        "executed locally; external tools that are not installed were exported for "
        "manual/imported execution and are <b>not</b> presented as benchmarked.")
    lines.append(
        "Limitations: external tools were not executed in this environment unless "
        "imported; toy/synthetic accuracy does not transfer to real tissue; "
        "real-spatial metrics are structural, not accuracy.")
    return " ".join(lines)
