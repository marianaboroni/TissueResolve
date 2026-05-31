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
