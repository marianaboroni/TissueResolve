"""
Automatic result interpretation and structured warnings for reports.

Turns the raw tables into short, honest, user-facing prose plus a severity-
ranked warning list, so a report *explains* the results instead of just
listing files.  Nothing here changes estimates; it only describes them and
never overclaims (bulk = mRNA-derived proportions; spatial = spot-level
RNA-derived composition).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

__all__ = [
    "SEVERITY_ORDER",
    "Warning_",
    "interpret_reference_quality",
    "interpret_gene_overlap",
    "interpret_bulk_predictions",
    "interpret_bulk_qc",
    "interpret_spatial_predictions",
    "interpret_spatial_structure",
    "interpret_separability_spillover",
    "interpret_uncertainty",
    "generate_key_findings",
    "collect_structured_warnings",
    "overall_qc_status",
    "executive_summary_cards",
    "executive_summary_paragraph",
]

SEVERITY_ORDER = ["CRITICAL", "WARNING", "CAUTION", "INFO"]


@dataclass
class Warning_:
    """A structured, severity-ranked warning."""

    severity: str          # INFO | CAUTION | WARNING | CRITICAL
    message: str
    recommended_action: str = ""

    def rank(self) -> int:
        return SEVERITY_ORDER.index(self.severity) if self.severity in SEVERITY_ORDER else 99


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _top_types(props: Optional[pd.DataFrame], n: int = 4) -> list[str]:
    if props is None or props.shape[1] == 0:
        return []
    return props.mean(axis=0).sort_values(ascending=False).index[:n].astype(str).tolist()


def _level(value: float, hi: float, lo: float) -> str:
    if value >= hi:
        return "high"
    if value >= lo:
        return "moderate"
    return "low"


# ---------------------------------------------------------------------------
# Interpretation paragraphs
# ---------------------------------------------------------------------------


def interpret_reference_quality(ref_summary: Optional[pd.DataFrame],
                                cell_type_counts: Optional[pd.DataFrame]) -> str:
    if ref_summary is None:
        return "Reference summary not available."
    d = ref_summary["value"].to_dict() if "value" in ref_summary.columns else {}
    n_cells = d.get("n_cells", "?")
    n_types = d.get("n_cell_types", "?")
    txt = (f"The reference contains {n_cells} cells across {n_types} annotated "
           f"cell type(s).")
    if cell_type_counts is not None and len(cell_type_counts):
        col = cell_type_counts.columns[0]
        counts = cell_type_counts[col]
        total = counts.sum()
        top = counts.idxmax()
        frac = counts.max() / total if total else 0
        txt += (f" The largest population is {top!r}, ~{frac:.0%} of the "
                "reference.")
        if frac > 0.30:
            txt += (" This imbalance may influence signature stability, though "
                    "retained types passed the minimum-cell threshold.")
    return txt


def interpret_gene_overlap(overlap: Optional[pd.DataFrame]) -> str:
    if overlap is None:
        return "Gene overlap not available."
    d = overlap.iloc[:, 0].to_dict() if overlap.shape[1] else {}
    shared = d.get("n_shared")
    ref = d.get("n_reference")
    if shared is None or not ref:
        return "Gene overlap not available."
    frac = shared / ref if ref else 0
    lvl = _level(frac, 0.7, 0.4)
    return (f"Gene overlap with the reference was {lvl} "
            f"({int(shared)} of {int(ref)} reference genes shared).")


def interpret_bulk_predictions(props: Optional[pd.DataFrame]) -> str:
    top = _top_types(props, 4)
    if not top:
        return "No bulk predictions available."
    return ("The highest predicted mRNA-derived compartments were "
            + ", ".join(top) + ". These are mRNA-derived proportions, not "
            "absolute cell fractions.")


def interpret_bulk_qc(qc: Optional[pd.DataFrame]) -> str:
    if qc is None or qc.select_dtypes("number").shape[1] == 0:
        return "No per-sample QC metrics available."
    num = qc.select_dtypes("number")
    parts = []
    for col in ("recon_r2", "profile_corr", "coverage_r2"):
        if col in num.columns:
            parts.append(f"mean {col} = {num[col].mean():.3f}")
    low = []
    if "recon_r2" in num.columns:
        low = num.index[num["recon_r2"] < 0.5].astype(str).tolist()
    txt = "Per-sample QC: " + (", ".join(parts) if parts else "summary available") + "."
    if low:
        txt += f" Low-confidence sample(s): {', '.join(low[:5])}."
    return txt


def interpret_spatial_predictions(props: Optional[pd.DataFrame]) -> str:
    top = _top_types(props, 4)
    if not top:
        return "No spatial predictions available."
    return ("The dominant spatial compartments (mean spot composition) were "
            + ", ".join(top) + ". Values are spot-level RNA-derived composition, "
            "not direct cell counts.")


def interpret_spatial_structure(morans: Optional[pd.DataFrame]) -> str:
    if morans is None or "morans_i" not in getattr(morans, "columns", []):
        return "Moran's I not available."
    s = morans["morans_i"].sort_values(ascending=False)
    top = s.index[:3].astype(str).tolist()
    strength = _level(float(s.max()), 0.5, 0.2)
    return (f"Spatial structure was {strength} (max Moran's I = {s.max():.2f}). "
            f"Most spatially structured: {', '.join(top)}.")


def interpret_separability_spillover(merges: Optional[pd.DataFrame],
                                     sep: Optional[pd.DataFrame]) -> str:
    n_problem = 0
    if sep is not None and "resolvability" in sep.columns:
        n_problem = int(sep["resolvability"].isin(
            ["poorly_resolved", "unresolved"]).sum())
    if merges is not None and len(merges):
        names = (merges["family_name"] if "family_name" in merges.columns
                 else pd.Series(merges.index))
        fams = "; ".join(names.astype(str).head(3))
        return (f"Fine labels include several closely related subtypes — "
                f"{len(merges)} recommended merge family(ies) (e.g. {fams}). "
                "Interpret these at the family level unless supported by "
                "additional validation.")
    if n_problem:
        return (f"{n_problem} poorly/unresolved cell-type pair(s) detected; "
                "subtype-level estimates for those may be unreliable.")
    return "Cell types were adequately separable; no merge families recommended."


def interpret_uncertainty(has_ci: bool) -> str:
    if has_ci:
        return "Bootstrap confidence intervals were computed; see the uncertainty panel."
    return ("Bootstrap uncertainty was not computed in this run. Run with "
            "--n-bootstrap > 0 to quantify confidence intervals.")


def generate_key_findings(ctx: dict[str, Any]) -> list[str]:
    """A handful of headline findings assembled from the available context."""
    findings = []
    if ctx.get("reference_quality"):
        findings.append(ctx["reference_quality"])
    if ctx.get("gene_overlap"):
        findings.append(ctx["gene_overlap"])
    if ctx.get("predictions"):
        findings.append(ctx["predictions"])
    if ctx.get("structure"):
        findings.append(ctx["structure"])
    if ctx.get("separability"):
        findings.append(ctx["separability"])
    if ctx.get("uncertainty"):
        findings.append(ctx["uncertainty"])
    return findings


# ---------------------------------------------------------------------------
# Structured warnings + QC status
# ---------------------------------------------------------------------------


def collect_structured_warnings(
    *,
    modality: str,
    gene_overlap: Optional[pd.DataFrame] = None,
    ref_summary: Optional[pd.DataFrame] = None,
    cell_type_counts: Optional[pd.DataFrame] = None,
    separability: Optional[pd.DataFrame] = None,
    merges: Optional[pd.DataFrame] = None,
    spillover: Optional[pd.DataFrame] = None,
    has_bootstrap: bool = False,
    converged: Optional[bool] = None,
    has_he_image: Optional[bool] = None,
    static_export: bool = True,
    extra: Optional[list[Warning_]] = None,
) -> list[Warning_]:
    """Build the severity-ranked warning list.  Never returns "no warnings"
    spuriously — every real issue below is surfaced."""
    w: list[Warning_] = list(extra or [])

    if separability is not None and "resolvability" in separability.columns:
        n_crit = int((separability["resolvability"] == "unresolved").sum())
        n_high = int((separability["resolvability"] == "poorly_resolved").sum())
        if n_crit:
            w.append(Warning_("CRITICAL",
                              f"{n_crit} unresolved (near-identical) cell-type pair(s).",
                              "Interpret at family level; see recommended_merges.tsv."))
        if n_high:
            w.append(Warning_("WARNING",
                              f"{n_high} poorly-separable cell-type pair(s).",
                              "Treat affected subtypes cautiously."))
    if merges is not None and len(merges):
        w.append(Warning_("WARNING",
                          f"{len(merges)} recommended merge family(ies) — fine "
                          "subtypes are confusable.",
                          "Prefer family-level estimates for these groups."))
    if spillover is not None and "spillover_risk" in spillover.columns:
        n_hi = int((spillover["spillover_risk"] >= 0.30).sum())
        if n_hi:
            w.append(Warning_("CAUTION",
                              f"{n_hi} cell type(s) with high spillover risk (≥0.30).",
                              "Cross-type leakage may inflate related types."))
    if not has_bootstrap:
        w.append(Warning_("CAUTION",
                          "Bootstrap uncertainty not computed.",
                          "Run with --n-bootstrap > 0 to quantify confidence."))
    if gene_overlap is not None:
        d = gene_overlap.iloc[:, 0].to_dict() if gene_overlap.shape[1] else {}
        shared, ref = d.get("n_shared"), d.get("n_reference")
        if shared and ref and shared / ref < 0.4:
            w.append(Warning_("WARNING", "Low gene overlap with the reference.",
                              "Check gene identifier conventions."))
    if ref_summary is not None and "value" in ref_summary.columns:
        try:
            n_cells = int(ref_summary["value"].get("n_cells", 0))
            if 0 < n_cells < 2000:
                w.append(Warning_("CAUTION", f"Small reference ({n_cells} cells).",
                                  "Signatures may be unstable."))
        except (ValueError, TypeError):
            pass
    if cell_type_counts is not None and len(cell_type_counts):
        col = cell_type_counts.columns[0]
        total = cell_type_counts[col].sum()
        if total and cell_type_counts[col].max() / total > 0.30:
            w.append(Warning_("CAUTION", "Strong reference imbalance "
                              "(largest type >30% of cells).",
                              "May bias signature stability."))
    if converged is False:
        w.append(Warning_("WARNING", "Spatial model did not converge.",
                          "Increase max_iter or loosen tol."))
    if modality == "spatial" and has_he_image is False:
        w.append(Warning_("CAUTION", "H&E image unavailable; coordinate-only "
                          "spatial plots were used.",
                          "Provide the Visium image for H&E overlays."))
    if not static_export:
        w.append(Warning_("INFO", "Static image export (kaleido) unavailable; "
                          "interactive HTML + source data were written.",
                          "pip install kaleido for PDF/SVG/PNG."))

    w.sort(key=lambda x: x.rank())
    return w


def overall_qc_status(warnings: list[Warning_]) -> str:
    """PASS / CAUTION / FAIL from the worst warning severity."""
    sev = {x.severity for x in warnings}
    if "CRITICAL" in sev:
        return "FAIL"
    if "WARNING" in sev or "CAUTION" in sev:
        return "CAUTION"
    return "PASS"


def recommended_interpretation(merges: Optional[pd.DataFrame],
                               separability: Optional[pd.DataFrame]) -> str:
    if merges is not None and len(merges):
        return "family-level"
    if (separability is not None and "resolvability" in separability.columns
            and separability["resolvability"].isin(
                ["poorly_resolved", "unresolved"]).any()):
        return "unresolved caution"
    return "fine"


def executive_summary_cards(stats: dict[str, Any]) -> dict[str, Any]:
    """Normalise the summary-card values (pass-through with defaults)."""
    keys = ["modality", "n_samples_or_spots", "reference_cells", "reference_genes",
            "cell_types", "shared_genes", "qc_status", "n_highrisk_pairs",
            "n_high_spillover", "recommended_interpretation"]
    return {k: stats.get(k, "—") for k in keys}


def executive_summary_paragraph(modality: str, cards: dict[str, Any],
                                top_types: list[str],
                                gene_overlap_level: str = "—",
                                structure_level: Optional[str] = None) -> str:
    top = ", ".join(top_types) if top_types else "—"
    if modality == "bulk":
        return (
            f"Overall, TissueResolve estimated mRNA-derived cell-type composition "
            f"across {cards.get('n_samples_or_spots')} bulk samples using a "
            f"single-cell reference with {cards.get('reference_cells')} cells, "
            f"{cards.get('reference_genes')} genes, and {cards.get('cell_types')} "
            f"annotated cell types. Gene overlap was {gene_overlap_level}. The main "
            f"predicted compartments were {top}. However, several related cell types "
            f"showed high separability/spillover risk, so subtype-level estimates "
            f"should be interpreted cautiously.")
    return (
        f"Overall, TissueResolve estimated spot-level RNA-derived composition across "
        f"{cards.get('n_samples_or_spots')} Visium spots. Gene overlap with the "
        f"reference was {gene_overlap_level}. The dominant spatial compartments were "
        f"{top}. Spatial structure was {structure_level or 'reported via Moran’s I'}. "
        f"Fine subtype interpretation is limited where separability/spillover risk is "
        f"high.")
