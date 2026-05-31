"""
Section builders for results-directory-driven HTML reports.

Each builder reads a results directory (tables/ + figures/ + metadata) and
returns the ordered list of ``(title, html)`` sections for the bulk or spatial
report.  The layout leads with an executive summary, key findings, the main
publication figure and interpretation; raw tables/matrices are placed in a
collapsible "Detailed outputs" section.  Missing pieces render as "not
available"; warnings and failed checks are surfaced, never hidden.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional
import json

import pandas as pd

from tissueresolve.report import assets, interpretation as I, templates as T
from tissueresolve.report.methods_text import estimate_type_statement

__all__ = ["bulk_sections", "spatial_sections", "build_sections_html"]

_BULK_ESTIMATE = ("These are <b>mRNA-derived proportions</b>, <b>not</b> "
                  "absolute cell fractions.")
_SPATIAL_ESTIMATE = ("These are <b>spot-level RNA-derived composition "
                     "estimates</b>, <b>not</b> direct cell counts.")


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------


def _fig(figures: dict, stem: str, title: str, caption: str = "") -> str:
    """Figure block referencing ``figures/<stem>.html`` if present."""
    if stem in figures:
        p = figures[stem]
        links = {}
        for fmt in ("pdf", "svg", "png"):
            sp = p.with_suffix(f".{fmt}")
            if sp.exists():
                links[fmt] = f"figures/{sp.name}"
        for suffix in ("data", "sample_order", "cell_type_order"):
            dp = p.with_name(f"{stem}.{suffix}.tsv")
            if dp.exists():
                links[f"{suffix}.tsv"] = f"figures/{dp.name}"
        cap = caption
        cf = p.with_name(f"{stem}.caption.txt")
        if cf.exists():
            cap = cf.read_text(encoding="utf-8").splitlines()[0]
        return T.figure_block(title, iframe_src=f"figures/{p.name}",
                              caption=cap, links=links)
    return f"<h3>{T.escape(title)}</h3><p class='caption'>figure not available</p>"


def _reference_quality(results_dir: Path) -> str:
    ref = assets.find_table(results_dir, "reference_summary.tsv")
    counts = assets.find_table(results_dir, "cell_type_counts.tsv")
    body = T.df_table(ref) if ref is not None else "<p>reference summary not available</p>"
    body += "<p>" + T.escape(I.interpret_reference_quality(ref, counts)) + "</p>"
    if counts is not None:
        body += T.collapsible("Cells per cell type (detailed)",
                              T.df_table(counts, max_rows=60))
    for base in (results_dir, results_dir / "tables"):
        f = base / "selected_gene_identifier_column.txt"
        if f.exists():
            body += (f"<p class='caption'>gene identifier source: "
                     f"{T.escape(f.read_text().strip())}</p>")
            break
    return body


# ---------------------------------------------------------------------------
# Shared context (interpretation, warnings, cards)
# ---------------------------------------------------------------------------


def _load(results_dir: Path):
    return {
        "props_bulk": assets.find_table(results_dir, "bulk_estimated_proportions.tsv"),
        "props_spatial": assets.find_table(results_dir, "spatial_spot_proportions.tsv"),
        "qc": assets.find_table(results_dir, "bulk_qc.tsv"),
        "overlap": assets.find_table(results_dir, "gene_overlap.tsv"),
        "sep": assets.find_table(results_dir, "separability_report.tsv",
                                 "pairwise_resolvability.tsv"),
        "spill": assets.find_table(results_dir, "spillover_report.tsv",
                                   "spillover_risk_by_celltype.tsv"),
        "merges": assets.find_table(results_dir, "recommended_merges.tsv"),
        "family_bulk": assets.find_table(results_dir, "bulk_family_proportions.tsv"),
        "family_spatial": assets.find_table(results_dir, "spatial_family_proportions.tsv"),
        "morans": assets.find_table(results_dir, "morans_i.tsv"),
        "ref_summary": assets.find_table(results_dir, "reference_summary.tsv"),
        "counts": assets.find_table(results_dir, "cell_type_counts.tsv"),
    }


def _ref_value(ref_summary, key, default="—"):
    if ref_summary is not None and "value" in getattr(ref_summary, "columns", []):
        try:
            return ref_summary["value"].get(key, default)
        except Exception:
            return default
    return default


def _has_bootstrap(results_dir: Path) -> bool:
    for n in ("bulk_uncertainty.tsv", "lower_ci.tsv"):
        if (results_dir / n).exists() or (results_dir / "tables" / n).exists():
            return True
    return False


def _build_cards_and_warnings(results_dir, modality, t, run_metadata):
    props = t["props_bulk"] if modality == "bulk" else t["props_spatial"]
    overlap = t["overlap"]
    n_shared = (overlap.iloc[:, 0].get("n_shared") if overlap is not None
                and overlap.shape[1] else "—")
    n_highrisk = (int(t["sep"]["resolvability"].isin(
        ["poorly_resolved", "unresolved"]).sum())
        if t["sep"] is not None and "resolvability" in t["sep"].columns else 0)
    n_high_spill = (int((t["spill"]["spillover_risk"] >= 0.30).sum())
                    if t["spill"] is not None
                    and "spillover_risk" in t["spill"].columns else 0)
    has_bs = _has_bootstrap(results_dir)
    converged = (run_metadata or {}).get("converged") if modality == "spatial" else None
    extra = [I.Warning_("CAUTION", m) for m in assets.collect_warnings(results_dir)]
    warns = I.collect_structured_warnings(
        modality=modality, gene_overlap=overlap, ref_summary=t["ref_summary"],
        cell_type_counts=t["counts"], separability=t["sep"], merges=t["merges"],
        spillover=t["spill"], has_bootstrap=has_bs, converged=converged,
        has_he_image=(run_metadata or {}).get("has_he_image"),
        static_export=(run_metadata or {}).get("static_export", True), extra=extra)
    status = I.overall_qc_status(warns)
    rec = I.recommended_interpretation(t["merges"], t["sep"])
    cards = I.executive_summary_cards({
        "modality": modality,
        "n_samples_or_spots": props.shape[0] if props is not None else "—",
        "reference_cells": _ref_value(t["ref_summary"], "n_cells"),
        "reference_genes": _ref_value(t["ref_summary"], "n_genes"),
        "cell_types": props.shape[1] if props is not None else "—",
        "shared_genes": n_shared,
        "qc_status": status,
        "n_highrisk_pairs": n_highrisk,
        "n_high_spillover": n_high_spill,
        "recommended_interpretation": rec,
    })
    return cards, warns, props


# ---------------------------------------------------------------------------
# Bulk
# ---------------------------------------------------------------------------


def bulk_sections(results_dir: Path, *, run_metadata: Optional[dict] = None,
                  warnings: Optional[list] = None) -> list[tuple[str, str]]:
    results_dir = Path(results_dir)
    figs = assets.list_figures(results_dir)
    t = _load(results_dir)
    cards, warns, props = _build_cards_and_warnings(
        results_dir, "bulk", t, run_metadata)

    top = (props.mean(0).sort_values(ascending=False).index[:4].astype(str).tolist()
           if props is not None else [])
    overlap_lvl = I.interpret_gene_overlap(t["overlap"]).split("was ")[-1].split(" (")[0] \
        if t["overlap"] is not None else "—"
    summary_para = I.executive_summary_paragraph("bulk", cards, top, overlap_lvl)
    findings = I.generate_key_findings({
        "reference_quality": I.interpret_reference_quality(t["ref_summary"], t["counts"]),
        "gene_overlap": I.interpret_gene_overlap(t["overlap"]),
        "predictions": I.interpret_bulk_predictions(props),
        "separability": I.interpret_separability_spillover(t["merges"], t["sep"]),
        "uncertainty": I.interpret_uncertainty(_has_bootstrap(results_dir)),
    })

    secs: list[tuple[str, str]] = []
    secs.append(("Executive summary",
                 T.estimate_box(_BULK_ESTIMATE) + T.summary_cards(cards)
                 + f"<p>{T.escape(summary_para)}</p>"))
    secs.append(("Key findings", T.key_findings(findings)))
    secs.append(("Main publication figure",
                 _fig(figs, "bulk_main_summary_figure",
                      "Bulk deconvolution summary (publication figure)")))
    secs.append(("Main results interpretation",
                 f"<p>{T.escape(I.interpret_bulk_predictions(props))}</p>"
                 f"<p>{T.escape(I.interpret_bulk_qc(t['qc']))}</p>"))
    secs.append(("Single-cell reference quality", _reference_quality(results_dir)))
    secs.append(("Input data summary",
                 T.kv_table({"samples": props.shape[0] if props is not None else "—",
                             "cell types": props.shape[1] if props is not None else "—",
                             "shared genes": cards["shared_genes"]})
                 + f"<p>{T.escape(I.interpret_gene_overlap(t['overlap']))}</p>"))
    secs.append(("Deconvolution predictions",
                 _fig(figs, "bulk_composition_clustered_barplot",
                      "Clustered composition (top types + Other)")
                 + f"<p>{T.escape(I.interpret_bulk_predictions(props))}</p>"))
    secs.append(("Prediction QC",
                 f"<p>{T.escape(I.interpret_bulk_qc(t['qc']))}</p>"
                 + (T.collapsible("Per-sample QC table", T.df_table(t["qc"]))
                    if t["qc"] is not None else "")
                 + _fig(figs, "bulk_qc_summary", "QC summary")))
    secs.append(("Separability and spillover", _separability_section(figs, t, "bulk")))
    secs.append(("Uncertainty", _uncertainty_section(results_dir, figs)))
    secs.append(("Warnings and limitations", T.severity_warning_box(warns)))
    secs.append(("Methods", _methods_html(results_dir, "bulk")))
    secs.append(("Detailed outputs", _detailed_outputs(t, props)))
    if run_metadata and run_metadata.get("analysis_plan") is not None:
        secs.append(("Analysis plan", _analysis_plan_section(run_metadata)))
    secs.append(("Output files", T.file_list(_output_files(results_dir))))
    return secs


# ---------------------------------------------------------------------------
# Spatial
# ---------------------------------------------------------------------------


def spatial_sections(results_dir: Path, *, run_metadata: Optional[dict] = None,
                     warnings: Optional[list] = None) -> list[tuple[str, str]]:
    results_dir = Path(results_dir)
    figs = assets.list_figures(results_dir)
    t = _load(results_dir)
    meta = run_metadata or {}
    cards, warns, props = _build_cards_and_warnings(
        results_dir, "spatial", t, meta)

    morans = t["morans"]
    top = (props.mean(0).sort_values(ascending=False).index[:4].astype(str).tolist()
           if props is not None else [])
    overlap_lvl = (I.interpret_gene_overlap(t["overlap"]).split("was ")[-1].split(" (")[0]
                   if t["overlap"] is not None else "—")
    struct_lvl = (I.interpret_spatial_structure(morans).split("was ")[-1].split(" (")[0]
                  if morans is not None else None)
    summary_para = I.executive_summary_paragraph("spatial", cards, top, overlap_lvl,
                                                 structure_level=struct_lvl)
    findings = I.generate_key_findings({
        "reference_quality": I.interpret_reference_quality(t["ref_summary"], t["counts"]),
        "gene_overlap": I.interpret_gene_overlap(t["overlap"]),
        "predictions": I.interpret_spatial_predictions(props),
        "structure": I.interpret_spatial_structure(morans),
        "separability": I.interpret_separability_spillover(t["merges"], t["sep"]),
    })

    secs: list[tuple[str, str]] = []
    secs.append(("Executive summary",
                 T.estimate_box(_SPATIAL_ESTIMATE) + T.summary_cards(cards)
                 + f"<p>{T.escape(summary_para)}</p>"))
    secs.append(("Key findings", T.key_findings(findings)))
    secs.append(("Main publication figure",
                 _fig(figs, "spatial_main_summary_figure",
                      "Spatial deconvolution summary (publication figure)")
                 + _fig(figs, "he_dominant_cell_type", "Dominant cell type on H&E")))
    secs.append(("Main results interpretation",
                 f"<p>{T.escape(I.interpret_spatial_predictions(props))}</p>"
                 f"<p>{T.escape(I.interpret_spatial_structure(morans))}</p>"))
    secs.append(("Single-cell reference quality", _reference_quality(results_dir)))
    secs.append(("Input data summary",
                 T.kv_table({"spots": props.shape[0] if props is not None else "—",
                             "cell types": props.shape[1] if props is not None else "—",
                             "shared genes": cards["shared_genes"]})
                 + f"<p>{T.escape(I.interpret_gene_overlap(t['overlap']))}</p>"))
    secs.append(("Spatial graph / model summary",
                 T.kv_table({"lambda_spatial": meta.get("lambda_spatial", "?"),
                             "alpha": meta.get("alpha", "?"),
                             "converged": meta.get("converged", "?")}, raw_values=True)))
    secs.append(("Spatial predictions",
                 _fig(figs, "spatial_mean_composition_barplot",
                      "Average spot-level composition")
                 + _fig(figs, "spatial_abundance_maps", "Abundance maps")
                 + _fig(figs, "spatial_dominant_cell_type_map", "Dominant cell type map")
                 + f"<p>{T.escape(I.interpret_spatial_predictions(props))}</p>"))
    secs.append(("Spatial structure",
                 _fig(figs, "spatial_morans_i_barplot", "Moran's I by cell type")
                 + f"<p>{T.escape(I.interpret_spatial_structure(morans))}</p>"))
    secs.append(("Separability and spillover", _separability_section(figs, t, "spatial")))
    secs.append(("Warnings and limitations", T.severity_warning_box(warns)))
    secs.append(("Methods", _methods_html(results_dir, "spatial")))
    secs.append(("Detailed outputs", _detailed_outputs(t, props)))
    if run_metadata and run_metadata.get("analysis_plan") is not None:
        secs.append(("Analysis plan", _analysis_plan_section(run_metadata)))
    secs.append(("Output files", T.file_list(_output_files(results_dir))))
    return secs


# ---------------------------------------------------------------------------
# Shared section bodies
# ---------------------------------------------------------------------------


def _separability_section(figs, t, modality) -> str:
    merges, sep = t["merges"], t["sep"]
    family = t["family_bulk"] if modality == "bulk" else t["family_spatial"]
    parts = [f"<p>{T.escape(I.interpret_separability_spillover(merges, sep))}</p>"]
    if merges is not None and len(merges):
        parts.append("<h3>Recommended merge families</h3>"
                     + T.df_table(merges, max_rows=40))
    # Top non-separable pairs up top (full matrix is collapsed in Detailed outputs).
    if sep is not None and {"type_a", "type_b"}.issubset(sep.columns):
        score_col = ("separability_score" if "separability_score" in sep.columns
                     else sep.columns[-1])
        top = sep.sort_values(score_col).head(10)
        parts.append("<h3>Top non-separable pairs</h3>" + T.df_table(top, max_rows=10))
    if family is not None and len(family):
        parts.append("<h3>Family-level estimates (safer interpretation)</h3>"
                     + T.df_table(family, max_rows=40))
    return "".join(parts)


def _uncertainty_section(results_dir: Path, figs) -> str:
    if _has_bootstrap(results_dir):
        return _fig(figs, "bulk_uncertainty_plot", "Bootstrap uncertainty")
    # Message card instead of a meaningless empty plot.
    return T.info_card(
        "Bootstrap uncertainty was not computed in this run. "
        "Run with <code>--n-bootstrap &gt; 0</code> to quantify confidence "
        "intervals (CI width by cell type, low-confidence flags).")


def _detailed_outputs(t, props) -> str:
    """Full tables/matrices, collapsed by default."""
    blocks = []
    if props is not None:
        blocks.append(T.collapsible("Full prediction table",
                                    T.df_table(props, max_rows=1000)))
    for label, key in (("Full separability table", "sep"),
                       ("Full spillover table", "spill")):
        if t.get(key) is not None:
            blocks.append(T.collapsible(label, T.df_table(t[key], max_rows=1000)))
    return "".join(blocks) or "<p>—</p>"


def _methods_html(results_dir: Path, modality: str) -> str:
    for base in (results_dir, results_dir / "tables"):
        mt = base / "methods.txt"
        if mt.exists():
            text = mt.read_text(encoding="utf-8")
            return f"<p>{T.escape(text)}</p>".replace("\n\n", "</p><p>")
    return f"<p>{T.escape(estimate_type_statement(modality))}</p>"


def _analysis_plan_section(run_metadata: dict) -> str:
    plan = run_metadata.get("analysis_plan")
    if not plan:
        return "<p>No analysis plan available.</p>"
    return f"<pre>{T.escape(json.dumps(plan, indent=2))}</pre>"


def _output_files(results_dir: Path) -> list[str]:
    results_dir = Path(results_dir)
    out = []
    for sub in ("tables", "figures"):
        d = results_dir / sub
        if d.is_dir():
            out += [f"{sub}/{p.name}" for p in sorted(d.iterdir()) if p.is_file()]
    return out


def build_sections_html(sections: list[tuple[str, str]]) -> tuple[str, list[str]]:
    """Render an ordered section list to HTML; return ``(html, toc_titles)``."""
    html = "".join(T.section(i, title, body) for i, (title, body) in enumerate(sections))
    return html, [t for t, _ in sections]
