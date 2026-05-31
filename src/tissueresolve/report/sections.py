"""
Section builders for results-directory-driven HTML reports.

Each builder reads a results directory (tables/ + figures/ + metadata) and
returns the ordered list of ``(title, html)`` sections for the bulk or spatial
report.  Missing pieces are rendered as "not available" — failed checks and
warnings are surfaced, never hidden.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from tissueresolve.report import assets, templates as T
from tissueresolve.report.methods_text import estimate_type_statement

__all__ = ["bulk_sections", "spatial_sections", "build_sections_html"]

_BULK_ESTIMATE = ("These are <b>mRNA-derived proportions</b>, <b>not</b> "
                  "absolute cell fractions.")
_SPATIAL_ESTIMATE = ("These are <b>spot-level RNA-derived composition "
                     "estimates</b>, <b>not</b> direct cell counts.")


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
        return T.figure_block(title, iframe_src=f"figures/{p.name}",
                              caption=caption, links=links)
    return f"<h3>{T.escape(title)}</h3><p class='caption'>figure not available</p>"


def _reference_quality(results_dir: Path) -> str:
    ref = assets.find_table(results_dir, "reference_summary.tsv")
    counts = assets.find_table(results_dir, "cell_type_counts.tsv")
    body = T.df_table(ref) if ref is not None else "<p>reference summary not available</p>"
    if counts is not None:
        body += "<h3>Cells per cell type</h3>" + T.df_table(counts, max_rows=40)
    src = results_dir / "selected_gene_identifier_column.txt"
    for base in (results_dir, results_dir / "tables"):
        f = base / "selected_gene_identifier_column.txt"
        if f.exists():
            body += f"<p class='caption'>gene identifier source: {T.escape(f.read_text().strip())}</p>"
            break
    return body


def bulk_sections(results_dir: Path, *, run_metadata: Optional[dict] = None,
                  warnings: Optional[list] = None) -> list[tuple[str, str]]:
    results_dir = Path(results_dir)
    figs = assets.list_figures(results_dir)
    props = assets.find_table(results_dir, "bulk_estimated_proportions.tsv")
    qc = assets.find_table(results_dir, "bulk_qc.tsv")
    overlap = assets.find_table(results_dir, "gene_overlap.tsv")
    selected = assets.find_table(results_dir, "selected_genes.tsv")
    risk = assets.find_table(results_dir, "protocol_risk.tsv")
    sep = assets.find_table(results_dir, "separability_report.tsv",
                            "pairwise_resolvability.tsv")
    spill = assets.find_table(results_dir, "spillover_report.tsv",
                              "spillover_risk_by_celltype.tsv")
    warns = list(warnings or []) + assets.collect_warnings(results_dir)

    n_samples = props.shape[0] if props is not None else "?"
    n_types = props.shape[1] if props is not None else "?"

    secs: list[tuple[str, str]] = []
    secs.append(("Run metadata & estimate type",
                 T.estimate_box(_BULK_ESTIMATE)
                 + T.kv_table(run_metadata or {"modality": "bulk"})))
    secs.append(("Input data summary",
                 T.kv_table({"samples": n_samples, "cell types": n_types,
                             "panel genes": (len(selected) if selected is not None
                                             else "?")})))
    secs.append(("Single-cell reference quality", _reference_quality(results_dir)))
    secs.append(("Gene overlap and filtering",
                 (T.df_table(overlap) if overlap is not None else "<p>not available</p>")
                 + (f"<h3>Selected genes ({len(selected)})</h3>"
                    + T.df_table(selected, max_rows=20) if selected is not None else "")
                 + (f"<h3>Protocol risk</h3>{T.df_table(risk)}" if risk is not None else "")))
    secs.append(("Deconvolution predictions",
                 (T.df_table(props) if props is not None else "<p>not available</p>")
                 + _fig(figs, "bulk_composition_clustered_barplot",
                        "Clustered composition barplot")
                 + _fig(figs, "bulk_composition_heatmap", "Composition heatmap")))
    secs.append(("Prediction QC",
                 (T.df_table(qc) if qc is not None else "<p>not available</p>")
                 + _fig(figs, "bulk_qc_summary", "QC summary")))
    secs.append(("Uncertainty",
                 _fig(figs, "bulk_uncertainty_plot", "Uncertainty (bootstrap CIs)",
                      "Wide intervals → low-confidence estimates.")))
    secs.append(("Separability and spillover",
                 (T.df_table(sep, max_rows=20) if sep is not None else "")
                 + _fig(figs, "bulk_separability_heatmap", "Separability heatmap")
                 + (T.df_table(spill, max_rows=20) if spill is not None else "")
                 + _fig(figs, "bulk_spillover_heatmap", "Spillover heatmap")
                 + _fig(figs, "spillover_network", "Spillover network")))
    secs.append(("Warnings and limitations", T.warning_box(warns)))
    secs.append(("Methods", _methods_html(results_dir, "bulk")))
    secs.append(("Output files", T.file_list(_output_files(results_dir))))
    return secs


def spatial_sections(results_dir: Path, *, run_metadata: Optional[dict] = None,
                     warnings: Optional[list] = None) -> list[tuple[str, str]]:
    results_dir = Path(results_dir)
    figs = assets.list_figures(results_dir)
    props = assets.find_table(results_dir, "spatial_spot_proportions.tsv")
    qc = assets.find_table(results_dir, "spatial_qc.tsv")
    morans = assets.find_table(results_dir, "morans_i.tsv")
    overlap = assets.find_table(results_dir, "gene_overlap.tsv")
    sep = assets.find_table(results_dir, "separability_report.tsv",
                            "pairwise_resolvability.tsv")
    spill = assets.find_table(results_dir, "spillover_report.tsv",
                              "spillover_risk_by_celltype.tsv")
    warns = list(warnings or []) + assets.collect_warnings(results_dir)

    n_spots = props.shape[0] if props is not None else "?"
    n_types = props.shape[1] if props is not None else "?"
    meta = run_metadata or {}

    secs: list[tuple[str, str]] = []
    secs.append(("Run metadata & estimate type",
                 T.estimate_box(_SPATIAL_ESTIMATE) + T.kv_table(meta or {"modality": "spatial"})))
    secs.append(("Input spatial data summary",
                 T.kv_table({"spots": n_spots, "cell types": n_types})))
    secs.append(("Single-cell reference quality", _reference_quality(results_dir)))
    secs.append(("Gene overlap",
                 T.df_table(overlap) if overlap is not None else "<p>not available</p>"))
    secs.append(("Spatial graph / model summary",
                 T.kv_table({"lambda_spatial": meta.get("lambda_spatial", "?"),
                             "alpha": meta.get("alpha", "?"),
                             "converged": meta.get("converged", "?")}, raw_values=True)))
    secs.append(("Spatial predictions",
                 _fig(figs, "spatial_mean_composition_barplot",
                      "Average spot-level composition")
                 + _fig(figs, "spatial_spot_pie_charts", "Per-spot composition (pies)")
                 + _fig(figs, "spatial_abundance_maps", "Abundance maps")
                 + _fig(figs, "spatial_dominant_cell_type_map", "Dominant cell type map")))
    secs.append(("Spatial QC",
                 (T.df_table(qc.describe().T, max_rows=40) if qc is not None else "")
                 + _fig(figs, "spatial_qc_maps", "Spatial QC maps")))
    secs.append(("Spatial structure",
                 (T.df_table(morans) if morans is not None else "")
                 + _fig(figs, "spatial_morans_i_barplot", "Moran's I by cell type")))
    secs.append(("Separability and spillover",
                 (T.df_table(sep, max_rows=20) if sep is not None else "")
                 + _fig(figs, "spatial_separability_heatmap", "Separability heatmap")
                 + (T.df_table(spill, max_rows=20) if spill is not None else "")
                 + _fig(figs, "spatial_spillover_heatmap", "Spillover heatmap")))
    secs.append(("Warnings and limitations", T.warning_box(warns)))
    secs.append(("Methods", _methods_html(results_dir, "spatial")))
    secs.append(("Output files", T.file_list(_output_files(results_dir))))
    return secs


def _methods_html(results_dir: Path, modality: str) -> str:
    for base in (results_dir, results_dir / "tables"):
        mt = base / "methods.txt"
        if mt.exists():
            text = mt.read_text(encoding="utf-8")
            return f"<p>{T.escape(text)}</p>".replace("\n\n", "</p><p>")
    return f"<p>{T.escape(estimate_type_statement(modality))}</p>"


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
