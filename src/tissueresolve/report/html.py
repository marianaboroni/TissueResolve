"""
Self-contained HTML report generation for TissueResolve (bulk + spatial).

The report is plain f-string HTML with embedded CSS (no Jinja2 required).
Design rules enforced here:

* **Warnings are surfaced, never hidden.**  A prominent warning box at the top
  lists the estimate-type caveat plus every QC recommendation, non-convergence,
  protocol-risk and separability problem.
* **Failed checks are shown**, not omitted (e.g. ``converged = False`` is
  rendered in red).
* Figures are embedded as base64 PNG with their captions, so the HTML is
  self-contained.
"""
from __future__ import annotations

import base64
import html as _html
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence, Union

import pandas as pd

from tissueresolve.report import methods_text as _mt

__all__ = ["generate_bulk_report", "generate_spatial_report"]

_CSS = """
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       margin: 0 auto; max-width: 980px; padding: 24px; color: #1a1a1a; }
h1 { border-bottom: 2px solid #4C72B0; padding-bottom: 6px; }
h2 { margin-top: 28px; color: #2a2a2a; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
table { border-collapse: collapse; margin: 8px 0; font-size: 13px; }
th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
th { background: #f0f3f8; }
td:first-child, th:first-child { text-align: left; }
.warnbox { background: #fff4e5; border: 1px solid #dd8452; border-radius: 6px;
           padding: 12px 16px; margin: 16px 0; }
.warnbox h2 { color: #b35900; border: none; margin-top: 0; }
.warn { color: #b35900; }
.fail { color: #c0392b; font-weight: bold; }
.ok { color: #1e7e34; font-weight: bold; }
.estimate { background: #eef3fb; border-left: 4px solid #4C72B0; padding: 8px 12px;
            margin: 12px 0; font-size: 14px; }
.caption { font-size: 12px; color: #555; margin: 4px 0 18px; }
.meta { font-family: monospace; font-size: 12px; }
img { max-width: 100%; border: 1px solid #eee; }
ul { font-size: 13px; }
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_bulk_report(
    result: Any,  # BulkPipelineResult
    output_path: Union[str, Path],
    *,
    separability: Any = None,
    figures: Optional[Sequence[Any]] = None,
    output_files: Optional[Iterable[Union[str, Path]]] = None,
    methods: Optional[str] = None,
    title: str = "TissueResolve — Bulk deconvolution report",
) -> Path:
    """Write a self-contained bulk HTML report.  Returns the written path."""
    deconv = result.deconv
    qc = getattr(result, "qc", None)
    meta = getattr(result, "run_metadata", {}) or {}
    risk = getattr(result, "protocol_risk", None)

    warnings = _collect_bulk_warnings(deconv, qc, risk, separability)

    parts: list[str] = []
    parts.append(_estimate_box(
        "These are <b>RNA-derived mRNA proportions</b>, <b>not</b> absolute "
        "cell fractions."
    ))
    parts.append(_warning_box(warnings))

    parts.append(_section("Input summary", _kv_table({
        "samples": deconv.n_samples,
        "cell types": deconv.n_cell_types,
        "panel genes": len(deconv.gene_panel),
        "reference genes": meta.get("n_reference_genes", "—"),
    })))
    parts.append(_section("Reference summary", _kv_table({
        "cell types": deconv.n_cell_types,
        "genome": meta.get("genome", "—"),
        "donor-aware": meta.get("donor_aware", "—"),
    })))
    if "n_shared_genes" in meta:
        parts.append(_section("Gene overlap", _kv_table({
            "shared genes": meta.get("n_shared_genes"),
        })))
    if risk is not None:
        parts.append(_section("Protocol risk", _kv_table({
            "risk level": getattr(risk, "risk_level", "unknown"),
            "mismatch types": ", ".join(getattr(risk, "mismatch_types", []) or []) or "—",
            "genes excluded": getattr(risk, "n_genes_excluded", 0),
            "genes down-weighted": len(getattr(risk, "downweighted_genes", []) or []),
        })))
    parts.append(_section("Selected genes", _gene_list_html(deconv.gene_panel)))
    parts.append(_section(
        "Deconvolution estimates (mRNA proportions)",
        _df_html(deconv.proportions.round(4))
        + "<h3>Summary</h3>" + _df_html(deconv.summary().round(4)),
    ))
    if qc is not None:
        parts.append(_section("QC metrics", _bulk_qc_html(qc)))
    parts.append(_section("Uncertainty", _bulk_uncertainty_html(deconv, qc)))
    if separability is not None:
        parts.append(_section("Separability", _separability_html(separability)))
    if figures:
        parts.append(_section("Figures", _figures_html(figures)))
    if output_files:
        parts.append(_section("Output files", _file_list_html(output_files)))
    parts.append(_section("Run metadata", _meta_html(meta)))

    methods = methods or _mt.compose_bulk_methods(result)
    parts.append(_section("Methods", f"<p>{_html.escape(methods)}</p>".replace("\n\n", "</p><p>")))

    return _write(output_path, title, parts)


def generate_spatial_report(
    result: Any,  # SpatialPipelineResult
    output_path: Union[str, Path],
    *,
    separability: Any = None,
    figures: Optional[Sequence[Any]] = None,
    output_files: Optional[Iterable[Union[str, Path]]] = None,
    methods: Optional[str] = None,
    title: str = "TissueResolve — Spatial deconvolution report",
) -> Path:
    """Write a self-contained spatial HTML report.  Returns the written path."""
    deconv = result.deconv
    meta = getattr(result, "run_metadata", {}) or {}
    morans = getattr(result, "morans_i", None)
    spot_qc = getattr(result, "spot_qc", None)

    warnings = _collect_spatial_warnings(result, deconv, separability)

    parts: list[str] = []
    parts.append(_estimate_box(
        "These are <b>spot-level RNA-derived composition estimates</b>, "
        "<b>not</b> direct single-cell counts."
    ))
    parts.append(_warning_box(warnings))

    parts.append(_section("Input summary", _kv_table({
        "spots": deconv.n_spots,
        "cell types": deconv.n_cell_types,
        "marker genes": len(deconv.marker_genes),
    })))
    parts.append(_section("Reference summary", _kv_table({
        "cell types": deconv.n_cell_types,
        "genome": meta.get("genome", "—"),
    })))
    if "n_genes_panel" in meta or "n_shared_genes" in meta:
        parts.append(_section("Gene overlap", _kv_table({
            "shared genes": meta.get("n_shared_genes", "—"),
            "panel genes": meta.get("n_genes_panel", len(deconv.marker_genes)),
        })))
    parts.append(_section("Spatial graph", _kv_table({
        "spots": deconv.n_spots,
        "alpha (mix weight)": meta.get("alpha", "—"),
    })))
    parts.append(_section("Smoothing parameters", _kv_table({
        "lambda_spatial": deconv.lambda_spatial,
        "alpha": meta.get("alpha", "—"),
        "n_smooth (niche)": deconv.n_smooth if deconv.n_smooth is not None else "not applied",
    })))
    conv_cls = "ok" if deconv.converged else "fail"
    parts.append(_section("Convergence", _kv_table({
        "converged": f'<span class="{conv_cls}">{deconv.converged}</span>',
        "iterations": deconv.n_iter,
        "final δΠ": round(deconv.convergence_trace[-1], 6) if deconv.convergence_trace else "—",
    }, raw_values=True)))
    parts.append(_section(
        "Spot-level estimates (RNA-derived composition)",
        f"<p>{deconv.n_spots} spots × {deconv.n_cell_types} cell types. "
        "Full table written to the output directory.</p>"
        + _df_html(deconv.proportions.head(10).round(4))
        + ("<p class='caption'>(first 10 spots shown)</p>"
           if deconv.n_spots > 10 else ""),
    ))
    if spot_qc is not None:
        parts.append(_section("Spatial QC", _df_html(spot_qc.describe().round(4))))
    if morans is not None:
        parts.append(_section("Moran's I", _df_html(morans.to_frame("morans_i").round(4))))
    if separability is not None:
        parts.append(_section("Separability", _separability_html(separability)))
    if figures:
        parts.append(_section("Figures", _figures_html(figures)))
    if output_files:
        parts.append(_section("Output files", _file_list_html(output_files)))
    parts.append(_section("Run metadata", _meta_html(meta)))

    methods = methods or _mt.compose_spatial_methods(result, separability=separability)
    parts.append(_section("Methods", f"<p>{_html.escape(methods)}</p>".replace("\n\n", "</p><p>")))

    return _write(output_path, title, parts)


# ---------------------------------------------------------------------------
# Warning collection
# ---------------------------------------------------------------------------


def _collect_bulk_warnings(deconv, qc, risk, separability) -> list[str]:
    w = [_mt.estimate_type_statement("bulk")]
    if qc is not None:
        w.extend(qc.recommendations or [])
        if qc.mismatch_flag is not None:
            high = [str(s) for s, f in qc.mismatch_flag.items() if str(f) == "high"]
            if high:
                w.append(f"High protocol-mismatch flag for sample(s): {', '.join(high)}.")
    if risk is not None and getattr(risk, "risk_level", "low") in ("medium", "high"):
        w.append(f"Protocol risk level is {risk.risk_level}; "
                 f"{getattr(risk, 'n_genes_excluded', 0)} gene(s) excluded.")
    if separability is not None and getattr(separability, "has_problems", False):
        w.append(f"{separability.n_critical} CRITICAL / {separability.n_high} HIGH "
                 "poorly-separable cell-type pair(s); their estimates are unreliable.")
    return w


def _collect_spatial_warnings(result, deconv, separability) -> list[str]:
    w = [_mt.estimate_type_statement("spatial")]
    if not deconv.converged:
        w.append(f"Model did NOT converge within {deconv.n_iter} iteration(s); "
                 "estimates may be suboptimal.")
    if deconv.lambda_spatial and deconv.lambda_spatial > 0:
        w.append(f"Spatial smoothing applied (λ_spatial = {deconv.lambda_spatial}); "
                 "parameter recorded, estimates are smoothed.")
    qc = getattr(result, "qc", None)
    if qc is not None and getattr(qc, "recommendations", None):
        w.extend(qc.recommendations)
    if separability is not None and getattr(separability, "has_problems", False):
        w.append(f"{separability.n_critical} CRITICAL / {separability.n_high} HIGH "
                 "poorly-separable cell-type pair(s); their estimates are unreliable.")
    return w


# ---------------------------------------------------------------------------
# HTML building blocks
# ---------------------------------------------------------------------------


def _write(output_path, title, parts) -> Path:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{_html.escape(title)}</title><style>{_CSS}</style></head><body>"
        f"<h1>{_html.escape(title)}</h1>"
        + "".join(parts)
        + "</body></html>"
    )
    out.write_text(doc, encoding="utf-8")
    return out


def _section(title: str, body_html: str) -> str:
    return f"<h2>{_html.escape(title)}</h2>{body_html}"


def _estimate_box(html_text: str) -> str:
    return f"<div class='estimate'>⚑ Estimate type: {html_text}</div>"


def _warning_box(warnings: list[str]) -> str:
    if not warnings:
        return "<div class='warnbox'><h2>Warnings</h2><p>No warnings.</p></div>"
    items = "".join(f"<li class='warn'>{_html.escape(w)}</li>" for w in warnings)
    return (
        "<div class='warnbox'><h2>⚠ Warnings &amp; caveats "
        f"({len(warnings)})</h2><ul>{items}</ul></div>"
    )


def _kv_table(d: dict, *, raw_values: bool = False) -> str:
    rows = []
    for k, v in d.items():
        val = str(v) if raw_values else _html.escape(str(v))
        rows.append(f"<tr><td>{_html.escape(str(k))}</td><td>{val}</td></tr>")
    return f"<table>{''.join(rows)}</table>"


def _df_html(df: pd.DataFrame) -> str:
    return df.to_html(border=0, na_rep="—")


def _gene_list_html(genes: Sequence[str], limit: int = 50) -> str:
    n = len(genes)
    shown = list(genes[:limit])
    more = f" … (+{n - limit} more)" if n > limit else ""
    return (f"<p>{n} gene(s) selected.</p>"
            f"<p class='meta'>{_html.escape(', '.join(shown))}{more}</p>")


def _bulk_qc_html(qc) -> str:
    blocks = []
    if qc.recon_r2 is not None:
        df = pd.DataFrame({"recon_r2": qc.recon_r2, "profile_corr": qc.profile_corr,
                           "mismatch_flag": qc.mismatch_flag})
        blocks.append("<h3>Per-sample</h3>" + _df_html(df.round(4)))
    if qc.marker_recall is not None:
        df = pd.DataFrame({"marker_recall": qc.marker_recall,
                           "spillover_risk": qc.spillover_risk})
        blocks.append("<h3>Per-cell-type</h3>" + _df_html(df.round(4)))
    if qc.condition_number is not None:
        blocks.append(_kv_table({"condition number κ(Φ)": round(qc.condition_number, 2)}))
    return "".join(blocks) or "<p>No QC metrics.</p>"


def _bulk_uncertainty_html(deconv, qc) -> str:
    if deconv.lower_ci is None:
        return "<p>No bootstrap confidence intervals were computed.</p>"
    note = ""
    if qc is not None and qc.mean_ci_width is not None:
        note = "<h3>Mean CI width</h3>" + _df_html(
            qc.mean_ci_width.to_frame("mean_ci_width").round(4))
    return "<p>Bootstrap confidence intervals available.</p>" + note


def _separability_html(report) -> str:
    cls = "fail" if report.has_problems else "ok"
    head = (f"<p class='{cls}'>CRITICAL: {report.n_critical}, HIGH: {report.n_high}, "
            f"MEDIUM: {report.n_medium}, OK: {report.n_ok}</p>")
    probs = [p for p in report.pairs if getattr(p, "is_problematic", False)]
    if not probs:
        return head + "<p>No poorly-separable pairs.</p>"
    rows = "".join(
        f"<tr><td>{_html.escape(p.type_a)} vs {_html.escape(p.type_b)}</td>"
        f"<td>{p.risk_level}</td><td>{p.bhattacharyya_coeff:.4f}</td></tr>"
        for p in probs
    )
    return head + ("<table><tr><th>pair</th><th>risk</th><th>BC</th></tr>"
                   f"{rows}</table>")


def _figures_html(figures: Sequence[Any]) -> str:
    blocks = []
    for i, fig in enumerate(figures):
        png = None
        caption = None
        # Accept PlotResult, a path, or a (path, caption) pair.
        if hasattr(fig, "figure_paths"):
            png = fig.figure_paths.get("png")
            caption = getattr(fig, "caption", None)
        elif isinstance(fig, (str, Path)):
            png = Path(fig)
        if png is not None and Path(png).exists():
            b64 = base64.b64encode(Path(png).read_bytes()).decode("ascii")
            blocks.append(f"<img src='data:image/png;base64,{b64}' alt='figure {i}'>")
        if caption:
            blocks.append(f"<p class='caption'>{_html.escape(caption)}</p>")
    return "".join(blocks) or "<p>No figures.</p>"


def _file_list_html(files: Iterable[Union[str, Path]]) -> str:
    items = "".join(f"<li class='meta'>{_html.escape(str(f))}</li>" for f in files)
    return f"<ul>{items}</ul>"


def _meta_html(meta: dict) -> str:
    if not meta:
        return "<p>No metadata.</p>"
    rows = "".join(
        f"<tr><td>{_html.escape(str(k))}</td><td class='meta'>{_html.escape(str(v))}</td></tr>"
        for k, v in meta.items()
    )
    return f"<table>{rows}</table>"
