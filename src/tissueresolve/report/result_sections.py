"""
In-memory result section builders for the canonical report path.

Builds the ordered ``(title, body_html)`` sections for a bulk or spatial
**in-memory pipeline result** (``BulkPipelineResult`` / ``SpatialPipelineResult``).
The canonical :mod:`tissueresolve.report.orchestration` layer renders these
through the unified single-page shell.  This is the in-memory counterpart of the
results-directory builders in :mod:`tissueresolve.report.sections`.

Design rules preserved: warnings are surfaced (never hidden), failed checks are
shown (e.g. ``converged = False``), and the estimate-type caveat leads every
report.
"""
from __future__ import annotations

import base64
import html as _html
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence, Union

import pandas as pd

from tissueresolve.report import methods_text as _mt

__all__ = [
    "bulk_result_sections", "spatial_result_sections", "read_run_metadata",
]


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def bulk_result_sections(
    result: Any,
    *,
    separability: Any = None,
    figures: Optional[Sequence[Any]] = None,
    output_files: Optional[Iterable[Union[str, Path]]] = None,
    methods: Optional[str] = None,
    fig_dir: Optional[Path] = None,
) -> list[tuple[str, str]]:
    """Ordered ``(title, body_html)`` sections for a bulk in-memory result."""
    deconv = result.deconv
    qc = getattr(result, "qc", None)
    meta = getattr(result, "run_metadata", {}) or {}
    risk = getattr(result, "protocol_risk", None)

    warnings = _collect_bulk_warnings(deconv, qc, risk, separability)

    secs: list[tuple[str, str]] = []
    secs.append(("Estimate type", _estimate_box(
        "These are <b>RNA-derived mRNA proportions</b>, <b>not</b> absolute "
        "cell fractions.")))
    secs.append(("Warnings", _warning_box(warnings)))

    res_html = _resolution_mode_html(result)
    if res_html is not None:
        secs.append(("Resolution mode", res_html))
    secs.append(("Input summary", _kv_table({
        "samples": deconv.n_samples,
        "cell types": deconv.n_cell_types,
        "panel genes": len(deconv.gene_panel),
        "reference genes": meta.get("n_reference_genes", "—"),
    })))
    secs.append(("Reference summary", _kv_table({
        "cell types": deconv.n_cell_types,
        "genome": meta.get("genome", "—"),
        "donor-aware": meta.get("donor_aware", "—"),
    })))
    if "n_shared_genes" in meta:
        secs.append(("Gene overlap", _kv_table({
            "shared genes": meta.get("n_shared_genes"),
        })))
    if risk is not None:
        secs.append(("Protocol risk", _kv_table({
            "risk level": getattr(risk, "risk_level", "unknown"),
            "mismatch types": ", ".join(getattr(risk, "mismatch_types", []) or []) or "—",
            "genes excluded": getattr(risk, "n_genes_excluded", 0),
            "genes down-weighted": len(getattr(risk, "downweighted_genes", []) or []),
        })))
    secs.append(("Selected genes", _gene_list_html(deconv.gene_panel)))
    secs.append((
        "Deconvolution estimates (mRNA proportions)",
        _df_html(deconv.proportions.round(4))
        + "<h3>Summary</h3>" + _df_html(deconv.summary().round(4)),
    ))
    hier_html = _hierarchical_html(result, "bulk")
    if hier_html is not None:
        secs.append(("Hierarchical resolution-aware deconvolution", hier_html))
    if qc is not None:
        secs.append(("QC metrics", _bulk_qc_html(qc)))
    secs.append(("Uncertainty", _bulk_uncertainty_html(deconv, qc)))
    if separability is not None:
        secs.append(("Separability", _separability_html(separability)))
    if figures:
        secs.append(("Figures", _figures_html(figures, fig_dir)))
    if output_files:
        secs.append(("Output files", _file_list_html(output_files)))
    secs.append(("Run metadata", _meta_html(meta)))

    methods = methods or _mt.compose_bulk_methods(result)
    secs.append(("Methods",
                 f"<p>{_html.escape(methods)}</p>".replace("\n\n", "</p><p>")))
    return secs


def spatial_result_sections(
    result: Any,
    *,
    separability: Any = None,
    figures: Optional[Sequence[Any]] = None,
    output_files: Optional[Iterable[Union[str, Path]]] = None,
    methods: Optional[str] = None,
    fig_dir: Optional[Path] = None,
) -> list[tuple[str, str]]:
    """Ordered ``(title, body_html)`` sections for a spatial in-memory result."""
    deconv = result.deconv
    meta = getattr(result, "run_metadata", {}) or {}
    morans = getattr(result, "morans_i", None)
    spot_qc = getattr(result, "spot_qc", None)

    warnings = _collect_spatial_warnings(result, deconv, separability)

    secs: list[tuple[str, str]] = []
    secs.append(("Estimate type", _estimate_box(
        "These are <b>spot-level RNA-derived composition estimates</b>, "
        "<b>not</b> direct single-cell counts.")))
    secs.append(("Warnings", _warning_box(warnings)))

    res_html = _resolution_mode_html(result)
    if res_html is not None:
        secs.append(("Resolution mode", res_html))
    secs.append(("Input summary", _kv_table({
        "spots": deconv.n_spots,
        "cell types": deconv.n_cell_types,
        "marker genes": len(deconv.marker_genes),
    })))
    secs.append(("Reference summary", _kv_table({
        "cell types": deconv.n_cell_types,
        "genome": meta.get("genome", "—"),
    })))
    if "n_genes_panel" in meta or "n_shared_genes" in meta:
        secs.append(("Gene overlap", _kv_table({
            "shared genes": meta.get("n_shared_genes", "—"),
            "panel genes": meta.get("n_genes_panel", len(deconv.marker_genes)),
        })))
    secs.append(("Spatial graph", _kv_table({
        "spots": deconv.n_spots,
        "alpha (mix weight)": meta.get("alpha", "—"),
    })))
    secs.append(("Smoothing parameters", _kv_table({
        "lambda_spatial": deconv.lambda_spatial,
        "alpha": meta.get("alpha", "—"),
        "n_smooth (niche)": deconv.n_smooth if deconv.n_smooth is not None else "not applied",
    })))
    conv_cls = "ok" if deconv.converged else "fail"
    secs.append(("Convergence", _kv_table({
        "converged": f'<span class="{conv_cls}">{deconv.converged}</span>',
        "iterations": deconv.n_iter,
        "final δΠ": round(deconv.convergence_trace[-1], 6) if deconv.convergence_trace else "—",
    }, raw_values=True)))
    secs.append((
        "Spot-level estimates (RNA-derived composition)",
        f"<p>{deconv.n_spots} spots × {deconv.n_cell_types} cell types. "
        "Full table written to the output directory.</p>"
        + _df_html(deconv.proportions.head(10).round(4))
        + ("<p class='caption'>(first 10 spots shown)</p>"
           if deconv.n_spots > 10 else ""),
    ))
    hier_html = _hierarchical_html(result, "spatial")
    if hier_html is not None:
        secs.append(("Hierarchical resolution-aware deconvolution", hier_html))
    if spot_qc is not None:
        secs.append(("Spatial QC", _df_html(spot_qc.describe().round(4))))
    if morans is not None:
        secs.append(("Moran's I", _df_html(morans.to_frame("morans_i").round(4))))
    if separability is not None:
        secs.append(("Separability", _separability_html(separability)))
    if figures:
        secs.append(("Figures", _figures_html(figures, fig_dir)))
    if output_files:
        secs.append(("Output files", _file_list_html(output_files)))
    secs.append(("Run metadata", _meta_html(meta)))

    methods = methods or _mt.compose_spatial_methods(result, separability=separability)
    secs.append(("Methods",
                 f"<p>{_html.escape(methods)}</p>".replace("\n\n", "</p><p>")))
    return secs


def read_run_metadata(results_dir: Path) -> dict:
    """Read ``run_metadata.json`` from a results directory (or ``tables/``)."""
    from tissueresolve.report import assets

    for base in (results_dir, results_dir / "tables"):
        m = assets.read_json(base / "run_metadata.json")
        if isinstance(m, dict):
            return m
    return {}


# ---------------------------------------------------------------------------
# Resolution mode + hierarchical (broad→fine) section
# ---------------------------------------------------------------------------


def _resolution_mode_html(result: Any) -> Optional[str]:
    """Render a short statement of whether the run was flat or hierarchical."""
    meta = getattr(result, "run_metadata", {}) or {}
    mode = meta.get("resolution_mode")
    if mode is None:
        deconv = getattr(result, "deconv", None)
        mode = getattr(deconv, "run_metadata", {}).get("resolution_mode") if deconv else None
    if mode is None:
        return None
    is_hier = (mode == "hierarchical") or (getattr(result, "estimates", None) is not None)
    label = "hierarchical (broad → fine)" if is_hier else f"flat ({mode})"
    reason = meta.get("resolution_mode_reason")
    body = (f"<p>This run used <b>{_html.escape(label)}</b> deconvolution.</p>")
    if reason:
        body += f"<p class='caption'>Reason: {_html.escape(str(reason))}</p>"
    if not is_hier:
        body += ("<p class='caption'>Flat mode estimates all fine cell types at "
                 "once.  When broad/fine annotations are available, hierarchical "
                 "broad→fine mode is recommended (it reduces spillover between "
                 "similar subpopulations).</p>")
    return body


def _hierarchical_html(result: Any, modality: str) -> Optional[str]:
    """Render the "Hierarchical resolution-aware deconvolution" section.

    Returns ``None`` when *result* is not a hierarchical result (no
    ``.estimates``), so flat reports are unaffected.
    """
    est = getattr(result, "estimates", None)
    if est is None or not hasattr(est, "family_proportions"):
        return None

    unit = "spots" if modality == "spatial" else "samples"
    intro = (
        "<p>Hierarchical mode first estimates <b>broad cell-type families</b> "
        "and then estimates <b>fine subpopulations within each family</b>.  "
        "This reduces spillover between unrelated compartments and yields more "
        "cautious subtype-level predictions: a family's mass is only split into "
        "subtypes when those subtypes are demonstrably separable within the "
        "family.  Families that are not separable are reported at the broad "
        f"level as <code>unresolved_&lt;family&gt;</code> (averaged over {unit}).</p>"
    )

    # mean family-level composition
    fam_mean = est.family_proportions.mean(axis=0).sort_values(ascending=False)
    fam_html = _df_html(fam_mean.to_frame("mean_proportion").round(4))

    # per-family resolvability + decision
    qc = est.qc.copy()
    if "resolvable" in qc.columns:
        qc["resolvable"] = qc["resolvable"].map(
            lambda b: '<span class="ok">yes</span>' if bool(b)
            else '<span class="warn">no</span>')
    qc_html = _df_html(qc, raw=True)

    # unresolved mass
    unresolved_html = ""
    if est.unresolved_mass is not None and est.unresolved_mass.shape[1] > 0:
        um = est.unresolved_mass.mean(axis=0).to_frame("mean_unresolved_mass").round(4)
        unresolved_html = (
            "<h3>Unresolved family mass</h3>"
            "<p>These families were detected but their subtypes could not be "
            "reliably separated; interpret them at the family level only.</p>"
            + _df_html(um)
        )
    else:
        unresolved_html = ("<p><i>All multi-subtype families were resolvable; "
                           "no unresolved mass was reported.</i></p>")

    # reliable vs family-level interpretation table
    rows = []
    for _, r in est.qc.iterrows():
        fam = r.get("broad_family", "")
        n = int(r.get("n_subtypes", 0)) if not pd.isna(r.get("n_subtypes", 0)) else 0
        reliable = bool(r.get("resolvable", False))
        rows.append({
            "broad_family": fam,
            "n_subtypes": n,
            "interpret_at": "subtype level" if reliable else "family level only",
        })
    interp_html = _df_html(pd.DataFrame(rows)) if rows else ""

    body = (
        intro
        + "<h3>Broad family composition (mean)</h3>" + fam_html
        + "<h3>Within-family resolvability</h3>" + qc_html
        + unresolved_html
        + "<h3>How to interpret each family</h3>" + interp_html
    )
    return body


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


def _df_html(df: pd.DataFrame, *, raw: bool = False) -> str:
    # raw=True renders embedded HTML (e.g. coloured status spans) unescaped.
    return df.to_html(border=0, na_rep="—", escape=not raw)


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


def _figures_html(figures: Sequence[Any], base_dir: Path | None = None) -> str:
    blocks = []
    for i, fig in enumerate(figures):
        png = None
        html_link = None
        caption = None
        # Accept PlotResult / FigureResult, a path, or a (path, caption) pair.
        if hasattr(fig, "figure_paths"):
            png = fig.figure_paths.get("png")
            caption = getattr(fig, "caption", None)
            html_link = getattr(fig, "html_path", None)
        elif isinstance(fig, (str, Path)):
            png = Path(fig)
        elif isinstance(fig, tuple) and len(fig) == 2:
            png, caption = fig
        if png is not None and Path(png).exists():
            b64 = base64.b64encode(Path(png).read_bytes()).decode("ascii")
            blocks.append(f"<img src='data:image/png;base64,{b64}' alt='figure {i}'>")
        elif html_link is not None:
            path = Path(html_link)
            if base_dir is not None:
                try:
                    rel = path.relative_to(base_dir)
                except ValueError:
                    rel = path
            else:
                rel = path
            blocks.append(
                f"<p><a href='{_html.escape(str(rel))}'>Open interactive figure {i}</a></p>"
            )
        if caption:
            blocks.append(f"<p class='caption'>{_html.escape(str(caption))}</p>")
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
