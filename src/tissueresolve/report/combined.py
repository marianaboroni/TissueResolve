"""
Combined bulk + spatial report.

Merges two *existing* ``tissueresolve run`` output directories — one bulk, one
spatial — into a single ``report.html`` (plus ``methods.txt``, ``warnings.json``
and ``run_metadata.json``).  Bulk and spatial sections (and any benchmark
summaries) are kept **separate**, and the report states explicitly that it
summarises two independent runs sharing a single-cell reference — **not** a
single joint bulk+spatial model.

This does not run any deconvolution; it reads what the two runs already wrote.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

from tissueresolve.report import assets
from tissueresolve.report import templates as T
from tissueresolve.report.orchestration import render_sections
from tissueresolve.report.result_sections import read_run_metadata

__all__ = ["generate_combined_report"]


# ---------------------------------------------------------------------------
# Run-directory readers (tolerant of missing pieces)
# ---------------------------------------------------------------------------


def _modality_of(run_dir: Path) -> Optional[str]:
    meta = read_run_metadata(run_dir)
    mode = meta.get("mode")
    if mode is None and isinstance(meta.get("analysis_plan"), dict):
        mode = meta["analysis_plan"].get("mode")
    return mode if mode in ("bulk", "spatial") else None


def _read_predictions(run_dir: Path, modality: str):
    """Predictions for a run dir: prefer ``deconv/proportions.tsv`` (what ``run``
    writes), fall back to the report-shaped table names."""
    p = run_dir / "deconv" / "proportions.tsv"
    if p.exists():
        return assets.read_tsv(p)
    name = ("bulk_estimated_proportions.tsv" if modality == "bulk"
            else "spatial_spot_proportions.tsv")
    return assets.find_table(run_dir, name)


def _read_warning_messages(run_dir: Path) -> list[str]:
    """Collect warning strings from a run dir's ``warnings.json`` (list-of-dicts
    or dict-of-lists) and ``qc/recommendations.txt``."""
    msgs: list[str] = []
    wj = run_dir / "warnings.json"
    if wj.exists():
        try:
            data = json.loads(wj.read_text(encoding="utf-8"))
        except Exception:
            data = None
        if isinstance(data, list):
            for w in data:
                if isinstance(w, dict):
                    msgs.append(str(w.get("message", w)))
                else:
                    msgs.append(str(w))
        elif isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, list) and "recommend" in k.lower():
                    msgs.extend(str(x) for x in v)
    rec = run_dir / "qc" / "recommendations.txt"
    if rec.exists():
        msgs.extend(line for line in rec.read_text(encoding="utf-8").splitlines()
                    if line.strip())
    # de-duplicate, preserve order
    seen: set[str] = set()
    out: list[str] = []
    for m in msgs:
        if m not in seen:
            seen.add(m)
            out.append(m)
    return out


def _read_methods(run_dir: Path) -> str:
    for cand in (run_dir / "methods.txt", run_dir / "tables" / "methods.txt"):
        if cand.exists():
            return cand.read_text(encoding="utf-8")
    return ""


def _embed_figures(run_dir: Path, out_dir: Path) -> str:
    """Embed a run directory's interpretive figures (``figures/*.png``, with
    interactive ``.html`` and ``.data.tsv`` linked) as figure cards. Returns an
    empty string when the run produced no figures."""
    fig_dir = run_dir / "figures"
    if not fig_dir.is_dir():
        return ""
    rel = lambda p: os.path.relpath(p, out_dir)  # noqa: E731
    pngs = sorted(fig_dir.glob("*.png"))
    if not pngs:
        return ""
    cards = []
    for png in pngs:
        title = png.stem.replace("_", " ").strip().capitalize()
        links = []
        for ext in ("html", "pdf", "svg"):
            alt = png.with_suffix("." + ext)
            if alt.exists():
                links.append((ext.upper() if ext != "html" else "interactive",
                              rel(alt)))
        data = png.with_suffix(".data.tsv")
        if data.exists():
            links.append(("source data (.tsv)", rel(data)))
        link_html = " ".join(f"<a href='{T.escape(h)}'>{T.escape(l)}</a>"
                             for l, h in links)
        cards.append(
            f"<figure style='margin:12px 0'>"
            f"<figcaption style='font-weight:600'>{T.escape(title)}</figcaption>"
            f"<img src='{T.escape(rel(png))}' alt='{T.escape(title)}' "
            "style='max-width:100%;height:auto;border:1px solid #e3e8ef;"
            "border-radius:6px'/>"
            + (f"<div class='links' style='font-size:12px'>{link_html}</div>"
               if link_html else "")
            + "</figure>")
    return "<h3>Figures</h3>" + "".join(cards)


def _find_benchmark(run_dir: Path) -> Optional[Path]:
    """A benchmark artifact inside the run dir, if any (benchmarks are normally
    produced separately under benchmarks/outputs/, so usually absent here)."""
    for pat in ("*benchmark*report*.html", "*benchmark*.tsv", "benchmark_metadata.json"):
        hits = sorted(run_dir.glob(pat)) + sorted((run_dir / "tables").glob(pat))
        if hits:
            return hits[0]
    return None


# ---------------------------------------------------------------------------
# Combined report
# ---------------------------------------------------------------------------


def generate_combined_report(bulk_dir, spatial_dir, out_dir) -> Path:
    """Build a combined bulk+spatial report from two existing run directories.

    Validates that *bulk_dir* is a bulk run and *spatial_dir* is a spatial run
    (when their metadata records a modality), writes ``report.html`` +
    ``methods.txt`` + ``warnings.json`` + ``run_metadata.json`` under *out_dir*,
    and returns the report path.
    """
    bulk_dir, spatial_dir, out_dir = Path(bulk_dir), Path(spatial_dir), Path(out_dir)
    for d, label in ((bulk_dir, "--bulk-dir"), (spatial_dir, "--spatial-dir")):
        if not d.is_dir():
            raise ValueError(f"{label} {str(d)!r} is not a directory.")

    b_mode, s_mode = _modality_of(bulk_dir), _modality_of(spatial_dir)
    if b_mode == "spatial":
        raise ValueError(f"--bulk-dir {str(bulk_dir)!r} looks like a SPATIAL run "
                         "(its metadata records mode='spatial'); pass it as "
                         "--spatial-dir.")
    if s_mode == "bulk":
        raise ValueError(f"--spatial-dir {str(spatial_dir)!r} looks like a BULK run "
                         "(its metadata records mode='bulk'); pass it as --bulk-dir.")

    b_meta, s_meta = read_run_metadata(bulk_dir), read_run_metadata(spatial_dir)
    b_props = _read_predictions(bulk_dir, "bulk")
    s_props = _read_predictions(spatial_dir, "spatial")
    b_warns = _read_warning_messages(bulk_dir)
    s_warns = _read_warning_messages(spatial_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    rel = lambda p: os.path.relpath(Path(p), out_dir)  # noqa: E731

    secs: list[tuple[str, str]] = []

    # 1. Executive summary -----------------------------------------------------
    b_n = b_props.shape[0] if b_props is not None else "—"
    b_k = b_props.shape[1] if b_props is not None else "—"
    s_n = s_props.shape[0] if s_props is not None else "—"
    s_k = s_props.shape[1] if s_props is not None else "—"
    exec_body = (
        "<div class='estimate'>This combined report summarises <b>two separate "
        "runs</b> — one bulk, one spatial — that share a single-cell reference. "
        "It is <b>not</b> a single joint bulk+spatial model; bulk and spatial "
        "estimates are produced independently and shown side by side.</div>"
        + T.summary_cards({
            "bulk samples": b_n, "bulk cell types": b_k,
            "spatial spots": s_n, "spatial cell types": s_k,
        })
        + "<p>Bulk estimates are RNA-derived mRNA proportions (not absolute cell "
        "fractions); spatial estimates are spot-level RNA-derived composition "
        "(not single-cell counts).  Read each modality's QC before its "
        "predictions.</p>")
    secs.append(("Executive summary", exec_body))

    # 2. Shared reference summary (if metadata available) ----------------------
    ref_rows = {}
    for label, meta in (("bulk run", b_meta), ("spatial run", s_meta)):
        plan = meta.get("analysis_plan") if isinstance(meta, dict) else None
        if isinstance(plan, dict):
            ref_rows[f"{label}: mode"] = plan.get("mode", "—")
            ref_rows[f"{label}: preset"] = plan.get("preset", "—")
            ref_rows[f"{label}: resolution mode"] = plan.get("resolution_mode", "—")
            ref_rows[f"{label}: reference detected"] = plan.get("detected_reference", "—")
    if ref_rows:
        note = ("<p>The two runs are expected to share one single-cell reference. "
                "If the reference detected/used differs between them, the "
                "combined view is not directly comparable.</p>")
        secs.append(("Shared reference summary", T.kv_table(ref_rows) + note))

    # 3. Bulk QC and predictions ----------------------------------------------
    bulk_body = [f"<p>Source: <a href='{T.escape(rel(bulk_dir))}'>"
                 f"{T.escape(str(bulk_dir))}</a>"]
    bp = bulk_dir / "report.html"
    if bp.exists():
        bulk_body.append(f" — <a href='{T.escape(rel(bp))}'>full bulk report</a>")
    bulk_body.append("</p>")
    if b_props is not None:
        bulk_body.append("<h3>Predicted composition (mRNA proportions)</h3>")
        bulk_body.append(T.df_table(b_props, max_rows=20))
    else:
        bulk_body.append("<p>No bulk prediction table found.</p>")
    bulk_body.append(_embed_figures(bulk_dir, out_dir))
    b_recon = assets.read_tsv(bulk_dir / "deconv" / "coverage_r2.tsv")
    if b_recon is not None:
        bulk_body.append("<h3>Reconstruction QC</h3>" + T.df_table(b_recon, max_rows=20))
    if b_warns:
        bulk_body.append("<h3>Bulk warnings</h3>" + T.warning_box(b_warns))
    secs.append(("Bulk QC and predictions", "".join(bulk_body)))

    # 4. Spatial QC and predictions -------------------------------------------
    sp_body = [f"<p>Source: <a href='{T.escape(rel(spatial_dir))}'>"
               f"{T.escape(str(spatial_dir))}</a>"]
    spp = spatial_dir / "report.html"
    if spp.exists():
        sp_body.append(f" — <a href='{T.escape(rel(spp))}'>full spatial report</a>")
    sp_body.append("</p>")
    if s_props is not None:
        sp_body.append("<h3>Predicted composition (spot-level RNA-derived)</h3>")
        sp_body.append(T.df_table(s_props, max_rows=20))
    else:
        sp_body.append("<p>No spatial prediction table found.</p>")
    sp_body.append(_embed_figures(spatial_dir, out_dir))
    morans = assets.find_table(spatial_dir, "morans_i.tsv") \
        or assets.read_tsv(spatial_dir / "qc" / "morans_i.tsv")
    if morans is not None:
        sp_body.append("<h3>Moran's I (spatial structure)</h3>"
                       + T.df_table(morans, max_rows=20))
    if s_warns:
        sp_body.append("<h3>Spatial warnings</h3>" + T.warning_box(s_warns))
    secs.append(("Spatial QC and predictions", "".join(sp_body)))

    # 5. Bulk benchmark summary (if available) --------------------------------
    b_bench = _find_benchmark(bulk_dir)
    if b_bench is not None:
        secs.append(("Bulk benchmark summary",
                     f"<p>Bulk benchmark artifact: "
                     f"<a href='{T.escape(rel(b_bench))}'>{T.escape(b_bench.name)}</a>.</p>"))
    else:
        secs.append(("Bulk benchmark summary",
                     "<p>No bulk benchmark found in this run directory. "
                     "Benchmarks are produced separately under "
                     "<code>benchmarks/outputs/</code>.</p>"))

    # 6. Spatial benchmark summary (if available) -----------------------------
    s_bench = _find_benchmark(spatial_dir)
    if s_bench is not None:
        secs.append(("Spatial benchmark summary",
                     f"<p>Spatial benchmark artifact: "
                     f"<a href='{T.escape(rel(s_bench))}'>{T.escape(s_bench.name)}</a>. "
                     "Real Visium has no spot-level ground truth, so any spatial "
                     "benchmark reports concordance / structure / runtime, never "
                     "accuracy.</p>"))
    else:
        secs.append(("Spatial benchmark summary",
                     "<p>No spatial benchmark found in this run directory. "
                     "Real Visium has no ground truth, so a spatial benchmark "
                     "reports concordance/structure, not accuracy.</p>"))

    # 7. Warnings and limitations (combined, attributable) --------------------
    warn_html = []
    if b_warns:
        warn_html.append("<h3>Bulk</h3>" + T.warning_box(b_warns))
    if s_warns:
        warn_html.append("<h3>Spatial</h3>" + T.warning_box(s_warns))
    if not warn_html:
        warn_html.append("<p>No warnings were recorded by either run.</p>")
    warn_html.append("<p class='caption'>This combined report does not run a "
                     "joint model; each modality's warnings apply only to that "
                     "run.</p>")
    secs.append(("Warnings and limitations", "".join(warn_html)))

    # 8. Methods and source data links ----------------------------------------
    b_methods, s_methods = _read_methods(bulk_dir), _read_methods(spatial_dir)
    methods_text = _combined_methods_text(b_methods, s_methods, bulk_dir, spatial_dir)
    methods_body = (f"<p>{T.escape(methods_text)}</p>".replace("\n\n", "</p><p>")
                    + "<h3>Source data</h3><ul>"
                    + f"<li>Bulk run directory: <a href='{T.escape(rel(bulk_dir))}'>"
                    f"{T.escape(str(bulk_dir))}</a></li>"
                    + f"<li>Spatial run directory: <a href='{T.escape(rel(spatial_dir))}'>"
                    f"{T.escape(str(spatial_dir))}</a></li>"
                    + "</ul>")
    secs.append(("Methods and source data", methods_body))

    # --- render + side files --------------------------------------------------
    out_html = render_sections(out_dir / "report.html",
                               "TissueResolve — Combined bulk + spatial report",
                               secs, subtitle="two runs sharing a reference")

    (out_dir / "methods.txt").write_text(methods_text, encoding="utf-8")

    combined_warns = (
        [{"modality": "bulk", "message": m} for m in b_warns]
        + [{"modality": "spatial", "message": m} for m in s_warns]
        + [{"modality": "combined", "severity": "info",
            "message": "Combined report of two separate runs sharing a "
                       "reference; not a single joint bulk+spatial model."}])
    (out_dir / "warnings.json").write_text(
        json.dumps(combined_warns, indent=2), encoding="utf-8")

    (out_dir / "run_metadata.json").write_text(json.dumps({
        "generated_by": "tissueresolve combine-report",
        "note": "Two separate runs sharing a single-cell reference; not a joint "
                "bulk+spatial model.",
        "bulk_dir": str(bulk_dir),
        "spatial_dir": str(spatial_dir),
        "bulk_run_metadata": b_meta or None,
        "spatial_run_metadata": s_meta or None,
    }, indent=2, default=str), encoding="utf-8")

    return out_html


def _combined_methods_text(bulk_methods: str, spatial_methods: str,
                           bulk_dir: Path, spatial_dir: Path) -> str:
    parts = [
        "Combined report. Bulk and spatial cell-type-level deconvolution were "
        "run separately on the same single-cell reference and summarised "
        "together. This is not a single joint bulk+spatial model.",
    ]
    parts.append("Bulk methods:\n" + (bulk_methods.strip()
                 or f"(see {bulk_dir}/methods.txt)"))
    parts.append("Spatial methods:\n" + (spatial_methods.strip()
                 or f"(see {spatial_dir}/methods.txt)"))
    return "\n\n".join(parts)
