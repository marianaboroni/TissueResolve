"""
Benchmark HTML reports: per-modality reports + a single summary report that
links to them (so the user opens one file).
"""
from __future__ import annotations

import html as _html
from pathlib import Path

import pandas as pd

from .io import OUTPUTS_DIR, ensure_dir

_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Arial,sans-serif;max-width:1000px;
margin:0 auto;padding:24px;color:#1a1a1a}
h1{border-bottom:2px solid #4C72B0;padding-bottom:6px}
h2{margin-top:26px;border-bottom:1px solid #ddd;padding-bottom:4px}
table{border-collapse:collapse;font-size:13px;margin:8px 0}
th,td{border:1px solid #ccc;padding:4px 8px;text-align:right}
td:first-child,th:first-child{text-align:left}th{background:#f0f3f8}
.note{background:#eef3fb;border-left:4px solid #4C72B0;padding:8px 12px;margin:12px 0}
.warn{background:#fff4e5;border:1px solid #dd8452;border-radius:6px;padding:10px 14px;margin:12px 0}
a{color:#2a6}
"""


def _page(title: str, body: str) -> str:
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title>"
            f"<style>{_CSS}</style></head><body><h1>{title}</h1>{body}</body></html>")


def _df(df: pd.DataFrame) -> str:
    return df.to_html(border=0, na_rep="—")


def _section(title: str, body: str) -> str:
    return f"<h2>{title}</h2>{body}"


def build_modality_report(modality: str, sections: dict[str, object],
                          out_path: Path) -> Path:
    """Sections: {title: DataFrame|html-str}.  Renders + writes the report."""
    out_path = Path(out_path)
    ensure_dir(out_path.parent)
    body = []
    intro = {
        "bulk": ("Bulk pseudobulk has known ground truth, so accuracy metrics are "
                 "valid here."),
        "spatial": ("Real spatial Visium has <b>no ground truth</b>: this report "
                    "evaluates concordance, spatial structure, stability, runtime "
                    "and failure modes, not absolute accuracy.  Synthetic spatial "
                    "mixtures (when used) do have ground truth."),
    }.get(modality, "")
    if intro:
        body.append(f"<div class='note'>{intro}</div>")
    for title, content in sections.items():
        html = _df(content) if isinstance(content, pd.DataFrame) else str(content)
        body.append(_section(title, html))
    out_path.write_text(_page(f"TissueResolve benchmark — {modality}", "".join(body)),
                        encoding="utf-8")
    return out_path


def build_summary_report(bulk_path: Path | None, spatial_path: Path | None,
                         highlights: dict, out_path: Path) -> Path:
    out_path = Path(out_path)
    ensure_dir(out_path.parent)
    body = ["<div class='note'>This is the single entry point — the detailed "
            "bulk and spatial reports are linked below.</div>"]
    links = []
    if bulk_path is not None:
        rel = Path(bulk_path).relative_to(out_path.parent) if Path(bulk_path).is_absolute() else bulk_path
        links.append(f"<li><a href='{rel}'>Bulk benchmark report</a></li>")
    if spatial_path is not None:
        rel = Path(spatial_path).relative_to(out_path.parent) if Path(spatial_path).is_absolute() else spatial_path
        links.append(f"<li><a href='{rel}'>Spatial benchmark report</a></li>")
    body.append(_section("Reports", "<ul>" + "".join(links) + "</ul>"))
    if highlights:
        rows = "".join(f"<tr><td>{_html.escape(str(k))}</td>"
                       f"<td>{_html.escape(str(v))}</td></tr>"
                       for k, v in highlights.items())
        body.append(_section("Highlights", f"<table><tr><th>metric</th>"
                             f"<th>value</th></tr>{rows}</table>"))
    body.append("<div class='warn'>Accuracy claims apply only where ground truth "
                "exists (bulk pseudobulk, synthetic spatial). For real Visium, "
                "concordance/structure/stability are reported, not accuracy. "
                "Batch-aware marker selection can improve robustness, but "
                "aggressive batch correction may remove biological signal. "
                "scRNA and snRNA references differ systematically; mixed "
                "references require library-type confounding checks.</div>")
    out_path.write_text(_page("TissueResolve benchmark — summary", "".join(body)),
                        encoding="utf-8")
    return out_path
