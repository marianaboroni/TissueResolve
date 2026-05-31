"""
HTML templates and reusable fragments for TissueResolve reports.

Plain f-string / string building with embedded CSS — no Jinja2 required (it is
an optional convenience only).  Figures are embedded as ``<iframe>`` references
to their interactive HTML files (which load plotly.js from CDN), keeping the
report light while staying interactive; static images and source data are
linked alongside.
"""
from __future__ import annotations

import html as _html
from typing import Iterable, Mapping, Optional

__all__ = [
    "CSS",
    "page",
    "section",
    "kv_table",
    "df_table",
    "warning_box",
    "estimate_box",
    "figure_block",
    "file_list",
    "escape",
]

escape = _html.escape

CSS = """
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       margin: 0 auto; max-width: 1080px; padding: 24px; color: #1a1a1a; }
h1 { border-bottom: 3px solid #4C72B0; padding-bottom: 8px; }
h2 { margin-top: 34px; color: #2a2a2a; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
h3 { margin-top: 18px; color: #333; }
table { border-collapse: collapse; margin: 8px 0; font-size: 13px; }
th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
th { background: #f0f3f8; } td:first-child, th:first-child { text-align: left; }
.warnbox { background: #fff4e5; border: 1px solid #dd8452; border-radius: 6px;
           padding: 12px 16px; margin: 16px 0; }
.warnbox h2 { color: #b35900; border: none; margin-top: 0; }
.warn { color: #b35900; } .fail { color: #c0392b; font-weight: bold; }
.ok { color: #1e7e34; font-weight: bold; }
.estimate { background: #eef3fb; border-left: 4px solid #4C72B0; padding: 10px 14px;
            margin: 12px 0; font-size: 14px; }
.caption { font-size: 12px; color: #555; margin: 4px 0 16px; }
.meta { font-family: monospace; font-size: 12px; }
iframe { width: 100%; height: 560px; border: 1px solid #eee; border-radius: 4px; }
.links a { margin-right: 12px; font-size: 12px; }
.toc a { display: inline-block; margin: 2px 10px 2px 0; font-size: 13px; }
"""


def page(title: str, body: str, *, toc: Optional[Iterable[str]] = None) -> str:
    toc_html = ""
    if toc:
        items = " ".join(
            f'<a href="#sec{i}">{escape(t)}</a>' for i, t in enumerate(toc))
        toc_html = f'<div class="toc"><b>Contents:</b> {items}</div>'
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{escape(title)}</title><style>{CSS}</style></head><body>"
        f"<h1>{escape(title)}</h1>{toc_html}{body}</body></html>"
    )


def section(index: int, title: str, body_html: str) -> str:
    return f'<h2 id="sec{index}">{index + 1}. {escape(title)}</h2>{body_html}'


def estimate_box(html_text: str) -> str:
    return f"<div class='estimate'>⚑ Estimate type: {html_text}</div>"


def warning_box(warnings: Iterable[str]) -> str:
    warnings = list(warnings)
    if not warnings:
        return "<div class='warnbox'><h2>Warnings &amp; limitations</h2><p>No warnings.</p></div>"
    items = "".join(f"<li class='warn'>{escape(str(w))}</li>" for w in warnings)
    return ("<div class='warnbox'><h2>⚠ Warnings &amp; limitations "
            f"({len(warnings)})</h2><ul>{items}</ul></div>")


def kv_table(d: Mapping, *, raw_values: bool = False) -> str:
    rows = []
    for k, v in d.items():
        val = str(v) if raw_values else escape(str(v))
        rows.append(f"<tr><td>{escape(str(k))}</td><td>{val}</td></tr>")
    return f"<table>{''.join(rows)}</table>" if rows else "<p>—</p>"


def df_table(df, *, max_rows: int = 30, float_fmt: str = "%.4f") -> str:
    if df is None or len(df) == 0:
        return "<p>—</p>"
    shown = df.head(max_rows)
    html = shown.to_html(border=0, na_rep="—", float_format=lambda x: float_fmt % x)
    if len(df) > max_rows:
        html += f"<p class='caption'>(first {max_rows} of {len(df)} rows shown)</p>"
    return html


def figure_block(title: str, *, iframe_src: Optional[str] = None,
                 caption: Optional[str] = None,
                 links: Optional[Mapping[str, str]] = None) -> str:
    parts = [f"<h3>{escape(title)}</h3>"]
    if iframe_src:
        parts.append(f"<iframe src='{escape(iframe_src)}'></iframe>")
    if links:
        link_html = " ".join(
            f"<a href='{escape(href)}'>{escape(label)}</a>"
            for label, href in links.items())
        parts.append(f"<div class='links'>{link_html}</div>")
    if caption:
        parts.append(f"<p class='caption'>{escape(caption)}</p>")
    return "".join(parts)


def file_list(files: Iterable[str]) -> str:
    items = "".join(f"<li class='meta'>{escape(str(f))}</li>" for f in files)
    return f"<ul>{items}</ul>" if items else "<p>—</p>"
