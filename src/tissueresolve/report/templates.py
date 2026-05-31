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
    "severity_warning_box",
    "estimate_box",
    "figure_block",
    "file_list",
    "summary_cards",
    "key_findings",
    "collapsible",
    "info_card",
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
.cards { display: flex; flex-wrap: wrap; gap: 10px; margin: 12px 0; }
.card { border: 1px solid #ddd; border-radius: 8px; padding: 10px 14px;
        min-width: 120px; background: #fafbfc; }
.card .v { font-size: 20px; font-weight: 700; color: #2a4d8f; }
.card .k { font-size: 11px; color: #666; text-transform: uppercase; }
.status-PASS .v { color: #1e7e34; } .status-CAUTION .v { color: #b35900; }
.status-FAIL .v { color: #c0392b; }
.sev { padding: 1px 6px; border-radius: 4px; font-size: 11px; font-weight: 700;
       color: #fff; }
.sev-CRITICAL { background: #c0392b; } .sev-WARNING { background: #d35400; }
.sev-CAUTION { background: #b8860b; } .sev-INFO { background: #5b8db8; }
.findings li { margin: 5px 0; font-size: 14px; }
.infocard { background: #eef3fb; border: 1px solid #b9cdec; border-radius: 6px;
            padding: 10px 14px; margin: 10px 0; font-size: 13px; }
details { margin: 8px 0; } details > summary { cursor: pointer; font-weight: 600;
          color: #2a4d8f; }
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


def summary_cards(cards: Mapping) -> str:
    """Executive-summary cards.  A ``qc_status`` key gets PASS/CAUTION/FAIL styling."""
    out = []
    for k, v in cards.items():
        cls = "card"
        if str(k).lower() in ("qc_status", "qc status", "overall qc status"):
            cls += f" status-{escape(str(v))}"
        out.append(f"<div class='{cls}'><div class='v'>{escape(str(v))}</div>"
                   f"<div class='k'>{escape(str(k).replace('_', ' '))}</div></div>")
    return f"<div class='cards'>{''.join(out)}</div>"


def key_findings(findings: Iterable[str]) -> str:
    items = "".join(f"<li>{escape(str(f))}</li>" for f in findings)
    return f"<ul class='findings'>{items}</ul>" if items else "<p>—</p>"


def severity_warning_box(warnings: Iterable, *, top_n: int = 5) -> str:
    """Warning box ranked by severity.  *warnings* are objects with .severity,
    .message, .recommended_action (e.g. interpretation.Warning_)."""
    warns = list(warnings)
    if not warns:
        return ("<div class='warnbox'><h2>Warnings &amp; limitations</h2>"
                "<p>No warnings.</p></div>")
    by_sev: dict[str, int] = {}
    for w in warns:
        by_sev[w.severity] = by_sev.get(w.severity, 0) + 1
    counts = " ".join(
        f"<span class='sev sev-{escape(s)}'>{escape(s)}: {n}</span>"
        for s, n in by_sev.items())
    items = "".join(
        f"<li><span class='sev sev-{escape(w.severity)}'>{escape(w.severity)}</span> "
        f"{escape(w.message)}"
        + (f" <i>→ {escape(w.recommended_action)}</i>" if w.recommended_action else "")
        + "</li>" for w in warns[:top_n])
    more = (f"<p class='caption'>… and {len(warns) - top_n} more.</p>"
            if len(warns) > top_n else "")
    return ("<div class='warnbox'><h2>⚠ Warnings &amp; limitations "
            f"({len(warns)})</h2><p>{counts}</p><ul>{items}</ul>{more}</div>")


def collapsible(title: str, body_html: str, *, open: bool = False) -> str:
    """A collapsible <details> block (raw tables/matrices go here, closed by default)."""
    attr = " open" if open else ""
    return f"<details{attr}><summary>{escape(title)}</summary>{body_html}</details>"


def info_card(html_text: str) -> str:
    """A neutral info card (e.g. 'bootstrap not computed' instead of an empty plot)."""
    return f"<div class='infocard'>{html_text}</div>"
