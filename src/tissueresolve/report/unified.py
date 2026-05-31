"""
Unified single-page report.

Builds one ``report.html`` with a sticky top navigation bar that links to every
section (executive summary, reference quality, input data, bulk, spatial,
hierarchical, resolution/spillover, benchmark, warnings, figures, methods,
output files).  Sections may embed content inline and/or link out to the
detailed sub-reports — so the user opens **one** file instead of hunting
through folders.
"""
from __future__ import annotations

import html as _html
from pathlib import Path
from typing import Optional, Sequence

__all__ = ["Section", "build_unified_report"]

_CSS = """
:root{--accent:#4C72B0}
*{box-sizing:border-box}
body{font-family:-apple-system,Segoe UI,Roboto,Arial,sans-serif;margin:0;
color:#1a1a1a;line-height:1.5}
nav{position:sticky;top:0;background:#fff;border-bottom:2px solid var(--accent);
padding:8px 16px;display:flex;flex-wrap:wrap;gap:6px;z-index:10;font-size:13px}
nav a{color:var(--accent);text-decoration:none;padding:3px 8px;border-radius:4px}
nav a:hover{background:#eef3fb}
main{max-width:1000px;margin:0 auto;padding:20px}
section{scroll-margin-top:54px;border-bottom:1px solid #eee;padding:18px 0}
h1{margin:16px 0 4px}h2{color:#222;border-bottom:1px solid #ddd;padding-bottom:4px}
table{border-collapse:collapse;font-size:13px;margin:8px 0}
th,td{border:1px solid #ccc;padding:4px 8px;text-align:right}
td:first-child,th:first-child{text-align:left}th{background:#f0f3f8}
.note{background:#eef3fb;border-left:4px solid var(--accent);padding:8px 12px;margin:10px 0}
.warn{background:#fff4e5;border:1px solid #dd8452;border-radius:6px;padding:10px 14px;margin:10px 0}
.links a{display:inline-block;margin:2px 10px 2px 0}
.muted{color:#777;font-size:12px}
"""


class Section:
    """One report section: a stable anchor id, a title, and HTML body."""

    def __init__(self, anchor: str, title: str, body_html: str,
                 links: Optional[Sequence[tuple]] = None):
        self.anchor = anchor
        self.title = title
        self.body_html = body_html or ""
        self.links = list(links or [])  # list of (label, href)

    def render(self) -> str:
        out = [f"<section id='{self.anchor}'><h2>{_html.escape(self.title)}</h2>"]
        if self.links:
            out.append("<div class='links'>" + " ".join(
                f"<a href='{href}'>↗ {_html.escape(label)}</a>"
                for label, href in self.links) + "</div>")
        out.append(self.body_html or "<p class='muted'>Not available for this run.</p>")
        out.append("</section>")
        return "".join(out)


def build_unified_report(out_path: Path | str, sections: Sequence[Section], *,
                         title: str = "TissueResolve — analysis report",
                         subtitle: str = "") -> Path:
    """Write the unified single-page report and return its path."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nav = " ".join(f"<a href='#{s.anchor}'>{_html.escape(s.title)}</a>" for s in sections)
    body = "".join(s.render() for s in sections)
    sub = f"<p class='muted'>{_html.escape(subtitle)}</p>" if subtitle else ""
    page = (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>{_html.escape(title)}</title><style>{_CSS}</style></head><body>"
            f"<nav>{nav}</nav><main><h1>{_html.escape(title)}</h1>{sub}{body}</main>"
            f"</body></html>")
    out_path.write_text(page, encoding="utf-8")
    return out_path
