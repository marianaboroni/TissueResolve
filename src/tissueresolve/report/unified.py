"""
Unified single-page report (design-system based).

One ``report.html`` with a sticky sidebar that links to every section
(executive summary, reference quality, input data, bulk, spatial, hierarchical,
resolution/spillover, benchmark, reconstruction, warnings, methods, outputs).
Sections embed metric cards, figure cards, methodology/how-to-read/limitation
boxes and collapsible tables, and may link out to detailed sub-reports — so the
user opens **one** clean, navigable file.
"""
from __future__ import annotations

import html as _html
from pathlib import Path
from typing import Optional, Sequence

from tissueresolve.report import components as C

__all__ = ["Section", "build_unified_report"]


class Section:
    """One report section: stable anchor id, title, HTML body, optional links."""

    def __init__(self, anchor: str, title: str, body_html: str,
                 links: Optional[Sequence] = None):
        self.anchor = anchor
        self.title = title
        self.body_html = body_html or ""
        self.links = list(links or [])  # list of (label, href)

    def render(self) -> str:
        body = self.body_html or "<p class='muted'>Not available for this run.</p>"
        if self.links:
            body += "<div style='margin-top:8px'>" + "".join(
                C.source_data_link(lbl, href) for lbl, href in self.links) + "</div>"
        return C.report_section(self.anchor, self.title, body)


def build_unified_report(out_path: Path | str, sections: Sequence[Section], *,
                         title: str = "TissueResolve — analysis report",
                         subtitle: str = "") -> Path:
    """Write the unified single-page report and return its path."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nav = C.navigation_sidebar(
        "TissueResolve", subtitle or "analysis report",
        [(s.anchor, s.title) for s in sections])
    main = "".join(s.render() for s in sections)
    html = C.page(title, nav, main, page_title=title, page_sub=subtitle)
    out_path.write_text(html, encoding="utf-8")
    return out_path
