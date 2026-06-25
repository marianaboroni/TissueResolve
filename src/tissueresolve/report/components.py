"""
Reusable HTML report components (design-system aware).

Pure functions returning HTML strings, used to assemble the unified report and
sub-reports with a consistent look: metric cards, methodology/warning/limitation
boxes, figure cards (with full captions + source-data links), collapsible
tables, sidebar navigation, status badges, and a variable dictionary.

These never invent biological conclusions — they present data and guidance.
"""
from __future__ import annotations

import html as _html
from typing import Iterable, Mapping, Optional, Sequence

from tissueresolve.report.style import REPORT_CSS

__all__ = [
    "esc", "metric_card", "metric_grid", "method_box", "warning_box",
    "limitation_box", "how_to_read_box", "interpretation_guide",
    "variable_dictionary", "figure_card", "collapsible_table",
    "source_data_link", "report_section", "navigation_sidebar", "status_badge",
    "methodology_summary", "estimate_note", "page",
]


def esc(x) -> str:
    return _html.escape(str(x))


# --- cards / boxes ----------------------------------------------------------

def metric_card(label: str, value, sub: str = "") -> str:
    sub_html = f"<div class='sub'>{esc(sub)}</div>" if sub else ""
    return (f"<div class='metric-card'><div class='label'>{esc(label)}</div>"
            f"<div class='value'>{esc(value)}</div>{sub_html}</div>")


def metric_grid(cards: Mapping) -> str:
    """cards: {label: value} or {label: (value, sub)}."""
    items = []
    for label, v in cards.items():
        if isinstance(v, (tuple, list)) and len(v) == 2:
            items.append(metric_card(label, v[0], v[1]))
        else:
            items.append(metric_card(label, v))
    return f"<div class='metric-grid'>{''.join(items)}</div>"


def decision_status_card(label: str, status: str, sub: str = "") -> str:
    """A single decision card: a label, a PASS/CAUTION/WARNING/FAIL (or
    yes/no/available) badge, and an optional one-line sub-note."""
    sub_html = f"<div class='sub'>{esc(sub)}</div>" if sub else ""
    return (f"<div class='metric-card status-card'>"
            f"<div class='label'>{esc(label)}</div>"
            f"<div class='value'>{status_badge(status)}</div>{sub_html}</div>")


def decision_status_grid(items: Sequence) -> str:
    """*items*: sequence of (label, status) or (label, status, sub)."""
    cards = []
    for it in items:
        label, status = it[0], it[1]
        sub = it[2] if len(it) > 2 else ""
        cards.append(decision_status_card(label, status, sub))
    return f"<div class='metric-grid'>{''.join(cards)}</div>"


def checklist(items: Sequence) -> str:
    """A "before interpreting results" checklist.

    *items*: sequence of (text, state) where state is True (done/✓),
    False (not done/✗) or None (not applicable/•).
    """
    rows = []
    for text, state in items:
        mark, cls = ("✓", "ok") if state is True else \
            (("✗", "bad") if state is False else ("•", "na"))
        rows.append(f"<li class='check-{cls}'><span class='check-mark'>{mark}</span> "
                    f"{esc(text)}</li>")
    return f"<ul class='checklist'>{''.join(rows)}</ul>"


def method_box(html_text: str) -> str:
    return f"<div class='methods-card'>{html_text}</div>"


def methodology_summary(points: Iterable[str]) -> str:
    lis = "".join(f"<li>{esc(p)}</li>" for p in points)
    return f"<div class='methods-card'><ul style='margin:4px 0 0;padding-left:18px'>{lis}</ul></div>"


def warning_box(warnings: Iterable[str]) -> str:
    """Render warnings.  Returns '' when there are none — callers should only
    add the box when warnings exist (never print a false 'No warnings')."""
    warns = [w for w in warnings if str(w).strip()]
    if not warns:
        return ""
    lis = "".join(f"<li>{esc(w)}</li>" for w in warns)
    return f"<div class='warning-card'><ul style='margin:4px 0 0;padding-left:18px'>{lis}</ul></div>"


def limitation_box(items: Iterable[str]) -> str:
    lis = "".join(f"<li>{esc(i)}</li>" for i in items if str(i).strip())
    return f"<div class='limitation-card'><ul style='margin:4px 0 0;padding-left:18px'>{lis}</ul></div>" if lis else ""


def how_to_read_box(text: str) -> str:
    return f"<div class='how-to-read'>{text}</div>"


def interpretation_guide(text: str) -> str:
    return f"<div class='interpretation-guide'>{text}</div>"


def variable_dictionary(defs: Mapping[str, str]) -> str:
    if not defs:
        return ""
    rows = "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in defs.items())
    return f"<div class='variable-card'><dl>{rows}</dl></div>"


def status_badge(status: str) -> str:
    s = str(status).lower()
    cls = {"pass": "status-pass", "caution": "status-caution",
           "warning": "status-warning", "warn": "status-warning",
           "fail": "status-critical", "critical": "status-critical"}.get(s, "status-caution")
    return f"<span class='badge {cls}'>{esc(status)}</span>"


def estimate_note(text: str) -> str:
    return f"<div class='estimate-note'>{text}</div>"


# --- figures / tables / links -----------------------------------------------

def source_data_link(label: str, href: str) -> str:
    return f"<a class='source-data-link' href='{esc(href)}'>⤓ {esc(label)}</a>"


def figure_card(*, title: str, subtitle: str = "", body_html: str = "",
                caption: str = "", legend: str = "",
                how_to_read: str = "", methodology: str = "",
                variables: Optional[Mapping[str, str]] = None,
                source_links: Optional[Sequence] = None) -> str:
    """A fully-labelled figure card.

    *source_links* is a sequence of (label, href).  Caption answers what is
    plotted / axes / colors / data / method / what to check / limitations.
    """
    parts = ["<div class='figure-card'>",
             "<div class='fig-head'>",
             f"<div class='fig-title'>{esc(title)}</div>"]
    if subtitle:
        parts.append(f"<div class='fig-subtitle'>{esc(subtitle)}</div>")
    parts.append("</div>")
    parts.append(f"<div class='fig-body'>{body_html}</div>")
    if caption:
        parts.append(f"<div class='figure-caption'><b>Figure.</b> {esc(caption)}</div>")
    if legend:
        parts.append(f"<div class='figure-legend'>{esc(legend)}</div>")
    if how_to_read:
        parts.append(how_to_read_box(esc(how_to_read)))
    if methodology:
        parts.append(method_box(esc(methodology)))
    if variables:
        parts.append(variable_dictionary(variables))
    if source_links:
        links = "".join(source_data_link(lbl, href) for lbl, href in source_links)
        parts.append(f"<div>{links}</div>")
    parts.append("</div>")
    return "".join(parts)


def collapsible_table(title: str, table_html: str, *, open: bool = False,
                      note: str = "") -> str:
    note_html = f"<div class='muted'>{esc(note)}</div>" if note else ""
    return (f"<details class='collapsible-table'{' open' if open else ''}>"
            f"<summary>{esc(title)}</summary>{note_html}{table_html}</details>")


# --- layout -----------------------------------------------------------------

def report_section(anchor: str, title: str, body_html: str) -> str:
    return (f"<section id='{esc(anchor)}'><div class='section-card'>"
            f"<h2>{esc(title)}</h2>{body_html}</div></section>")


def navigation_sidebar(title: str, subtitle: str, items: Sequence) -> str:
    """items: sequence of (anchor, label)."""
    links = "".join(
        f"<a href='#{esc(a)}'><span class='navnum'>{i+1:02d}</span>{esc(l)}</a>"
        for i, (a, l) in enumerate(items))
    return (f"<nav class='sidebar-nav'><h1>{esc(title)}</h1>"
            f"<div class='tagline'>{esc(subtitle)}</div>{links}</nav>")


def page(title: str, sidebar_html: str, main_html: str, *,
         page_title: str = "", page_sub: str = "") -> str:
    pt = esc(page_title or title)
    sub = f"<p class='page-sub'>{esc(page_sub)}</p>" if page_sub else ""
    return (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width, initial-scale=1'>"
            f"<title>{esc(title)}</title><style>{REPORT_CSS}</style></head><body>"
            f"<div class='report-container'>{sidebar_html}"
            f"<main><h1 class='page-title'>{pt}</h1>{sub}{main_html}</main>"
            f"</div></body></html>")
