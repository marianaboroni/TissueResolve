"""
Minimal, dependency-light benchmark plotting.

Uses Plotly when available (interactive HTML); otherwise writes a small static
HTML table.  **Every figure also writes its source data** as ``<name>.data.tsv``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io import save_plot_data

__all__ = ["bar_figure", "heatmap_figure", "grouped_bar_figure"]


def _have_plotly() -> bool:
    try:
        import plotly.express  # noqa: F401
        return True
    except Exception:
        return False


def _fallback_html(df: pd.DataFrame, title: str, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    html = (f"<h2>{title}</h2><p><i>Plotly not installed — static table; "
            f"source data in {path.with_suffix('.data.tsv').name}.</i></p>"
            + df.to_html(border=0))
    path.write_text(html, encoding="utf-8")
    return path


def bar_figure(df: pd.DataFrame, *, x: str, y: str, title: str, path: Path,
               color: str | None = None) -> Path:
    save_plot_data(df, path)
    if not _have_plotly():
        return _fallback_html(df, title, path)
    import plotly.express as px
    fig = px.bar(df, x=x, y=y, color=color, title=title)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn")
    return Path(path)


def grouped_bar_figure(df: pd.DataFrame, *, x: str, y: str, color: str,
                       title: str, path: Path) -> Path:
    save_plot_data(df, path)
    if not _have_plotly():
        return _fallback_html(df, title, path)
    import plotly.express as px
    fig = px.bar(df, x=x, y=y, color=color, barmode="group", title=title)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn")
    return Path(path)


def heatmap_figure(df: pd.DataFrame, *, title: str, path: Path) -> Path:
    save_plot_data(df, path)
    if not _have_plotly():
        return _fallback_html(df, title, path)
    import plotly.express as px
    fig = px.imshow(df, title=title, aspect="auto", color_continuous_scale="Viridis")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn")
    return Path(path)
