"""
Figure export for the publication layer.

Every figure goes out as **interactive HTML** plus **source data TSV**; static
PDF/SVG/PNG are written when ``kaleido`` is available, otherwise a clear
warning is recorded and HTML + data are still saved (never a silent failure,
never a figure without its data).

All plot functions return a :class:`FigureResult`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence, Union

import pandas as pd

__all__ = [
    "DEFAULT_STATIC_FORMATS",
    "FigureResult",
    "kaleido_available",
    "require_plotly",
    "export_figure",
    "save_source_data",
]

DEFAULT_STATIC_FORMATS: tuple[str, ...] = ("pdf", "svg", "png")


@dataclass
class FigureResult:
    """Outputs of one publication figure.

    Attributes
    ----------
    figure:
        The Plotly ``Figure`` (still open; caller may edit/show it).
    html_path:
        Path to the interactive HTML (always written).
    static_paths:
        ``format -> Path`` for each static image actually written (may be empty
        when kaleido is unavailable).
    data_paths:
        ``name -> Path`` of source-data TSVs.  Always non-empty.
    warnings:
        Non-fatal issues (e.g. "kaleido not installed; static export skipped").
    caption:
        Auto-generated caption string.
    """

    figure: Any
    html_path: Path | None = None
    static_paths: dict[str, Path] = field(default_factory=dict)
    data_paths: dict[str, Path] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    caption: str | None = None


def require_plotly():
    """Import and return ``plotly.graph_objects`` or raise a clear error."""
    try:
        import plotly.graph_objects as go

        return go
    except Exception as exc:  # pragma: no cover - only without plotly
        raise ImportError(
            "Publication figures require plotly.  Install the report extra:\n"
            "    pip install 'tissueresolve[report]'"
        ) from exc


def kaleido_available() -> bool:
    """Whether static image export (PDF/SVG/PNG) is available."""
    try:
        import kaleido  # noqa: F401

        return True
    except Exception:
        return False


def _write_tsv(df: pd.DataFrame, path: Path, comment: Sequence[str] | None) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for line in comment or []:
            fh.write(f"# {line}\n")
        df.to_csv(fh, sep="\t")


def save_source_data(
    data: Union[pd.DataFrame, pd.Series, Mapping[str, Union[pd.DataFrame, pd.Series]]],
    out_dir: Union[str, Path],
    name: str,
    *,
    comment: Sequence[str] | None = None,
) -> dict[str, Path]:
    """Write a figure's source data as ``<name>.data.tsv`` (or per-suffix TSVs)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tables = dict(data) if isinstance(data, Mapping) else {"data": data}
    paths: dict[str, Path] = {}
    for suffix, table in tables.items():
        df = table.to_frame() if isinstance(table, pd.Series) else table
        p = out / f"{name}.{suffix}.tsv"
        _write_tsv(df, p, comment)
        paths[suffix] = p
    return paths


def export_figure(
    fig: Any,
    out_dir: Union[str, Path],
    name: str,
    *,
    data: Union[pd.DataFrame, pd.Series, Mapping[str, Any]],
    formats: Sequence[str] = DEFAULT_STATIC_FORMATS,
    caption: str | None = None,
    data_comment: Sequence[str] | None = None,
    write_html: bool = True,
) -> FigureResult:
    """Write *fig* (HTML + static when possible) and its *data* (TSV).

    HTML and source data are always written.  Static formats are written only
    when kaleido is importable; otherwise a warning is recorded.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    html_path = None
    if write_html:
        html_path = out / f"{name}.html"
        try:
            fig.write_html(str(html_path), include_plotlyjs="cdn",
                           full_html=True)
        except Exception as exc:  # pragma: no cover
            warnings.append(f"HTML export failed: {exc}")
            html_path = None

    static_paths: dict[str, Path] = {}
    if formats:
        if kaleido_available():
            for fmt in formats:
                p = out / f"{name}.{fmt}"
                try:
                    fig.write_image(str(p), format=fmt)
                    static_paths[fmt] = p
                except Exception as exc:  # pragma: no cover
                    warnings.append(f"static {fmt} export failed: {exc}")
        else:
            warnings.append(
                "kaleido not installed; static export (PDF/SVG/PNG) skipped. "
                "Interactive HTML and source data were written. "
                "Install with: pip install kaleido"
            )

    data_paths = save_source_data(data, out, name, comment=data_comment)

    return FigureResult(
        figure=fig, html_path=html_path, static_paths=static_paths,
        data_paths=data_paths, warnings=warnings, caption=caption,
    )
