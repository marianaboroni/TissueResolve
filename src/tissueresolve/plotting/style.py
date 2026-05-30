"""
Unified plotting style and figure/data persistence for TissueResolve.

Design rules (enforced here so every plot module inherits them):

* **matplotlib only** — seaborn is never required (it may be used elsewhere
  but this package must work without it).
* **Every figure has saved data.**  :func:`save_outputs` writes both the
  figure (PNG/PDF/SVG) *and* the underlying plotting data (TSV) and returns
  the paths.  A figure without saved data violates the DESIGN_SPEC.
* **Every plotting function returns a** :class:`PlotResult` carrying the
  figure object plus the figure and data paths.

The matplotlib ``Agg`` backend is selected on import so plotting works
headless (CI, servers) without a display.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence, Union

import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")  # headless / deterministic
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    _HAS_MPL = True
except Exception:  # pragma: no cover - exercised only without matplotlib
    _HAS_MPL = False
    Figure = Any  # type: ignore[assignment,misc]

__all__ = [
    "DEFAULT_FORMATS",
    "PlotResult",
    "require_matplotlib",
    "apply_style",
    "get_palette",
    "new_figure",
    "save_outputs",
]

#: Figure formats saved by default.  All three are publication-usable.
DEFAULT_FORMATS: tuple[str, ...] = ("png", "pdf", "svg")

#: Base qualitative palette (colour-blind-friendly, seaborn-free).
_BASE_PALETTE: tuple[str, ...] = (
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
)

_STYLE: dict[str, Any] = {
    "figure.figsize": (6.0, 4.0),
    "figure.dpi": 110,
    # Headless batch plotting opens many figures; this is matplotlib's memory
    # advisory, not a scientific warning. Disable the count cap, not warnings.
    "figure.max_open_warning": 0,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
}


@dataclass
class PlotResult:
    """Outputs of a single plotting function.

    Attributes
    ----------
    figure:
        The matplotlib ``Figure`` object (still open; caller may further edit
        or close it).
    figure_paths:
        Mapping ``format -> Path`` of every saved figure file.
    data_paths:
        Mapping ``name -> Path`` of every saved source-data TSV.  Always
        non-empty: a figure without underlying data is not produced.
    caption:
        Optional auto-generated caption string.
    """

    figure: Figure
    figure_paths: dict[str, Path] = field(default_factory=dict)
    data_paths: dict[str, Path] = field(default_factory=dict)
    caption: str | None = None


def require_matplotlib() -> None:
    """Raise a clear error when matplotlib is not installed."""
    if not _HAS_MPL:
        raise ImportError(
            "Plotting requires matplotlib.  Install the plotting extra:\n"
            "    pip install 'tissueresolve[plots]'"
        )


def apply_style() -> None:
    """Apply the shared rcParams.  Idempotent; safe to call repeatedly."""
    require_matplotlib()
    plt.rcParams.update(_STYLE)


def get_palette(n: int) -> list[str]:
    """Return *n* distinct hex colours (no seaborn dependency).

    Uses the base qualitative palette for ``n ≤ 10``, matplotlib ``tab20``
    for ``n ≤ 20``, and evenly-spaced HSV samples beyond that.
    """
    require_matplotlib()
    if n <= len(_BASE_PALETTE):
        return list(_BASE_PALETTE[:n])
    if n <= 20:
        cmap = plt.get_cmap("tab20")
        return [matplotlib.colors.to_hex(cmap(i / 20)) for i in range(n)]
    cmap = plt.get_cmap("hsv")
    return [matplotlib.colors.to_hex(cmap(i / max(n, 1))) for i in range(n)]


def new_figure(
    *, figsize: tuple[float, float] | None = None, nrows: int = 1, ncols: int = 1
):
    """Create a styled figure and axes.  Returns ``(fig, ax_or_axes)``."""
    apply_style()
    fig, ax = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=figsize or plt.rcParams["figure.figsize"],
    )
    return fig, ax


def save_outputs(
    fig: Figure,
    data: Union[pd.DataFrame, pd.Series, Mapping[str, Union[pd.DataFrame, pd.Series]]],
    output_dir: Union[str, Path],
    name: str,
    *,
    formats: Sequence[str] = DEFAULT_FORMATS,
    data_comment: Sequence[str] | None = None,
) -> tuple[dict[str, Path], dict[str, Path]]:
    """Save *fig* in each requested format and *data* as TSV.

    Parameters
    ----------
    fig:
        Figure to save.
    data:
        A DataFrame/Series, or a mapping of suffix → DataFrame/Series for
        plots backed by several tables.  **Never omitted** — every figure
        must ship its source data.
    output_dir:
        Destination directory (created if needed).
    name:
        Base filename (without extension).
    formats:
        Figure formats to write (subset of png/pdf/svg).
    data_comment:
        Optional ``#``-prefixed comment lines prepended to each TSV (e.g.
        the estimate type, smoothing parameters).

    Returns
    -------
    (figure_paths, data_paths)
    """
    require_matplotlib()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    figure_paths: dict[str, Path] = {}
    for fmt in formats:
        p = out / f"{name}.{fmt}"
        fig.savefig(p, format=fmt)
        figure_paths[fmt] = p

    if isinstance(data, Mapping):
        tables = dict(data)
    else:
        tables = {"": data}

    data_paths: dict[str, Path] = {}
    for suffix, table in tables.items():
        df = table.to_frame() if isinstance(table, pd.Series) else table
        fname = f"{name}{('_' + suffix) if suffix else ''}.tsv"
        p = out / fname
        _write_tsv(df, p, data_comment)
        data_paths[suffix or name] = p

    return figure_paths, data_paths


def _write_tsv(df: pd.DataFrame, path: Path, comment: Sequence[str] | None) -> None:
    with path.open("w", encoding="utf-8") as fh:
        if comment:
            for line in comment:
                fh.write(f"# {line}\n")
        df.to_csv(fh, sep="\t")
