"""
Spatial (10x Visium) deconvolution plots.

Every function:

* takes a :class:`~tissueresolve.results.SpatialDeconvResult` (or its
  proportion table / QC frame) plus the **array coordinates**,
* saves the figure and the underlying data **including the coordinates used**,
* returns a :class:`~tissueresolve.plotting.style.PlotResult`,
* never silently transforms or smooths the estimates — what is plotted is the
  recorded proportion matrix, and λ_spatial is surfaced in the caption.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence, Union

import numpy as np
import pandas as pd

from tissueresolve.plotting import captions
from tissueresolve.plotting.style import (
    PlotResult,
    apply_style,
    get_palette,
    new_figure,
    save_outputs,
)
from tissueresolve.results import SpatialDeconvResult

__all__ = [
    "abundance_map",
    "abundance_map_multi",
    "dominant_type_map",
    "spatial_qc_map",
    "morans_i_barplot",
]

_COMP_COMMENT = [
    "estimate_type: spot_rna_composition",
    "WARNING: RNA-derived composition estimates, NOT single-cell counts.",
]


def _coords_frame(array_row, array_col) -> pd.DataFrame:
    return pd.DataFrame({
        "array_row": np.asarray(array_row),
        "array_col": np.asarray(array_col),
    })


def _scatter(ax, array_col, array_row, values, *, cmap="viridis", vmin=None, vmax=None):
    # Visium array_row increases downward; invert y so the map reads like tissue.
    sc = ax.scatter(
        np.asarray(array_col), np.asarray(array_row),
        c=values, cmap=cmap, s=18, marker="h",
        vmin=vmin, vmax=vmax, linewidths=0,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("array_col")
    ax.set_ylabel("array_row")
    ax.grid(visible=False)
    return sc


def abundance_map(
    result: SpatialDeconvResult,
    array_row: Sequence[float],
    array_col: Sequence[float],
    cell_type: str,
    output_dir: Union[str, Path],
    *,
    name: str | None = None,
) -> PlotResult:
    """Spatial abundance map of one cell type."""
    if cell_type not in result.proportions.columns:
        raise ValueError(
            f"cell_type {cell_type!r} not in proportions {list(result.proportions.columns)}."
        )
    vals = result.proportions[cell_type].to_numpy()
    name = name or f"spatial_abundance_{cell_type}"

    fig, ax = new_figure(figsize=(5.0, 4.5))
    sc = _scatter(ax, array_col, array_row, vals, vmin=0, vmax=1)
    ax.set_title(f"{cell_type} — RNA-derived composition")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("spot composition")

    data = _coords_frame(array_row, array_col)
    data[cell_type] = vals
    data.index = result.proportions.index
    fig_paths, data_paths = save_outputs(
        fig, data, output_dir, name,
        data_comment=_COMP_COMMENT + [f"lambda_spatial: {result.lambda_spatial}"],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.spatial_abundance_caption(
            cell_type, lambda_spatial=result.lambda_spatial, n_spots=result.n_spots,
        ),
    )


def abundance_map_multi(
    result: SpatialDeconvResult,
    array_row: Sequence[float],
    array_col: Sequence[float],
    cell_types: Sequence[str],
    output_dir: Union[str, Path],
    *,
    name: str = "spatial_abundance_multi",
    ncols: int = 3,
) -> PlotResult:
    """Multi-panel spatial abundance map for selected cell types."""
    cts = list(cell_types)
    missing = [c for c in cts if c not in result.proportions.columns]
    if missing:
        raise ValueError(f"cell types not in proportions: {missing}")

    n = len(cts)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    apply_style()
    fig, axes = new_figure(figsize=(4.0 * ncols, 3.6 * nrows), nrows=nrows, ncols=ncols)
    axes = np.atleast_1d(axes).ravel()

    for i, ct in enumerate(cts):
        sc = _scatter(axes[i], array_col, array_row,
                      result.proportions[ct].to_numpy(), vmin=0, vmax=1)
        axes[i].set_title(str(ct))
        fig.colorbar(sc, ax=axes[i], fraction=0.046, pad=0.04)
    for j in range(n, len(axes)):
        axes[j].axis("off")
    fig.suptitle("Spatial RNA-derived composition")

    data = _coords_frame(array_row, array_col)
    for ct in cts:
        data[ct] = result.proportions[ct].to_numpy()
    data.index = result.proportions.index
    fig_paths, data_paths = save_outputs(
        fig, data, output_dir, name,
        data_comment=_COMP_COMMENT + [f"lambda_spatial: {result.lambda_spatial}"],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.spatial_abundance_caption(
            cts, lambda_spatial=result.lambda_spatial, n_spots=result.n_spots,
        ),
    )


def dominant_type_map(
    result: SpatialDeconvResult,
    array_row: Sequence[float],
    array_col: Sequence[float],
    output_dir: Union[str, Path],
    *,
    name: str = "spatial_dominant_type",
) -> PlotResult:
    """Map of the dominant (argmax) cell type per spot."""
    props = result.proportions
    dom_idx = props.to_numpy().argmax(axis=1)
    dom_frac = props.to_numpy().max(axis=1)
    cell_types = list(props.columns)
    colours = get_palette(len(cell_types))

    fig, ax = new_figure(figsize=(5.5, 4.5))
    ax.set_aspect("equal")
    ac = np.asarray(array_col)
    ar = np.asarray(array_row)
    for k, ct in enumerate(cell_types):
        mask = dom_idx == k
        if mask.any():
            ax.scatter(ac[mask], ar[mask], s=18, marker="h",
                       color=colours[k], label=str(ct), linewidths=0)
    ax.invert_yaxis()
    ax.set_xlabel("array_col")
    ax.set_ylabel("array_row")
    ax.set_title("Dominant cell type per spot")
    ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", title="Cell type")
    ax.grid(visible=False)

    data = _coords_frame(array_row, array_col)
    data["dominant_type"] = [cell_types[i] for i in dom_idx]
    data["dominant_fraction"] = dom_frac
    data.index = props.index
    fig_paths, data_paths = save_outputs(
        fig, data, output_dir, name,
        data_comment=_COMP_COMMENT + [f"lambda_spatial: {result.lambda_spatial}"],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.dominant_type_caption(
            lambda_spatial=result.lambda_spatial, n_spots=result.n_spots,
        ),
    )


def spatial_qc_map(
    spot_qc: pd.DataFrame,
    array_row: Sequence[float],
    array_col: Sequence[float],
    metric: str,
    output_dir: Union[str, Path],
    *,
    lambda_spatial: float | None = None,
    name: str | None = None,
) -> PlotResult:
    """Map a per-spot QC metric (e.g. ``nb_loglik``, ``spatial_resid``,
    ``prop_entropy``) at the array coordinates."""
    if metric not in spot_qc.columns:
        raise ValueError(
            f"metric {metric!r} not in spot_qc columns {list(spot_qc.columns)}."
        )
    vals = spot_qc[metric].to_numpy()
    name = name or f"spatial_qc_{metric}"

    fig, ax = new_figure(figsize=(5.0, 4.5))
    sc = _scatter(ax, array_col, array_row, vals, cmap="magma")
    ax.set_title(f"Spatial QC: {metric}")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(metric)

    data = _coords_frame(array_row, array_col)
    data[metric] = vals
    comment = ["Per-spot spatial QC metric at recorded coordinates."]
    if lambda_spatial is not None:
        comment.append(f"lambda_spatial: {lambda_spatial}")
    fig_paths, data_paths = save_outputs(fig, data, output_dir, name, data_comment=comment)
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.spatial_qc_caption(metric, lambda_spatial=lambda_spatial),
    )


def morans_i_barplot(
    morans_i: pd.Series,
    output_dir: Union[str, Path],
    *,
    name: str = "spatial_morans_i",
) -> PlotResult:
    """Bar plot of Moran's I spatial autocorrelation per cell type."""
    s = morans_i.sort_values(ascending=False)
    fig, ax = new_figure(figsize=(max(5.0, 0.5 * len(s) + 2), 4.0))
    x = np.arange(len(s))
    ax.bar(x, s.to_numpy(), color="#55A868")
    ax.axhline(0, color="#8C8C8C", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(s.index.astype(str), rotation=90)
    ax.set_ylabel("Moran's I")
    ax.set_title("Spatial autocorrelation (Moran's I) per cell type")
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(
        fig, s.to_frame("morans_i"), output_dir, name,
        data_comment=["Moran's I per cell-type composition map."],
    )
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)
