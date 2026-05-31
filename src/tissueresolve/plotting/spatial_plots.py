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


# ===========================================================================
# Publication layer (Plotly) — interactive HTML + static + source data
# ===========================================================================

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly  # noqa: E402
from tissueresolve.plotting.style import (  # noqa: E402
    ESTIMATE_SUBTITLE,
    color_sequence,
    plotly_layout,
)

_SP_SUB = ESTIMATE_SUBTITLE["spatial"]


def _collapse_spatial_top_n(props, top_n):
    if not top_n or top_n >= props.shape[1]:
        return props
    order = props.mean(axis=0).sort_values(ascending=False).index.tolist()
    keep = order[:top_n]
    other = [c for c in props.columns if c not in keep]
    out = props[keep].copy()
    if other:
        out["Other"] = props[other].sum(axis=1)
    return out


def _coords_df(array_row, array_col, index=None):
    return pd.DataFrame({"array_row": np.asarray(array_row),
                         "array_col": np.asarray(array_col)}, index=index)


def plot_spatial_mean_composition_barplot(
    proportions, output_dir, *, sample_labels=None,
    name="spatial_mean_composition_barplot",
):
    """Stacked bar of *average* spot-level composition (one bar per sample)."""
    go = require_plotly()
    props = proportions.fillna(0.0)
    if sample_labels is not None:
        labels = pd.Series(list(sample_labels), index=props.index)
        means = props.groupby(labels).mean()
    else:
        means = props.mean(axis=0).to_frame("all_spots").T
    colours = color_sequence(means.shape[1])
    fig = go.Figure()
    for j, ct in enumerate(means.columns):
        fig.add_bar(x=[str(s) for s in means.index], y=means[ct].to_numpy(),
                    name=str(ct), marker_color=colours[j])
    fig.update_layout(barmode="stack",
                      **plotly_layout("Average spot-level composition",
                                      subtitle="average " + _SP_SUB, height=480))
    fig.update_yaxes(title_text="mean spot composition", range=[0, 1])
    fig.update_xaxes(title_text="sample / section", tickangle=45)
    return export_figure(fig, output_dir, name, data={"data": means},
                         caption="Average spot-level RNA-derived composition. " + _SP_SUB,
                         data_comment=["estimate_type: spot_rna_composition (mean over spots)"])


def plot_spatial_spot_pie_charts(
    proportions, array_row, array_col, output_dir, *, top_n=5, max_spots=200,
    min_prop=0.0, pie_size=0.045, name="spatial_spot_pie_charts",
):
    """One pie chart per spot, placed at its array coordinates.

    Collapses to *top_n* types + 'Other'; samples down to *max_spots* (with a
    warning) to keep the interactive figure responsive.
    """
    go = require_plotly()
    props = _collapse_spatial_top_n(proportions.fillna(0.0), top_n)
    ar = np.asarray(array_row, dtype=float)
    ac = np.asarray(array_col, dtype=float)
    n = len(props)
    warns: list[str] = []
    idx = np.arange(n)
    if n > max_spots:
        step = int(np.ceil(n / max_spots))
        idx = idx[::step]
        warns.append(f"{n} spots exceeds max_spots={max_spots}; showing a "
                     f"sampled subset of {len(idx)} spots for performance.")

    def _norm(v):
        lo, hi = float(np.min(v)), float(np.max(v))
        return (v - lo) / (hi - lo) if hi > lo else np.full_like(v, 0.5)

    xn, yn = _norm(ac), 1.0 - _norm(ar)  # invert row so map reads like tissue
    cell_types = [str(c) for c in props.columns]
    colours = color_sequence(len(cell_types))
    fig = go.Figure()
    half = pie_size / 2
    for s in idx:
        vals = props.iloc[int(s)].to_numpy().astype(float)
        vals = np.where(vals >= min_prop, vals, 0.0)
        if vals.sum() <= 0:
            continue
        fig.add_trace(go.Pie(
            labels=cell_types, values=vals, sort=False, textinfo="none",
            marker={"colors": colours}, showlegend=bool(s == idx[0]),
            domain={"x": [max(0, xn[s] - half), min(1, xn[s] + half)],
                    "y": [max(0, yn[s] - half), min(1, yn[s] + half)]}))
    fig.update_layout(**plotly_layout("Per-spot composition (pie charts)",
                                      subtitle=_SP_SUB, height=620))

    long = []
    for s in idx:
        for ct in cell_types:
            long.append({"spot": str(props.index[int(s)]),
                         "array_row": ar[s], "array_col": ac[s],
                         "cell_type": ct, "fraction": float(props.iloc[int(s)][ct])})
    res = export_figure(fig, output_dir, name,
                        data={"data": pd.DataFrame(long)},
                        caption="Per-spot RNA-derived composition. " + _SP_SUB,
                        data_comment=["estimate_type: spot_rna_composition"],
                        formats=())  # static export of many pies is not meaningful
    res.warnings.extend(warns)
    return res


def plot_spatial_abundance_maps(
    proportions, array_row, array_col, output_dir, *, top_n=6,
    name="spatial_abundance_maps",
):
    """Per-cell-type abundance scatter with a dropdown to switch cell type."""
    go = require_plotly()
    props = proportions.fillna(0.0)
    top = props.mean(axis=0).sort_values(ascending=False).index[:top_n].tolist()
    ar = np.asarray(array_row, dtype=float)
    ac = np.asarray(array_col, dtype=float)
    fig = go.Figure()
    for i, ct in enumerate(top):
        fig.add_trace(go.Scatter(
            x=ac, y=ar, mode="markers", name=str(ct), visible=(i == 0),
            marker={"color": props[ct].to_numpy(), "colorscale": "Viridis",
                    "cmin": 0, "cmax": 1, "size": 7, "symbol": "hexagon",
                    "colorbar": {"title": "composition"}}))
    buttons = [{"label": str(ct), "method": "update",
                "args": [{"visible": [j == i for j in range(len(top))]},
                         {"title": f"Abundance: {ct}<br><sup>{_SP_SUB}</sup>"}]}
               for i, ct in enumerate(top)]
    layout = plotly_layout(f"Abundance: {top[0]}", subtitle=_SP_SUB, height=560)
    layout["updatemenus"] = [{"buttons": buttons, "x": 1.02, "y": 1.0}]
    layout["yaxis"] = {"autorange": "reversed"}
    fig.update_layout(**layout)

    data = _coords_df(ar, ac, index=props.index)
    for ct in top:
        data[str(ct)] = props[ct].to_numpy()
    return export_figure(fig, output_dir, name, data={"data": data},
                         caption="Spatial abundance maps. " + _SP_SUB,
                         data_comment=["estimate_type: spot_rna_composition"])


def plot_spatial_dominant_cell_type_map(
    proportions, array_row, array_col, output_dir, *,
    low_conf_threshold=0.5, name="spatial_dominant_cell_type_map",
):
    """Map each spot by dominant cell type; low-confidence spots de-emphasised."""
    go = require_plotly()
    props = proportions.fillna(0.0)
    dom_idx = props.to_numpy().argmax(axis=1)
    dom_frac = props.to_numpy().max(axis=1)
    cell_types = list(props.columns)
    colours = color_sequence(len(cell_types))
    ar = np.asarray(array_row, dtype=float)
    ac = np.asarray(array_col, dtype=float)
    fig = go.Figure()
    for k, ct in enumerate(cell_types):
        mask = dom_idx == k
        if not mask.any():
            continue
        conf = dom_frac[mask] >= low_conf_threshold
        fig.add_trace(go.Scatter(
            x=ac[mask], y=ar[mask], mode="markers", name=str(ct),
            marker={"color": colours[k], "size": 8, "symbol": "hexagon",
                    "opacity": np.where(conf, 1.0, 0.35).tolist()},
            text=[f"{ct}: {f:.2f}" for f in dom_frac[mask]]))
    layout = plotly_layout("Dominant cell type per spot", subtitle=_SP_SUB, height=560)
    layout["yaxis"] = {"autorange": "reversed"}
    fig.update_layout(**layout)
    data = _coords_df(ar, ac, index=props.index)
    data["dominant_type"] = [cell_types[i] for i in dom_idx]
    data["dominant_fraction"] = dom_frac
    data["low_confidence"] = dom_frac < low_conf_threshold
    return export_figure(fig, output_dir, name, data={"data": data},
                         caption="Dominant cell type per spot (low-confidence "
                                 "spots faded). " + _SP_SUB,
                         data_comment=["estimate_type: spot_rna_composition",
                                       f"low_conf_threshold: {low_conf_threshold}"])


def plot_spatial_qc_maps(
    spot_qc, array_row, array_col, output_dir, *, metrics=None,
    name="spatial_qc_maps",
):
    """Per-spot QC metric maps with a dropdown to switch metric."""
    go = require_plotly()
    default = ["prop_entropy", "dominant_frac", "n_types_gt05", "nb_loglik",
               "spatial_resid"]
    cols = [m for m in (metrics or default) if m in spot_qc.columns]
    if not cols:
        cols = list(spot_qc.select_dtypes("number").columns)
    ar = np.asarray(array_row, dtype=float)
    ac = np.asarray(array_col, dtype=float)
    fig = go.Figure()
    for i, m in enumerate(cols):
        fig.add_trace(go.Scatter(
            x=ac, y=ar, mode="markers", name=m, visible=(i == 0),
            marker={"color": spot_qc[m].to_numpy(), "colorscale": "Magma",
                    "size": 7, "symbol": "hexagon",
                    "colorbar": {"title": m}}))
    buttons = [{"label": m, "method": "update",
                "args": [{"visible": [j == i for j in range(len(cols))]},
                         {"title": f"Spatial QC: {m}"}]}
               for i, m in enumerate(cols)]
    layout = plotly_layout(f"Spatial QC: {cols[0]}", subtitle=_SP_SUB, height=560)
    layout["updatemenus"] = [{"buttons": buttons, "x": 1.02, "y": 1.0}]
    layout["yaxis"] = {"autorange": "reversed"}
    fig.update_layout(**layout)
    data = _coords_df(ar, ac)
    for m in cols:
        data[m] = spot_qc[m].to_numpy()
    return export_figure(fig, output_dir, name, data={"data": data},
                         caption="Per-spot spatial QC metrics. " + _SP_SUB,
                         data_comment=["per-spot QC at recorded coordinates"])


def plot_spatial_morans_i_barplot(
    morans_i, output_dir, *, name="spatial_morans_i_barplot",
):
    """Moran's I per cell type, sorted descending (spatial structure)."""
    go = require_plotly()
    s = morans_i.sort_values(ascending=False)
    colours = ["#55A868" if v > 0 else "#C44E52" for v in s.to_numpy()]
    fig = go.Figure(go.Bar(x=[str(c) for c in s.index], y=s.to_numpy(),
                           marker_color=colours))
    fig.update_layout(**plotly_layout("Spatial autocorrelation (Moran's I)",
                                      subtitle="higher = more spatially structured",
                                      height=460))
    fig.update_xaxes(tickangle=90, title_text="cell type")
    fig.update_yaxes(title_text="Moran's I")
    return export_figure(fig, output_dir, name,
                         data={"data": s.to_frame("morans_i")},
                         caption="Moran's I per cell-type composition map.",
                         data_comment=["Moran's I per cell type"])
