"""
Publication-ready bulk RNA-seq deconvolution plots.

Every function:

* takes a result/QC container discovered from :mod:`tissueresolve.results`,
* saves the figure (PNG/PDF/SVG) **and** the underlying data (TSV),
* returns a :class:`~tissueresolve.plotting.style.PlotResult`,
* labels outputs as **mRNA proportions**, never absolute cell fractions.
"""
from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

from tissueresolve.plotting import captions
from tissueresolve.plotting.style import (
    PlotResult,
    get_palette,
    new_figure,
    save_outputs,
)
from tissueresolve.results import BulkDeconvResult, QCReport

__all__ = [
    "composition_barplot",
    "proportion_heatmap",
    "qc_barplot",
    "bootstrap_ci_plot",
    "marker_recall_plot",
]

_PROP_COMMENT = [
    "estimate_type: mRNA_proportion",
    "WARNING: mRNA proportions, NOT absolute cell fractions.",
]


def composition_barplot(
    result: BulkDeconvResult,
    output_dir: Union[str, Path],
    *,
    name: str = "bulk_composition",
) -> PlotResult:
    """Stacked bar of mRNA proportions per sample."""
    props = result.proportions
    colours = get_palette(props.shape[1])

    fig, ax = new_figure(figsize=(max(6.0, 0.5 * len(props) + 2), 4.5))
    bottom = np.zeros(len(props))
    x = np.arange(len(props))
    for j, ct in enumerate(props.columns):
        ax.bar(x, props[ct].to_numpy(), bottom=bottom, label=str(ct), color=colours[j])
        bottom += props[ct].to_numpy()
    ax.set_xticks(x)
    ax.set_xticklabels(props.index.astype(str), rotation=90)
    ax.set_ylabel("mRNA proportion")
    ax.set_ylim(0, 1)
    ax.set_title("Bulk cell-type composition (mRNA proportions)")
    ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", title="Cell type")
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(
        fig, props, output_dir, name, data_comment=_PROP_COMMENT
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.bulk_composition_caption(
            result.n_samples, result.n_cell_types,
            gene_panel_size=len(result.gene_panel),
        ),
    )


def proportion_heatmap(
    result: BulkDeconvResult,
    output_dir: Union[str, Path],
    *,
    name: str = "bulk_proportion_heatmap",
) -> PlotResult:
    """Heatmap of mRNA proportions (samples × cell types)."""
    props = result.proportions
    fig, ax = new_figure(figsize=(max(5.0, 0.4 * props.shape[1] + 3),
                                  max(4.0, 0.3 * props.shape[0] + 2)))
    im = ax.imshow(props.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(props.shape[1]))
    ax.set_xticklabels(props.columns.astype(str), rotation=90)
    ax.set_yticks(np.arange(props.shape[0]))
    ax.set_yticklabels(props.index.astype(str))
    ax.set_title("Bulk mRNA proportions")
    ax.grid(visible=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("mRNA proportion")

    fig_paths, data_paths = save_outputs(
        fig, props, output_dir, name, data_comment=_PROP_COMMENT
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.bulk_composition_caption(
            result.n_samples, result.n_cell_types,
            gene_panel_size=len(result.gene_panel),
        ),
    )


def qc_barplot(
    qc: QCReport,
    output_dir: Union[str, Path],
    *,
    deconv: BulkDeconvResult | None = None,
    name: str = "bulk_qc",
    r2_warn: float | None = None,
) -> PlotResult:
    """Grouped bar of per-sample reconstruction R² and profile correlation.

    Uses ``qc.recon_r2`` / ``qc.profile_corr`` when present, falling back to
    ``deconv.coverage_r2`` for R².
    """
    r2 = qc.recon_r2
    if r2 is None and deconv is not None:
        r2 = deconv.coverage_r2
    pc = qc.profile_corr

    cols: dict[str, pd.Series] = {}
    if r2 is not None:
        cols["recon_r2"] = r2
    if pc is not None:
        cols["profile_corr"] = pc
    if not cols:
        raise ValueError(
            "qc_barplot needs at least one of recon_r2/profile_corr (or a "
            "deconv with coverage_r2)."
        )
    table = pd.DataFrame(cols)

    fig, ax = new_figure(figsize=(max(6.0, 0.45 * len(table) + 2), 4.0))
    x = np.arange(len(table))
    width = 0.8 / max(1, table.shape[1])
    colours = get_palette(table.shape[1])
    for j, col in enumerate(table.columns):
        ax.bar(x + j * width, table[col].to_numpy(), width, label=col, color=colours[j])
    if r2_warn is not None:
        ax.axhline(r2_warn, ls="--", lw=1, color="firebrick",
                   label=f"R² warn = {r2_warn}")
    ax.set_xticks(x + width * (table.shape[1] - 1) / 2)
    ax.set_xticklabels(table.index.astype(str), rotation=90)
    ax.set_ylabel("metric value")
    ax.set_ylim(0, 1.02)
    ax.set_title("Bulk per-sample QC (heuristic thresholds)")
    ax.legend()
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(
        fig, table, output_dir, name,
        data_comment=["Per-sample bulk QC. Thresholds are heuristic."],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.bulk_qc_caption(r2_warn=r2_warn),
    )


def bootstrap_ci_plot(
    result: BulkDeconvResult,
    output_dir: Union[str, Path],
    *,
    sample: str | None = None,
    ci_level: float = 0.95,
    name: str = "bulk_bootstrap_ci",
) -> PlotResult:
    """Point estimate with bootstrap CI per cell type for one sample.

    Raises if the result carries no bootstrap CIs (never invents intervals).
    """
    if result.lower_ci is None or result.upper_ci is None:
        raise ValueError(
            "bootstrap_ci_plot requires bootstrap CIs; run the pipeline with "
            "n_bootstrap > 0 first."
        )
    if sample is None:
        sample = str(result.proportions.index[0])

    est = result.proportions.loc[sample]
    lo = result.lower_ci.loc[sample]
    hi = result.upper_ci.loc[sample]
    table = pd.DataFrame({"estimate": est, "ci_lower": lo, "ci_upper": hi})

    y = np.arange(len(table))
    err = np.vstack([
        (est - lo).to_numpy().clip(min=0),
        (hi - est).to_numpy().clip(min=0),
    ])
    fig, ax = new_figure(figsize=(6.0, max(3.0, 0.4 * len(table) + 1)))
    ax.errorbar(est.to_numpy(), y, xerr=err, fmt="o", color="#4C72B0",
                ecolor="#8C8C8C", capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels(table.index.astype(str))
    ax.set_xlabel("mRNA proportion")
    ax.set_xlim(0, 1)
    ax.set_title(f"Bootstrap {int(ci_level * 100)}% CI — sample {sample}")

    n_boot = result.run_metadata.get("n_bootstrap")
    fig_paths, data_paths = save_outputs(
        fig, table, output_dir, name,
        data_comment=_PROP_COMMENT + [f"sample: {sample}", f"ci_level: {ci_level}"],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.bulk_ci_caption(sample, ci_level=ci_level, n_bootstrap=n_boot),
    )


def marker_recall_plot(
    qc: QCReport,
    output_dir: Union[str, Path],
    *,
    name: str = "bulk_marker_recall",
) -> PlotResult:
    """Per-cell-type marker recall and spillover risk bars.

    Raises if neither metric is present (these come from the bulk QC step).
    """
    if qc.marker_recall is None and qc.spillover_risk is None:
        raise ValueError("marker_recall_plot needs marker_recall/spillover_risk.")

    cols: dict[str, pd.Series] = {}
    if qc.marker_recall is not None:
        cols["marker_recall"] = qc.marker_recall
    if qc.spillover_risk is not None:
        cols["spillover_risk"] = qc.spillover_risk
    table = pd.DataFrame(cols)

    fig, ax = new_figure(figsize=(max(5.0, 0.5 * len(table) + 2), 4.0))
    x = np.arange(len(table))
    width = 0.8 / max(1, table.shape[1])
    colours = get_palette(table.shape[1])
    for j, col in enumerate(table.columns):
        ax.bar(x + j * width, table[col].to_numpy(), width, label=col, color=colours[j])
    ax.set_xticks(x + width * (table.shape[1] - 1) / 2)
    ax.set_xticklabels(table.index.astype(str), rotation=90)
    ax.set_ylabel("metric value")
    ax.set_ylim(0, 1.02)
    ax.set_title("Per-cell-type marker recall & spillover risk")
    ax.legend()
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(fig, table, output_dir, name)
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)


# ===========================================================================
# Publication layer (Plotly) — interactive HTML + static + source data
# ===========================================================================

from tissueresolve.plotting import captions as _captions  # noqa: E402
from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly  # noqa: E402
from tissueresolve.plotting.style import (  # noqa: E402
    ESTIMATE_SUBTITLE,
    color_sequence,
    plotly_layout,
)

_BULK_SUB = ESTIMATE_SUBTITLE["bulk"]


def _cluster_order(matrix: np.ndarray) -> list[int]:
    """Hierarchical-clustering leaf order for rows of *matrix* (offline, scipy)."""
    n = matrix.shape[0]
    if n < 3:
        return list(range(n))
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage

        Z = linkage(matrix, method="average", metric="euclidean")
        return list(leaves_list(Z))
    except Exception:
        return list(range(n))


def _collapse_to_top_n(props: pd.DataFrame, top_n: Optional[int]) -> pd.DataFrame:
    """Keep the *top_n* cell types by mean proportion; sum the rest into 'Other'."""
    if not top_n or top_n >= props.shape[1]:
        return props
    order = props.mean(axis=0).sort_values(ascending=False).index.tolist()
    keep = order[:top_n]
    other = [c for c in props.columns if c not in keep]
    out = props[keep].copy()
    if other:
        out["Other"] = props[other].sum(axis=1)
    return out


def plot_bulk_composition_clustered_barplot(
    proportions: pd.DataFrame,
    output_dir,
    *,
    name: str = "bulk_composition_clustered_barplot",
    top_n: Optional[int] = None,
    cluster_samples: bool = True,
    cluster_cell_types: bool = True,
    sample_metadata: Optional[pd.DataFrame] = None,
) -> FigureResult:
    """Clustered stacked bar of bulk mRNA-derived composition (samples × types).

    Samples are ordered by hierarchical clustering of their composition; cell
    types optionally clustered too; rare types optionally collapsed into
    'Other'.  Saves the plotted matrix plus sample and cell-type orders.
    """
    go = require_plotly()
    props = _collapse_to_top_n(proportions.fillna(0.0), top_n)

    sample_order = (_cluster_order(props.to_numpy()) if cluster_samples
                    else list(range(props.shape[0])))
    samples = [props.index[i] for i in sample_order]
    if cluster_cell_types and props.shape[1] >= 3:
        ct_order = _cluster_order(props.to_numpy().T)
        cell_types = [props.columns[i] for i in ct_order]
    else:
        cell_types = list(props.columns)

    ordered = props.loc[samples, cell_types]
    colours = color_sequence(len(cell_types))

    fig = go.Figure()
    for j, ct in enumerate(cell_types):
        fig.add_bar(x=[str(s) for s in samples], y=ordered[ct].to_numpy(),
                    name=str(ct), marker_color=colours[j])
    fig.update_layout(barmode="stack",
                      **plotly_layout("Bulk cell-type composition (clustered)",
                                      subtitle=_BULK_SUB,
                                      height=520))
    fig.update_yaxes(title_text="mRNA-derived proportion", range=[0, 1])
    fig.update_xaxes(title_text="sample (clustered order)", tickangle=90)

    data = {
        "data": ordered,
        "sample_order": pd.DataFrame({"order": range(len(samples)),
                                      "sample": [str(s) for s in samples]}),
        "cell_type_order": pd.DataFrame({"order": range(len(cell_types)),
                                         "cell_type": [str(c) for c in cell_types]}),
    }
    caption = _captions.bulk_composition_caption(props.shape[0], proportions.shape[1])
    return export_figure(fig, output_dir, name, data=data, caption=caption,
                         data_comment=["estimate_type: mRNA_proportion (NOT cell fractions)"])


def plot_bulk_composition_heatmap(
    proportions: pd.DataFrame,
    output_dir,
    *,
    name: str = "bulk_composition_heatmap",
    cluster_samples: bool = True,
    cluster_cell_types: bool = True,
) -> FigureResult:
    """Clustered heatmap of bulk composition — useful for many cell types."""
    go = require_plotly()
    props = proportions.fillna(0.0)
    s_order = (_cluster_order(props.to_numpy()) if cluster_samples
               else list(range(props.shape[0])))
    c_order = (_cluster_order(props.to_numpy().T)
               if cluster_cell_types and props.shape[1] >= 3
               else list(range(props.shape[1])))
    samples = [props.index[i] for i in s_order]
    cell_types = [props.columns[i] for i in c_order]
    ordered = props.loc[samples, cell_types]

    fig = go.Figure(go.Heatmap(
        z=ordered.to_numpy(), x=[str(c) for c in cell_types],
        y=[str(s) for s in samples], colorscale="Viridis", zmin=0, zmax=1,
        colorbar={"title": "proportion"}))
    fig.update_layout(**plotly_layout("Bulk composition heatmap (clustered)",
                                      subtitle=_BULK_SUB, height=560))
    fig.update_xaxes(tickangle=90)
    return export_figure(fig, output_dir, name, data={"data": ordered},
                         caption=_captions.bulk_composition_caption(
                             props.shape[0], props.shape[1]),
                         data_comment=["estimate_type: mRNA_proportion"])


def plot_bulk_qc_summary(
    qc_table: pd.DataFrame,
    output_dir,
    *,
    name: str = "bulk_qc_summary",
    r2_warn: Optional[float] = None,
) -> FigureResult:
    """Grouped bar of per-sample QC metrics (R², profile correlation, …)."""
    go = require_plotly()
    table = qc_table.select_dtypes("number")
    fig = go.Figure()
    colours = color_sequence(table.shape[1])
    for j, col in enumerate(table.columns):
        fig.add_bar(x=[str(s) for s in table.index], y=table[col].to_numpy(),
                    name=str(col), marker_color=colours[j])
    if r2_warn is not None:
        fig.add_hline(y=r2_warn, line_dash="dash", line_color="firebrick",
                      annotation_text=f"R² warn = {r2_warn}")
    fig.update_layout(barmode="group",
                      **plotly_layout("Bulk per-sample QC (heuristic thresholds)",
                                      subtitle=_BULK_SUB, height=460))
    fig.update_xaxes(title_text="sample", tickangle=90)
    fig.update_yaxes(title_text="metric value")
    return export_figure(fig, output_dir, name, data={"data": table},
                         caption=_captions.bulk_qc_caption(r2_warn=r2_warn),
                         data_comment=["Per-sample bulk QC; thresholds heuristic."])


def plot_bulk_uncertainty(
    proportions: pd.DataFrame,
    lower_ci: Optional[pd.DataFrame],
    upper_ci: Optional[pd.DataFrame],
    output_dir,
    *,
    name: str = "bulk_uncertainty_plot",
) -> FigureResult:
    """Mean estimate and bootstrap CI width per cell type (uncertainty surfaced)."""
    go = require_plotly()
    mean_est = proportions.mean(axis=0)
    fig = go.Figure()
    warnings: list[str] = []
    if lower_ci is not None and upper_ci is not None:
        width = (upper_ci - lower_ci).mean(axis=0).reindex(mean_est.index)
        table = pd.DataFrame({"mean_estimate": mean_est, "mean_ci_width": width})
        fig.add_bar(x=table.index.astype(str), y=table["mean_estimate"].to_numpy(),
                    name="mean estimate",
                    error_y={"type": "data", "array": (width / 2).to_numpy(),
                             "visible": True},
                    marker_color="#4C72B0")
        sub = "mean mRNA-derived proportion ± half mean bootstrap CI width"
    else:
        table = pd.DataFrame({"mean_estimate": mean_est})
        fig.add_bar(x=table.index.astype(str), y=table["mean_estimate"].to_numpy(),
                    marker_color="#8C8C8C")
        sub = "no bootstrap CIs available — uncertainty not quantified"
        warnings.append("No bootstrap CIs available; uncertainty not shown.")
    fig.update_layout(**plotly_layout("Bulk uncertainty", subtitle=sub, height=460))
    fig.update_xaxes(title_text="cell type", tickangle=90)
    fig.update_yaxes(title_text="mRNA-derived proportion", range=[0, 1])
    res = export_figure(fig, output_dir, name, data={"data": table},
                        caption="Bulk estimate uncertainty. " + _BULK_SUB,
                        data_comment=["estimate_type: mRNA_proportion"])
    res.warnings.extend(warnings)
    return res
