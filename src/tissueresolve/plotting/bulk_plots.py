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
