"""
Lightweight benchmark comparison plots (not full dashboards).

* estimated-vs-true proportion scatter
* per-cell-type RMSE bar plot
* method-comparison bar plot across a list of
  :class:`~tissueresolve.results.BenchmarkResult`

These operate on synthetic-benchmark outputs only and are explicitly labelled
as such; they are not real-data validation.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence, Union

import numpy as np
import pandas as pd

from tissueresolve.plotting import captions
from tissueresolve.plotting.style import (
    PlotResult,
    get_palette,
    new_figure,
    save_outputs,
)
from tissueresolve.results import BenchmarkResult

__all__ = [
    "estimated_vs_true_scatter",
    "per_celltype_rmse_barplot",
    "method_comparison_plot",
]


def estimated_vs_true_scatter(
    Pi_pred: np.ndarray,
    Pi_true: np.ndarray,
    cell_types: Sequence[str],
    output_dir: Union[str, Path],
    *,
    method: str = "SpatCAR",
    scenario: str = "basic",
    name: str = "benchmark_estimated_vs_true",
) -> PlotResult:
    """Scatter of estimated vs ground-truth proportions, coloured by cell type."""
    Pi_pred = np.asarray(Pi_pred)
    Pi_true = np.asarray(Pi_true)
    if Pi_pred.shape != Pi_true.shape:
        raise ValueError(
            f"shape mismatch: Pi_pred {Pi_pred.shape} vs Pi_true {Pi_true.shape}."
        )
    cts = list(cell_types)
    colours = get_palette(len(cts))

    fig, ax = new_figure(figsize=(5.0, 5.0))
    for k, ct in enumerate(cts):
        ax.scatter(Pi_true[:, k], Pi_pred[:, k], s=10, alpha=0.6,
                   color=colours[k], label=str(ct), linewidths=0)
    ax.plot([0, 1], [0, 1], "--", color="#8C8C8C", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("true proportion")
    ax.set_ylabel("estimated proportion")
    ax.set_title(f"Estimated vs true — {method} / {scenario}")
    ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", fontsize=7)

    # Source data: long-form (spot, cell_type, true, pred).
    rows = []
    for k, ct in enumerate(cts):
        for s in range(Pi_pred.shape[0]):
            rows.append({"spot": s, "cell_type": ct,
                         "true": float(Pi_true[s, k]), "pred": float(Pi_pred[s, k])})
    data = pd.DataFrame(rows)
    fig_paths, data_paths = save_outputs(
        fig, data, output_dir, name,
        data_comment=["Synthetic benchmark — not real-data validation.",
                      f"method: {method}", f"scenario: {scenario}"],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.benchmark_estimated_vs_true_caption(method, scenario),
    )


def per_celltype_rmse_barplot(
    per_celltype_metrics: pd.DataFrame,
    output_dir: Union[str, Path],
    *,
    column: str = "rmse",
    name: str = "benchmark_per_celltype_rmse",
) -> PlotResult:
    """Bar plot of a per-cell-type benchmark metric (RMSE by default)."""
    if column not in per_celltype_metrics.columns:
        raise ValueError(
            f"column {column!r} not in metrics {list(per_celltype_metrics.columns)}."
        )
    s = per_celltype_metrics[column]
    fig, ax = new_figure(figsize=(max(5.0, 0.5 * len(s) + 2), 4.0))
    x = np.arange(len(s))
    ax.bar(x, s.to_numpy(), color="#C44E52")
    ax.set_xticks(x)
    ax.set_xticklabels(s.index.astype(str), rotation=90)
    ax.set_ylabel(column)
    ax.set_title(f"Per-cell-type {column}")
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(
        fig, per_celltype_metrics, output_dir, name,
        data_comment=["Synthetic benchmark — not real-data validation."],
    )
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)


def method_comparison_plot(
    results: Sequence[BenchmarkResult],
    output_dir: Union[str, Path],
    *,
    metric: str = "rmse",
    name: str = "benchmark_method_comparison",
) -> PlotResult:
    """Grouped bar comparing one metric across methods, faceted by scenario."""
    if not results:
        raise ValueError("method_comparison_plot needs at least one BenchmarkResult.")
    rows = []
    for r in results:
        if metric not in r.metrics:
            raise ValueError(
                f"metric {metric!r} not in BenchmarkResult.metrics for "
                f"{r.method}/{r.scenario}: {list(r.metrics)}."
            )
        rows.append({"scenario": r.scenario, "method": r.method,
                     metric: float(r.metrics[metric])})
    long = pd.DataFrame(rows)
    table = long.pivot(index="scenario", columns="method", values=metric)

    methods = list(table.columns)
    scenarios = list(table.index)
    colours = get_palette(len(methods))
    fig, ax = new_figure(figsize=(max(6.0, 1.2 * len(scenarios) + 2), 4.0))
    x = np.arange(len(scenarios))
    width = 0.8 / max(1, len(methods))
    for j, m in enumerate(methods):
        ax.bar(x + j * width, table[m].to_numpy(), width, label=str(m), color=colours[j])
    ax.set_xticks(x + width * (len(methods) - 1) / 2)
    ax.set_xticklabels(scenarios, rotation=20)
    ax.set_ylabel(metric)
    ax.set_title(f"Method comparison — {metric}")
    ax.legend(title="method")
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(
        fig, table, output_dir, name,
        data_comment=["Synthetic benchmark — not real-data validation."],
    )
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)
