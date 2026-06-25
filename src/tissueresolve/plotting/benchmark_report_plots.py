"""
Benchmark publication figures (report section 9).

Simple, honest benchmark visuals built from the benchmark harness TSVs: method
status, measured bulk accuracy (executed/imported only), runtime, and the
weighted composite *scorecard* (explicitly labelled as a scorecard, not
objective accuracy).  Spatial methods are never ranked in the bulk accuracy
leaderboard (real Visium has no ground truth).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.style import color_sequence, plotly_layout

__all__ = [
    "plot_benchmark_method_status",
    "plot_bulk_accuracy_leaderboard",
    "plot_runtime_comparison",
    "plot_composite_scorecard",
]

_STATUS_COLOR = {
    "executed": "#3a9b62", "success": "#3a9b62", "imported": "#4e79a7",
    "executed_imported": "#4e79a7", "exported_only": "#d8a93a",
    "exported_not_run": "#d8a93a", "skipped": "#b0b0b0", "failed": "#cc4b3e",
}
_SCORED = ("executed", "success", "imported", "executed_imported")


def _ensure_method_column(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with a 'method' column (from a named index if needed)."""
    df = df.copy()
    if "method" in df.columns:
        return df
    df = df.reset_index()
    if "method" not in df.columns:
        df = df.rename(columns={df.columns[0]: "method"})
    return df


def plot_benchmark_method_status(
    status: pd.DataFrame,
    output_dir,
    *,
    name: str = "benchmark_method_status_summary",
) -> FigureResult:
    """Bar of method counts per status (executed/imported/skipped/…)."""
    go = require_plotly()
    df = status.copy()
    if "status" not in df.columns:
        raise ValueError("benchmark status needs a 'status' column.")
    counts = df["status"].astype(str).value_counts()
    colors = [_STATUS_COLOR.get(s, "#999999") for s in counts.index]
    fig = go.Figure(go.Bar(x=[str(s) for s in counts.index], y=counts.values,
                           marker_color=colors, text=counts.values,
                           textposition="auto"))
    fig.update_layout(**plotly_layout("Benchmark method status",
                                      subtitle="how many tools ran, were imported, or skipped",
                                      height=420))
    fig.update_xaxes(title_text="status")
    fig.update_yaxes(title_text="number of methods")
    fig.update_layout(showlegend=False)
    table = counts.rename("n_methods").reset_index().rename(columns={"index": "status"})
    return export_figure(
        fig, output_dir, name, data={"data": table},
        caption=("Number of benchmark methods per status; only executed/imported "
                 "tools carry measured metrics, the rest are listed but not "
                 "scored."),
        data_comment=["Benchmark method status counts."])


def plot_bulk_accuracy_leaderboard(
    status: pd.DataFrame,
    output_dir,
    *,
    name: str = "bulk_accuracy_leaderboard",
) -> FigureResult:
    """Measured accuracy for executed/imported **bulk** methods only."""
    go = require_plotly()
    df = _ensure_method_column(status)
    if "accuracy" not in df.columns:
        raise ValueError("leaderboard needs an 'accuracy' column.")
    keep = df.copy()
    if "modality" in keep.columns:
        keep = keep[keep["modality"].astype(str) == "bulk"]
    if "status" in keep.columns:
        keep = keep[keep["status"].astype(str).isin(_SCORED)]
    keep = keep[pd.to_numeric(keep["accuracy"], errors="coerce").notna()]
    if keep.empty:
        raise ValueError("no executed/imported bulk methods with accuracy.")
    keep = keep.sort_values("accuracy", ascending=True)
    fig = go.Figure(go.Bar(
        x=pd.to_numeric(keep["accuracy"]), y=keep["method"].astype(str),
        orientation="h", marker_color="#4e79a7",
        text=[f"{a:.3f}" for a in pd.to_numeric(keep["accuracy"])],
        textposition="auto"))
    fig.update_layout(**plotly_layout(
        "Bulk benchmark accuracy",
        subtitle="executed/imported bulk methods only (pseudobulk ground truth)",
        height=max(360, 28 * len(keep) + 130)))
    fig.update_xaxes(title_text="accuracy (e.g. Pearson vs true proportions)")
    fig.update_yaxes(title_text="method")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": keep[["method", "accuracy"]]},
        caption=("Measured bulk accuracy for executed/imported methods against "
                 "pseudobulk ground truth; skipped/exported/failed tools are "
                 "excluded."),
        data_comment=["Bulk accuracy leaderboard (executed/imported only)."])


def plot_runtime_comparison(
    status: pd.DataFrame,
    output_dir,
    *,
    name: str = "runtime_comparison",
) -> FigureResult:
    """Runtime (seconds) for executed/imported methods."""
    go = require_plotly()
    df = _ensure_method_column(status)
    if "runtime_seconds" not in df.columns:
        raise ValueError("runtime comparison needs 'runtime_seconds'.")
    keep = df.copy()
    if "status" in keep.columns:
        keep = keep[keep["status"].astype(str).isin(_SCORED)]
    keep = keep[pd.to_numeric(keep["runtime_seconds"], errors="coerce").fillna(0) > 0]
    if keep.empty:
        raise ValueError("no executed/imported methods with runtime > 0.")
    keep = keep.sort_values("runtime_seconds", ascending=True)
    fig = go.Figure(go.Bar(
        x=pd.to_numeric(keep["runtime_seconds"]), y=keep["method"].astype(str),
        orientation="h", marker_color="#9c755f"))
    fig.update_layout(**plotly_layout("Benchmark runtime",
                                      subtitle="executed/imported methods",
                                      height=max(360, 26 * len(keep) + 130)))
    fig.update_xaxes(title_text="runtime (seconds)")
    fig.update_yaxes(title_text="method")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": keep[["method", "runtime_seconds"]]},
        caption="Runtime per executed/imported method (seconds).",
        data_comment=["Runtime comparison (executed/imported)."])


def plot_composite_scorecard(
    composite: pd.DataFrame,
    output_dir,
    *,
    name: str = "composite_scorecard",
) -> FigureResult:
    """Stacked score-dimension bars per method (a weighted *scorecard*).

    Explicitly a scorecard, not objective accuracy; only scored methods shown.
    """
    go = require_plotly()
    df = _ensure_method_column(composite)
    if "final_score" in df.columns:
        df = df[pd.to_numeric(df["final_score"], errors="coerce").notna()]
    if df.empty:
        raise ValueError("no scored methods in composite table.")
    dims = [c for c in df.columns if c.endswith("_score") and c != "final_score"]
    if not dims:
        raise ValueError("composite table has no *_score dimension columns.")
    df = df.sort_values("final_score" if "final_score" in df.columns else dims[0],
                        ascending=True)
    colors = color_sequence(len(dims))
    fig = go.Figure()
    for j, d in enumerate(dims):
        fig.add_bar(y=df["method"].astype(str), x=pd.to_numeric(df[d], errors="coerce"),
                    name=d.replace("_score", "").replace("_", " "),
                    orientation="h", marker_color=colors[j])
    fig.update_layout(barmode="stack",
                      **plotly_layout(
                          "Composite benchmark scorecard",
                          subtitle="weighted QC scorecard — NOT objective accuracy",
                          height=max(380, 30 * len(df) + 150)))
    fig.update_xaxes(title_text="weighted score contribution")
    fig.update_yaxes(title_text="method")
    fig.update_layout(legend={"title": {"text": "dimension"}})
    keep = ["method"] + dims + (["final_score"] if "final_score" in df.columns else [])
    return export_figure(
        fig, output_dir, name, data={"data": df[keep]},
        caption=("Weighted composite *scorecard* by dimension (accuracy, "
                 "robustness, usability, interpretability, resolution-awareness, "
                 "runtime). This is a scorecard, not an objective accuracy "
                 "measure; bulk and spatial are scored separately."),
        data_comment=["Composite scorecard (scored methods only)."])
