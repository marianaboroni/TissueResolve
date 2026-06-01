"""
Bulk-benchmark publication figures (report section "Bulk benchmark").

Bulk pseudobulk mixtures have **known ground truth**, so accuracy metrics
(Pearson, RMSE, MAE) are valid here.  Only executed/imported **bulk** methods are
ranked; spatial methods never appear in these figures.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.style import plotly_layout

__all__ = [
    "plot_bulk_method_status_summary",
    "plot_bulk_fine_accuracy_leaderboard",
    "plot_bulk_family_accuracy_leaderboard",
    "plot_bulk_rmse_mae_comparison",
    "plot_bulk_runtime_comparison",
]

_SCORED = ("executed", "success", "imported", "executed_imported")
_STATUS_COLOR = {"executed": "#3a9b62", "success": "#3a9b62",
                 "imported": "#4e79a7", "executed_imported": "#4e79a7",
                 "exported_only": "#d8a93a", "exported_not_run": "#d8a93a",
                 "skipped": "#b0b0b0", "failed": "#cc4b3e"}


def _ensure_method(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "method" in df.columns:
        return df
    df = df.reset_index()
    if "method" not in df.columns:
        df = df.rename(columns={df.columns[0]: "method"})
    return df


def _bulk_only(df: pd.DataFrame) -> pd.DataFrame:
    df = _ensure_method(df)
    if "modality" in df.columns:
        df = df[df["modality"].astype(str) == "bulk"]
    return df


def plot_bulk_method_status_summary(status: pd.DataFrame, output_dir, *,
                                    name: str = "bulk_method_status_summary") -> FigureResult:
    """Status counts for **bulk** benchmark methods only."""
    go = require_plotly()
    df = _bulk_only(status)
    if "status" not in df.columns or df.empty:
        raise ValueError("bulk status summary needs bulk rows with a 'status' column.")
    counts = df["status"].astype(str).value_counts()
    colors = [_STATUS_COLOR.get(s, "#999999") for s in counts.index]
    fig = go.Figure(go.Bar(x=[str(s) for s in counts.index], y=counts.values,
                           marker_color=colors, text=counts.values, textposition="auto"))
    fig.update_layout(**plotly_layout("Bulk benchmark — method status",
                                      subtitle="bulk methods only", height=400))
    fig.update_xaxes(title_text="status"); fig.update_yaxes(title_text="n methods")
    fig.update_layout(showlegend=False)
    return export_figure(fig, output_dir, name,
                         data={"data": counts.rename("n_methods").reset_index().rename(
                             columns={"index": "status"})},
                         caption="Bulk benchmark method status (bulk methods only). "
                                 "Only executed/imported methods carry accuracy.",
                         data_comment=["Bulk benchmark method status counts."])


def _accuracy_leaderboard(status, output_dir, name, *, metric_col, title, level):
    go = require_plotly()
    df = _bulk_only(status)
    if metric_col not in df.columns:
        raise ValueError(f"leaderboard needs a '{metric_col}' column.")
    if "status" in df.columns:
        df = df[df["status"].astype(str).isin(_SCORED)]
    df = df[pd.to_numeric(df[metric_col], errors="coerce").notna()]
    if df.empty:
        raise ValueError(f"no executed/imported bulk methods with {metric_col}.")
    df = df.sort_values(metric_col, ascending=True)
    fig = go.Figure(go.Bar(x=pd.to_numeric(df[metric_col]), y=df["method"].astype(str),
                           orientation="h", marker_color="#4e79a7",
                           text=[f"{v:.3f}" for v in pd.to_numeric(df[metric_col])],
                           textposition="auto"))
    fig.update_layout(**plotly_layout(
        title, subtitle="executed/imported bulk methods; pseudobulk ground truth",
        height=max(340, 28 * len(df) + 130)))
    fig.update_xaxes(title_text=f"{level} accuracy ({metric_col})")
    fig.update_yaxes(title_text="method"); fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": df[["method", metric_col]]},
        caption=(f"Bulk {level}-level accuracy leaderboard. Bars show {metric_col} "
                 "between predicted and known pseudobulk proportions for "
                 "executed/imported bulk methods only. Higher is better; valid "
                 "because pseudobulk ground truth is available."),
        data_comment=[f"Bulk {level} accuracy leaderboard (executed/imported)."])


def plot_bulk_fine_accuracy_leaderboard(status: pd.DataFrame, output_dir, *,
                                        metric_col: str = "accuracy",
                                        name: str = "bulk_fine_accuracy_leaderboard") -> FigureResult:
    """Fine-level accuracy (Pearson vs pseudobulk truth) for bulk methods."""
    return _accuracy_leaderboard(status, output_dir, name, metric_col=metric_col,
                                 title="Bulk fine-level benchmark accuracy", level="fine")


def plot_bulk_family_accuracy_leaderboard(status: pd.DataFrame, output_dir, *,
                                          metric_col: str = "family_accuracy",
                                          name: str = "bulk_family_accuracy_leaderboard") -> FigureResult:
    """Family-level accuracy leaderboard (requires a family-accuracy column)."""
    return _accuracy_leaderboard(status, output_dir, name, metric_col=metric_col,
                                 title="Bulk family-level benchmark accuracy", level="family")


def plot_bulk_rmse_mae_comparison(metrics: pd.DataFrame, output_dir, *,
                                  name: str = "bulk_rmse_mae_comparison") -> FigureResult:
    """Grouped RMSE / MAE bars per bulk method (lower is better)."""
    go = require_plotly()
    df = _ensure_method(metrics)
    cols = [c for c in ("rmse", "mae", "RMSE", "MAE") if c in df.columns]
    if not cols:
        raise ValueError("rmse/mae comparison needs an 'rmse' and/or 'mae' column.")
    fig = go.Figure()
    palette = {"rmse": "#e15759", "mae": "#f28e2b", "RMSE": "#e15759", "MAE": "#f28e2b"}
    for c in cols:
        fig.add_bar(x=df["method"].astype(str), y=pd.to_numeric(df[c], errors="coerce"),
                    name=c.upper(), marker_color=palette.get(c, "#999999"))
    fig.update_layout(barmode="group",
                      **plotly_layout("Bulk benchmark error (RMSE / MAE)",
                                      subtitle="lower is better; pseudobulk ground truth",
                                      height=420))
    fig.update_xaxes(title_text="method", tickangle=45)
    fig.update_yaxes(title_text="error")
    return export_figure(fig, output_dir, name, data={"data": df[["method"] + cols]},
                         caption="Bulk benchmark RMSE/MAE vs pseudobulk ground truth "
                                 "(lower is better). Valid because ground truth exists.",
                         data_comment=["Bulk RMSE/MAE comparison."])


def plot_bulk_runtime_comparison(status: pd.DataFrame, output_dir, *,
                                 name: str = "bulk_runtime_comparison") -> FigureResult:
    """Runtime (seconds) for executed/imported bulk methods."""
    go = require_plotly()
    df = _bulk_only(status)
    if "runtime_seconds" not in df.columns:
        raise ValueError("bulk runtime needs 'runtime_seconds'.")
    if "status" in df.columns:
        df = df[df["status"].astype(str).isin(_SCORED)]
    df = df[pd.to_numeric(df["runtime_seconds"], errors="coerce").fillna(0) > 0]
    if df.empty:
        raise ValueError("no executed/imported bulk methods with runtime.")
    df = df.sort_values("runtime_seconds", ascending=True)
    fig = go.Figure(go.Bar(x=pd.to_numeric(df["runtime_seconds"]),
                           y=df["method"].astype(str), orientation="h",
                           marker_color="#9c755f"))
    fig.update_layout(**plotly_layout("Bulk benchmark runtime",
                                      subtitle="executed/imported bulk methods",
                                      height=max(320, 26 * len(df) + 120)))
    fig.update_xaxes(title_text="runtime (seconds)"); fig.update_yaxes(title_text="method")
    fig.update_layout(showlegend=False)
    return export_figure(fig, output_dir, name,
                         data={"data": df[["method", "runtime_seconds"]]},
                         caption="Runtime per executed/imported bulk method (seconds). "
                                 "A cost (runtime) metric, not accuracy.",
                         data_comment=["Bulk runtime comparison."])
