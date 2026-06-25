"""
Spatial-benchmark publication figures (report section "Spatial benchmark").

Real Visium data have **no spot-level ground truth** in this benchmark, so these
figures evaluate **concordance and spatial structure**, never accuracy.  No
Pearson/RMSE-vs-truth is shown unless a synthetic spatial benchmark with known
simulated truth is explicitly provided.  Only **spatial** methods appear here.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.style import plotly_layout

__all__ = [
    "plot_spatial_method_status_summary",
    "plot_spatial_concordance_heatmap",
    "plot_spatial_structure_metrics_summary",
    "plot_spatial_runtime_comparison",
    "plot_spatial_output_completeness_summary",
]

_SCORED = ("executed", "success", "imported", "executed_imported")
_STATUS_COLOR = {"executed": "#3a9b62", "success": "#3a9b62",
                 "imported": "#4e79a7", "executed_imported": "#4e79a7",
                 "exported_only": "#d8a93a", "exported_not_run": "#d8a93a",
                 "skipped": "#b0b0b0", "failed": "#cc4b3e"}
_NO_GT = ("Real Visium has no spot-level ground truth here; this is a "
          "concordance/structure metric, NOT accuracy.")


def _ensure_method(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "method" in df.columns:
        return df
    df = df.reset_index()
    if "method" not in df.columns:
        df = df.rename(columns={df.columns[0]: "method"})
    return df


def _spatial_only(df: pd.DataFrame) -> pd.DataFrame:
    df = _ensure_method(df)
    if "modality" in df.columns:
        df = df[df["modality"].astype(str) == "spatial"]
    return df


def plot_spatial_method_status_summary(status: pd.DataFrame, output_dir, *,
                                       name: str = "spatial_method_status_summary") -> FigureResult:
    """Status counts for **spatial** benchmark methods only."""
    go = require_plotly()
    df = _spatial_only(status)
    if "status" not in df.columns or df.empty:
        raise ValueError("spatial status summary needs spatial rows with 'status'.")
    counts = df["status"].astype(str).value_counts()
    colors = [_STATUS_COLOR.get(s, "#999999") for s in counts.index]
    fig = go.Figure(go.Bar(x=[str(s) for s in counts.index], y=counts.values,
                           marker_color=colors, text=counts.values, textposition="auto"))
    fig.update_layout(**plotly_layout("Spatial benchmark — method status",
                                      subtitle="spatial methods only", height=400))
    fig.update_xaxes(title_text="status"); fig.update_yaxes(title_text="n methods")
    fig.update_layout(showlegend=False)
    return export_figure(fig, output_dir, name,
                         data={"data": counts.rename("n_methods").reset_index().rename(
                             columns={"index": "status"})},
                         caption="Spatial benchmark method status (spatial methods only). "
                                 + _NO_GT,
                         data_comment=["Spatial benchmark method status counts."])


def plot_spatial_concordance_heatmap(concordance: pd.DataFrame, output_dir, *,
                                     name: str = "spatial_concordance_heatmap") -> FigureResult:
    """Method×method concordance heatmap (agreement, NOT accuracy).

    *concordance* is a square method×method matrix (e.g. mean per-spot
    correlation).  Requires ≥2 methods.
    """
    go = require_plotly()
    M = concordance.copy()
    if M.shape[0] < 2 or M.shape[1] < 2:
        raise ValueError("concordance heatmap needs >=2 spatial methods.")
    fig = go.Figure(go.Heatmap(z=M.to_numpy(), x=[str(c) for c in M.columns],
                               y=[str(i) for i in M.index], colorscale="Viridis",
                               zmin=0, zmax=1, colorbar={"title": "agreement"}))
    fig.update_layout(**plotly_layout("Spatial method concordance",
                                      subtitle="pairwise agreement (NOT accuracy)",
                                      height=460))
    fig.update_xaxes(tickangle=45)
    return export_figure(fig, output_dir, name, data={"data": M},
                         caption="Spatial method concordance heatmap. Cells show "
                                 "pairwise agreement between spatial methods. " + _NO_GT,
                         data_comment=["Spatial method-method concordance."])


def plot_spatial_structure_metrics_summary(metrics: pd.DataFrame, output_dir, *,
                                           name: str = "spatial_structure_metrics_summary") -> FigureResult:
    """Spatial-structure metrics per method/population (Moran's I, entropy, …).

    *metrics*: long form with columns ``metric``, ``value`` and a label column
    (``population`` or ``method``).
    """
    go = require_plotly()
    df = metrics.copy()
    if not {"metric", "value"} <= set(df.columns):
        raise ValueError("structure summary needs 'metric' and 'value' columns.")
    label_col = next((c for c in ("population", "method", "label")
                      if c in df.columns), None)
    fig = go.Figure()
    if label_col is None:
        agg = df.groupby("metric")["value"].mean()
        fig.add_bar(x=[str(m) for m in agg.index], y=agg.values, marker_color="#59a14f")
    else:
        for met, grp in df.groupby("metric"):
            fig.add_bar(x=grp[label_col].astype(str), y=pd.to_numeric(grp["value"]),
                        name=str(met))
        fig.update_layout(barmode="group")
    fig.update_layout(**plotly_layout(
        "Spatial structure metrics",
        subtitle="Moran's I / entropy / dominant fraction (structure, not accuracy)",
        height=440))
    fig.update_xaxes(title_text=label_col or "metric", tickangle=45)
    fig.update_yaxes(title_text="value")
    return export_figure(fig, output_dir, name, data={"data": df},
                         caption="Spatial structure metrics (Moran's I, entropy, "
                                 "dominant/near-zero fraction). " + _NO_GT,
                         data_comment=["Spatial structure metrics summary."])


def plot_spatial_runtime_comparison(status: pd.DataFrame, output_dir, *,
                                    name: str = "spatial_runtime_comparison") -> FigureResult:
    """Runtime (seconds) for executed/imported spatial methods."""
    go = require_plotly()
    df = _spatial_only(status)
    if "runtime_seconds" not in df.columns:
        raise ValueError("spatial runtime needs 'runtime_seconds'.")
    if "status" in df.columns:
        df = df[df["status"].astype(str).isin(_SCORED)]
    df = df[pd.to_numeric(df["runtime_seconds"], errors="coerce").fillna(0) > 0]
    if df.empty:
        raise ValueError("no executed/imported spatial methods with runtime.")
    df = df.sort_values("runtime_seconds", ascending=True)
    fig = go.Figure(go.Bar(x=pd.to_numeric(df["runtime_seconds"]),
                           y=df["method"].astype(str), orientation="h",
                           marker_color="#9c755f"))
    fig.update_layout(**plotly_layout("Spatial benchmark runtime",
                                      subtitle="executed/imported spatial methods",
                                      height=max(320, 26 * len(df) + 120)))
    fig.update_xaxes(title_text="runtime (seconds)"); fig.update_yaxes(title_text="method")
    fig.update_layout(showlegend=False)
    return export_figure(fig, output_dir, name,
                         data={"data": df[["method", "runtime_seconds"]]},
                         caption="Runtime per executed/imported spatial method (seconds). "
                                 "A cost metric, not accuracy.",
                         data_comment=["Spatial runtime comparison."])


def plot_spatial_output_completeness_summary(status: pd.DataFrame, output_dir, *,
                                             name: str = "spatial_output_completeness_summary") -> FigureResult:
    """Per-method output completeness (executed / imported / exported / skipped / failed)."""
    go = require_plotly()
    df = _spatial_only(status)
    if "status" not in df.columns or df.empty:
        raise ValueError("completeness summary needs spatial rows with 'status'.")
    df = df.sort_values("status")
    colors = [_STATUS_COLOR.get(str(s), "#999999") for s in df["status"]]
    fig = go.Figure(go.Bar(x=df["method"].astype(str),
                           y=np.ones(len(df)), marker_color=colors,
                           text=df["status"].astype(str), textposition="auto"))
    fig.update_layout(**plotly_layout("Spatial benchmark — output completeness",
                                      subtitle="did each spatial method produce usable output?",
                                      height=400))
    fig.update_xaxes(title_text="method", tickangle=45)
    fig.update_yaxes(title_text="", showticklabels=False)
    fig.update_layout(showlegend=False)
    return export_figure(fig, output_dir, name,
                         data={"data": df[["method", "status"]]},
                         caption="Spatial method output completeness (status per method). "
                                 + _NO_GT,
                         data_comment=["Spatial output completeness."])
