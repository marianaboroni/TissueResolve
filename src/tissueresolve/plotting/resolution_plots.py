"""
Resolution / spillover / unresolved-mass publication figures (report section 8).

These summarise *what resolution the user can trust*: the global separability
distribution, unresolved family mass, and a per-family recommended-interpretation
table.  All save source data and use the hierarchical palette.
"""
from __future__ import annotations

from typing import Mapping, Optional

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.palette import (
    OTHER_COLOR,
    build_hierarchical_color_map,
    color_dicts_from_map,
)
from tissueresolve.plotting.style import plotly_layout

__all__ = [
    "plot_separability_distribution",
    "plot_unresolved_mass_by_family",
    "plot_trusted_resolution_summary",
]


def plot_separability_distribution(
    separability: pd.DataFrame,
    output_dir,
    *,
    name: str = "separability_distribution",
    high_risk_threshold: float = 0.25,
) -> FigureResult:
    """Histogram of pairwise separability with a high-risk threshold line."""
    go = require_plotly()
    df = separability.copy()
    if "separability_score" not in df.columns and "bhattacharyya" in df.columns:
        df["separability_score"] = 1.0 - df["bhattacharyya"]
    if "separability_score" not in df.columns:
        raise ValueError("needs 'separability_score' or 'bhattacharyya'.")
    vals = df["separability_score"].astype(float)
    n_risk = int((vals < high_risk_threshold).sum())

    fig = go.Figure(go.Histogram(x=vals, nbinsx=30, marker_color="#4e79a7"))
    fig.add_vline(x=high_risk_threshold, line_dash="dash", line_color="firebrick",
                  annotation_text=f"high-risk < {high_risk_threshold}")
    fig.update_layout(**plotly_layout(
        "Distribution of pairwise cell-type separability",
        subtitle=f"{n_risk} of {len(vals)} pairs below the high-risk threshold",
        height=440))
    fig.update_xaxes(title_text="separability score (1 − Bhattacharyya)", range=[0, 1])
    fig.update_yaxes(title_text="number of cell-type pairs")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name,
        data={"data": df[[c for c in ("type_a", "type_b", "separability_score")
                          if c in df.columns]]},
        caption=(f"Histogram of pairwise separability; {n_risk} of {len(vals)} "
                 f"pairs fall below the high-risk threshold ({high_risk_threshold}) "
                 "and are unreliable at the fine level."),
        data_comment=[f"high_risk_threshold = {high_risk_threshold}"])


def _strip_unresolved_prefix(col: str) -> str:
    c = str(col)
    return c[len("unresolved_"):] if c.lower().startswith("unresolved_") else c


def plot_unresolved_mass_by_family(
    unresolved_mass: pd.DataFrame,
    output_dir,
    *,
    name: str = "unresolved_mass_by_family",
    mapping: Optional[Mapping[str, str]] = None,
) -> FigureResult:
    """Mean unresolved family mass per broad family (horizontal bar).

    *unresolved_mass*: samples × ``unresolved_<family>`` columns.
    """
    go = require_plotly()
    df = unresolved_mass.copy()
    # drop a leading non-numeric index/sample column if present
    num = df.select_dtypes("number")
    if num.empty:
        raise ValueError("unresolved-mass table has no numeric columns.")
    means = num.mean(axis=0)
    means.index = [_strip_unresolved_prefix(c) for c in means.index]
    means = means.sort_values(ascending=True)
    fams = list(means.index)
    cm = build_hierarchical_color_map([], {f: f for f in fams})  # broad rows only
    _, broad_colors = color_dicts_from_map(cm)
    colors = [broad_colors.get(f, OTHER_COLOR) for f in fams]

    fig = go.Figure(go.Bar(x=means.values, y=[str(f) for f in fams],
                           orientation="h", marker_color=colors))
    fig.update_layout(**plotly_layout(
        "Unresolved mass by broad family",
        subtitle="mean RNA-derived mass not split into subtypes",
        height=max(340, 32 * len(fams) + 120)))
    fig.update_xaxes(title_text="mean unresolved mass (proportion)")
    fig.update_yaxes(title_text="broad family")
    fig.update_layout(showlegend=False)
    table = means.rename("mean_unresolved_mass").reset_index().rename(
        columns={"index": "broad_family"})
    return export_figure(
        fig, output_dir, name, data={"data": table},
        caption=("Mean unresolved (family-level) mass per broad family; higher "
                 "values mean more mass could not be split into reliable "
                 "subtypes."),
        data_comment=["Mean unresolved mass per family across samples."])


def _level_for(sep_status: str, unresolved: bool) -> str:
    if unresolved:
        return "broad / family-level"
    if str(sep_status).upper() in ("FAIL", "WARNING"):
        return "caution"
    if str(sep_status).upper() == "CAUTION":
        return "selected fine (caution)"
    return "fine"


_LEVEL_COLOR = {"fine": "#3a9b62", "selected fine (caution)": "#d8a93a",
                "caution": "#dd7d34", "broad / family-level": "#636363"}


def plot_trusted_resolution_summary(
    mapping: Mapping[str, str],
    separability: pd.DataFrame,
    output_dir,
    *,
    name: str = "trusted_resolution_summary",
    unresolved_families: Optional[list] = None,
    high_risk_threshold: float = 0.25,
) -> FigureResult:
    """Per-family recommended interpretation level (broad/fine/caution).

    Combines within-family separability and the unresolved-family list into a
    single recommended-level bar per family.
    """
    go = require_plotly()
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    if not mp:
        raise ValueError("trusted_resolution_summary needs a fine→broad mapping.")
    df = separability.copy()
    if "separability_score" not in df.columns and "bhattacharyya" in df.columns:
        df["separability_score"] = 1.0 - df["bhattacharyya"]
    unres = set(map(str, unresolved_families or []))

    fams = sorted(set(mp.values()))
    n_fine = {f: sum(1 for v in mp.values() if v == f) for f in fams}
    # within-family min separability per family (worst pair)
    rows = []
    if {"type_a", "type_b", "separability_score"} <= set(df.columns):
        fa = df["type_a"].astype(str).map(lambda c: mp.get(c, c))
        fb = df["type_b"].astype(str).map(lambda c: mp.get(c, c))
        within = df[fa == fb].copy()
        within["family"] = fa[fa == fb]
        min_sep = within.groupby("family")["separability_score"].min()
    else:
        min_sep = pd.Series(dtype=float)
    for f in fams:
        ms = float(min_sep.get(f, np.nan))
        if np.isnan(ms):
            sep_status = "UNKNOWN" if n_fine[f] > 1 else "PASS"
        elif ms < 0.1:
            sep_status = "FAIL"
        elif ms < high_risk_threshold:
            sep_status = "WARNING"
        elif ms < 0.5:
            sep_status = "CAUTION"
        else:
            sep_status = "PASS"
        level = _level_for(sep_status, f in unres or (n_fine[f] <= 1 and False))
        if f in unres:
            level = "broad / family-level"
        rows.append({"broad_family": f, "n_fine_labels": n_fine[f],
                     "min_within_separability": (None if np.isnan(ms) else round(ms, 3)),
                     "separability_status": sep_status,
                     "unresolved": f in unres,
                     "recommended_level": level})
    table = pd.DataFrame(rows).sort_values("broad_family")
    # encode level as an ordinal bar for a readable visual
    order = {"broad / family-level": 1, "caution": 2,
             "selected fine (caution)": 3, "fine": 4}
    table["_rank"] = table["recommended_level"].map(order).fillna(2)
    colors = [_LEVEL_COLOR.get(l, "#999999") for l in table["recommended_level"]]

    fig = go.Figure(go.Bar(
        x=table["_rank"], y=table["broad_family"].astype(str), orientation="h",
        marker_color=colors, text=table["recommended_level"], textposition="auto",
        customdata=np.stack([table["n_fine_labels"],
                             table["separability_status"]], axis=-1),
        hovertemplate="%{y}<br>fine labels=%{customdata[0]}"
        "<br>separability=%{customdata[1]}<br>level=%{text}<extra></extra>"))
    fig.update_layout(**plotly_layout(
        "Recommended interpretation level by family",
        subtitle="broad-only → caution → selected fine → fine",
        height=max(360, 30 * len(table) + 140)))
    fig.update_xaxes(title_text="recommended resolution (1=broad … 4=fine)",
                     range=[0, 4.4])
    fig.update_yaxes(title_text="broad family")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": table.drop(columns=["_rank"])},
        caption=("Per-family recommended interpretation level, combining "
                 "within-family separability and unresolved-family status; "
                 "families flagged broad/family-level should not be read as "
                 "confident subtypes."),
        data_comment=[f"high_risk_threshold = {high_risk_threshold}"])
