"""
Signature-quality & hierarchy publication figures (report section 3).

These answer "are the signatures good enough to tell cell types apart, and is
the broad→fine hierarchy usable?" — signature heatmap, most-confusable pairs,
within- vs between-family separability, the annotation hierarchy, and per-family
marker support.

All figures use the deterministic hierarchical palette and save source data.
Inputs come from the reference signature matrix and the
``resolution/pairwise_separability.tsv`` diagnostics.
"""
from __future__ import annotations

from typing import Mapping, Optional, Union

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
    "plot_signature_matrix_heatmap",
    "plot_top_confusable_pairs",
    "plot_within_vs_between_family_separability",
    "plot_hierarchy_map",
    "plot_marker_support_by_family",
]


def _palettes(fine_types, mapping):
    cm = build_hierarchical_color_map(
        list(map(str, fine_types)),
        {str(k): str(v) for k, v in (mapping or {}).items()})
    return color_dicts_from_map(cm)


def plot_signature_matrix_heatmap(
    matrix: np.ndarray,
    gene_names,
    cell_types,
    mapping: Mapping[str, str],
    output_dir,
    *,
    name: str = "signature_matrix_heatmap",
    top_n_genes: int = 40,
) -> FigureResult:
    """Heatmap of the most informative signature genes × cell types.

    Columns are grouped by broad family; only the *top_n_genes* most variable
    genes are shown (full matrix stays as source data).  Values are row z-scored
    for visualization only.
    """
    go = require_plotly()
    M = np.asarray(matrix, dtype=float)
    genes = [str(g) for g in gene_names]
    cts = [str(c) for c in cell_types]
    # Accept either orientation: (genes × cell_types) or (cell_types × genes).
    if M.shape == (len(cts), len(genes)) and len(cts) != len(genes):
        M = M.T
    if M.shape != (len(genes), len(cts)):
        raise ValueError(f"matrix shape {M.shape} not compatible with "
                         f"(genes={len(genes)}, cell_types={len(cts)})")
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    # order columns by family then name
    col_order = sorted(range(len(cts)), key=lambda j: (mp.get(cts[j], "ZZ"), cts[j]))
    cts_o = [cts[j] for j in col_order]
    Mo = M[:, col_order]
    # top genes by variance across cell types
    var = Mo.var(axis=1)
    top = np.argsort(var)[::-1][:max(1, min(top_n_genes, len(genes)))]
    top = sorted(top)
    sub = Mo[top, :]
    sub_genes = [genes[i] for i in top]
    # row z-score for display (guard zero variance)
    mu = sub.mean(axis=1, keepdims=True)
    sd = sub.std(axis=1, keepdims=True)
    z = np.divide(sub - mu, sd, out=np.zeros_like(sub), where=sd > 0)

    fig = go.Figure(go.Heatmap(
        z=z, x=cts_o, y=sub_genes, colorscale="RdBu_r", zmid=0,
        colorbar={"title": "row z-score"}))
    fig.update_layout(**plotly_layout(
        "Reference signature matrix",
        subtitle=f"top {len(sub_genes)} most variable genes; columns grouped by family",
        height=max(460, 12 * len(sub_genes) + 160)))
    fig.update_xaxes(title_text="cell type (grouped by broad family)", tickangle=90)
    fig.update_yaxes(title_text="signature gene")
    display = pd.DataFrame(z, index=sub_genes, columns=cts_o)
    return export_figure(
        fig, output_dir, name, data={"data": display},
        caption=(f"Top {len(sub_genes)} most variable signature genes (rows) "
                 "across cell types (columns, grouped by broad family); colour "
                 "is the row z-scored signature value (display only)."),
        data_comment=["Row z-scored for display; full matrix is the reference."])


def _severity(score: float) -> str:
    if score < 0.1:
        return "FAIL"
    if score < 0.25:
        return "WARNING"
    if score < 0.5:
        return "CAUTION"
    return "PASS"


_SEV_COLOR = {"FAIL": "#cc4b3e", "WARNING": "#dd7d34", "CAUTION": "#d8a93a",
              "PASS": "#3a9b62"}


def plot_top_confusable_pairs(
    separability: pd.DataFrame,
    output_dir,
    *,
    name: str = "top_confusable_pairs",
    top_n: int = 20,
) -> FigureResult:
    """Lollipop of the *top_n* least-separable (most confusable) cell-type pairs."""
    go = require_plotly()
    df = separability.copy()
    if not {"type_a", "type_b"} <= set(df.columns):
        raise ValueError("separability needs 'type_a' and 'type_b' columns.")
    score_col = "separability_score" if "separability_score" in df.columns else None
    if score_col is None and "bhattacharyya" in df.columns:
        df["separability_score"] = 1.0 - df["bhattacharyya"]
        score_col = "separability_score"
    if score_col is None:
        raise ValueError("separability needs 'separability_score' or 'bhattacharyya'.")
    df = df.sort_values(score_col, ascending=True).head(top_n).copy()
    df["pair"] = df["type_a"].astype(str) + " ↔ " + df["type_b"].astype(str)
    df["severity"] = df[score_col].apply(_severity)
    colors = [_SEV_COLOR[s] for s in df["severity"]]
    # plot in ascending so the most confusable is at top
    df = df.iloc[::-1]
    colors = colors[::-1]

    fig = go.Figure()
    for i, (_, r) in enumerate(df.iterrows()):
        fig.add_trace(go.Scatter(
            x=[0, r[score_col]], y=[r["pair"], r["pair"]], mode="lines",
            line={"color": "#cccccc", "width": 2}, showlegend=False,
            hoverinfo="skip"))
    fig.add_trace(go.Scatter(
        x=df[score_col], y=df["pair"], mode="markers",
        marker={"color": colors, "size": 11},
        text=df["severity"], hovertemplate="%{y}<br>separability=%{x:.3f}"
        "<br>%{text}<extra></extra>", showlegend=False))
    fig.update_layout(**plotly_layout(
        "Most difficult-to-distinguish cell-type pairs",
        subtitle=f"{len(df)} least-separable pairs (low = confusable)",
        height=max(420, 22 * len(df) + 140)))
    fig.update_xaxes(title_text="separability score (1 − Bhattacharyya)", range=[0, 1])
    fig.update_yaxes(title_text="cell-type pair")
    keep = [c for c in ("type_a", "type_b", score_col, "bhattacharyya",
                        "n_discriminating_genes", "resolvability", "severity")
            if c in df.columns]
    return export_figure(
        fig, output_dir, name, data={"data": df[keep][::-1]},
        caption=(f"The {len(df)} least-separable cell-type pairs (1 − "
                 "Bhattacharyya); low scores (red) mark pairs whose fine labels "
                 "are unreliable — interpret them at the family level."),
        data_comment=["Top confusable pairs (ascending separability)."])


def plot_within_vs_between_family_separability(
    separability: pd.DataFrame,
    mapping: Mapping[str, str],
    output_dir,
    *,
    name: str = "within_vs_between_family_separability",
) -> FigureResult:
    """Box+points of separability for within-family vs between-family pairs."""
    go = require_plotly()
    df = separability.copy()
    if not {"type_a", "type_b"} <= set(df.columns):
        raise ValueError("separability needs 'type_a' and 'type_b' columns.")
    if "separability_score" not in df.columns and "bhattacharyya" in df.columns:
        df["separability_score"] = 1.0 - df["bhattacharyya"]
    if "separability_score" not in df.columns:
        raise ValueError("separability needs 'separability_score' or 'bhattacharyya'.")
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    fam_a = df["type_a"].astype(str).map(lambda c: mp.get(c, c))
    fam_b = df["type_b"].astype(str).map(lambda c: mp.get(c, c))
    df["group"] = np.where(fam_a == fam_b, "within-family", "between-family")

    fig = go.Figure()
    palette = {"within-family": "#e15759", "between-family": "#4e79a7"}
    for grp in ("within-family", "between-family"):
        vals = df.loc[df["group"] == grp, "separability_score"]
        if len(vals):
            fig.add_trace(go.Box(y=vals, name=grp, boxpoints="all",
                                 jitter=0.4, pointpos=0,
                                 marker_color=palette[grp]))
    fig.update_layout(**plotly_layout(
        "Separability within and between broad families",
        subtitle="within-family pairs are typically harder to separate",
        height=460))
    fig.update_yaxes(title_text="separability score (1 − Bhattacharyya)", range=[0, 1])
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": df[["type_a", "type_b", "group",
                                                  "separability_score"]]},
        caption=("Separability of within-family vs between-family cell-type "
                 "pairs; lower within-family separability means subtypes of a "
                 "family are hard to resolve."),
        data_comment=["Within vs between family separability."])


def plot_hierarchy_map(
    counts: Union[Mapping[str, int], pd.Series],
    mapping: Mapping[str, str],
    output_dir,
    *,
    name: str = "hierarchy_map",
    unresolved_families: Optional[list] = None,
) -> FigureResult:
    """Treemap of the broad→fine annotation hierarchy, sized by reference cells.

    Unresolved families are flagged in the label.
    """
    go = require_plotly()
    s = (pd.Series(dict(counts)) if not isinstance(counts, pd.Series) else counts).astype(float)
    s.index = s.index.astype(str)
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    if not mp:
        raise ValueError("hierarchy_map needs a fine→broad mapping.")
    fine_colors, broad_colors = _palettes(s.index, mp)
    unres = set(map(str, unresolved_families or []))

    labels, parents, values, colors, rows = [], [], [], [], []
    fam_total: dict[str, float] = {}
    for ct, n in s.items():
        fam_total[mp.get(ct, "Other")] = fam_total.get(mp.get(ct, "Other"), 0.0) + n
    for fam, tot in fam_total.items():
        flag = " ⚠ unresolved" if fam in unres else ""
        labels.append(f"{fam}{flag}")
        parents.append("")
        values.append(tot)
        colors.append(broad_colors.get(fam, OTHER_COLOR))
        rows.append({"node": fam, "parent": "", "level": "broad", "n_cells": tot,
                     "unresolved": fam in unres})
    for ct, n in s.items():
        fam = mp.get(ct, "Other")
        flag = " ⚠ unresolved" if fam in unres else ""
        labels.append(str(ct))
        parents.append(f"{fam}{flag}")
        values.append(float(n))
        colors.append(fine_colors.get(ct, OTHER_COLOR))
        rows.append({"node": ct, "parent": fam, "level": "fine", "n_cells": float(n),
                     "unresolved": fam in unres})

    fig = go.Figure(go.Treemap(
        labels=labels, parents=parents, values=values, branchvalues="total",
        marker_colors=colors, hovertemplate="%{label}<br>cells=%{value}<extra></extra>"))
    fig.update_layout(**plotly_layout(
        "Broad-to-fine annotation hierarchy",
        subtitle="box size = reference cells; ⚠ = family kept at family level",
        height=520))
    return export_figure(
        fig, output_dir, name, data={"data": pd.DataFrame(rows)},
        caption=("Annotation hierarchy: broad families contain fine "
                 "subpopulations, box size = number of reference cells; ⚠ marks "
                 "families reported only at the family level (unresolved)."),
        data_comment=["Hierarchy map (broad and fine nodes)."])


def plot_marker_support_by_family(
    separability: pd.DataFrame,
    mapping: Mapping[str, str],
    output_dir,
    *,
    name: str = "marker_support_by_family",
    warn_threshold: int = 25,
) -> FigureResult:
    """Mean number of discriminating markers for within-family subtype pairs."""
    go = require_plotly()
    df = separability.copy()
    if "n_discriminating_genes" not in df.columns:
        raise ValueError("marker support needs 'n_discriminating_genes'.")
    if not {"type_a", "type_b"} <= set(df.columns):
        raise ValueError("separability needs 'type_a' and 'type_b' columns.")
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    fam_a = df["type_a"].astype(str).map(lambda c: mp.get(c, c))
    fam_b = df["type_b"].astype(str).map(lambda c: mp.get(c, c))
    within = df[fam_a == fam_b].copy()
    within["family"] = fam_a[fam_a == fam_b]
    if within.empty:
        raise ValueError("no within-family pairs to summarise marker support.")
    agg = (within.groupby("family")["n_discriminating_genes"]
           .mean().sort_values(ascending=True))
    _, broad_colors = _palettes(df["type_a"].astype(str).tolist(), mp)
    colors = [broad_colors.get(f, OTHER_COLOR) for f in agg.index]

    fig = go.Figure(go.Bar(
        x=agg.values, y=[str(f) for f in agg.index], orientation="h",
        marker_color=colors))
    fig.add_vline(x=warn_threshold, line_dash="dash", line_color="firebrick",
                  annotation_text=f"low support < {warn_threshold}")
    fig.update_layout(**plotly_layout(
        "Subtype marker support within each broad family",
        subtitle="mean discriminating genes per within-family subtype pair",
        height=max(360, 30 * len(agg) + 140)))
    fig.update_xaxes(title_text="mean number of discriminating markers")
    fig.update_yaxes(title_text="broad family")
    fig.update_layout(showlegend=False)
    table = agg.rename("mean_discriminating_genes").reset_index()
    return export_figure(
        fig, output_dir, name, data={"data": table},
        caption=("Mean discriminating markers for within-family subtype pairs; "
                 f"families below {warn_threshold} (dashed) have weak subtype "
                 "marker support."),
        data_comment=[f"warn_threshold = {warn_threshold} markers"])
