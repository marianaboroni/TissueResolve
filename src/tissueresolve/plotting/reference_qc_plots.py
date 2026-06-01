"""
Reference-QC publication figures (report section 2).

These answer "is the single-cell/nucleus reference good enough?" — composition
by broad family, fine-subpopulation support, cell-type imbalance, gene overlap
with the query, and the suitability-score components.

Every function returns a :class:`~tissueresolve.plotting.export.FigureResult`
(interactive HTML + static when kaleido is present + ``<name>.data.tsv``) and
uses the deterministic hierarchical colour palette so colours match the rest of
the report.  Functions never invent data: when an input is missing they raise a
clear error (callers in the report harness catch it and record a
``missing_data`` manifest row).
"""
from __future__ import annotations

from pathlib import Path
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
    "plot_reference_broad_family_composition",
    "plot_reference_fine_subpopulation_support",
    "plot_reference_celltype_imbalance",
    "plot_gene_overlap_by_modality",
    "plot_reference_suitability_components",
]

# traffic-light colours for suitability statuses
_STATUS_COLOR = {
    "PASS": "#3a9b62", "CAUTION": "#d8a93a", "WARNING": "#dd7d34",
    "FAIL": "#cc4b3e", "UNKNOWN": "#b0b0b0",
}


def _counts_series(counts: Union[Mapping[str, int], pd.Series]) -> pd.Series:
    s = pd.Series(dict(counts)) if not isinstance(counts, pd.Series) else counts.copy()
    s.index = s.index.astype(str)
    return s.astype(float)


def _palettes(fine_types, mapping):
    cm = build_hierarchical_color_map(list(map(str, fine_types)),
                                      {str(k): str(v) for k, v in (mapping or {}).items()})
    return color_dicts_from_map(cm)


def plot_reference_broad_family_composition(
    counts: Union[Mapping[str, int], pd.Series],
    mapping: Mapping[str, str],
    output_dir,
    *,
    name: str = "reference_broad_family_composition",
) -> FigureResult:
    """Horizontal bar of reference cells per broad family (sorted by abundance)."""
    go = require_plotly()
    s = _counts_series(counts)
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    if not mp:
        raise ValueError("broad-family composition needs a fine→broad mapping.")
    fam = {}
    for ct, n in s.items():
        fam[mp.get(ct, "Other")] = fam.get(mp.get(ct, "Other"), 0.0) + n
    table = (pd.DataFrame({"n_cells": pd.Series(fam)})
             .sort_values("n_cells", ascending=True))
    total = float(table["n_cells"].sum()) or 1.0
    table["percent"] = 100.0 * table["n_cells"] / total
    _, broad_colors = _palettes(s.index, mp)
    colors = [broad_colors.get(f, OTHER_COLOR) for f in table.index]

    fig = go.Figure(go.Bar(
        x=table["n_cells"], y=[str(f) for f in table.index], orientation="h",
        marker_color=colors,
        text=[f"{int(n):,} ({p:.0f}%)" for n, p in zip(table["n_cells"], table["percent"])],
        textposition="auto"))
    fig.update_layout(**plotly_layout("Reference composition by broad cell family",
                                      subtitle="number of reference profiles per family",
                                      height=420))
    fig.update_xaxes(title_text="number of reference cells/nuclei")
    fig.update_yaxes(title_text="broad family")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": table.reset_index().rename(
            columns={"index": "broad_family"})},
        caption="Reference cells per broad family (colours = broad-family palette).",
        data_comment=["Reference composition by broad family."])


def plot_reference_fine_subpopulation_support(
    counts: Union[Mapping[str, int], pd.Series],
    mapping: Mapping[str, str],
    output_dir,
    *,
    name: str = "reference_fine_subpopulation_support",
    min_cells: int = 50,
) -> FigureResult:
    """Per-fine-label cell counts, grouped by broad family, with a min-cells line.

    Fine labels are coloured with the subtone of their broad family.
    """
    go = require_plotly()
    s = _counts_series(counts)
    mp = {str(k): str(v) for k, v in (mapping or {}).items()}
    if not mp:
        raise ValueError("fine-subpopulation support needs a fine→broad mapping.")
    fine_colors, _ = _palettes(s.index, mp)
    rows = [{"fine_cell_type": ct, "broad_family": mp.get(ct, "Other"),
             "n_cells": n} for ct, n in s.items()]
    table = pd.DataFrame(rows).sort_values(["broad_family", "n_cells"],
                                           ascending=[True, True])
    colors = [fine_colors.get(ct, OTHER_COLOR) for ct in table["fine_cell_type"]]

    fig = go.Figure(go.Bar(
        x=table["n_cells"], y=table["fine_cell_type"].astype(str),
        orientation="h", marker_color=colors,
        customdata=table["broad_family"],
        hovertemplate="%{y}<br>family=%{customdata}<br>n=%{x}<extra></extra>"))
    fig.add_vline(x=min_cells, line_dash="dash", line_color="firebrick",
                  annotation_text=f"min cells = {min_cells}")
    fig.update_layout(**plotly_layout("Reference support for fine subpopulations",
                                      subtitle="fine labels coloured by broad family",
                                      height=max(420, 18 * len(table) + 120)))
    fig.update_xaxes(title_text="number of reference cells/nuclei")
    fig.update_yaxes(title_text="fine subpopulation")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": table},
        caption=("Reference cells per fine subpopulation, grouped by broad "
                 f"family; dashed line marks the {min_cells}-cell support "
                 "threshold."),
        data_comment=[f"min_cells threshold = {min_cells}"])


def _gini(values: np.ndarray) -> float:
    v = np.sort(np.asarray(values, dtype=float))
    n = v.size
    if n == 0 or v.sum() == 0:
        return float("nan")
    cum = np.cumsum(v)
    return float((n + 1 - 2 * (cum.sum() / cum[-1])) / n)


def plot_reference_celltype_imbalance(
    counts: Union[Mapping[str, int], pd.Series],
    output_dir,
    *,
    name: str = "reference_celltype_imbalance",
    top_label: int = 5,
) -> FigureResult:
    """Ranked cell-type counts + Lorenz curve; reports a Gini imbalance score."""
    go = require_plotly()
    s = _counts_series(counts).sort_values(ascending=False)
    total = float(s.sum()) or 1.0
    table = pd.DataFrame({"cell_type": s.index, "n_cells": s.values})
    table["fraction"] = table["n_cells"] / total
    table["cumulative_fraction"] = table["fraction"].cumsum()
    gini = _gini(s.values)

    fig = go.Figure(go.Bar(
        x=table["cell_type"].astype(str), y=table["n_cells"],
        marker_color="#4e79a7",
        text=[f"{int(n):,}" if i < top_label else "" for i, n in enumerate(table["n_cells"])],
        textposition="outside"))
    fig.update_layout(**plotly_layout(
        "Reference imbalance across cell types",
        subtitle=f"Gini imbalance score = {gini:.2f} (0 = even, 1 = one type dominates)",
        height=460))
    fig.update_xaxes(title_text="cell type (ranked by abundance)", tickangle=90)
    fig.update_yaxes(title_text="number of reference cells/nuclei")
    fig.update_layout(showlegend=False)
    return export_figure(
        fig, output_dir, name, data={"data": table},
        caption=(f"Cell types ranked by abundance; Gini imbalance = {gini:.2f}. "
                 "A few dominant types can destabilise rare-type signatures."),
        data_comment=[f"gini_imbalance = {gini:.4f}"])


def plot_gene_overlap_by_modality(
    overlap: Mapping[str, Mapping[str, int]],
    output_dir,
    *,
    name: str = "gene_overlap_by_modality",
) -> FigureResult:
    """Grouped bar of reference / query / shared gene counts per input modality.

    *overlap*: ``{modality: {"n_reference":, "n_query":, "n_shared":}}``.
    """
    go = require_plotly()
    rows = []
    for modality, d in (overlap or {}).items():
        rows.append({"modality": str(modality),
                     "reference_genes": int(d.get("n_reference", 0)),
                     "query_genes": int(d.get("n_query", 0)),
                     "shared_genes": int(d.get("n_shared", 0))})
    if not rows:
        raise ValueError("gene_overlap_by_modality needs at least one modality.")
    table = pd.DataFrame(rows).set_index("modality")
    cats = ["reference_genes", "query_genes", "shared_genes"]
    palette = {"reference_genes": "#9c755f", "query_genes": "#4e79a7",
               "shared_genes": "#59a14f"}
    fig = go.Figure()
    for cat in cats:
        fig.add_bar(x=table.index.astype(str), y=table[cat].to_numpy(),
                    name=cat.replace("_", " "), marker_color=palette[cat])
    fig.update_layout(barmode="group",
                      **plotly_layout("Gene overlap between reference and query data",
                                      subtitle="shared genes are usable for deconvolution",
                                      height=420))
    fig.update_xaxes(title_text="input modality")
    fig.update_yaxes(title_text="number of genes")
    fig.update_layout(legend={"title": {"text": "gene set"}})
    return export_figure(
        fig, output_dir, name, data={"data": table.reset_index()},
        caption=("Reference, query and shared gene counts per modality; only "
                 "shared genes inform the deconvolution."),
        data_comment=["Gene overlap by modality."])


def plot_reference_suitability_components(
    components: pd.DataFrame,
    output_dir,
    *,
    name: str = "reference_suitability_components",
) -> FigureResult:
    """Horizontal traffic-light bar of suitability components (status-coloured)."""
    go = require_plotly()
    df = components.copy()
    if "component" not in df.columns:
        df = df.reset_index().rename(columns={df.index.name or "index": "component"})
    if "component" not in df.columns or "status" not in df.columns:
        raise ValueError("suitability components need 'component' and 'status' columns.")
    df["status"] = df["status"].astype(str).str.upper()
    # score may be missing/NaN for UNKNOWN components; show a unit bar so the
    # status colour is still visible.
    score = pd.to_numeric(df.get("score"), errors="coerce")
    df["_bar"] = score.fillna(1.0).clip(0, 1)
    colors = [_STATUS_COLOR.get(s, _STATUS_COLOR["UNKNOWN"]) for s in df["status"]]

    fig = go.Figure(go.Bar(
        x=df["_bar"], y=df["component"].astype(str), orientation="h",
        marker_color=colors, text=df["status"], textposition="auto",
        customdata=score, hovertemplate="%{y}<br>status=%{text}<br>score=%{customdata}<extra></extra>"))
    fig.update_layout(**plotly_layout("Reference suitability components",
                                      subtitle="traffic-light QC per component",
                                      height=max(360, 34 * len(df) + 120)))
    fig.update_xaxes(title_text="component score (0–1; full bar when score N/A)",
                     range=[0, 1.02])
    fig.update_yaxes(title_text="suitability component")
    fig.update_layout(showlegend=False)
    keep = [c for c in ("component", "score", "status", "detail") if c in df.columns]
    return export_figure(
        fig, output_dir, name, data={"data": df[keep]},
        caption=("Suitability components coloured PASS/CAUTION/WARNING/FAIL/"
                 "UNKNOWN; WARNING/FAIL components should be addressed before "
                 "trusting fine predictions."),
        data_comment=["Reference suitability components (traffic light)."])
