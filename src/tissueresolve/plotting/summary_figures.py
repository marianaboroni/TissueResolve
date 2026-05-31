"""
Main multi-panel publication summary figures for bulk and spatial runs.

One clean, white-background, family-coloured figure per modality that answers
the headline questions: reference quality, predicted composition, reliability,
and which types are confusable.  Main panels show top cell types + "Other";
full matrices live in the report's collapsible detailed sections.

Each figure writes ``<name>.html`` (+ static when kaleido is present),
``<name>.data.tsv`` and ``<name>.caption.txt``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.palette import (
    assign_family_palette, shorten_cell_type_label,
)

__all__ = ["bulk_main_summary_figure", "spatial_main_summary_figure"]

_BULK_SUB = "mRNA-derived proportions, not absolute cell fractions"
_SP_SUB = "spot-level RNA-derived composition, not cell counts"


def _top_n_other(props: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, list[str]]:
    order = props.mean(axis=0).sort_values(ascending=False).index.tolist()
    keep = order[:top_n]
    out = props[keep].copy()
    if len(order) > top_n:
        out["Other"] = props[order[top_n:]].sum(axis=1)
    return out, list(out.columns)


def _schematic(fig, row, col, steps: list[str]):
    """A simple left→right boxes-and-arrows workflow schematic."""
    n = len(steps)
    xs = np.linspace(0.06, 0.94, n)
    for i, (x, label) in enumerate(zip(xs, steps)):
        fig.add_annotation(x=x, y=0.5, xref=f"x{_axn(row, col)}",
                           yref=f"y{_axn(row, col)}", text=f"<b>{label}</b>",
                           showarrow=False, font={"size": 10},
                           bordercolor="#4C72B0", borderwidth=1, borderpad=4,
                           bgcolor="#eef3fb")
        if i < n - 1:
            fig.add_annotation(x=(xs[i] + xs[i + 1]) / 2, y=0.5,
                               xref=f"x{_axn(row, col)}", yref=f"y{_axn(row, col)}",
                               text="→", showarrow=False, font={"size": 18})
    fig.update_xaxes(visible=False, range=[0, 1], row=row, col=col)
    fig.update_yaxes(visible=False, range=[0, 1], row=row, col=col)


def _axn(row, col, ncols=3):
    n = (row - 1) * ncols + col
    return "" if n == 1 else str(n)


def _write_caption(out_dir: Path, name: str, caption: str, panels: dict[str, str]):
    body = caption + "\n\nPanels:\n" + "\n".join(
        f"  {k}: {v}" for k, v in panels.items())
    (Path(out_dir) / f"{name}.caption.txt").write_text(body, encoding="utf-8")


# ---------------------------------------------------------------------------
# Bulk
# ---------------------------------------------------------------------------


def bulk_main_summary_figure(
    proportions: pd.DataFrame, output_dir, *,
    reference_family_counts: Optional[pd.Series] = None,
    validation: Optional[pd.DataFrame] = None,   # columns true,pred (long) optional
    top_pairs: Optional[pd.DataFrame] = None,     # type_a,type_b,separability_score
    color_map: Optional[dict] = None,
    top_n: int = 10, name: str = "bulk_main_summary_figure",
) -> FigureResult:
    """Multi-panel bulk summary: workflow, reference, composition, QC, separability."""
    go = require_plotly()
    from plotly.subplots import make_subplots

    comp, cts = _top_n_other(proportions.fillna(0.0), top_n)
    cmap = color_map or assign_family_palette(list(proportions.columns) + ["Other"])

    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=("A. Workflow", "B. Reference composition",
                        "C. Top non-separable pairs",
                        "D. Bulk composition (clustered)", "E. Validation / QC",
                        "F. Notes"),
        specs=[[{"type": "xy"}, {"type": "xy"}, {"type": "xy"}],
               [{"type": "xy"}, {"type": "xy"}, {"type": "xy"}]])

    _schematic(fig, 1, 1, ["scRNA ref", "signature", "bulk", "deconv", "mRNA prop."])

    if reference_family_counts is not None and len(reference_family_counts):
        rc = reference_family_counts.sort_values()
        fig.add_trace(go.Bar(y=[str(i) for i in rc.index], x=rc.to_numpy(),
                             orientation="h", marker_color="#4C72B0",
                             showlegend=False), row=1, col=2)

    if top_pairs is not None and len(top_pairs):
        tp = top_pairs.head(10)
        labels = [f"{shorten_cell_type_label(a)} / {shorten_cell_type_label(b)}"
                  for a, b in zip(tp["type_a"], tp["type_b"])]
        score = tp["separability_score"] if "separability_score" in tp else tp.iloc[:, -1]
        fig.add_trace(go.Bar(y=labels, x=score.to_numpy(), orientation="h",
                             marker_color="#C44E52", showlegend=False), row=1, col=3)
        fig.update_xaxes(title_text="separability score", row=1, col=3)

    samples = comp.index.astype(str).tolist()
    for ct in cts:
        fig.add_trace(go.Bar(x=samples, y=comp[ct].to_numpy(),
                             name=shorten_cell_type_label(ct),
                             legendgroup=str(ct),
                             marker_color=cmap.get(ct, "#999999")), row=2, col=1)
    fig.update_yaxes(title_text="mRNA-derived prop.", range=[0, 1], row=2, col=1)

    if validation is not None and {"true", "pred"} <= set(validation.columns):
        fig.add_trace(go.Scatter(x=validation["true"], y=validation["pred"],
                                 mode="markers", marker={"size": 5, "opacity": 0.5,
                                 "color": "#4C72B0"}, showlegend=False), row=2, col=2)
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                                 line={"dash": "dash", "color": "#888"},
                                 showlegend=False), row=2, col=2)
        fig.update_xaxes(title_text="true", row=2, col=2)
        fig.update_yaxes(title_text="predicted", row=2, col=2)
    else:
        fig.add_annotation(x=0.5, y=0.5, xref=f"x{_axn(2, 2)}", yref=f"y{_axn(2, 2)}",
                           text="No ground truth;<br>see QC panel/tables",
                           showarrow=False, font={"size": 10})
        fig.update_xaxes(visible=False, row=2, col=2)
        fig.update_yaxes(visible=False, row=2, col=2)

    fig.add_annotation(
        x=0.5, y=0.5, xref=f"x{_axn(2, 3)}", yref=f"y{_axn(2, 3)}", showarrow=False,
        font={"size": 9}, align="left",
        text=("Top types shown; rest in 'Other'.<br>"
              "Confusable types → interpret at<br>family level (see report)."))
    fig.update_xaxes(visible=False, row=2, col=3)
    fig.update_yaxes(visible=False, row=2, col=3)

    fig.update_layout(
        barmode="stack", template="plotly_white", height=760, width=1100,
        title={"text": "TissueResolve — bulk deconvolution summary"
               f"<br><sup>{_BULK_SUB}</sup>", "x": 0.5},
        legend={"font": {"size": 8}, "orientation": "h", "y": -0.08})

    caption = ("Main bulk summary. Composition shows the top cell types "
               "(rest in 'Other') with family-aware colours; samples can be "
               f"clustered by composition. {_BULK_SUB}.")
    _write_caption(output_dir, name, caption, {
        "A": "workflow schematic", "B": "reference composition by family",
        "C": "top non-separable pairs", "D": "bulk composition (top types + Other)",
        "E": "validation scatter or QC", "F": "interpretation notes"})
    return export_figure(fig, output_dir, name,
                         data={"data": comp,
                               "reference_family": (reference_family_counts.to_frame("n_cells")
                                                    if reference_family_counts is not None
                                                    else pd.DataFrame())},
                         caption=caption,
                         data_comment=["estimate_type: mRNA_proportion (NOT cell fractions)"])


# ---------------------------------------------------------------------------
# Spatial
# ---------------------------------------------------------------------------


def spatial_main_summary_figure(
    proportions: pd.DataFrame, output_dir, *,
    coords: Optional[pd.DataFrame] = None,        # array_row, array_col
    morans: Optional[pd.Series] = None,
    reference_family_counts: Optional[pd.Series] = None,
    color_map: Optional[dict] = None,
    top_n: int = 10, name: str = "spatial_main_summary_figure",
) -> FigureResult:
    """Multi-panel spatial summary: workflow, reference, dominant map, composition,
    Moran's I."""
    go = require_plotly()
    from plotly.subplots import make_subplots

    comp, cts = _top_n_other(proportions.fillna(0.0), top_n)
    cmap = color_map or assign_family_palette(list(proportions.columns) + ["Other"])
    mean_comp = comp.mean(axis=0)

    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=("A. Workflow", "B. Reference composition",
                        "C. Dominant cell type (spatial)",
                        "D. Mean composition", "E. Moran's I (spatial structure)",
                        "F. Notes"),
        specs=[[{"type": "xy"}, {"type": "xy"}, {"type": "xy"}],
               [{"type": "xy"}, {"type": "xy"}, {"type": "xy"}]])

    _schematic(fig, 1, 1, ["scRNA ref", "Visium", "graph/model", "spot comp.", "maps"])

    if reference_family_counts is not None and len(reference_family_counts):
        rc = reference_family_counts.sort_values()
        fig.add_trace(go.Bar(y=[str(i) for i in rc.index], x=rc.to_numpy(),
                             orientation="h", marker_color="#4C72B0",
                             showlegend=False), row=1, col=2)

    if coords is not None and {"array_row", "array_col"}.issubset(coords.columns):
        dom = proportions.idxmax(axis=1).astype(str)
        ar = coords["array_row"].to_numpy(dtype=float)
        ac = coords["array_col"].to_numpy(dtype=float)
        for ct in dom.unique():
            mask = (dom == ct).to_numpy()
            fig.add_trace(go.Scatter(x=ac[mask], y=ar[mask], mode="markers",
                                     marker={"color": cmap.get(ct, "#999"), "size": 5},
                                     name=shorten_cell_type_label(ct),
                                     legendgroup=str(ct), showlegend=False),
                          row=1, col=3)
        fig.update_yaxes(autorange="reversed", row=1, col=3)
    else:
        fig.add_annotation(x=0.5, y=0.5, xref=f"x{_axn(1, 3)}", yref=f"y{_axn(1, 3)}",
                           text="No coordinates available", showarrow=False,
                           font={"size": 10})
        fig.update_xaxes(visible=False, row=1, col=3)
        fig.update_yaxes(visible=False, row=1, col=3)

    fig.add_trace(go.Bar(x=[shorten_cell_type_label(c) for c in mean_comp.index],
                         y=mean_comp.to_numpy(),
                         marker_color=[cmap.get(c, "#999") for c in mean_comp.index],
                         showlegend=False), row=2, col=1)
    fig.update_yaxes(title_text="mean spot comp.", row=2, col=1)
    fig.update_xaxes(tickangle=90, row=2, col=1)

    if morans is not None and len(morans):
        s = morans.sort_values(ascending=False).head(12)
        fig.add_trace(go.Bar(x=[shorten_cell_type_label(c) for c in s.index],
                             y=s.to_numpy(),
                             marker_color=[cmap.get(c, "#55A868") for c in s.index],
                             showlegend=False), row=2, col=2)
        fig.update_yaxes(title_text="Moran's I", row=2, col=2)
        fig.update_xaxes(tickangle=90, row=2, col=2)

    fig.add_annotation(
        x=0.5, y=0.5, xref=f"x{_axn(2, 3)}", yref=f"y{_axn(2, 3)}", showarrow=False,
        font={"size": 9}, align="left",
        text=("Spot-level RNA-derived composition,<br>not cell counts. "
              "Top types shown;<br>confusable types → family level."))
    fig.update_xaxes(visible=False, row=2, col=3)
    fig.update_yaxes(visible=False, row=2, col=3)

    fig.update_layout(
        barmode="stack", template="plotly_white", height=760, width=1100,
        title={"text": "TissueResolve — spatial deconvolution summary"
               f"<br><sup>{_SP_SUB}</sup>", "x": 0.5})

    caption = ("Main spatial summary: reference composition, dominant cell-type "
               f"map, mean composition and Moran's I spatial structure. {_SP_SUB}.")
    _write_caption(output_dir, name, caption, {
        "A": "workflow schematic", "B": "reference composition by family",
        "C": "dominant cell type per spot", "D": "mean spot composition (top + Other)",
        "E": "Moran's I per cell type", "F": "interpretation notes"})
    return export_figure(fig, output_dir, name,
                         data={"data": comp, "mean_composition": mean_comp.to_frame("mean")},
                         caption=caption,
                         data_comment=["estimate_type: spot_rna_composition (NOT cell counts)"])
