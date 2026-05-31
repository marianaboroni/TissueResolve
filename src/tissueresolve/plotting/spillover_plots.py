"""
Plotly spillover figures (shared by bulk and spatial reports).

Visualises the spillover matrix from
:func:`tissueresolve.benchmark.spillover.estimate_spillover_matrix`
(row = true type, column = predicted type, rows sum to 1) as a heatmap and as a
directed leakage network.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.style import color_sequence, plotly_layout

__all__ = ["plot_spillover_heatmap", "plot_spillover_network"]


def plot_spillover_heatmap(
    spillover_matrix: pd.DataFrame, output_dir, *,
    name: str = "spillover_heatmap",
    title: str = "Spillover matrix (true → predicted)",
) -> FigureResult:
    """Heatmap of estimated cross-type leakage (rows = true type, sum to 1)."""
    go = require_plotly()
    M = spillover_matrix.to_numpy(dtype=float)
    rows = [str(r) for r in spillover_matrix.index]
    cols = [str(c) for c in spillover_matrix.columns]
    fig = go.Figure(go.Heatmap(
        z=M, x=cols, y=rows, colorscale="Magma", zmin=0, zmax=1,
        colorbar={"title": "fraction"}))
    fig.update_layout(**plotly_layout(
        title, subtitle="diagonal = self-retention; off-diagonal = leakage",
        height=560))
    fig.update_xaxes(title_text="predicted type", tickangle=90)
    fig.update_yaxes(title_text="true type")
    return export_figure(
        fig, output_dir, name, data={"data": spillover_matrix},
        caption="Spillover: deconvolving a pure type and seeing mass land on "
                "another. Off-diagonal mass marks confusable pairs.",
        data_comment=["row = true type, column = predicted; rows sum to 1"])


def plot_spillover_network(
    spillover_matrix: pd.DataFrame, output_dir, *,
    threshold: float = 0.10, name: str = "spillover_network",
) -> FigureResult:
    """Directed leakage network: an edge a→b when type *a* leaks ≥ *threshold* into *b*.

    Nodes are laid out on a circle (no networkx dependency).  Edge source data
    are saved alongside the figure.
    """
    go = require_plotly()
    cts = [str(c) for c in spillover_matrix.index]
    n = len(cts)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xs, ys = np.cos(angles), np.sin(angles)
    pos = {ct: (xs[i], ys[i]) for i, ct in enumerate(cts)}
    colours = color_sequence(n)

    edges = []
    edge_x, edge_y = [], []
    M = spillover_matrix.to_numpy(dtype=float)
    cols = [str(c) for c in spillover_matrix.columns]
    for i, a in enumerate(cts):
        for j, b in enumerate(cols):
            if a == b:
                continue
            frac = float(M[i, j])
            if frac >= threshold:
                edges.append({"from": a, "into": b, "fraction": round(frac, 4)})
                (x0, y0), (x1, y1) = pos[a], pos.get(b, (0, 0))
                edge_x += [x0, x1, None]
                edge_y += [y0, y1, None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                             line={"color": "rgba(120,120,120,0.4)", "width": 1},
                             hoverinfo="none", showlegend=False))
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="markers+text", text=cts, textposition="top center",
        marker={"size": 14, "color": colours}, showlegend=False,
        hovertext=cts, hoverinfo="text"))
    layout = plotly_layout(f"Spillover network (≥ {threshold:.0%} leakage)",
                           subtitle="edge a→b: type a leaks into type b", height=600)
    layout["xaxis"] = {"visible": False}
    layout["yaxis"] = {"visible": False}
    fig.update_layout(**layout)

    edges_df = pd.DataFrame(edges, columns=["from", "into", "fraction"])
    return export_figure(
        fig, output_dir, name, data={"data": edges_df},
        caption=f"Directed spillover network (edges ≥ {threshold:.0%}). "
                "Dense neighbourhoods are confusable families.",
        data_comment=[f"edges with leakage >= {threshold}"])
