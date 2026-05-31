"""
Plotly separability figures (shared by bulk and spatial reports).

Shows which cell-type pairs are confusable, using the Bhattacharyya-coefficient
matrix from :func:`tissueresolve.reference.separability.separability_heatmap_data`.
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from tissueresolve.plotting.export import FigureResult, export_figure, require_plotly
from tissueresolve.plotting.style import plotly_layout

__all__ = ["plot_separability_heatmap"]


def plot_separability_heatmap(
    matrix: Union["pd.DataFrame", "np.ndarray"],
    output_dir,
    *,
    cell_types: Optional[Sequence[str]] = None,
    name: str = "separability_heatmap",
    title: str = "Pairwise cell-type separability",
) -> FigureResult:
    """Heatmap of the Bhattacharyya coefficient (1 = indistinguishable).

    Accepts a precomputed K×K DataFrame/array.  A
    :class:`~tissueresolve.results.SeparabilityReport` can be converted first
    with ``separability_heatmap_data(report, cell_types)``.
    """
    go = require_plotly()
    if isinstance(matrix, pd.DataFrame):
        M = matrix.to_numpy(dtype=float)
        labels = [str(c) for c in matrix.index]
    else:
        M = np.asarray(matrix, dtype=float)
        labels = [str(c) for c in (cell_types or range(M.shape[0]))]
    df = pd.DataFrame(M, index=labels, columns=labels)

    fig = go.Figure(go.Heatmap(
        z=M, x=labels, y=labels, colorscale="Inferno", zmin=0, zmax=1,
        colorbar={"title": "Bhattacharyya"}))
    fig.update_layout(**plotly_layout(
        title, subtitle="1 = indistinguishable (poorly separable), 0 = orthogonal",
        height=560))
    fig.update_xaxes(tickangle=90)
    return export_figure(
        fig, output_dir, name, data={"data": df},
        caption="Pairwise separability (Bhattacharyya coefficient). "
                "High values mark poorly separable pairs whose estimates are unreliable.",
        data_comment=["Bhattacharyya coefficient; >0.90 = poorly separable"])
