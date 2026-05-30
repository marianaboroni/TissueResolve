"""
Shared QC plots used by both bulk and spatial workflows.

* separability heatmap (Bhattacharyya coefficient, K × K)
* protocol-risk summary (excluded / down-weighted gene counts + risk level)
* gene-overlap summary (shared / query-only / reference-only)
* warning summary (counts of surfaced warnings by category)

Every function saves the figure and its source data and returns a
:class:`~tissueresolve.plotting.style.PlotResult`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence, Union

import numpy as np
import pandas as pd

from tissueresolve.plotting import captions
from tissueresolve.plotting.style import PlotResult, new_figure, save_outputs
from tissueresolve.reference.separability import separability_heatmap_data
from tissueresolve.results import SeparabilityReport

__all__ = [
    "separability_heatmap",
    "protocol_risk_summary",
    "gene_overlap_plot",
    "warning_summary_plot",
]


def separability_heatmap(
    report: SeparabilityReport,
    cell_types: Sequence[str],
    output_dir: Union[str, Path],
    *,
    warn_threshold: float = 0.90,
    name: str = "separability_heatmap",
) -> PlotResult:
    """K × K Bhattacharyya-coefficient heatmap (worse pairs are brighter)."""
    cts = list(cell_types)
    M = separability_heatmap_data(report, cts)

    fig, ax = new_figure(figsize=(max(4.5, 0.5 * len(cts) + 2),
                                  max(4.0, 0.5 * len(cts) + 1.5)))
    im = ax.imshow(M, cmap="rocket" if _has_cmap("rocket") else "inferno",
                   vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(cts)))
    ax.set_xticklabels(cts, rotation=90)
    ax.set_yticks(np.arange(len(cts)))
    ax.set_yticklabels(cts)
    ax.set_title("Pairwise separability (Bhattacharyya coeff.)")
    ax.grid(visible=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Bhattacharyya coeff. (1 = indistinguishable)")
    # Annotate the problematic pairs so warnings are visible on the figure.
    for p in report.pairs:
        if p.is_problematic and p.type_a in cts and p.type_b in cts:
            i, j = cts.index(p.type_a), cts.index(p.type_b)
            ax.text(j, i, "!", ha="center", va="center", color="white", fontsize=9)

    data = pd.DataFrame(M, index=cts, columns=cts)
    fig_paths, data_paths = save_outputs(
        fig, data, output_dir, name,
        data_comment=[
            f"warn_threshold(BC): {warn_threshold}",
            f"n_critical: {report.n_critical}", f"n_high: {report.n_high}",
        ],
    )
    return PlotResult(
        figure=fig, figure_paths=fig_paths, data_paths=data_paths,
        caption=captions.separability_heatmap_caption(
            n_critical=report.n_critical, n_high=report.n_high,
            warn_threshold=warn_threshold,
        ),
    )


def protocol_risk_summary(
    risk,  # ProtocolRiskReport
    output_dir: Union[str, Path],
    *,
    name: str = "protocol_risk_summary",
) -> PlotResult:
    """Bar of excluded vs down-weighted gene counts, annotated with risk level."""
    counts = pd.Series({
        "excluded": int(getattr(risk, "n_genes_excluded", len(getattr(risk, "excluded_genes", [])))),
        "down_weighted": int(len(getattr(risk, "downweighted_genes", []))),
    })
    fig, ax = new_figure(figsize=(4.5, 4.0))
    ax.bar(np.arange(len(counts)), counts.to_numpy(),
           color=["#C44E52", "#DD8452"])
    ax.set_xticks(np.arange(len(counts)))
    ax.set_xticklabels(counts.index)
    ax.set_ylabel("number of genes")
    ax.set_title(f"Protocol risk: {getattr(risk, 'risk_level', 'unknown')}")
    ax.grid(axis="x", visible=False)

    data = counts.to_frame("n_genes")
    fig_paths, data_paths = save_outputs(
        fig, data, output_dir, name,
        data_comment=[
            f"risk_level: {getattr(risk, 'risk_level', 'unknown')}",
            f"mismatch_types: {getattr(risk, 'mismatch_types', [])}",
        ],
    )
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)


def gene_overlap_plot(
    query_genes: Sequence[str],
    reference_genes: Sequence[str],
    output_dir: Union[str, Path],
    *,
    name: str = "gene_overlap",
) -> PlotResult:
    """Bar of shared / query-only / reference-only gene counts."""
    q = set(query_genes)
    r = set(reference_genes)
    shared = q & r
    counts = pd.Series({
        "shared": len(shared),
        "query_only": len(q - r),
        "reference_only": len(r - q),
    })
    fig, ax = new_figure(figsize=(4.5, 4.0))
    ax.bar(np.arange(len(counts)), counts.to_numpy(),
           color=["#55A868", "#4C72B0", "#8172B3"])
    ax.set_xticks(np.arange(len(counts)))
    ax.set_xticklabels(counts.index, rotation=20)
    ax.set_ylabel("number of genes")
    ax.set_title("Gene overlap: query vs reference")
    ax.grid(axis="x", visible=False)

    fig_paths, data_paths = save_outputs(
        fig, counts.to_frame("n_genes"), output_dir, name,
        data_comment=[f"n_query: {len(q)}", f"n_reference: {len(r)}"],
    )
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)


def warning_summary_plot(
    warnings_list: Sequence[str],
    output_dir: Union[str, Path],
    *,
    name: str = "warning_summary",
) -> PlotResult:
    """Horizontal bar of surfaced-warning counts (never hides warnings)."""
    table = pd.DataFrame({"warning": list(warnings_list)})
    fig, ax = new_figure(figsize=(7.0, max(2.0, 0.4 * max(1, len(table)) + 1)))
    if len(table):
        y = np.arange(len(table))
        ax.barh(y, np.ones(len(table)), color="#CCB974")
        ax.set_yticks(y)
        ax.set_yticklabels([w[:70] + ("…" if len(w) > 70 else "")
                            for w in table["warning"]])
        ax.invert_yaxis()
        ax.set_xticks([])
    else:
        ax.text(0.5, 0.5, "No warnings", ha="center", va="center")
        ax.axis("off")
    ax.set_title(f"Surfaced warnings ({len(table)})")
    ax.grid(visible=False)

    fig_paths, data_paths = save_outputs(fig, table, output_dir, name)
    return PlotResult(figure=fig, figure_paths=fig_paths, data_paths=data_paths)


def _has_cmap(name: str) -> bool:
    try:
        import matplotlib.pyplot as plt

        plt.get_cmap(name)
        return True
    except Exception:
        return False
