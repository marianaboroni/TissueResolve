#!/usr/bin/env python
"""Figures for the Phase-2 held-out-donor bulk benchmark (reads the metric TSVs).

Each figure saves PNG + SVG + .data.tsv + .caption.txt (CLAUDE.md plotting rules).
Usage:  python benchmarks/bulk/plot_holdout_results.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "benchmarks" / "outputs" / "holdout_bulk"
FIG = OUT / "figures"
METHOD_ORDER = ["flat_nnls", "flat_weighted_nnls", "flat_ridge_nnls", "flat_auto", "hierarchical"]


def _save(fig, name, data, caption):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=150)
    fig.savefig(FIG / f"{name}.svg")
    plt.close(fig)
    data.to_csv(FIG / f"{name}.data.tsv", sep="\t")
    (FIG / f"{name}.caption.txt").write_text(caption.strip() + "\n", encoding="utf-8")
    print(f"  wrote {name}")


def main():
    fine = pd.read_csv(OUT / "bulk_fine_metrics.tsv", sep="\t", comment="#")
    broad = pd.read_csv(OUT / "bulk_broad_metrics.tsv", sep="\t", comment="#")
    cx = pd.read_csv(OUT / "bulk_complexity_metrics.tsv", sep="\t", comment="#")

    # A. predicted vs true effective-N (richness preservation) — PART-16 fig E
    fig, ax = plt.subplots(figsize=(6.5, 6))
    markers = {"flat_nnls": "o", "flat_weighted_nnls": "s", "flat_ridge_nnls": "^",
               "flat_auto": "D", "hierarchical": "x"}
    for m in METHOD_ORDER:
        sub = cx[cx.method == m]
        ax.scatter(sub["true_effective_n"], sub["pred_effective_n"],
                   marker=markers.get(m, "o"), s=70, label=m, alpha=0.8)
    lim = max(cx["true_effective_n"].max(), cx["pred_effective_n"].max()) + 2
    ax.plot([0, lim], [0, lim], "k--", lw=1, alpha=0.6, label="perfect (y=x)")
    ax.set_xlabel("true effective # populations  exp(H)")
    ax.set_ylabel("predicted effective # populations")
    ax.set_title("Composition-complexity preservation (bulk, held-out donors)")
    ax.legend(fontsize=8)
    _save(fig, "figA_richness_preservation", cx.set_index(["scenario", "method"]),
          "Predicted vs true effective number of populations per scenario. Points on "
          "y=x preserve complexity. ridge tracks the diffuse truth (sometimes "
          "over-disperses); nnls under-detects rich mixtures; hierarchical collapses "
          "to ~2-3 regardless of truth (over-aggressive within-family gating).")

    # B. fine & broad Pearson by method (mean over scenarios)
    fp = fine.groupby("method")["pearson"].mean().reindex(METHOD_ORDER)
    bp = broad.groupby("method")["pearson"].mean().reindex(METHOD_ORDER)
    d = pd.DataFrame({"fine_pearson": fp, "broad_pearson": bp})
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    x = np.arange(len(d)); w = 0.38
    ax.bar(x - w / 2, d["fine_pearson"], w, label="fine-level")
    ax.bar(x + w / 2, d["broad_pearson"], w, label="broad-family")
    ax.set_xticks(x); ax.set_xticklabels(d.index, rotation=20, ha="right")
    ax.set_ylabel("mean Pearson r vs mRNA-proportion truth"); ax.set_ylim(0, 1.05)
    ax.set_title("Accuracy by method (mean over 5 scenarios)")
    ax.legend()
    _save(fig, "figB_accuracy_by_method", d,
          "Mean Pearson correlation to the mRNA-proportion ground truth, fine-level "
          "and broad-family. Broad recovery is strong for all flat solvers; fine is "
          "led by nnls. Hierarchical trails at both levels (note: fine strips "
          "abstained/unresolved columns, but broad is fair and also lower).")

    # C. complexity-preservation error heatmap (scenario × method)
    piv = cx.pivot(index="scenario", columns="method", values="effective_n_abs_error")
    piv = piv.reindex(columns=[m for m in METHOD_ORDER if m in piv.columns])
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.imshow(piv.to_numpy(), aspect="auto", cmap="magma_r")
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=20, ha="right")
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.iloc[i, j]:.1f}", ha="center", va="center", fontsize=8,
                    color="white" if piv.iloc[i, j] > piv.to_numpy().mean() else "black")
    fig.colorbar(im, ax=ax, label="|pred − true| effective-N (lower = better)")
    ax.set_title("Complexity-preservation error by scenario × method")
    _save(fig, "figC_complexity_error_heatmap", piv,
          "Absolute error in effective number of populations. Lower is better. "
          "hierarchical has the largest errors everywhere (collapse); ridge errs by "
          "over-dispersion on concentrated mixtures (imbalanced/similar).")

    print(f"\nFigures -> {FIG}/")


if __name__ == "__main__":
    main()
