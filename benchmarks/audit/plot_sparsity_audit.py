#!/usr/bin/env python
"""Render the 8 PART-1 sparsity-audit diagnostic figures.

Reads the already-computed audit tables (no deconvolution re-run):
- benchmarks/outputs/tcga_prediction_complexity_audit.tsv
- benchmarks/outputs/raw/<method>_proportions.tsv

Every figure is saved as PNG + SVG with a sidecar ``.data.tsv`` (the exact
plotted values) and ``.caption.txt`` — per CLAUDE.md plotting rules ("Never
generate a figure that cannot be reproduced from saved data and metadata").

Usage:  python benchmarks/audit/plot_sparsity_audit.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "benchmarks" / "outputs"
RAW = OUT / "raw"
FIG = OUT / "figures" / "part1_audit"

FLAT_ORDER = ["flat_nnls", "flat_weighted_nnls", "flat_marker_nnls",
              "flat_ridge_nnls", "flat_auto"]
THRESH_COLS = ["n_gt_0", "n_gt_0.0001", "n_gt_0.001", "n_gt_0.005", "n_gt_0.01"]
THRESH_LABELS = [">0", ">1e-4", ">1e-3", ">5e-3", ">1e-2"]


def _save(fig, name, data: pd.DataFrame, caption: str):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=150)
    fig.savefig(FIG / f"{name}.svg")
    plt.close(fig)
    data.to_csv(FIG / f"{name}.data.tsv", sep="\t")
    (FIG / f"{name}.caption.txt").write_text(caption.strip() + "\n", encoding="utf-8")
    print(f"  wrote {name} (.png/.svg/.data.tsv/.caption.txt)")


def main() -> int:
    audit = pd.read_csv(OUT / "tcga_prediction_complexity_audit.tsv", sep="\t", comment="#")

    # 1. detected populations per sample across abundance thresholds (flat solvers)
    rows = []
    for m in FLAT_ORDER:
        sub = audit[audit.method == m]
        rows.append({"method": m, **{lbl: sub[c].mean()
                     for c, lbl in zip(THRESH_COLS, THRESH_LABELS)}})
    d1 = pd.DataFrame(rows).set_index("method")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(THRESH_LABELS))
    for m in d1.index:
        ax.plot(x, d1.loc[m], marker="o", label=m)
    ax.set_xticks(x); ax.set_xticklabels(THRESH_LABELS)
    ax.set_xlabel("abundance threshold"); ax.set_ylabel("mean # populations > threshold")
    ax.set_title("1. Detected populations per sample vs abundance threshold")
    ax.legend(fontsize=8)
    _save(fig, "fig1_detected_vs_threshold", d1,
          "Mean number of TCGA populations above each abundance threshold, by flat "
          "solver. nnls/weighted are sparse (~3-4); ridge/auto dense (~20). Source: "
          "tcga_prediction_complexity_audit.tsv.")

    # 2. raw (flat solver) vs post-processed (hierarchical combined) population count
    flat = audit[audit.method == "flat_auto"].set_index("sample")["n_gt_0"]
    hier = audit[audit.method == "hierarchical_hier_combined"].set_index("sample")["n_gt_0"]
    nnls = audit[audit.method == "flat_nnls"].set_index("sample")["n_gt_0"]
    d2 = pd.DataFrame({"flat_auto(ridge)": flat, "flat_nnls": nnls,
                       "hierarchical_combined": hier}).dropna()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.boxplot([d2["flat_auto(ridge)"], d2["flat_nnls"], d2["hierarchical_combined"]],
               tick_labels=["flat auto\n(ridge)", "flat\nnnls", "hierarchical\ncombined"])
    ax.set_ylabel("# non-zero populations per sample")
    ax.set_title("2. Population count: dense solver vs sparse solver vs gated hierarchy")
    _save(fig, "fig2_raw_vs_postprocessed", d2,
          "Per-sample non-zero population counts. The flat pipeline applies no "
          "data-level zeroing (ridge keeps ~23); sparsity is the nnls solver and the "
          "hierarchical within-family gating, not post-processing.")

    # 3. composition entropy by sample (flat solvers)
    piv = audit[audit.method.isin(FLAT_ORDER)].pivot_table(
        index="sample", columns="method", values="shannon_entropy_nats")
    d3 = piv[[m for m in FLAT_ORDER if m in piv.columns]]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.boxplot([d3[m].dropna() for m in d3.columns], tick_labels=[m.replace("flat_", "") for m in d3.columns])
    ax.set_ylabel("Shannon entropy (nats)")
    ax.set_title("3. Composition entropy by sample (flat solvers)")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    _save(fig, "fig3_entropy_by_sample", d3,
          "Per-sample Shannon entropy of the predicted composition. Higher = more "
          "diffuse. nnls collapses entropy; ridge/auto retain it.")

    # 4. dominant fraction & top-5 cumulative (mean per method)
    d4 = audit.groupby("method")[["dominant_fraction", "top5_cumulative"]].mean()
    d4 = d4.reindex([m for m in FLAT_ORDER + ["hierarchical_hier_combined"] if m in d4.index])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(d4)); w = 0.38
    ax.bar(x - w / 2, d4["dominant_fraction"], w, label="dominant fraction")
    ax.bar(x + w / 2, d4["top5_cumulative"], w, label="top-5 cumulative")
    ax.set_xticks(x); ax.set_xticklabels([m.replace("flat_", "").replace("hierarchical_", "")
                                          for m in d4.index], rotation=20, ha="right")
    ax.set_ylabel("fraction of mass"); ax.set_ylim(0, 1.05)
    ax.set_title("4. Dominant fraction and top-5 cumulative abundance")
    ax.legend()
    _save(fig, "fig4_dominant_top5", d4,
          "Mean dominant-population fraction and top-5 cumulative mass per method. "
          "Sparse methods concentrate ~100% in the top 5.")

    # 5. fraction removed by thresholding / aggregation
    rows = []
    for m in FLAT_ORDER + ["hierarchical_hier_combined"]:
        sub = audit[audit.method == m]
        # mass below 1% threshold = 1 - sum of mass in pops>1% (approx via mass_sum & top counts)
        # we report unresolved mass (aggregation) and mass not in top-5 (long tail)
        removed_agg = sub["unresolved_mass"].mean() if "unresolved_mass" in sub else 0.0
        tail = 1.0 - sub["top5_cumulative"].mean()
        rows.append({"method": m, "moved_to_unresolved": removed_agg,
                     "mass_outside_top5": tail})
    d5 = pd.DataFrame(rows).set_index("method")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(d5)); w = 0.38
    ax.bar(x - w / 2, d5["moved_to_unresolved"], w, label="moved to unresolved (aggregation)")
    ax.bar(x + w / 2, d5["mass_outside_top5"], w, label="mass outside top-5 (tail)")
    ax.set_xticks(x); ax.set_xticklabels([m.replace("flat_", "").replace("hierarchical_", "")
                                          for m in d5.index], rotation=20, ha="right")
    ax.set_ylabel("fraction of mass")
    ax.set_title("5. Mass moved by aggregation vs retained in long tail")
    ax.legend(fontsize=8)
    _save(fig, "fig5_fraction_removed", d5,
          "Flat solvers move ~0 mass to aggregation (no 'Other'/threshold zeroing); "
          "only the hierarchical path parks mass in unresolved_* (~0.27). The dense "
          "solvers hold substantial mass outside the top-5 (a real long tail).")

    # 6. broad vs fine detected populations (hierarchical)
    fam = audit[audit.method == "hierarchical_hier_family"].set_index("sample")["n_gt_0"]
    fine = audit[audit.method == "hierarchical_hier_fine"].set_index("sample")["n_gt_0"]
    d6 = pd.DataFrame({"broad_families_detected": fam, "fine_subtypes_detected": fine}).dropna()
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(d6["broad_families_detected"], d6["fine_subtypes_detected"], alpha=0.7)
    lim = max(d6.max().max(), 8) + 1
    ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.5)
    ax.set_xlabel("# broad families detected (of 8)")
    ax.set_ylabel("# fine subtypes detected (of 32)")
    ax.set_title("6. Broad vs fine detected populations (hierarchical)")
    _save(fig, "fig6_broad_vs_fine", d6,
          "Hierarchical mode resolves ~3 broad families but collapses fine subtypes "
          "to ~3 (of 32) because the within-family gate marks all multi-member "
          "families unresolved.")

    # 7. population prevalence across samples (flat_auto = shipped, threshold 1e-3)
    props = pd.read_csv(RAW / "flat_auto_proportions.tsv", sep="\t", index_col=0)
    prevalence = (props > 1e-3).mean(axis=0).sort_values(ascending=False)
    d7 = prevalence.to_frame("prevalence_gt_1e-3")
    fig, ax = plt.subplots(figsize=(7, 8))
    ax.barh(range(len(prevalence)), prevalence.values[::-1])
    ax.set_yticks(range(len(prevalence)))
    ax.set_yticklabels(prevalence.index[::-1], fontsize=6)
    ax.set_xlabel("fraction of samples with proportion > 1e-3")
    ax.set_title("7. Population prevalence across TCGA samples (auto/ridge)")
    _save(fig, "fig7_population_prevalence", d7,
          "Fraction of the 40 TCGA samples in which each population exceeds 1e-3, "
          "under the shipped auto(ridge) output.")

    # 8. rare-population detection frequency: mean abundance vs detection rate
    mean_ab = props.mean(axis=0)
    det_rate = (props > 1e-3).mean(axis=0)
    d8 = pd.DataFrame({"mean_abundance": mean_ab, "detection_rate_gt_1e-3": det_rate})
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.scatter(d8["mean_abundance"], d8["detection_rate_gt_1e-3"], alpha=0.7)
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set_xlabel("mean predicted abundance (symlog)")
    ax.set_ylabel("detection rate (fraction of samples > 1e-3)")
    ax.set_title("8. Rare-population detection frequency (auto/ridge)")
    _save(fig, "fig8_rare_detection", d8,
          "Detection rate vs mean predicted abundance per population. Low-abundance "
          "populations are detected in fewer samples — the empirical detection floor "
          "for the dense solver. NB: no ground truth, so this is descriptive only.")

    print(f"\nAll 8 figures -> {FIG}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
