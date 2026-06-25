#!/usr/bin/env python
"""Score external (MuSiC, BisqueRNA) held-out predictions and merge with TissueResolve.

Reads predictions written by run_external_holdout.R, scores them against the SAME
mRNA-proportion truth with the SAME metrics used for TissueResolve, and produces a
combined comparison table + status table + figure.  Answers PART 17: do external
methods also concentrate composition, or is TissueResolve uniquely sparse?

Usage:  python benchmarks/bulk/score_external_holdout.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from benchmarks.shared import metrics as M  # noqa: E402

HB = REPO / "benchmarks" / "outputs" / "holdout_bulk"
EXT = HB / "external_inputs"
SCENARIOS = ["balanced", "imbalanced", "rare", "similar_subtypes", "missing_population"]
EXT_METHODS = ["MuSiC", "BisqueRNA"]


def _score(truth: pd.DataFrame, pred: pd.DataFrame) -> dict:
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    comp = M.compositional_metrics(t, p)
    cxt = M.composition_complexity(t)["effective_n_populations"].mean()
    cxp = M.composition_complexity(p)["effective_n_populations"].mean()
    return {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
            "fine_jsd": comp["jsd_mean"], "fine_ccc": comp["ccc"],
            "pred_effective_n": float(cxp), "true_effective_n": float(cxt),
            "complexity_abs_error": float(abs(cxp - cxt))}


def main() -> int:
    rows = []
    for scen in SCENARIOS:
        truth = pd.read_csv(EXT / scen / "truth_mrna.tsv", sep="\t", index_col=0)
        for method in EXT_METHODS:
            pf = EXT / scen / f"{method}_pred.tsv"
            if not pf.exists():
                continue
            pred = pd.read_csv(pf, sep="\t", index_col=0)
            rows.append({"scenario": scen, "method": method, **_score(truth, pred)})
    ext = pd.DataFrame(rows)
    ext.to_csv(HB / "external_fine_metrics.tsv", sep="\t", index=False)

    # merge with TissueResolve single-split fine metrics (same scenarios/truth)
    tr = pd.read_csv(HB / "bulk_fine_metrics.tsv", sep="\t", comment="#")
    tr_keep = tr[["scenario", "method", "pearson", "rmse", "jsd", "ccc"]].rename(
        columns={"pearson": "fine_pearson", "rmse": "fine_rmse",
                 "jsd": "fine_jsd", "ccc": "fine_ccc"})
    tr_keep["family"] = "TissueResolve"
    ext_keep = ext[["scenario", "method", "fine_pearson", "fine_rmse",
                    "fine_jsd", "fine_ccc", "pred_effective_n", "true_effective_n"]].copy()
    ext_keep["family"] = "external"
    combined = pd.concat([tr_keep, ext_keep], ignore_index=True, sort=False)
    combined.to_csv(HB / "combined_fine_metrics.tsv", sep="\t", index=False)

    # mean over scenarios per method
    summ = combined.groupby("method").agg(
        family=("family", "first"),
        mean_fine_pearson=("fine_pearson", "mean"),
        mean_fine_rmse=("fine_rmse", "mean"),
        mean_fine_jsd=("fine_jsd", "mean"),
        mean_fine_ccc=("fine_ccc", "mean")).sort_values("mean_fine_pearson", ascending=False)
    summ.to_csv(HB / "combined_method_summary.tsv", sep="\t")

    # complexity comparison (external vs true)
    print("=== external method complexity (pred vs true effective-N) ===")
    print(ext.groupby("method")[["pred_effective_n", "true_effective_n",
                                 "complexity_abs_error"]].mean().round(2).to_string())
    print("\n=== combined mean fine accuracy (TissueResolve + external) ===")
    print(summ.round(3).to_string())

    # figure: mean fine Pearson by method, colored by family
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors = {"TissueResolve": "#2c6fbb", "external": "#b5562a"}
    s = summ.sort_values("mean_fine_pearson")
    ax.barh(range(len(s)), s["mean_fine_pearson"],
            color=[colors[f] for f in s["family"]])
    ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index)
    ax.set_xlabel("mean fine Pearson r vs mRNA truth (5 scenarios)")
    ax.set_title("Bulk fine accuracy: TissueResolve vs external (MuSiC, BisqueRNA)")
    import matplotlib.patches as mp
    ax.legend(handles=[mp.Patch(color=c, label=l) for l, c in colors.items()], fontsize=8)
    fig.tight_layout()
    FIGD = HB / "figures"; FIGD.mkdir(exist_ok=True)
    fig.savefig(FIGD / "figF_external_vs_tissueresolve.png", dpi=150)
    fig.savefig(FIGD / "figF_external_vs_tissueresolve.svg")
    plt.close(fig)
    s.to_csv(FIGD / "figF_external_vs_tissueresolve.data.tsv", sep="\t")
    (FIGD / "figF_external_vs_tissueresolve.caption.txt").write_text(
        "Mean fine-level Pearson to mRNA-proportion truth over 5 held-out-donor "
        "scenarios. External methods (MuSiC v1.0.0, BisqueRNA v1.0.5) executed on "
        "the SAME reference + bulk as TissueResolve. Addresses whether composition "
        "concentration is unique to TissueResolve or shared across methods.\n",
        encoding="utf-8")
    print(f"\nWrote combined tables + figF -> {HB}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
