#!/usr/bin/env python
"""Score multi-split external predictions; compute CIs + paired Wilcoxon vs TissueResolve.

Reads predictions from run_external_multisplit.R, scores each (split, scenario)
against its truth with the SAME metrics as the TissueResolve multi-split run,
merges into one long table, and reports bootstrap 95% CIs and paired Wilcoxon
(external methods vs the best TissueResolve flat solver) over the 25 paired
(split × scenario) replicates.

Usage:  python benchmarks/bulk/score_external_multisplit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402

HB = REPO / "benchmarks" / "outputs" / "holdout_bulk"
EXT = HB / "external_inputs_multisplit"
MULTI = REPO / "benchmarks" / "outputs" / "holdout_bulk_multisplit"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
SCENARIOS = ["balanced", "imbalanced", "rare", "similar_subtypes", "missing_population"]
EXT_METHODS = ["MuSiC", "BisqueRNA"]
COMPARE_TO = "flat_nnls"   # best TissueResolve flat solver from the multi-split run


def main() -> int:
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)

    rows = []
    for sdir in sorted(EXT.glob("split*")):
        s = int(sdir.name.replace("split", ""))
        for scen in SCENARIOS:
            truth_f = sdir / scen / "truth_mrna.tsv"
            if not truth_f.exists():
                continue
            truth = pd.read_csv(truth_f, sep="\t", index_col=0)
            cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
            t = truth[cols]
            mapping = build_cell_type_hierarchy(list(t.columns), raw_map)
            for method in EXT_METHODS:
                pf = sdir / scen / f"{method}_pred.tsv"
                if not pf.exists():
                    continue
                pred = pd.read_csv(pf, sep="\t", index_col=0)
                p = pred.reindex(index=t.index, columns=t.columns).fillna(0.0)
                acc = M.accuracy_metrics(t, p)
                comp = M.compositional_metrics(t, p)
                tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
                cxt = M.composition_complexity(t)["effective_n_populations"].mean()
                cxp = M.composition_complexity(p)["effective_n_populations"].mean()
                vals = {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
                        "fine_jsd": comp["jsd_mean"], "fine_ccc": comp["ccc"],
                        "broad_pearson": M.accuracy_metrics(tfam, pfam)["pearson"],
                        "complexity_abs_error": abs(cxp - cxt)}
                for metric, val in vals.items():
                    rows.append({"split": s, "scenario": scen, "method": method,
                                 "metric": metric, "value": float(val)})
    ext_long = pd.DataFrame(rows)
    if ext_long.empty:
        print("No external predictions found — run run_external_multisplit.R first.")
        return 1

    # merge with TissueResolve multi-split long-form (shared metric names)
    tr_long = pd.read_csv(MULTI / "metrics_long.tsv", sep="\t", comment="#")
    shared_metrics = ["fine_pearson", "fine_rmse", "fine_jsd", "fine_ccc",
                      "broad_pearson", "complexity_abs_error"]
    tr_long = tr_long[tr_long.metric.isin(shared_metrics)]
    combined = pd.concat([tr_long, ext_long], ignore_index=True)
    HB.mkdir(parents=True, exist_ok=True)
    out = HB / "multisplit_combined"
    out.mkdir(exist_ok=True)
    combined.to_csv(out / "metrics_long.tsv", sep="\t", index=False)

    # bootstrap CIs per (method, metric)
    ci = []
    for (method, metric), g in combined.groupby(["method", "metric"]):
        ci.append({"method": method, "metric": metric, **ST.bootstrap_ci(g["value"].to_numpy(), seed=0)})
    ci = pd.DataFrame(ci)
    ci.to_csv(out / "metrics_ci.tsv", sep="\t", index=False)

    # paired Wilcoxon: each method vs COMPARE_TO on fine_pearson (aligned by split×scenario)
    wil = []
    for metric in ["fine_pearson", "complexity_abs_error"]:
        piv = combined[combined.metric == metric].pivot_table(
            index=["split", "scenario"], columns="method", values="value")
        if COMPARE_TO not in piv.columns:
            continue
        for m in piv.columns:
            if m == COMPARE_TO:
                continue
            r = ST.paired_wilcoxon(piv[COMPARE_TO].to_numpy(), piv[m].to_numpy())
            wil.append({"metric": metric, "reference": COMPARE_TO, "other": m, **r})
    pd.DataFrame(wil).to_csv(out / "pairwise_wilcoxon_vs_nnls.tsv", sep="\t", index=False)

    # report
    print("=== fine Pearson, mean [95% CI] (TissueResolve + external, n=25) ===")
    fp = ci[ci.metric == "fine_pearson"].sort_values("point", ascending=False)
    for _, r in fp.iterrows():
        print(f"  {r['method']:24s} {r['point']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}]  n={int(r['n'])}")
    print(f"\n=== paired Wilcoxon vs {COMPARE_TO} (fine Pearson) ===")
    w = pd.DataFrame(wil)
    wf = w[w.metric == "fine_pearson"]
    for _, r in wf.iterrows():
        sig = "ns" if (not np.isfinite(r["p_value"]) or r["p_value"] >= 0.05) else "*"
        print(f"  {COMPARE_TO} vs {r['other']:22s} p={r['p_value']:.4f} effect={r['effect_size']:+.2f} {sig}")

    # figure: fine-Pearson forest, external highlighted
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fp2 = fp.sort_values("point")
    ext_set = set(EXT_METHODS)
    colors = ["#b5562a" if m in ext_set else "#2c6fbb" for m in fp2["method"]]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    y = np.arange(len(fp2))
    ax.errorbar(fp2["point"], y, xerr=[fp2["point"] - fp2["lo"], fp2["hi"] - fp2["point"]],
                fmt="o", capsize=4, ecolor="grey", linestyle="none")
    ax.scatter(fp2["point"], y, c=colors, s=60, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(fp2["method"])
    ax.set_xlabel("fine Pearson r vs mRNA truth (mean, 95% bootstrap CI, n=25)")
    ax.set_title("Bulk fine accuracy with CIs: TissueResolve vs MuSiC/BisqueRNA")
    import matplotlib.patches as mp
    ax.legend(handles=[mp.Patch(color="#2c6fbb", label="TissueResolve"),
                       mp.Patch(color="#b5562a", label="external")], fontsize=8)
    fig.tight_layout()
    figd = out / "figures"; figd.mkdir(exist_ok=True)
    fig.savefig(figd / "figG_multisplit_external_ci.png", dpi=150)
    fig.savefig(figd / "figG_multisplit_external_ci.svg")
    plt.close(fig)
    fp2.to_csv(figd / "figG_multisplit_external_ci.data.tsv", sep="\t", index=False)
    (figd / "figG_multisplit_external_ci.caption.txt").write_text(
        "Fine-level Pearson with bootstrap 95% CIs over 25 split×scenario "
        "replicates. External MuSiC/BisqueRNA executed on the same held-out "
        "inputs as TissueResolve. hierarchical_resolvable is a resolved-subset "
        "score (not full-panel comparable).\n", encoding="utf-8")
    print(f"\nWrote -> {out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
