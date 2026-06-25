#!/usr/bin/env python
"""Score multi-seed spatial (internal + external) → bootstrap CIs + paired Wilcoxon.

Combines internal_metrics_long.tsv (TissueResolve, NNLS) with per-seed RCTD /
cell2location predictions scored against each seed's truth, then reports
bootstrap 95% CIs and paired Wilcoxon (vs TissueResolve_spatial_flat) over seeds.

Usage:  python benchmarks/spatial/score_spatial_multiseed.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import spatial_metrics as SM  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402

BASE = REPO / "benchmarks" / "outputs" / "spatial_multiseed"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
EXTERNAL = ["RCTD", "cell2location"]
REF_METHOD = "TissueResolve_spatial_flat"


def main() -> int:
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)

    long = pd.read_csv(BASE / "internal_metrics_long.tsv", sep="\t")
    ext_rows = []
    for sdir in sorted(BASE.glob("seed*")):
        s = int(sdir.name.replace("seed", ""))
        truth = pd.read_csv(sdir / "truth.tsv", sep="\t", index_col=0)
        dl = pd.read_csv(sdir / "domain_labels.tsv", sep="\t")
        coords = dl[["row", "col"]].to_numpy(float); domains = dl["domain"].to_numpy()
        cols = list(truth.columns)
        mapping = build_cell_type_hierarchy(cols, raw_map)
        for m in EXTERNAL:
            pf = sdir / f"{m}_pred.tsv"
            if not pf.exists():
                continue
            pred = pd.read_csv(pf, sep="\t", index_col=0); pred.index = [str(i) for i in pred.index]
            t = truth.loc[[i for i in truth.index if i in set(pred.index)]]
            p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
            acc = M.accuracy_metrics(t, p); comp = M.compositional_metrics(t, p)
            tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
            fid = SM.spatial_fidelity_metrics(t, p, coords[:len(t)], domain_labels=domains[:len(t)], k=6)
            vals = {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
                    "fine_jsd": comp["jsd_mean"], "broad_pearson": M.accuracy_metrics(tfam, pfam)["pearson"],
                    "dominant_accuracy": M.dominant_accuracy(t, p), "local_rmse": fid["local_rmse"],
                    "oversmoothing_score": fid["oversmoothing_score"],
                    "domain_recovery_ari": fid.get("domain_recovery_ari", float("nan"))}
            for metric, val in vals.items():
                ext_rows.append({"seed": s, "method": m, "metric": metric, "value": float(val)})

    combined = pd.concat([long, pd.DataFrame(ext_rows)], ignore_index=True)
    combined.to_csv(BASE / "combined_metrics_long.tsv", sep="\t", index=False)

    ci = []
    for (method, metric), g in combined.groupby(["method", "metric"]):
        ci.append({"method": method, "metric": metric, **ST.bootstrap_ci(g["value"].to_numpy(), seed=0)})
    ci = pd.DataFrame(ci); ci.to_csv(BASE / "metrics_ci.tsv", sep="\t", index=False)

    wil = []
    for metric in ["fine_pearson", "broad_pearson", "oversmoothing_score"]:
        piv = combined[combined.metric == metric].pivot_table(index="seed", columns="method", values="value")
        if REF_METHOD not in piv.columns:
            continue
        for m in piv.columns:
            if m == REF_METHOD:
                continue
            wil.append({"metric": metric, "reference": REF_METHOD, "other": m,
                        **ST.paired_wilcoxon(piv[REF_METHOD].to_numpy(), piv[m].to_numpy())})
    pd.DataFrame(wil).to_csv(BASE / "pairwise_wilcoxon.tsv", sep="\t", index=False)

    print("=== fine Pearson, mean [95% CI] over seeds ===")
    fp = ci[ci.metric == "fine_pearson"].sort_values("point", ascending=False)
    for _, r in fp.iterrows():
        print(f"  {r['method']:36s} {r['point']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}]  n={int(r['n'])}")
    print(f"\n=== paired Wilcoxon vs {REF_METHOD} (fine Pearson) ===")
    for _, r in pd.DataFrame(wil)[lambda d: d.metric == "fine_pearson"].iterrows():
        print(f"  vs {r['other']:34s} p={r['p_value']:.4f} effect={r['effect_size']:+.2f} n={int(r['n'])}")

    # forest figure
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ext_set = set(EXTERNAL)
    fp2 = fp.sort_values("point")
    colors = ["#b5562a" if m in ext_set else "#2c6fbb" for m in fp2["method"]]
    fig, ax = plt.subplots(figsize=(8.5, 4.5)); y = np.arange(len(fp2))
    ax.errorbar(fp2["point"], y, xerr=[fp2["point"]-fp2["lo"], fp2["hi"]-fp2["point"]],
                fmt="o", capsize=4, ecolor="grey", linestyle="none")
    ax.scatter(fp2["point"], y, c=colors, s=60, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(fp2["method"]); ax.set_xlabel("fine Pearson (mean, 95% CI)")
    ax.set_title(f"Spatial fine accuracy across {int(fp['n'].max())} seeds (144-spot grids)")
    import matplotlib.patches as mp
    ax.legend(handles=[mp.Patch(color="#2c6fbb", label="TissueResolve/baseline"),
                       mp.Patch(color="#b5562a", label="external")], fontsize=8)
    fig.tight_layout()
    figd = BASE / "figures"; figd.mkdir(exist_ok=True)
    fig.savefig(figd / "figK_spatial_multiseed_ci.png", dpi=150)
    fig.savefig(figd / "figK_spatial_multiseed_ci.svg"); plt.close(fig)
    fp2.to_csv(figd / "figK_spatial_multiseed_ci.data.tsv", sep="\t", index=False)
    (figd / "figK_spatial_multiseed_ci.caption.txt").write_text(
        "Fine-level Pearson with bootstrap 95% CIs across seeds (each seed = a new "
        "donor split + spatial layout, 144-spot grid). External RCTD/cell2location "
        "on the same per-seed inputs. n = number of seeds (small).\n", encoding="utf-8")
    print(f"\nWrote CIs + Wilcoxon + figK -> {BASE}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
