#!/usr/bin/env python
"""Score external spatial predictions (RCTD, cell2location) + merge with TissueResolve.

Scores whatever external methods produced predictions against the SAME synthetic
spatial truth, with the SAME accuracy + spatial-fidelity metrics, and merges with
the TissueResolve spatial results.  Honest status table for all external methods
(executed / failed / skipped), never fabricated.

Usage:  python benchmarks/spatial/score_external_synthetic.py
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

SS = REPO / "benchmarks" / "outputs" / "spatial_synthetic"
EXT = SS / "external_inputs"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
EXTERNAL = ["RCTD", "cell2location", "CARD"]
NOT_INSTALLABLE = {"SPOTlight": "not installed (Bioconductor; Seurat absent)",
                   "Tangram": "not installed",
                   "DestVI": "not installed"}


def main() -> int:
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    truth = pd.read_csv(EXT / "spatial" / "truth_mrna.tsv", sep="\t", index_col=0)
    coords = pd.read_csv(EXT / "spatial" / "spot_coords.tsv", sep="\t")[["row", "col"]].to_numpy(float)
    domains = pd.read_csv(EXT / "spatial" / "domain_labels.tsv", sep="\t")["domain"].to_numpy()
    cols = list(truth.columns)
    mapping = build_cell_type_hierarchy(cols, load_hierarchy_mapping(HIERARCHY_TSV))

    rows, status = [], []
    for m in EXTERNAL:
        sf = EXT / f"{m}_status.tsv"
        st = pd.read_csv(sf, sep="\t").iloc[0].to_dict() if sf.exists() else \
            {"method": m, "status": "not_run", "version": "", "runtime_seconds": np.nan}
        status.append(st)
        pf = EXT / f"{m}_pred.tsv"
        if st.get("status") != "executed" or not pf.exists():
            print(f"{m}: {st.get('status')} — not scored")
            continue
        pred = pd.read_csv(pf, sep="\t", index_col=0)
        pred.index = [str(i) for i in pred.index]
        t = truth.loc[[i for i in truth.index if i in set(pred.index)]]
        p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
        acc = M.accuracy_metrics(t, p); comp = M.compositional_metrics(t, p)
        tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
        fid = SM.spatial_fidelity_metrics(t, p, coords[:len(t)], domain_labels=domains[:len(t)], k=6)
        rows.append({"method": m, "family": "external",
                     "runtime_seconds": st.get("runtime_seconds"),
                     "fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
                     "fine_jsd": comp["jsd_mean"], "fine_ccc": comp["ccc"],
                     "broad_pearson": M.accuracy_metrics(tfam, pfam)["pearson"],
                     "dominant_accuracy": M.dominant_accuracy(t, p), **fid})
        print(f"{m}: scored fine_r={acc['pearson']:.3f} oversmooth={fid['oversmoothing_score']:.2f}")

    for m, why in NOT_INSTALLABLE.items():
        status.append({"method": m, "status": "failed" if "load" in why else "skipped",
                       "version": "", "runtime_seconds": np.nan, "error": why})

    # merge with TissueResolve spatial metrics
    tr = pd.read_csv(SS / "spatial_synthetic_metrics.tsv", sep="\t", comment="#")
    tr["family"] = "TissueResolve"
    keep = ["method", "family", "fine_pearson", "broad_pearson", "dominant_accuracy",
            "local_rmse", "oversmoothing_score", "domain_recovery_ari"]
    combined = pd.concat([tr, pd.DataFrame(rows)], ignore_index=True, sort=False)
    combined = combined[[c for c in keep if c in combined.columns]]
    combined.to_csv(SS / "combined_spatial_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(status).to_csv(SS / "external_spatial_status.tsv", sep="\t", index=False)

    print("\n=== combined spatial metrics ===")
    print(combined.round(3).to_string(index=False))
    print("\n=== external status ===")
    print(pd.DataFrame(status)[["method", "status", "version"]].to_string(index=False))

    # figure
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    c = combined.dropna(subset=["fine_pearson"]).sort_values("fine_pearson")
    colors = {"TissueResolve": "#2c6fbb", "external": "#b5562a"}
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.barh(range(len(c)), c["fine_pearson"],
            color=[colors.get(f, "#888") for f in c["family"]])
    ax.set_yticks(range(len(c))); ax.set_yticklabels(c["method"])
    ax.set_xlabel("fine Pearson r vs spot mRNA truth")
    ax.set_title("Synthetic spatial: TissueResolve vs external (RCTD, cell2location)")
    import matplotlib.patches as mp
    ax.legend(handles=[mp.Patch(color=v, label=k) for k, v in colors.items()], fontsize=8)
    fig.tight_layout()
    figd = SS / "figures"; figd.mkdir(exist_ok=True)
    fig.savefig(figd / "figJ_spatial_external.png", dpi=150)
    fig.savefig(figd / "figJ_spatial_external.svg"); plt.close(fig)
    c.to_csv(figd / "figJ_spatial_external.data.tsv", sep="\t", index=False)
    (figd / "figJ_spatial_external.caption.txt").write_text(
        "Fine-level Pearson on the synthetic spatial benchmark — TissueResolve "
        "vs external methods that executed (RCTD, cell2location). CARD failed to "
        "build; SPOTlight not attempted (see external_spatial_status.tsv).\n",
        encoding="utf-8")
    print(f"\nWrote -> {SS}/combined_spatial_metrics.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
