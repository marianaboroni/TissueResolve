#!/usr/bin/env python
"""Synthetic spatial benchmark with ground truth (PART 12) — TissueResolve only.

Donor-disjoint, count-level, structured (sharp border + gradient + rare niche).
Builds the reference from reference donors, realises spots from held-out query
donors, runs TissueResolve spatial (flat + hierarchical, spatial-aware NB-CAR)
and an NNLS-per-spot (spatial-naive) baseline, and scores per-spot accuracy +
spatial-fidelity metrics (Moran's I preservation, oversmoothing, local RMSE,
domain recovery).  No core algorithm modified.

Outputs (benchmarks/outputs/spatial_synthetic/)
- spatial_synthetic_metrics.tsv, manifest.json
- raw/<method>_proportions.tsv, raw/truth.tsv, raw/domain_labels.tsv
- figures/ (spatial maps + metric bars + Moran scatter, each with .data.tsv)

Usage:  python benchmarks/spatial/run_synthetic_spatial.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import spatial_metrics as SM  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "spatial_synthetic"
RAW = OUT / "raw"
FIG = OUT / "figures"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
SEED = 0
N_SIDE = 20
CELLS_PER_SPOT = 40
MIN_CELLS = 30


def _score(truth, pred, coords, domains, mapping):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    comp = M.compositional_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
    fid = SM.spatial_fidelity_metrics(t, p, coords, domain_labels=domains, k=6)
    return {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
            "fine_jsd": comp["jsd_mean"], "fine_ccc": comp["ccc"],
            "broad_pearson": M.accuracy_metrics(tfam, pfam)["pearson"],
            "dominant_accuracy": M.dominant_accuracy(t, p), **fid}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--n-side", type=int, default=N_SIDE)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import anndata as ad
    import tissueresolve as tr
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

    RAW.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    ref_donors, query_donors = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=SEED)
    ref_mask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[ref_mask].copy(), min_cells=MIN_CELLS,
                                  estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    mapping = build_cell_type_hierarchy(ref_types, load_hierarchy_mapping(HIERARCHY_TSV))
    rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                      if t in ref_types), ref_types[-1])
    print(f"reference: {len(ref_types)} types; donors ref={len(ref_donors)} query={len(query_donors)}")

    sc = SSP.generate_spatial_scenario(adata, ref_types, query_donors, mapping,
                                       n_side=args.n_side, cells_per_spot=CELLS_PER_SPOT,
                                       seed=SEED, rare_type=rare_type)
    truth = sc["truth"]; coords = sc["coords"]; domains = sc["domain_labels"]
    truth.to_csv(RAW / "truth.tsv", sep="\t")
    pd.DataFrame({"spot": truth.index, "domain": domains,
                  "row": sc["array_row"], "col": sc["array_col"]}).to_csv(
        RAW / "domain_labels.tsv", sep="\t", index=False)
    print(f"scenario: {len(truth)} spots, true mean Moran's I="
          f"{SM.morans_i_per_column(truth, coords, 6).mean():.2f}")

    bulk_df = pd.DataFrame(sc["Y"].T, index=sc["gene_names"], columns=sc["spot_ids"])
    rows = []

    def _run(name, fn):
        t0 = time.perf_counter()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pred = fn()
        except Exception as exc:  # noqa: BLE001
            print(f"  {name:34s} FAILED: {exc}"); return
        rt = time.perf_counter() - t0
        pred.to_csv(RAW / f"{name}_proportions.tsv", sep="\t")
        m = _score(truth, pred, coords, domains, mapping)
        rows.append({"method": name, "runtime_seconds": round(rt, 1), **m})
        print(f"  {name:34s} {rt:5.1f}s fine_r={m['fine_pearson']:.3f} "
              f"local_rmse={m['local_rmse']:.3f} oversmooth={m['oversmoothing_score']:.2f} "
              f"domain_ARI={m.get('domain_recovery_ari', float('nan')):.2f}")

    def _spatial(mode, mp=None):
        res = tr.deconv_spatial(sc["Y"], ref, sc["array_row"], sc["array_col"],
                                sc["lib_sizes"], sc["gene_names"], spot_ids=sc["spot_ids"],
                                resolution_mode=mode, hierarchy_mapping=mp,
                                run_neighbourhood=False)
        return res.deconv.proportions if mode != "hierarchical" else res.estimates.combined_fine

    _run("NNLS_per_spot", lambda: deconv_bulk(bulk_df, ref, solver="nnls",
                                              resolution_mode="none").deconv.proportions)
    _run("TissueResolve_spatial_flat", lambda: _spatial("none"))
    _run("TissueResolve_spatial_hierarchical", lambda: _spatial("hierarchical", mapping))

    metrics = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "spatial_synthetic_metrics.tsv", "w") as fh:
        fh.write("# Synthetic spatial benchmark (ground truth). Accuracy + spatial fidelity. "
                 "TissueResolve spatial vs NNLS-per-spot (spatial-naive) baseline.\n")
        metrics.to_csv(fh, sep="\t", index=False)
    (OUT / "manifest.json").write_text(json.dumps({
        "n_spots": int(len(truth)), "n_side": args.n_side, "cells_per_spot": CELLS_PER_SPOT,
        "n_ref_types": len(ref_types), "reference_donors": ref_donors,
        "query_donors": query_donors, "donor_disjoint": True, "rare_type": rare_type,
        "seed": SEED, "structure": "sharp vertical border + top-bottom gradient + rare niche",
        "ground_truth": "mRNA-proportion per spot",
    }, indent=2), encoding="utf-8")

    _figures(metrics, truth, coords, domains, rare_type)
    print(f"\nWrote -> {OUT}/")
    print(metrics[["method", "fine_pearson", "broad_pearson", "local_rmse",
                   "oversmoothing_score", "domain_recovery_ari"]].round(3).to_string(index=False))
    return 0


def _figures(metrics, truth, coords, domains, rare_type):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def save(fig, name, data, cap):
        fig.tight_layout(); fig.savefig(FIG / f"{name}.png", dpi=150)
        fig.savefig(FIG / f"{name}.svg"); plt.close(fig)
        data.to_csv(FIG / f"{name}.data.tsv", sep="\t")
        (FIG / f"{name}.caption.txt").write_text(cap.strip() + "\n", encoding="utf-8")

    # spatial map: true vs each method for the rare cell type
    methods = [m for m in metrics["method"]]
    panels = [("truth", truth)] + [(m, pd.read_csv(RAW / f"{m}_proportions.tsv",
                                    sep="\t", index_col=0)) for m in methods]
    fig, axes = plt.subplots(1, len(panels), figsize=(4 * len(panels), 3.6))
    for ax, (name, df) in zip(np.atleast_1d(axes), panels):
        v = df.reindex(columns=[rare_type]).fillna(0.0).to_numpy().ravel() \
            if rare_type in df.columns else np.zeros(len(coords))
        scat = ax.scatter(coords[:, 1], -coords[:, 0], c=v, cmap="viridis", s=18)
        ax.set_title(f"{name}\n{rare_type[:18]}", fontsize=8); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(scat, ax=ax, fraction=0.046)
    save(fig, "figH_spatial_map_rare", truth[[rare_type]] if rare_type in truth else truth.iloc[:, :1],
         f"Spatial maps of the rare cell type ({rare_type}) — true vs each method. "
         "Tests whether the rare niche (top-right block) is recovered.")

    # metric bars
    sub = metrics.set_index("method")[["fine_pearson", "broad_pearson", "domain_recovery_ari"]]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    sub.plot.bar(ax=ax); ax.set_ylabel("score"); ax.set_ylim(0, 1.05)
    ax.set_title("Synthetic spatial: accuracy + domain recovery")
    plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
    save(fig, "figI_spatial_metrics", sub,
         "Per-method fine/broad Pearson and domain-recovery ARI on the synthetic "
         "spatial benchmark (ground truth). Spatial-aware vs NNLS-per-spot.")


if __name__ == "__main__":
    raise SystemExit(main())
