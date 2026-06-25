#!/usr/bin/env python
"""TASK D — level-specific spatial smoothing benchmark (experimental).

10 seeds × 6 modes on structured synthetic spatial (sharp border + gradient + rare
niche), donor-disjoint reference. Modes vary (λ_broad, λ_fine) via the experimental
level-specific wrapper (existing solver, no NB-CAR VI). Tests whether λ_broad>λ_fine
reduces fine over-smoothing/edge blurring without hurting broad accuracy.

Outputs: benchmarks/outputs/level_specific_spatial_smoothing.tsv,
spatial_boundary_metrics.tsv, taskD_promotion_gates.tsv; figure.
Usage:  python benchmarks/diagnostics/phase_taskD_level_spatial.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
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
from benchmarks.shared import stats as ST  # noqa: E402
from tissueresolve.experimental.spatial_smoothing import fit_level_specific_spatial  # noqa: E402

OUT = REPO / "benchmarks" / "outputs"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
SEEDS = list(range(10))
N_SIDE, CELLS_PER_SPOT, MIN_CELLS = 14, 40, 30
MODES = {  # (lambda_broad, lambda_fine)
    "none": (0.0, 0.0), "default_equal": (0.1, 0.1), "weak_equal": (0.02, 0.02),
    "broad_only": (0.1, 0.0), "fine_only": (0.0, 0.1), "level_specific": (0.1, 0.02),
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--n-seeds", type=int, default=len(SEEDS))
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIER)

    rows = []
    for s in range(args.n_seeds):
        ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=s)
        mapped = adata.obs["cell_type"].astype(str).isin(set(raw_map)).to_numpy()
        rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy() & mapped
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[rmask].copy(), min_cells=MIN_CELLS,
                                      estimate_overdispersion=True).reference
        ref_types = [str(c) for c in ref.cell_types]
        mapping = build_cell_type_hierarchy(ref_types, raw_map)
        rare = next((t for t in ["regulatory T cell", "natural killer cell"] if t in ref_types), ref_types[-1])
        sc = SSP.generate_spatial_scenario(adata, ref_types, qry, mapping, n_side=N_SIDE,
                                           cells_per_spot=CELLS_PER_SPOT, seed=s, rare_type=rare)
        truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
        for mode, (lb, lf) in MODES.items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = fit_level_specific_spatial(sc["Y"], ref, sc["array_row"], sc["array_col"],
                                               sc["lib_sizes"], sc["gene_names"], mapping,
                                               lambda_broad=lb, lambda_fine=lf, spot_ids=sc["spot_ids"])
            pred = r.fine_proportions
            cols = list(truth.columns)
            p = pred.reindex(index=truth.index, columns=cols).fillna(0.0)
            acc = M.accuracy_metrics(truth, p); comp = M.compositional_metrics(truth, p)
            tfam = M.aggregate_to_families(truth, mapping); pfam = M.aggregate_to_families(p, mapping)
            fid = SM.spatial_fidelity_metrics(truth, p, coords, domain_labels=domains, k=6)
            bnd = SM.boundary_metrics(truth, p, coords, domains, mapping=mapping, k=6)
            rns = SM.rare_niche_sensitivity(truth, p, domains, rare)
            cx = M.composition_complexity(p)
            rows.append({"seed": s, "mode": mode, "lambda_broad": lb, "lambda_fine": lf,
                         "runtime_s": round(r.runtime, 2),
                         "broad_rmse": M.accuracy_metrics(tfam, pfam)["rmse"],
                         "fine_rmse": acc["rmse"], "fine_pearson": acc["pearson"],
                         "jsd": comp["jsd_mean"], "local_rmse": fid["local_rmse"],
                         "oversmoothing": fid["oversmoothing_score"],
                         "morans_i_mae": fid["morans_i_mae"],
                         "boundary_f1": bnd["boundary_f1"], "edge_blurring": bnd["edge_blurring"],
                         "rare_niche_sensitivity": rns,
                         "richness": float((p > 0.01).sum(1).mean()),
                         "entropy": float(cx["shannon_entropy_nats"].mean())})
        print(f"  seed{s} done", flush=True)

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "level_specific_spatial_smoothing.tsv", "w") as fh:
        fh.write("# Task D level-specific spatial smoothing. 10 seeds × 6 modes (structured synthetic). Experimental.\n")
        df.to_csv(fh, sep="\t", index=False)
    bnd_cols = ["seed", "mode", "boundary_f1", "edge_blurring", "local_rmse", "oversmoothing", "morans_i_mae"]
    df[bnd_cols].to_csv(OUT / "spatial_boundary_metrics.tsv", sep="\t", index=False)

    g = df.groupby("mode").mean(numeric_only=True)
    ref_mode = "default_equal"
    ls, de = g.loc["level_specific"], g.loc[ref_mode]
    def piv(metric):
        return df.pivot_table(index="seed", columns="mode", values=metric)
    # paired Wilcoxon level_specific vs default on fine_rmse + edge_blurring
    wil = {}
    for met in ["fine_rmse", "edge_blurring", "broad_rmse", "rare_niche_sensitivity", "boundary_f1"]:
        pv = piv(met)
        wil[met] = ST.paired_wilcoxon(pv["level_specific"].to_numpy(), pv[ref_mode].to_numpy())
    rel = lambda a, b: (a - b) / abs(b) if b else float("inf")
    gates = [
        ("fine RMSE improves vs default", ls["fine_rmse"], de["fine_rmse"], ls["fine_rmse"] < de["fine_rmse"]),
        ("broad RMSE preserved (≤+2%)", ls["broad_rmse"], de["broad_rmse"], rel(ls["broad_rmse"], de["broad_rmse"]) <= 0.02),
        ("edge blurring reduced", ls["edge_blurring"], de["edge_blurring"], ls["edge_blurring"] < de["edge_blurring"]),
        ("boundary F1 preserved (≥ default −0.02)", ls["boundary_f1"], de["boundary_f1"], ls["boundary_f1"] >= de["boundary_f1"] - 0.02),
        ("rare-niche sensitivity not reduced", ls["rare_niche_sensitivity"], de["rare_niche_sensitivity"],
         ls["rare_niche_sensitivity"] >= de["rare_niche_sensitivity"] - 0.02),
        ("oversmoothing closer to 1.0", abs(ls["oversmoothing"] - 1), abs(de["oversmoothing"] - 1),
         abs(ls["oversmoothing"] - 1) < abs(de["oversmoothing"] - 1)),
    ]
    gate_df = pd.DataFrame([{"gate": n, "level_specific": round(float(a), 4),
                             "default_equal": round(float(b), 4), "pass": bool(p)} for n, a, b, p in gates])
    gate_df["fine_rmse_wilcoxon_p"] = round(wil["fine_rmse"]["p_value"], 4)
    with open(OUT / "taskD_promotion_gates.tsv", "w") as fh:
        fh.write("# Task D promotion gates: level_specific vs current default (equal λ=0.1). Experimental.\n")
        gate_df.to_csv(fh, sep="\t", index=False)

    _figure(g)
    print("\n=== Task D mode means (10 seeds) ===")
    print(g[["fine_rmse", "broad_rmse", "edge_blurring", "boundary_f1", "oversmoothing",
             "rare_niche_sensitivity", "local_rmse", "runtime_s"]].round(3).to_string())
    print("\n=== PROMOTION GATES (level_specific vs default_equal) ===")
    print(gate_df.to_string(index=False))
    print(f"\nWilcoxon fine_rmse p={wil['fine_rmse']['p_value']:.4f}; "
          f"GATES PASSED: {int(gate_df['pass'].sum())}/{len(gate_df)}")
    return 0


def _figure(g):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG = OUT / "figures" / "taskD"; FIG.mkdir(parents=True, exist_ok=True)
    order = ["none", "fine_only", "broad_only", "weak_equal", "default_equal", "level_specific"]
    gg = g.reindex(order)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    gg[["fine_rmse", "broad_rmse"]].plot.bar(ax=ax[0]); ax[0].set_title("RMSE by smoothing mode")
    ax[0].set_ylabel("RMSE"); plt.setp(ax[0].get_xticklabels(), rotation=25, ha="right")
    gg[["edge_blurring", "boundary_f1"]].plot.bar(ax=ax[1]); ax[1].set_title("Edge blurring & boundary F1")
    plt.setp(ax[1].get_xticklabels(), rotation=25, ha="right")
    fig.tight_layout(); fig.savefig(FIG / "figT_level_specific.png", dpi=150)
    fig.savefig(FIG / "figT_level_specific.svg"); plt.close(fig)
    gg.to_csv(FIG / "figT_level_specific.data.tsv", sep="\t")
    (FIG / "figT_level_specific.caption.txt").write_text(
        "Level-specific vs equal/single-level spatial smoothing across 10 seeds. "
        "Hypothesis: λ_broad>λ_fine lowers fine RMSE + edge blurring while preserving "
        "broad RMSE and boundaries.\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
