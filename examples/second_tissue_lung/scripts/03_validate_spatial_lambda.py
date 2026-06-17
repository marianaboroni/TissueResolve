#!/usr/bin/env python
"""STAGE 5 — validate the weaker-λ spatial finding on the second tissue (lung).

HLCA is scRNA (no Visium), so we test the λ finding on HLCA-derived SYNTHETIC
spatial (structured: sharp border + gradient + rare niche), 5 seeds, comparing the
current default (λ=0.1) vs the candidate weaker λ=0.02 via the existing solver
(level-specific wrapper at equal λ). Question: does weak-λ beat the default on a
SECOND tissue (does the breast finding generalize)?

Outputs: benchmarks/outputs/lung_spatial_lambda.tsv.
Usage:  python 03_validate_spatial_lambda.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import spatial_metrics as SM  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402
from tissueresolve.experimental.spatial_smoothing import fit_level_specific_spatial  # noqa: E402

SUBSET = Path(__file__).resolve().parents[1] / "data" / "derived" / "hlca_subset.h5ad"
HMAP = Path(__file__).resolve().parents[1] / "data" / "derived" / "hlca_hierarchy.tsv"
OUT = REPO / "benchmarks" / "outputs"
SEEDS = [0, 1, 2, 3, 4]
N_SIDE, CELLS_PER_SPOT, MIN_CELLS = 12, 40, 20
MODES = {"default_0.1": 0.1, "weak_0.02": 0.02}   # equal λ at both levels


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    adata = ad.read_h5ad(SUBSET)
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).to_numpy(); adata.var_names_make_unique()
    hmap = pd.read_csv(HMAP, sep="\t")
    mapping = dict(zip(hmap["fine_cell_type"].astype(str), hmap["broad_cell_type"].astype(str)))

    rows = []
    for s in SEEDS:
        ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=s)
        rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[rmask].copy(), cell_type_col="cell_type_fine",
                                      min_cells=MIN_CELLS, estimate_overdispersion=True).reference
        ref_types = [str(c) for c in ref.cell_types]
        mp = build_cell_type_hierarchy(ref_types, mapping)
        rare = ref_types[-1]
        sc = SSP.generate_spatial_scenario(adata, ref_types, qry, mp, n_side=N_SIDE,
                                           cells_per_spot=CELLS_PER_SPOT, seed=s, rare_type=rare,
                                           celltype_col="cell_type_fine")
        truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
        for mode, lam in MODES.items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = fit_level_specific_spatial(sc["Y"], ref, sc["array_row"], sc["array_col"],
                                               sc["lib_sizes"], sc["gene_names"], mp,
                                               lambda_broad=lam, lambda_fine=lam, spot_ids=sc["spot_ids"])
            p = r.fine_proportions.reindex(index=truth.index, columns=list(truth.columns)).fillna(0.0)
            acc = M.accuracy_metrics(truth, p)
            tfam, pfam = M.aggregate_to_families(truth, mp), M.aggregate_to_families(p, mp)
            fid = SM.spatial_fidelity_metrics(truth, p, coords, domain_labels=domains, k=6)
            bnd = SM.boundary_metrics(truth, p, coords, domains, mapping=mp, k=6)
            rows.append({"seed": s, "mode": mode, "lambda": lam, "fine_rmse": acc["rmse"],
                         "fine_pearson": acc["pearson"], "broad_rmse": M.accuracy_metrics(tfam, pfam)["rmse"],
                         "oversmoothing": fid["oversmoothing_score"], "local_rmse": fid["local_rmse"],
                         "edge_blurring": bnd["edge_blurring"], "boundary_f1": bnd["boundary_f1"]})
        print(f"  seed{s} done", flush=True)

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "lung_spatial_lambda.tsv", "w") as fh:
        fh.write("# Second-tissue (HLCA synthetic spatial) weak-λ vs default-λ. 5 seeds. Experimental.\n")
        df.to_csv(fh, sep="\t", index=False)
    g = df.groupby("mode").mean(numeric_only=True)
    piv = df.pivot_table(index="seed", columns="mode", values="fine_rmse")
    w = ST.paired_wilcoxon(piv["weak_0.02"].to_numpy(), piv["default_0.1"].to_numpy())
    print("\n=== lung synthetic-spatial λ comparison (5 seeds) ===")
    print(g[["fine_rmse", "broad_rmse", "oversmoothing", "local_rmse", "edge_blurring", "boundary_f1"]].round(3).to_string())
    we, de = g.loc["weak_0.02"], g.loc["default_0.1"]
    print(f"\nweak λ vs default: fine RMSE {we['fine_rmse']:.4f} vs {de['fine_rmse']:.4f} "
          f"(Wilcoxon p={w['p_value']:.3f}); oversmoothing {we['oversmoothing']:.2f} vs {de['oversmoothing']:.2f}")
    better = we["fine_rmse"] <= de["fine_rmse"] and abs(we["oversmoothing"] - 1) < abs(de["oversmoothing"] - 1)
    print(f"weaker-λ generalizes to lung: {better}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
