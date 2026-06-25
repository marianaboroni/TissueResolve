#!/usr/bin/env python
"""Export the LUNG synthetic-spatial seed-0 scenario for external tools (CARD).

Reuses run_weak_smoothing_grid._load_dataset so the reference, donor split, and
scenario are IDENTICAL to the lung weak-smoothing grid (seed 0). Writes the same
external_inputs layout that run_card_synthetic.R consumes, under a lung-specific
directory.

Usage:  PYTHONPATH=src:. python benchmarks/spatial/export_lung_spatial_external.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "spatial_synthetic_lung" / "external_inputs"
SEED = 0
N_SIDE = 20
CELLS_PER_SPOT = 40
MAX_CELLS_PER_TYPE = 150


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    D = _load_dataset("lung")
    adata, ref_types = D["adata"], set(D["ref_types"])
    donor_col, ct_col = D["donor_col"], D["celltype_col"]
    ref_donors = set(D["ref_donors"])

    # subsample reference-donor cells of the reference types for the external tool
    rng = np.random.default_rng(SEED)
    donors = adata.obs[donor_col].astype(str).to_numpy()
    ct = adata.obs[ct_col].astype(str).to_numpy()
    in_ref = np.isin(donors, list(ref_donors)) & np.isin(ct, list(ref_types))
    idx_all = np.where(in_ref)[0]
    keep = []
    for t in np.unique(ct[idx_all]):
        pool = idx_all[ct[idx_all] == t]
        if pool.size > MAX_CELLS_PER_TYPE:
            pool = rng.choice(pool, MAX_CELLS_PER_TYPE, replace=False)
        keep.extend(pool.tolist())
    keep = np.sort(np.array(keep))
    X = adata.X[keep]
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    cell_ids = [f"c{i}" for i in keep]
    refd = OUT / "reference"; refd.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=list(map(str, adata.var_names)),
                 columns=cell_ids).to_csv(refd / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cell_ids, "cellType": ct[keep]}).to_csv(
        refd / "reference_cell_metadata.tsv", sep="\t", index=False)

    # identical seed-0 lung scenario
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc = SSP.generate_spatial_scenario(
            adata, D["ref_types"], D["query_donors"], D["mapping"],
            celltype_col=ct_col, donor_col=donor_col, n_side=N_SIDE,
            cells_per_spot=CELLS_PER_SPOT, seed=SEED, rare_type=D["rare_type"],
            domain_families=D["domain_families"])
    spd = OUT / "spatial"; spd.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sc["Y"].T, index=sc["gene_names"], columns=sc["spot_ids"]).to_csv(
        spd / "spot_counts_genes_by_spots.tsv", sep="\t")
    pd.DataFrame({"spot": sc["spot_ids"], "row": sc["array_row"], "col": sc["array_col"]}).to_csv(
        spd / "spot_coords.tsv", sep="\t", index=False)
    sc["truth"].to_csv(spd / "truth_mrna.tsv", sep="\t")
    pd.DataFrame({"spot": sc["truth"].index, "domain": sc["domain_labels"]}).to_csv(
        spd / "domain_labels.tsv", sep="\t", index=False)
    (OUT / "manifest.json").write_text(json.dumps({
        "dataset": "lung", "seed": SEED, "n_side": N_SIDE,
        "n_spots": int(len(sc["truth"])), "n_ref_cells": int(keep.size),
        "rare_type": D["rare_type"], "ground_truth": "mRNA-proportion per spot",
    }, indent=2), encoding="utf-8")
    print(f"Exported {keep.size} ref cells + {len(sc['truth'])} spots -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
