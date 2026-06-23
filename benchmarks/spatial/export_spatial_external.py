#!/usr/bin/env python
"""Export the synthetic-spatial scenario for external methods (RCTD, cell2location).

Regenerates the EXACT seed-0 scenario used by run_synthetic_spatial.py and writes:
- reference/reference_counts_genes_by_cells.tsv + reference_cell_metadata.tsv
  (held-out reference-donor cells, subsampled — same donor split as TissueResolve)
- spatial/spot_counts_genes_by_spots.tsv, spot_coords.tsv, truth_mrna.tsv

Usage:  python benchmarks/spatial/export_spatial_external.py --run-real-data
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
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import synthetic_spatial as SSP  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "spatial_synthetic" / "external_inputs"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
SEED = 0
N_SIDE = 20
CELLS_PER_SPOT = 40
MIN_CELLS = 30
MAX_CELLS_PER_TYPE = 120


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

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
                                  estimate_overdispersion=False).reference
    ref_types = [str(c) for c in ref.cell_types]
    mapping = build_cell_type_hierarchy(ref_types, load_hierarchy_mapping(HIERARCHY_TSV))
    rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                      if t in ref_types), ref_types[-1])

    # --- sc reference: subsample ref-donor cells of the reference types ---
    rng = np.random.default_rng(SEED)
    type_ok = adata.obs["cell_type"].astype(str).isin(set(ref_types)).to_numpy()
    idx_all = np.where(ref_mask & type_ok)[0]
    ct = adata.obs["cell_type"].astype(str).to_numpy()
    keep = []
    for t in np.unique(ct[idx_all]):
        pool = idx_all[ct[idx_all] == t]
        if pool.size > MAX_CELLS_PER_TYPE:
            pool = rng.choice(pool, MAX_CELLS_PER_TYPE, replace=False)
        keep.extend(pool.tolist())
    keep = np.sort(np.array(keep))
    X = adata.X[keep]; X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    cell_ids = [f"c{i}" for i in keep]
    refd = OUT / "reference"; refd.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=list(map(str, adata.var_names)),
                 columns=cell_ids).to_csv(refd / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cell_ids, "cellType": ct[keep]}).to_csv(
        refd / "reference_cell_metadata.tsv", sep="\t", index=False)

    # --- spatial scenario (identical to run_synthetic_spatial.py seed 0) ---
    sc = SSP.generate_spatial_scenario(adata, ref_types, query_donors, mapping,
                                       n_side=N_SIDE, cells_per_spot=CELLS_PER_SPOT,
                                       seed=SEED, rare_type=rare_type)
    spd = OUT / "spatial"; spd.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sc["Y"].T, index=sc["gene_names"], columns=sc["spot_ids"]).to_csv(
        spd / "spot_counts_genes_by_spots.tsv", sep="\t")
    pd.DataFrame({"spot": sc["spot_ids"], "row": sc["array_row"], "col": sc["array_col"]}).to_csv(
        spd / "spot_coords.tsv", sep="\t", index=False)
    sc["truth"].to_csv(spd / "truth_mrna.tsv", sep="\t")
    pd.DataFrame({"spot": sc["truth"].index, "domain": sc["domain_labels"]}).to_csv(
        spd / "domain_labels.tsv", sep="\t", index=False)

    (OUT / "manifest.json").write_text(json.dumps({
        "seed": SEED, "n_side": N_SIDE, "n_spots": int(len(sc["truth"])),
        "cells_per_spot": CELLS_PER_SPOT, "n_ref_cells": len(keep), "n_ref_types": len(ref_types),
        "reference_donors": ref_donors, "query_donors": query_donors, "donor_disjoint": True,
        "rare_type": rare_type, "max_cells_per_type": MAX_CELLS_PER_TYPE,
    }, indent=2), encoding="utf-8")
    print(f"Exported {len(keep)} ref cells ({len(ref_types)} types) + "
          f"{len(sc['truth'])} spots -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
