#!/usr/bin/env python
"""Multi-seed synthetic spatial — internal methods + export inputs for externals.

For each seed: donor split → reference (reference donors) → structured spatial
scenario realised from held-out query donors → run TissueResolve spatial (flat +
hierarchical) and NNLS-per-spot, score (accuracy + spatial fidelity), and export
the sc reference + spot counts so RCTD / cell2location can run on the SAME data.

Smaller grid (12×12 = 144 spots) than the single-run benchmark to bound the
~30-min-per-run external cost across seeds.

Outputs (benchmarks/outputs/spatial_multiseed/)
- internal_metrics_long.tsv  (seed, method, metric, value)
- seed<s>/external_inputs/{reference,spatial}/...   (for the external runners)
- seed<s>/truth.tsv, domain_labels.tsv

Usage:  python benchmarks/spatial/run_spatial_multiseed.py --run-real-data
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

OUT = REPO / "benchmarks" / "outputs" / "spatial_multiseed"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
N_SEEDS = 5
N_SIDE = 12
CELLS_PER_SPOT = 40
MIN_CELLS = 30
MAX_CELLS_PER_TYPE = 120


def _score(truth, pred, coords, domains, mapping):
    cols = [c for c in truth.columns if not str(c).startswith("unresolved_")]
    t = truth[cols]; p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p); comp = M.compositional_metrics(t, p)
    tfam, pfam = M.aggregate_to_families(t, mapping), M.aggregate_to_families(p, mapping)
    fid = SM.spatial_fidelity_metrics(t, p, coords, domain_labels=domains, k=6)
    return {"fine_pearson": acc["pearson"], "fine_rmse": acc["rmse"],
            "fine_jsd": comp["jsd_mean"], "broad_pearson": M.accuracy_metrics(tfam, pfam)["pearson"],
            "dominant_accuracy": M.dominant_accuracy(t, p),
            "local_rmse": fid["local_rmse"], "oversmoothing_score": fid["oversmoothing_score"],
            "domain_recovery_ari": fid.get("domain_recovery_ari", float("nan"))}


def _export_inputs(adata, ref_mask, ref_types, sc_dir, rng):
    ct = adata.obs["cell_type"].astype(str).to_numpy()
    type_ok = np.isin(ct, list(ref_types))
    idx_all = np.where(ref_mask & type_ok)[0]
    keep = []
    for t in np.unique(ct[idx_all]):
        pool = idx_all[ct[idx_all] == t]
        if pool.size > MAX_CELLS_PER_TYPE:
            pool = rng.choice(pool, MAX_CELLS_PER_TYPE, replace=False)
        keep.extend(pool.tolist())
    keep = np.sort(np.array(keep))
    X = adata.X[keep]; X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    cids = [f"c{i}" for i in keep]
    sc_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=list(map(str, adata.var_names)),
                 columns=cids).to_csv(sc_dir / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cids, "cellType": ct[keep]}).to_csv(
        sc_dir / "reference_cell_metadata.tsv", sep="\t", index=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--n-seeds", type=int, default=N_SEEDS)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    import anndata as ad
    import tissueresolve as tr
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)
    rows = []

    for s in range(args.n_seeds):
        ref_donors, query_donors = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=s)
        ref_mask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[ref_mask].copy(), min_cells=MIN_CELLS,
                                      estimate_overdispersion=True).reference
        ref_types = [str(c) for c in ref.cell_types]
        mapping = build_cell_type_hierarchy(ref_types, raw_map)
        rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                          if t in ref_types), ref_types[-1])
        sc = SSP.generate_spatial_scenario(adata, ref_types, query_donors, mapping,
                                           n_side=N_SIDE, cells_per_spot=CELLS_PER_SPOT,
                                           seed=s, rare_type=rare_type)
        truth, coords, domains = sc["truth"], sc["coords"], sc["domain_labels"]
        sdir = OUT / f"seed{s}"
        (sdir).mkdir(parents=True, exist_ok=True)
        truth.to_csv(sdir / "truth.tsv", sep="\t")
        pd.DataFrame({"spot": truth.index, "domain": domains,
                      "row": sc["array_row"], "col": sc["array_col"]}).to_csv(
            sdir / "domain_labels.tsv", sep="\t", index=False)
        # export for external methods
        ei = sdir / "external_inputs"
        _export_inputs(adata, ref_mask, ref_types, ei / "reference", np.random.default_rng(s))
        pd.DataFrame(sc["Y"].T, index=sc["gene_names"], columns=sc["spot_ids"]).to_csv(
            ei / "spatial_spot_counts.tsv", sep="\t")
        pd.DataFrame({"spot": sc["spot_ids"], "row": sc["array_row"],
                      "col": sc["array_col"]}).to_csv(ei / "spatial_coords.tsv", sep="\t", index=False)

        bulk_df = pd.DataFrame(sc["Y"].T, index=sc["gene_names"], columns=sc["spot_ids"])

        def emit(name, pred):
            for metric, val in _score(truth, pred, coords, domains, mapping).items():
                rows.append({"seed": s, "method": name, "metric": metric, "value": float(val)})

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t0 = time.perf_counter()
            emit("NNLS_per_spot", deconv_bulk(bulk_df, ref, solver="nnls",
                                              resolution_mode="none").deconv.proportions)
            emit("TissueResolve_spatial_flat", tr.deconv_spatial(
                sc["Y"], ref, sc["array_row"], sc["array_col"], sc["lib_sizes"],
                sc["gene_names"], spot_ids=sc["spot_ids"], resolution_mode="none",
                run_neighbourhood=False).deconv.proportions)
            emit("TissueResolve_spatial_hierarchical", tr.deconv_spatial(
                sc["Y"], ref, sc["array_row"], sc["array_col"], sc["lib_sizes"],
                sc["gene_names"], spot_ids=sc["spot_ids"], resolution_mode="hierarchical",
                hierarchy_mapping=mapping, run_neighbourhood=False).estimates.combined_fine)
        print(f"seed{s}: {len(ref_types)} types, {len(truth)} spots, "
              f"internal done {time.perf_counter()-t0:.0f}s", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / "internal_metrics_long.tsv", sep="\t", index=False)
    (OUT / "manifest.json").write_text(json.dumps({
        "n_seeds": args.n_seeds, "n_side": N_SIDE, "cells_per_spot": CELLS_PER_SPOT,
        "replicate_unit": "seed (donor split + spatial layout)",
        "ground_truth": "per-spot mRNA-proportion"}, indent=2), encoding="utf-8")
    print(f"\nInternal metrics + per-seed external inputs -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
