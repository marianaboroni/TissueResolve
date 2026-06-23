#!/usr/bin/env python
"""Export ALL (split, scenario) held-out inputs for external methods (multi-split).

Replicates the EXACT data of `run_holdout_bulk_multisplit.py` (5 donor splits ×
5 scenarios, identical seeds) so MuSiC/BisqueRNA can be scored on the same 25
paired replicates as TissueResolve → enables bootstrap CIs + paired Wilcoxon.

Seeds (must match run_holdout_bulk_multisplit.py exactly):
  donor split seed      = split index s
  reference min_cells   = 30
  target seed           = 10 + SCEN_OFFSET[scen] + s*1000   (SCEN_OFFSET=i*100)
  realize seed          = s*7 + 1

Layout: benchmarks/outputs/holdout_bulk/external_inputs_multisplit/
  split<s>/reference/{reference_counts_genes_by_cells,reference_cell_metadata}.tsv
  split<s>/reference_missing/...
  split<s>/<scenario>/{bulk_counts_genes_by_samples,truth_mrna}.tsv

Usage:  python benchmarks/bulk/export_external_multisplit.py --run-real-data
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

OUT = REPO / "benchmarks" / "outputs" / "holdout_bulk" / "external_inputs_multisplit"
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"
SCENARIOS = ["balanced", "imbalanced", "rare", "similar_subtypes", "missing_population"]
SCEN_OFFSET = {s: i * 100 for i, s in enumerate(SCENARIOS)}
N_SPLITS = 5
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
MAX_CELLS_PER_TYPE = 120


def _export_reference(adata, mask, outdir, rng):
    obs = adata.obs
    ct = obs["cell_type"].astype(str).to_numpy()
    donor = obs["donor_id"].astype(str).to_numpy()
    idx_all = np.where(mask)[0]
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
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=list(map(str, adata.var_names)),
                 columns=cell_ids).to_csv(outdir / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cell_ids, "cellType": ct[keep],
                  "SubjectName": donor[keep]}).to_csv(
        outdir / "reference_cell_metadata.tsv", sep="\t", index=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--n-splits", type=int, default=N_SPLITS)
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
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)

    manifest = {"n_splits": args.n_splits, "scenarios": SCENARIOS, "splits": {}}
    for s in range(args.n_splits):
        ref_donors, query_donors = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=s)
        ref_mask = adata.obs["donor_id"].astype(str).isin(set(ref_donors)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[ref_mask].copy(), min_cells=30,
                                      estimate_overdispersion=False).reference
        ref_types = [str(c) for c in ref.cell_types]
        mapping = build_cell_type_hierarchy(ref_types, raw_map)
        fam_members: dict = {}
        for t in ref_types:
            fam_members.setdefault(mapping.get(t, t), []).append(t)
        sim_fam = max(fam_members, key=lambda f: len(fam_members[f]))
        similar_pair = tuple(fam_members[sim_fam][:2])
        rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                          if t in ref_types), ref_types[-1])
        missing_type = next((t for t in ref_types if mapping.get(t) == "Epithelial"),
                            ref_types[0])
        type_ok = adata.obs["cell_type"].astype(str).isin(set(ref_types)).to_numpy()
        rmask = ref_mask & type_ok
        mmask = rmask & (adata.obs["cell_type"].astype(str) != missing_type).to_numpy()

        sdir = OUT / f"split{s}"
        rng = np.random.default_rng(s)
        _export_reference(adata, rmask, sdir / "reference", rng)
        _export_reference(adata, mmask, sdir / "reference_missing", rng)

        for scen in SCENARIOS:
            kw = {}
            if scen == "rare":
                kw = {"rare_type": rare_type, "rare_level": 0.01}
            elif scen == "similar_subtypes":
                kw = {"similar_pair": similar_pair}
            elif scen == "missing_population":
                kw = {"missing_type": missing_type}
            tseed = 10 + SCEN_OFFSET[scen] + s * 1000
            targets = SH.build_target_proportions(ref_types, N_SAMPLES, scen, seed=tseed, **kw)
            ds = SH.realize_pseudobulk(adata, targets, celltype_col="cell_type",
                                       donor_col="donor_id", query_donors=query_donors,
                                       seed=s * 7 + 1, cells_per_sample=CELLS_PER_SAMPLE)
            d = sdir / scen
            d.mkdir(parents=True, exist_ok=True)
            ds.counts.to_csv(d / "bulk_counts_genes_by_samples.tsv", sep="\t")
            ds.true_mrna_proportions.to_csv(d / "truth_mrna.tsv", sep="\t")
        manifest["splits"][s] = {"reference_donors": ref_donors, "query_donors": query_donors,
                                 "n_ref_types": len(ref_types), "missing_type": missing_type}
        print(f"split{s}: {len(ref_types)} types, exported {len(SCENARIOS)} scenarios", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "export_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nExported -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
