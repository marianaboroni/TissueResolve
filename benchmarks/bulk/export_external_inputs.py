#!/usr/bin/env python
"""Export held-out-donor inputs for external bulk methods (MuSiC / BisqueRNA).

Regenerates the SAME seed-0 donor-disjoint scenarios that
`run_holdout_bulk_benchmark.py` used, and writes them in the format the existing
R wrappers expect, so MuSiC/Bisque are scored on **identical** reference + bulk
as TissueResolve (fair-input contract, PART 4).

Outputs (benchmarks/outputs/holdout_bulk/external_inputs/)
- reference/reference_counts_genes_by_cells.tsv   (ref-donor cells, subsampled)
- reference/reference_cell_metadata.tsv           (cell_id, cellType, SubjectName)
- reference_missing/<...>                          (same, excluding the missing type)
- <scenario>/bulk_counts_genes_by_samples.tsv
- <scenario>/truth_mrna.tsv
- export_manifest.json

Usage:  python benchmarks/bulk/export_external_inputs.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "holdout_bulk" / "external_inputs"
SEED = 0
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
MAX_CELLS_PER_TYPE = 120          # subsample sc reference for tractable R runtime
SCEN_SEED = {"balanced": 11, "imbalanced": 22, "rare": 33,
             "similar_subtypes": 44, "missing_population": 55}
HIERARCHY_TSV = REPO / "examples" / "real_breast_cancer" / "config" / \
    "breast_cancer_cell_type_hierarchy.tsv"


def _export_reference(adata, mask, name, rng):
    """Write a subsampled genes×cells sc reference + metadata for cells in *mask*."""
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
    genes = list(map(str, adata.var_names))
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=genes, columns=cell_ids).to_csv(
        d / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cell_ids, "cellType": ct[keep],
                  "SubjectName": donor[keep]}).to_csv(
        d / "reference_cell_metadata.tsv", sep="\t", index=False)
    return len(keep), len(np.unique(ct[keep]))


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
    # Match the TissueResolve runner exactly: reference cell types = those that
    # survive min_cells=30 on the reference donors (not all 41 atlas types).
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        _ref = H.prepare_reference(adata[ref_mask].copy(), min_cells=30,
                                   estimate_overdispersion=False).reference
    all_types = [str(c) for c in _ref.cell_types]
    raw_map = load_hierarchy_mapping(HIERARCHY_TSV)
    mapping = build_cell_type_hierarchy(all_types, raw_map)
    # restrict exported single cells to the reference types
    type_ok = adata.obs["cell_type"].astype(str).isin(set(all_types)).to_numpy()
    ref_mask = ref_mask & type_ok
    missing_type = next((t for t in all_types if mapping.get(t) == "Epithelial"), all_types[0])
    # match the runner's reference types (min_cells=30 on ref donors) for rare/similar params
    fam_members: dict = {}
    for t in all_types:
        fam_members.setdefault(mapping.get(t, t), []).append(t)
    sim_fam = max(fam_members, key=lambda f: len(fam_members[f]))
    similar_pair = tuple(fam_members[sim_fam][:2])
    rare_type = next((t for t in ["regulatory T cell", "natural killer cell"]
                      if t in all_types), all_types[-1])

    rng = np.random.default_rng(SEED)
    n_ref, n_ref_types = _export_reference(adata, ref_mask, "reference", rng)
    miss_mask = ref_mask & (adata.obs["cell_type"].astype(str) != missing_type).to_numpy()
    n_miss, _ = _export_reference(adata, miss_mask, "reference_missing", rng)
    print(f"reference: {n_ref} cells / {n_ref_types} types; missing-ref: {n_miss} cells")

    scen_info = {}
    for scen, sseed in SCEN_SEED.items():
        kw = {}
        if scen == "rare":
            kw = {"rare_type": rare_type, "rare_level": 0.01}
        elif scen == "similar_subtypes":
            kw = {"similar_pair": similar_pair}
        elif scen == "missing_population":
            kw = {"missing_type": missing_type}
        targets = SH.build_target_proportions(all_types, N_SAMPLES, scen, seed=sseed, **kw)
        ds = SH.realize_pseudobulk(adata, targets, celltype_col="cell_type",
                                   donor_col="donor_id", query_donors=query_donors,
                                   seed=SEED, cells_per_sample=CELLS_PER_SAMPLE)
        d = OUT / scen
        d.mkdir(parents=True, exist_ok=True)
        ds.counts.to_csv(d / "bulk_counts_genes_by_samples.tsv", sep="\t")
        ds.true_mrna_proportions.to_csv(d / "truth_mrna.tsv", sep="\t")
        scen_info[scen] = {"reference": "reference_missing" if scen == "missing_population" else "reference",
                           "n_samples": ds.counts.shape[1]}
        print(f"  exported {scen}: {ds.counts.shape}")

    (OUT / "export_manifest.json").write_text(json.dumps({
        "seed": SEED, "n_samples": N_SAMPLES, "cells_per_sample": CELLS_PER_SAMPLE,
        "max_cells_per_type": MAX_CELLS_PER_TYPE,
        "reference_donors": ref_donors, "query_donors": query_donors,
        "donor_disjoint": True, "missing_type": missing_type,
        "rare_type": rare_type, "similar_pair": list(similar_pair),
        "scenarios": scen_info,
    }, indent=2), encoding="utf-8")
    print(f"\nExported -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
