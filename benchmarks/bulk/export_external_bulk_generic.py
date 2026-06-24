#!/usr/bin/env python
"""Export held-out-donor bulk inputs for external tools (MuSiC/Bisque), generic over
dataset (breast / HLCA-lung) and seeds. Writes the SAME format the R wrapper expects.

For each dataset it writes one subsampled single-cell reference (ref-donor cells) and,
per (scenario, seed), a pseudobulk count matrix + mRNA-proportion truth — all from the
identical donor-disjoint split used elsewhere. Reuses the spatial benchmark's
`_load_dataset` so the reference/donor split matches.

Outputs (benchmarks/outputs/external_bulk/<dataset>/):
  reference/reference_counts_genes_by_cells.tsv, reference_cell_metadata.tsv
  <scenario>_seed<k>/bulk_counts_genes_by_samples.tsv, truth_mrna.tsv
  manifest.json  (mapping, rare_type, scenarios, seeds)

Usage:
  PYTHONPATH=src:. python benchmarks/bulk/export_external_bulk_generic.py \
      --run-real-data --dataset lung --seeds 0 1 2 3 4
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
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "external_bulk"
N_SAMPLES = 12
CELLS_PER_SAMPLE = 600
MAX_CELLS_PER_TYPE = 120
SCENARIOS = ["imbalanced", "rare", "similar_subtypes"]


def _export_reference(adata, mask, ct_col, donor_col, ref_types, out_dir, rng):
    obs = adata.obs
    ct = obs[ct_col].astype(str).to_numpy()
    donor = obs[donor_col].astype(str).to_numpy()
    idx_all = np.where(mask & np.isin(ct, list(ref_types)))[0]
    keep = []
    for t in np.unique(ct[idx_all]):
        pool = idx_all[ct[idx_all] == t]
        if pool.size > MAX_CELLS_PER_TYPE:
            pool = rng.choice(pool, MAX_CELLS_PER_TYPE, replace=False)
        keep.extend(pool.tolist())
    keep = np.sort(np.array(keep))
    X = adata.X[keep]
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    genes = list(map(str, adata.var_names))
    cell_ids = [f"c{i}" for i in keep]
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(X.T.astype(int), index=genes, columns=cell_ids).to_csv(
        out_dir / "reference_counts_genes_by_cells.tsv", sep="\t")
    pd.DataFrame({"cell_id": cell_ids, "cellType": ct[keep], "SubjectName": donor[keep]}).to_csv(
        out_dir / "reference_cell_metadata.tsv", sep="\t", index=False)
    return len(keep), int(len(np.unique(ct[keep])))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--dataset", required=True, choices=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    D = _load_dataset(args.dataset)
    adata, mapping = D["adata"], D["mapping"]
    ct_col, donor_col = D["celltype_col"], D["donor_col"]
    ref_types = set(D["ref_types"])
    rare = D["rare_type"]
    ref_donors = set(D["ref_donors"])
    q_donors = D["query_donors"]
    ref_mask = adata.obs[donor_col].astype(str).isin(ref_donors).to_numpy()

    base = OUT / args.dataset
    rng = np.random.default_rng(0)
    n_cells, n_types = _export_reference(adata, ref_mask, ct_col, donor_col, ref_types,
                                         base / "reference", rng)
    print(f"[{args.dataset}] reference: {n_cells} cells / {n_types} types; rare={rare}")

    # collinear within-family pair for similar_subtypes
    fam: dict = {}
    for t in sorted(ref_types):
        fam.setdefault(mapping.get(t, t), []).append(t)
    sim_pair = next((m[:2] for m in fam.values() if len(m) >= 2), None)

    scen_dirs = []
    for seed in args.seeds:
        for scen in SCENARIOS:
            kw = {}
            if scen == "rare":
                kw = {"rare_type": rare, "rare_level": 0.01}
            elif scen == "similar_subtypes" and sim_pair:
                kw = {"similar_pair": tuple(sim_pair)}
            targets = SH.build_target_proportions(sorted(ref_types), N_SAMPLES, scen,
                                                  seed=seed, **kw)
            ds = SH.realize_pseudobulk(adata, targets, celltype_col=ct_col,
                                       donor_col=donor_col, query_donors=q_donors,
                                       seed=seed, cells_per_sample=CELLS_PER_SAMPLE)
            name = f"{scen}_seed{seed}"
            d = base / name
            d.mkdir(parents=True, exist_ok=True)
            ds.counts.to_csv(d / "bulk_counts_genes_by_samples.tsv", sep="\t")
            ds.true_mrna_proportions.to_csv(d / "truth_mrna.tsv", sep="\t")
            scen_dirs.append(name)
            print(f"  exported {name}: {ds.counts.shape}")

    (base / "manifest.json").write_text(json.dumps({
        "dataset": args.dataset, "seeds": args.seeds, "scenarios": SCENARIOS,
        "scenario_dirs": scen_dirs, "rare_type": rare,
        "mapping": {str(k): str(v) for k, v in mapping.items()},
        "n_ref_cells": n_cells, "n_ref_types": n_types,
        "max_cells_per_type": MAX_CELLS_PER_TYPE, "n_samples": N_SAMPLES,
        "cells_per_sample": CELLS_PER_SAMPLE, "donor_disjoint": True,
    }, indent=2), encoding="utf-8")
    print(f"Exported -> {base}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
