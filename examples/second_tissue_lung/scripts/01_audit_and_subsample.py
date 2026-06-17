#!/usr/bin/env python
"""STAGE 3 — audit HLCA-core labels + subsample to a tractable validation atlas.

Audits obs (broad/fine annotation levels, donors, raw-count availability) WITHOUT
loading the full 5.9 GB matrix (anndata backed mode), then subsamples to a
donor-rich, fine-labelled atlas with RAW COUNTS and writes a small h5ad +
fine→broad hierarchy for the validation. Documents the broad/fine column choice
(rule 9) and flags low-support labels (rule 8 — flagged, not silently dropped).

Outputs (data/derived/): hlca_subset.h5ad, hlca_hierarchy.tsv,
hlca_label_support.tsv, hlca_audit.tsv.
Usage:  python 01_audit_and_subsample.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
H5AD = DATA / "hlca_core.h5ad"
DERIVED = DATA / "derived"
MAX_CELLS_PER_FINE = 200       # subsample cap per fine subtype
MIN_CELLS_FINE, MIN_DONORS_FINE = 80, 4
MAX_DONORS = 40                # cap donors for tractable CV
SEED = 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    if not H5AD.exists():
        print(f"Missing {H5AD}. Run 00_download_hlca.py first.", file=sys.stderr); return 2
    import anndata as ad
    DERIVED.mkdir(parents=True, exist_ok=True)
    print("opening (backed) ...", flush=True)
    A = ad.read_h5ad(H5AD, backed="r")
    obs = A.obs
    print(f"shape {A.shape}; obs cols: {list(obs.columns)}", flush=True)

    # discover annotation-level columns + donor column
    ann_cols = [c for c in obs.columns if c.lower().startswith("ann_level") or c == "ann_finest_level"]
    if "cell_type" in obs.columns and "cell_type" not in ann_cols:
        ann_cols.append("cell_type")
    donor_col = next((c for c in ["donor_id", "donor", "subject_ID", "patient"] if c in obs.columns), None)
    audit = []
    for c in ann_cols + ([donor_col] if donor_col else []):
        audit.append({"column": c, "n_unique": int(obs[c].nunique()),
                      "examples": "; ".join(map(str, pd.Series(obs[c].unique()).astype(str)[:6]))})
    adf = pd.DataFrame(audit); adf.to_csv(DERIVED / "hlca_audit.tsv", sep="\t", index=False)
    print("=== annotation audit ===\n" + adf.to_string(index=False), flush=True)

    # choose broad (cardinality ~6..15 → a mid level like ann_level_2, NOT the
    # 4-class ann_level_1 which collapses all immune lineages) and a finer fine level
    cand = [(c, obs[c].nunique()) for c in ann_cols]
    broad_col = next((c for c, n in sorted(cand, key=lambda x: x[1]) if 6 <= n <= 15), None)
    fine_col = next((c for c, n in sorted(cand, key=lambda x: -x[1]) if 18 <= n <= 70), None)
    if fine_col is None:  # fall back to the finest available
        fine_col = max(cand, key=lambda x: x[1])[0]
    if broad_col is None:
        broad_col = min(cand, key=lambda x: x[1])[0]
    print(f"\nchosen: broad='{broad_col}' ({obs[broad_col].nunique()}), "
          f"fine='{fine_col}' ({obs[fine_col].nunique()}), donor='{donor_col}'", flush=True)
    if not donor_col:
        print("ERROR: no donor column — cannot do donor-disjoint validation.", file=sys.stderr); return 3

    fine = obs[fine_col].astype(str); broad = obs[broad_col].astype(str); don = obs[donor_col].astype(str)
    # label support
    sup = pd.DataFrame({"broad": fine.groupby(fine).apply(lambda s: broad[s.index].mode().iloc[0]),
                        "n_cells": fine.groupby(fine).size(),
                        "n_donors": don.groupby(fine).nunique()})
    sup["low_support"] = (sup["n_cells"] < MIN_CELLS_FINE) | (sup["n_donors"] < MIN_DONORS_FINE)
    sup = sup.sort_values("n_cells"); sup.to_csv(DERIVED / "hlca_label_support.tsv", sep="\t")
    keep_fine = sup.index[~sup["low_support"]].tolist()
    print(f"\nfine subtypes: {len(sup)} total, {len(keep_fine)} pass support "
          f"(≥{MIN_CELLS_FINE} cells & ≥{MIN_DONORS_FINE} donors). {int(sup['low_support'].sum())} low-support FLAGGED.", flush=True)

    # hierarchy mapping (fine → broad), documented
    hmap = sup.loc[keep_fine, ["broad"]].reset_index()
    hmap.columns = ["fine_cell_type", "broad_cell_type"]
    hmap.to_csv(DERIVED / "hlca_hierarchy.tsv", sep="\t", index=False)
    nfam = hmap["broad_cell_type"].nunique()
    multi = hmap.groupby("broad_cell_type").size()
    print(f"hierarchy: {nfam} broad families, {len(hmap)} fine; "
          f"{int((multi>=2).sum())} families with ≥2 subtypes (testable).", flush=True)

    # choose donors (cap MAX_DONORS, prefer those covering many kept fine types)
    rng = np.random.default_rng(SEED)
    mask_keep = fine.isin(keep_fine).to_numpy()
    donors_all = sorted(don[mask_keep].unique())
    if len(donors_all) > MAX_DONORS:
        # rank donors by # kept fine types present, take top MAX_DONORS
        cov = don[mask_keep].groupby(don[mask_keep]).apply(lambda s: fine[s.index].nunique())
        donors_sel = cov.sort_values(ascending=False).index[:MAX_DONORS].tolist()
    else:
        donors_sel = donors_all
    sel = mask_keep & don.isin(donors_sel).to_numpy()
    idx_sel = np.where(sel)[0]

    # subsample cells per fine type
    keep_idx = []
    for ft in keep_fine:
        pool = idx_sel[(fine.iloc[idx_sel] == ft).to_numpy()]
        if pool.size > MAX_CELLS_PER_FINE:
            pool = rng.choice(pool, MAX_CELLS_PER_FINE, replace=False)
        keep_idx.extend(pool.tolist())
    keep_idx = np.sort(np.array(keep_idx))
    print(f"\nsubsampling {len(keep_idx)} cells from {len(donors_sel)} donors ...", flush=True)

    sub = A[keep_idx].to_memory()
    # RAW counts: prefer .raw.X (CELLxGENE schema); harmonise var to symbols
    if sub.raw is not None:
        raw = sub.raw.to_adata()
    else:
        raw = sub
    import scipy.sparse as sp
    X = raw.X
    is_int = bool(np.allclose((X[:50].toarray() if sp.issparse(X) else X[:50]),
                              np.round(X[:50].toarray() if sp.issparse(X) else X[:50])))
    # var symbols
    if "feature_name" in raw.var.columns:
        raw.var_names = raw.var["feature_name"].astype(str).to_numpy()
    raw.var_names_make_unique()
    out = raw  # genes × ... no, cells × genes AnnData
    out.obs = pd.DataFrame({
        "cell_type_fine": fine.iloc[keep_idx].to_numpy(),
        "cell_type_broad": broad.iloc[keep_idx].to_numpy(),
        "donor_id": don.iloc[keep_idx].to_numpy()}, index=sub.obs_names)
    out.write_h5ad(DERIVED / "hlca_subset.h5ad")
    print(f"raw-count integer check: {is_int}; wrote subset {out.shape} -> {DERIVED/'hlca_subset.h5ad'}")
    print(f"  donors={out.obs.donor_id.nunique()} broad={out.obs.cell_type_broad.nunique()} "
          f"fine={out.obs.cell_type_fine.nunique()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
