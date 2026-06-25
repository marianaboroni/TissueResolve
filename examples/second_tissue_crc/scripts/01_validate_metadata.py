#!/usr/bin/env python
"""STAGE 3 — audit GSE200997 metadata & labels before building a reference.

Audits (does NOT assume the published annotations are usable; rule 6): raw-count
integrity, donor identifiers, broad/fine label columns, cells-per-label,
donors-per-label, low-support labels (flagged, not silently dropped). Proposes a
fine→broad mapping for review (documented, not auto-merged; rule 9).

Outputs (examples/second_tissue_crc/data/derived/):
  crc_metadata_audit.tsv, crc_label_support.tsv, crc_proposed_hierarchy.tsv
  + docs note appended by the caller.
Usage:  python 01_validate_metadata.py --run-real-data
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data"
DERIVED = DATA / "derived"
ANNOT = DATA / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
COUNTS = DATA / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
MIN_CELLS, MIN_DONORS = 50, 3


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    if not ANNOT.exists():
        print(f"Missing {ANNOT}. Run 00_download_data.py first.", file=sys.stderr); return 2
    DERIVED.mkdir(parents=True, exist_ok=True)

    annot = pd.read_csv(ANNOT)
    print(f"annotation: {annot.shape} | columns: {list(annot.columns)}", flush=True)
    # report candidate columns for donor / cell-type (broad/fine)
    audit = []
    for c in annot.columns:
        nun = annot[c].nunique()
        kind = ("donor?" if any(k in c.lower() for k in ("patient", "sample", "donor", "subject"))
                else "label?" if any(k in c.lower() for k in ("type", "cluster", "cell", "anno", "annot", "ident", "lineage"))
                else "other")
        audit.append({"column": c, "n_unique": int(nun), "kind_guess": kind,
                      "examples": "; ".join(map(str, annot[c].dropna().unique()[:6]))})
    adf = pd.DataFrame(audit)
    adf.to_csv(DERIVED / "crc_metadata_audit.tsv", sep="\t", index=False)
    print("=== column audit ===\n" + adf.to_string(index=False), flush=True)

    # heuristic pick: donor col = a sample/patient-like column with 5..100 uniques;
    # fine label = the cell-type-like column with the MOST uniques; broad = fewer
    donor_cands = adf[(adf.kind_guess == "donor?") & (adf.n_unique.between(2, 200))]
    label_cands = adf[adf.kind_guess == "label?"].sort_values("n_unique", ascending=False)
    donor_col = donor_cands.sort_values("n_unique").iloc[0]["column"] if len(donor_cands) else None
    fine_col = label_cands.iloc[0]["column"] if len(label_cands) else None
    broad_col = label_cands.iloc[-1]["column"] if len(label_cands) >= 2 else None
    print(f"\nheuristic picks → donor={donor_col} fine={fine_col} broad={broad_col}")

    # label support (cells + donors per label)
    if fine_col and donor_col:
        sup = pd.DataFrame({
            "n_cells": annot.groupby(fine_col).size(),
            "n_donors": annot.groupby(fine_col)[donor_col].nunique()})
        sup["low_support"] = (sup["n_cells"] < MIN_CELLS) | (sup["n_donors"] < MIN_DONORS)
        sup = sup.sort_values("n_cells")
        sup.to_csv(DERIVED / "crc_label_support.tsv", sep="\t")
        print(f"\n=== label support ({fine_col}) — {int(sup['low_support'].sum())} low-support of {len(sup)} ===")
        print(sup.to_string())
        # proposed hierarchy (broad if a broad col exists, else fine==broad placeholder for review)
        if broad_col and broad_col != fine_col:
            hmap = annot[[fine_col, broad_col]].drop_duplicates().rename(
                columns={fine_col: "fine_cell_type", broad_col: "broad_cell_type"})
        else:
            hmap = pd.DataFrame({"fine_cell_type": sup.index,
                                 "broad_cell_type": "REVIEW_REQUIRED"})
        hmap.to_csv(DERIVED / "crc_proposed_hierarchy.tsv", sep="\t", index=False)
        print(f"\nproposed hierarchy -> {DERIVED/'crc_proposed_hierarchy.tsv'} "
              f"({'from ' + broad_col if broad_col and broad_col!=fine_col else 'REVIEW REQUIRED'})")

    # raw-count integrity (peek header + a few rows without loading all)
    with gzip.open(COUNTS, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        first = fh.readline().rstrip("\n").split(",")
    n_cols = len(header)
    vals = first[1:21]
    is_int = all(v.strip().lstrip("-").isdigit() or v.strip() == "" for v in vals if v != "")
    print(f"\n=== counts integrity ===\nmatrix columns (cells+1): {n_cols}; first gene row "
          f"'{first[0]}' looks integer-count: {is_int} (sample vals: {vals[:6]})")
    print(f"\nWrote audit tables -> {DERIVED}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
