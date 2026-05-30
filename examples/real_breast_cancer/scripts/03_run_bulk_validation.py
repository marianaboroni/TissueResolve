#!/usr/bin/env python
"""
03 — Bulk validation: deconvolve pseudobulk mixtures and score vs ground truth.

Loads the pseudobulk counts and the saved TissueResolve reference, explicitly
intersects genes (reporting overlap), runs the real bulk pipeline, and compares
estimated mRNA proportions against the known mRNA-proportion ground truth.

Usage
-----
    python scripts/03_run_bulk_validation.py
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

import _harness as H


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-bootstrap", type=int, default=0,
                        help="Bootstrap iterations (0 = skip for speed).")
    args = parser.parse_args(argv)

    from pathlib import Path

    for required in (H.PSEUDOBULK_COUNTS, H.PSEUDOBULK_TRUE_PROPS):
        if not Path(required).exists():
            print(f"ERROR: {required} not found.  Run 02_make_pseudobulk.py first.",
                  file=sys.stderr)
            return 1
    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"ERROR: reference not found at {H.SAVED_REFERENCE_DIR}.  "
              "Run 01_prepare_reference.py first.", file=sys.stderr)
        return 1

    from tissueresolve.io.bulk import read_bulk_counts
    from tissueresolve.results import ReferenceSignature
    import tissueresolve as tr

    bulk = read_bulk_counts(H.PSEUDOBULK_COUNTS)            # genes × samples
    true_props = pd.read_csv(H.PSEUDOBULK_TRUE_PROPS, sep="\t",
                             comment="#", index_col=0)
    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)

    overlap = H.gene_overlap(bulk.index, ref.gene_names)
    print(f"Gene overlap (bulk vs reference): {overlap}")
    shared = [g for g in bulk.index if g in set(ref.gene_names)]
    if not shared:
        print("ERROR: no shared genes between pseudobulk and reference.",
              file=sys.stderr)
        return 1
    bulk = bulk.loc[shared]

    print(f"Running bulk deconvolution on {bulk.shape[1]} samples …")
    result = tr.deconv_bulk(bulk, ref, n_bootstrap=args.n_bootstrap)
    est = result.deconv.proportions                         # samples × cell types

    metrics = H.compute_bulk_metrics(true_props, est)
    per_ct = H.per_celltype_metrics(true_props, est)
    print(f"Overall: pearson={metrics['pearson']:.3f} "
          f"spearman={metrics['spearman']:.3f} rmse={metrics['rmse']:.3f} "
          f"mae={metrics['mae']:.3f}")

    # QC table
    qc = result.qc
    qc_rows = {}
    if getattr(qc, "recon_r2", None) is not None:
        qc_rows["recon_r2"] = qc.recon_r2
    if getattr(qc, "profile_corr", None) is not None:
        qc_rows["profile_corr"] = qc.profile_corr
    qc_df = pd.DataFrame(qc_rows) if qc_rows else pd.DataFrame(
        {"coverage_r2": result.deconv.coverage_r2})

    warnings = {
        "estimate_type": result.deconv.ESTIMATE_TYPE,
        "estimate_type_note": "Bulk values are mRNA proportions, not cell fractions.",
        "gene_overlap": overlap,
        "qc_recommendations": list(getattr(qc, "recommendations", []) or []),
    }

    H.ensure_dirs()
    H.write_tsv(est, H.OUT_BULK_DIR / "bulk_estimated_proportions.tsv",
                comment=["estimate_type: mRNA_proportion (NOT cell fractions)"])
    H.write_tsv(pd.DataFrame([metrics]).T.rename(columns={0: "value"}),
                H.OUT_BULK_DIR / "bulk_validation_metrics.tsv")
    H.write_tsv(per_ct, H.OUT_BULK_DIR / "bulk_per_celltype_metrics.tsv")
    H.write_tsv(qc_df, H.OUT_BULK_DIR / "bulk_qc.tsv")
    H.write_json(warnings, H.OUT_BULK_DIR / "bulk_warnings.json")
    print(f"Saved bulk validation outputs -> {H.OUT_BULK_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
