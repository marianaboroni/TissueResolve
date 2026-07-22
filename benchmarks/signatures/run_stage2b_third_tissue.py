#!/usr/bin/env python
"""Stage 2B — third-tissue slot for the calibrated identifiability certificate.

Status: the calibrated certificate is tissue-agnostic (it takes any ReferenceSignature; the synthetic
tests exercise arbitrary references). A *real* third tissue could not be executed this stage because
no donor-annotated, cell-typed single-cell reference is available offline:

  * CRC (GSE200997) — the only candidate present — annotates cells by CONSENSUS MOLECULAR SUBTYPE
    (CMS1-4), NOT cell type, and ~40% of cells are unlabelled. This script verifies that fact.
  * scanpy's bundled pbmc68k_reduced is a single donor (no donor-disjoint split possible).

This script (a) verifies the CRC finding, (b) if the user drops a prepared third-tissue reference at
the documented path, builds it and runs the calibrated certificate (proving the pipeline generalises),
otherwise (c) prints exact, reproducible instructions and exits cleanly with NOT EXECUTED. It never
downloads silently and never fails silently.

Usage:
  PYTHONPATH=src:. python benchmarks/signatures/run_stage2b_third_tissue.py [--run-real-data]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

CRC_ANNOT = REPO / "examples" / "second_tissue_crc" / "data" / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
THIRD_REF = REPO / "examples" / "third_tissue" / "data" / "reference.h5ad"   # documented drop-in path
DONOR_CANDIDATES = ("donor_id", "donor", "samples", "subject", "patient", "individual")
CELLTYPE_CANDIDATES = ("cell_type", "cell_type_fine", "celltype", "annotation", "cell_ontology_class")


def audit_crc():
    if not CRC_ANNOT.exists():
        print(f"[crc] annotation not found at {CRC_ANNOT}")
        return
    a = pd.read_csv(CRC_ANNOT, index_col=0)
    pred = a["prediction"].astype(str) if "prediction" in a else pd.Series(dtype=str)
    vals = sorted(v for v in pred.unique() if v != "nan")
    cms = all(str(v).upper().startswith("CMS") for v in vals) if vals else False
    print(f"[crc] prediction labels = {vals}  (all CMS molecular subtypes = {cms}); "
          f"nan fraction = {(pred == 'nan').mean():.2f}")
    print("[crc] VERDICT: molecular-subtype annotation, NOT cell type -> unusable as a deconvolution "
          "reference. (This is the precise reason the third tissue is NOT EXECUTED.)")


def try_third_reference(run_real):
    if not THIRD_REF.exists():
        print(f"\n[third-tissue] no prepared reference at {THIRD_REF}")
        print("[third-tissue] TO ENABLE (reproducible, no code change needed):")
        print("  1. Obtain a multi-donor, cell-typed scRNA-seq reference (stable identifier), e.g.")
        print("     Tabula Sapiens (a tissue with >=3 donors) or a multi-donor PBMC atlas.")
        print(f"  2. Save it as an AnnData .h5ad at: {THIRD_REF}")
        print(f"     with obs columns for donor (one of {DONOR_CANDIDATES}) and")
        print(f"     cell type (one of {CELLTYPE_CANDIDATES}); raw counts in .X or .layers['counts'].")
        print("  3. Re-run this script; it will build the reference and run the calibrated certificate.")
        print("  4. Record a download_manifest.json (dataset, source, ID, date, versions, filters).")
        print("\n[third-tissue] NOT EXECUTED — real third tissue requires the data above.")
        return 0

    if not run_real:
        print(f"[third-tissue] reference present at {THIRD_REF}; pass --run-real-data to run.")
        return 0

    import anndata as adread
    import numpy as np
    import warnings
    import tissueresolve as tr
    from tissueresolve.reference.identifiability_calibration import calibrated_identifiability_certificate as ccert

    ad = adread.read_h5ad(THIRD_REF)
    donor = next((c for c in DONOR_CANDIDATES if c in ad.obs), None)
    ctc = next((c for c in CELLTYPE_CANDIDATES if c in ad.obs), None)
    if donor is None or ctc is None:
        print(f"[third-tissue] reference lacks a donor/cell-type column "
              f"(have {list(ad.obs.columns)}) — cannot run. NOT EXECUTED.")
        return 0
    print(f"[third-tissue] building reference (donor={donor}, cell_type={ctc}, "
          f"{ad.n_obs} cells, {ad.obs[donor].nunique()} donors, {ad.obs[ctc].nunique()} types)")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = tr.build_reference(ad, cell_type_col=ctc)
    genes = list(ref.gene_names)
    panel = list(np.random.default_rng(12345).choice(genes, size=min(1000, len(genes)), replace=False))
    cert = ccert(ref, query_detectable_genes=panel, library_size=1e6, n_sim=200, seed=0)
    out = REPO / "benchmarks" / "results" / "signatures" / "stage2b" / "third_tissue"
    cert.save(out)
    print("[third-tissue] class counts:", cert.per_type.recoverability.value_counts().to_dict())
    print("[third-tissue] confusable clusters:", len(cert.clusters))
    print(f"[third-tissue] EXECUTED — wrote {out}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    audit_crc()
    return try_third_reference(args.run_real_data)


if __name__ == "__main__":
    raise SystemExit(main())
