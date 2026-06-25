#!/usr/bin/env python
"""
Download real **bulk TCGA-BRCA TNBC** RNA-seq (STAR gene counts) from the GDC.

This is real-data validation tooling and runs **only** when explicitly invoked
(``--run-real-data`` or ``TISSUERESOLVE_RUN_REAL_DATA=1``).  It never runs as
part of the default offline test suite, and downloaded data is git-ignored.

TNBC selection (reproducible, from IHC receptor status)
-------------------------------------------------------
1. Download the TCGA-BRCA clinical Biotab supplement
   (``nationwidechildrens.org_clinical_patient_brca.txt``) via the GDC API and
   select **triple-negative** patients: ``er_status_by_ihc == 'Negative'`` and
   ``pr_status_by_ihc == 'Negative'`` and ``her2_status_by_ihc == 'Negative'``.
2. Query the GDC files endpoint for those patients' **STAR - Counts**
   ``Gene Expression Quantification`` files on **Primary Tumor** samples
   (one file per patient).
3. Download each file, take the ``unstranded`` counts column keyed by
   ``gene_name`` (gene symbols, matching the single-cell reference), aggregate
   duplicate symbols by summing, and assemble a genes × samples matrix.

Outputs (under ``data/bulk_tcga_tnbc/``)
----------------------------------------
- ``tcga_tnbc_counts.tsv``         genes × samples raw counts
- ``tcga_tnbc_metadata.tsv``       per-sample patient barcode, file id, sample type
- ``tcga_tnbc_download.json``      provenance (selection, file list, versions)
and a ``bulk_tcga_tnbc`` entry is added to ``data/download_manifest.json``.

**No ground truth.** Real bulk tumours have no known per-cell-type proportions,
so this dataset supports a *concordance / plausibility / robustness* benchmark,
**not** an accuracy benchmark (use the pseudobulk harness for ground-truth
accuracy).

Usage
-----
    python scripts/00b_download_tcga_tnbc_bulk.py --dry-run
    python scripts/00b_download_tcga_tnbc_bulk.py --run-real-data            # all 116 TNBC samples
    python scripts/00b_download_tcga_tnbc_bulk.py --run-real-data --max-samples 40
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import tarfile  # noqa: F401  (kept for optional batch-tar support)
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _harness as H  # noqa: E402

GDC_API = "https://api.gdc.cancer.gov"
CLINICAL_FILE_NAME = "nationwidechildrens.org_clinical_patient_brca.txt"

BULK_TNBC_DIR = H.DATA_DIR / "bulk_tcga_tnbc"
COUNTS_TSV = BULK_TNBC_DIR / "tcga_tnbc_counts.tsv"
META_TSV = BULK_TNBC_DIR / "tcga_tnbc_metadata.tsv"
PROVENANCE_JSON = BULK_TNBC_DIR / "tcga_tnbc_download.json"


def _gdc_post(endpoint: str, payload: dict, timeout: int = 60) -> dict:
    req = urllib.request.Request(
        f"{GDC_API}/{endpoint}", data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "tissueresolve"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def _gdc_get_bytes(path: str, timeout: int = 120, retries: int = 3) -> bytes:
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(f"{GDC_API}/{path}",
                                         headers={"User-Agent": "tissueresolve"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * attempt)
    raise RuntimeError(f"GDC download failed for {path!r} after {retries} tries: {last}")


def find_tnbc_barcodes() -> tuple[list[str], dict]:
    """Return (TNBC patient barcodes, selection summary) from the clinical Biotab."""
    clin = _gdc_post("files", {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Clinical Supplement"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["BCR Biotab"]}},
        ]},
        "fields": "file_id,file_name", "size": "50", "format": "JSON"})
    hits = [h for h in clin["data"]["hits"] if h["file_name"] == CLINICAL_FILE_NAME]
    if not hits:
        raise RuntimeError(
            "Could not locate the TCGA-BRCA clinical patient Biotab "
            f"({CLINICAL_FILE_NAME}) via the GDC API.  The GDC schema may have "
            "changed; check https://portal.gdc.cancer.gov/.")
    fid = hits[0]["file_id"]
    raw = _gdc_get_bytes(f"data/{fid}")
    df = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str).iloc[2:]  # drop 2 description rows
    needed = {"bcr_patient_barcode", "er_status_by_ihc",
              "pr_status_by_ihc", "her2_status_by_ihc"}
    if not needed <= set(df.columns):
        raise RuntimeError(f"clinical file missing receptor columns; has {list(df.columns)[:20]}")
    tnbc = df[(df.er_status_by_ihc == "Negative")
              & (df.pr_status_by_ihc == "Negative")
              & (df.her2_status_by_ihc == "Negative")]
    barcodes = sorted(tnbc["bcr_patient_barcode"].dropna().unique().tolist())
    summary = {
        "clinical_file_id": fid,
        "n_patients_total": int(df["bcr_patient_barcode"].nunique()),
        "n_tnbc": len(barcodes),
        "criteria": "er_status_by_ihc==Negative AND pr_status_by_ihc==Negative "
                    "AND her2_status_by_ihc==Negative",
    }
    return barcodes, summary


def find_rnaseq_files(barcodes: list[str]) -> pd.DataFrame:
    """One STAR-Counts primary-tumor RNA-seq file per TNBC patient."""
    res = _gdc_post("files", {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.submitter_id", "value": barcodes}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
        ]},
        "fields": "file_id,file_name,file_size,cases.submitter_id,"
                  "cases.samples.sample_type",
        "size": "2000", "format": "JSON"})
    rows = []
    for h in res["data"]["hits"]:
        case = (h.get("cases") or [{}])[0]
        patient = case.get("submitter_id", "")
        sample_type = ((case.get("samples") or [{}])[0]).get("sample_type", "")
        rows.append({"patient_barcode": patient, "file_id": h["file_id"],
                     "file_name": h["file_name"], "file_size": int(h.get("file_size", 0)),
                     "sample_type": sample_type})
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No TNBC RNA-seq STAR-Counts files returned by the GDC.")
    # one file per patient (deterministic: smallest file_id)
    df = df.sort_values(["patient_barcode", "file_id"]).drop_duplicates("patient_barcode")
    return df.reset_index(drop=True)


def _parse_star_counts(raw: bytes) -> "pd.Series":
    """gene_name → unstranded counts (duplicate symbols summed).

    The GDC ``augmented_star_gene_counts.tsv`` begins with a ``# gene-model:``
    comment line, then the header (``gene_id  gene_name  gene_type  unstranded
    …``), then four ``N_*`` summary rows before the per-gene rows.
    """
    df = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, comment="#")
    # drop the N_* summary rows (gene_id like N_unmapped)
    df = df[~df["gene_id"].astype(str).str.startswith("N_")]
    df = df[df["gene_name"].notna() & (df["gene_name"].astype(str) != "")]
    counts = pd.to_numeric(df["unstranded"], errors="coerce").fillna(0.0)
    s = pd.Series(counts.to_numpy(), index=df["gene_name"].astype(str))
    return s.groupby(level=0).sum()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-real-data", action="store_true",
                    help="Required to actually download (else dry-run plan only).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Show the TNBC selection + file count, download nothing.")
    ap.add_argument("--max-samples", type=int, default=None,
                    help="Cap the number of TNBC samples downloaded (logged).")
    args = ap.parse_args(argv)

    if not (H.real_data_enabled(args.run_real_data) or args.dry_run):
        print("Refusing to download without --run-real-data (or "
              f"{H.REAL_DATA_ENV}=1).  Use --dry-run to preview the plan.",
              file=sys.stderr)
        return 2

    print("Selecting TNBC patients from TCGA-BRCA clinical IHC status …")
    barcodes, sel = find_tnbc_barcodes()
    print(f"  {sel['n_tnbc']} TNBC patients of {sel['n_patients_total']} "
          f"({sel['criteria']}).")
    files = find_rnaseq_files(barcodes)
    total_mb = files["file_size"].sum() / 1e6
    print(f"  {len(files)} STAR-Counts primary-tumor RNA-seq files "
          f"(~{total_mb:.0f} MB total).")

    n_avail = len(files)
    if args.max_samples is not None and args.max_samples < n_avail:
        print(f"  NOTE: --max-samples={args.max_samples} → downloading "
              f"{args.max_samples} of {n_avail} TNBC samples; "
              f"{n_avail - args.max_samples} dropped.  Omit --max-samples for "
              "the full cohort.")
        files = files.head(args.max_samples).reset_index(drop=True)

    if args.dry_run:
        print(f"Dry run: would download {len(files)} files into {BULK_TNBC_DIR}/.")
        return 0

    BULK_TNBC_DIR.mkdir(parents=True, exist_ok=True)
    series_by_sample: dict[str, pd.Series] = {}
    meta_rows = []
    for i, row in files.iterrows():
        sid = row["patient_barcode"]
        print(f"  [{i+1}/{len(files)}] {sid}  {row['file_name'][:40]}…", flush=True)
        raw = _gdc_get_bytes(f"data/{row['file_id']}")
        series_by_sample[sid] = _parse_star_counts(raw)
        meta_rows.append({"sample": sid, "patient_barcode": sid,
                          "file_id": row["file_id"], "file_name": row["file_name"],
                          "sample_type": row["sample_type"]})

    counts = pd.DataFrame(series_by_sample).fillna(0.0)
    counts.index.name = "gene"
    counts = counts.round().astype("int64")
    H.write_tsv(counts, COUNTS_TSV, comment=[
        "TCGA-BRCA TNBC bulk RNA-seq (STAR - Counts, 'unstranded'); genes × samples.",
        "Gene identifiers are HGNC symbols (gene_name); duplicates summed.",
        "REAL bulk tumour — NO ground-truth cell proportions (concordance/QC only).",
    ])
    meta = pd.DataFrame(meta_rows).set_index("sample")
    H.write_tsv(meta, META_TSV)

    provenance = {
        "dataset_name": "TCGA-BRCA TNBC bulk RNA-seq (STAR gene counts)",
        "source": "NCI Genomic Data Commons (GDC) API",
        "selection": sel,
        "n_samples_downloaded": int(counts.shape[1]),
        "n_samples_available": int(n_avail),
        "max_samples": args.max_samples,
        "n_genes": int(counts.shape[0]),
        "files": meta_rows,
        "package_versions": H.package_versions(),
        "access_date": datetime.now(timezone.utc).isoformat(),
        "ground_truth": False,
    }
    H.write_json(provenance, PROVENANCE_JSON)

    _update_manifest(provenance)
    print(f"\nWrote {counts.shape[0]} genes × {counts.shape[1]} samples → {COUNTS_TSV}")
    print(f"Provenance → {PROVENANCE_JSON}")
    return 0


def _update_manifest(provenance: dict) -> None:
    """Add/refresh the ``bulk_tcga_tnbc`` dataset entry in the harness manifest."""
    try:
        manifest = H.read_manifest()
    except Exception:
        manifest = H.default_manifest()
    entry = {
        "dataset_name": provenance["dataset_name"],
        "source": provenance["source"],
        "url_or_id": "GDC project TCGA-BRCA; STAR - Counts; Primary Tumor; "
                     "TNBC by IHC (ER-/PR-/HER2-)",
        "access_date": provenance["access_date"],
        "package_versions": provenance["package_versions"],
        "filters_applied": [provenance["selection"]["criteria"],
                            "STAR - Counts", "Primary Tumor"],
        "downsampling": {"max_samples": provenance["max_samples"],
                         "n_downloaded": provenance["n_samples_downloaded"],
                         "n_available": provenance["n_samples_available"],
                         "applied": provenance["max_samples"] is not None},
        "file": str(COUNTS_TSV.relative_to(H.HARNESS_DIR)),
        "downloaded": True,
        "ground_truth": False,
    }
    manifest.setdefault("datasets", {})["bulk_tcga_tnbc"] = entry
    # write directly (validate_manifest only requires reference+spatial; this
    # extra dataset entry is additive provenance).
    H.MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
