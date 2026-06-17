#!/usr/bin/env python
"""Download the OPEN colorectal-cancer second-tissue dataset (GSE200997).

Committee-approved open substitute for the credential-gated Pelka atlas (see
docs/PELKA_CRC_DATASET_VERIFICATION.md). Open access (CC BY 4.0), raw UMI counts +
cell annotations + donor IDs. Resumable, checksummed, no credentials, --dry-run /
--force. Writes a download_manifest.json. Data dir is git-ignored (rule 11).

Usage:
    python 00_download_data.py --run-real-data        # download
    python 00_download_data.py --dry-run              # show plan only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/suppl"
FILES = {
    "counts": "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz",
    "annotation": "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz",
}
DATA = Path(__file__).resolve().parents[1] / "data"
MANIFEST = DATA / "download_manifest.json"
DATASET = {
    "name": "Colorectal cancer single-cell atlas (open substitute for Pelka)",
    "accession": "GSE200997", "repository": "NCBI GEO",
    "publication": "Refining CRC Classification through a Single-Cell Atlas (HCA project 4d9d56e4)",
    "license": "CC BY 4.0", "protocol": "10x Genomics scRNA-seq", "counts": "raw UMI",
    "approx_cells": 49589, "approx_samples": 23,
    "note": "Open substitute (committee-approved). NOT the 62-donor Pelka Cell 2021 atlas; "
            "smaller cohort; fine-annotation depth audited in 01_validate_metadata.py.",
}


def _sha256(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(buf):
            h.update(chunk)
    return h.hexdigest()


def _download(url, dest, force):
    if dest.exists() and not force:
        print(f"  exists, skip: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return False
    tmp = dest.with_suffix(dest.suffix + ".part")
    # resumable
    pos = tmp.stat().st_size if tmp.exists() else 0
    req = urllib.request.Request(url, headers={"Range": f"bytes={pos}-"} if pos else {})
    print(f"  downloading {dest.name} (from byte {pos}) ...", flush=True)
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "ab") as out:
        while chunk := r.read(1 << 20):
            out.write(chunk)
    tmp.rename(dest)
    print(f"  done: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    plan = {k: f"{BASE}/{v}" for k, v in FILES.items()}
    if args.dry_run or not args.run_real_data:
        print("DRY RUN — would download:")
        for k, u in plan.items():
            print(f"  {k}: {u}")
        print(f"into {DATA}/  (git-ignored). Re-run with --run-real-data to download.")
        return 0
    DATA.mkdir(parents=True, exist_ok=True)
    rec = {"dataset": DATASET, "download_date": time.strftime("%Y-%m-%d"),
           "dataset_version": "GEO 2022-04-18", "files": {}}
    for key, name in FILES.items():
        dest = DATA / name
        for attempt in range(3):
            try:
                _download(f"{BASE}/{name}", dest, args.force)
                break
            except Exception as exc:  # noqa: BLE001
                print(f"  attempt {attempt+1} failed: {exc}", file=sys.stderr); time.sleep(3)
        else:
            print(f"FAILED to download {name} after retries.", file=sys.stderr); return 1
        rec["files"][key] = {"file_name": name, "url": f"{BASE}/{name}",
                             "local_path": str(dest), "size_bytes": dest.stat().st_size,
                             "sha256": _sha256(dest), "processing_status": "downloaded"}
    MANIFEST.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(f"\nManifest -> {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
