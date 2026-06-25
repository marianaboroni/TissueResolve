#!/usr/bin/env python
"""Download the OPEN + LABELED second tissue: Human Lung Cell Atlas (HLCA) core.

Open (CZ CELLxGENE Discover), raw counts at adata.raw.X, manual consensus broad→
fine annotations (ann_level_1..5), ~107 donors. Distinct tissue (healthy lung) with
epithelial/immune/stromal/endothelial compartments — for testing whether the breast
identifiability/shared-lineage pattern generalises.

Resumable, checksummed, no credentials, --dry-run/--force. Data git-ignored (rule 11).

Source: CELLxGENE collection 6f6d381a-7701-4781-935c-db10d30de293,
dataset 066943a2-fdac-4b29-b348-40cede398e4e ("...(core)").
Usage:  python 00_download_hlca.py --run-real-data
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

URL = "https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"
DATA = Path(__file__).resolve().parents[1] / "data"
DEST = DATA / "hlca_core.h5ad"
MANIFEST = DATA / "download_manifest.json"
META = {
    "name": "Human Lung Cell Atlas (HLCA) core", "tissue": "lung (healthy)",
    "source": "CZ CELLxGENE Discover", "collection": "6f6d381a-7701-4781-935c-db10d30de293",
    "dataset_id": "066943a2-fdac-4b29-b348-40cede398e4e",
    "publication": "Sikkema et al. 2023, Nat Med (HLCA)",
    "approx_cells": 584944, "approx_donors": 107,
    "counts": "raw at adata.raw.X", "labels": "ann_level_1..5 (broad→fine) + cell_type ontology",
    "license": "CC BY 4.0 (CELLxGENE; verify on collection page)",
}


def _sha256(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while c := fh.read(buf):
            h.update(c)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    if args.dry_run or not args.run_real_data:
        print(f"DRY RUN — would download HLCA core (~5.9 GB):\n  {URL}\n  -> {DEST} (git-ignored)")
        return 0
    DATA.mkdir(parents=True, exist_ok=True)
    # expected total size (for completeness verification — the file was truncated once)
    with urllib.request.urlopen(urllib.request.Request(URL, method="HEAD"), timeout=60) as h:
        expected = int(h.headers.get("Content-Length", "0"))
    tmp = DEST.with_suffix(".h5ad.part")
    if DEST.exists() and not args.force:
        if DEST.stat().st_size == expected:
            print(f"exists & complete, skip: {DEST} ({DEST.stat().st_size/1e9:.2f} GB)")
        else:  # previously-renamed truncated file → resume
            print(f"incomplete prior file ({DEST.stat().st_size/1e9:.2f} GB ≠ "
                  f"{expected/1e9:.2f} GB) — resuming."); DEST.rename(tmp)
    if not (DEST.exists() and DEST.stat().st_size == expected):
        t0 = time.time()
        for attempt in range(20):                          # resume loop across dropped connections
            pos = tmp.stat().st_size if tmp.exists() else 0
            if expected and pos >= expected:
                break
            try:
                req = urllib.request.Request(URL, headers={"Range": f"bytes={pos}-"} if pos else {})
                with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "ab") as out:
                    while chunk := r.read(1 << 22):
                        out.write(chunk)
            except Exception as exc:  # noqa: BLE001
                print(f"  drop at {tmp.stat().st_size/1e9:.2f} GB ({exc}); retry {attempt+1}", flush=True)
                time.sleep(3)
        size = tmp.stat().st_size if tmp.exists() else 0
        if expected and size != expected:
            print(f"FAILED: got {size} of {expected} bytes after retries.", file=sys.stderr); return 1
        tmp.rename(DEST)
        print(f"done: {DEST} ({DEST.stat().st_size/1e9:.2f} GB, {(time.time()-t0):.0f}s)")
    rec = {"dataset": META, "download_date": time.strftime("%Y-%m-%d"), "url": URL,
           "local_path": str(DEST), "size_bytes": DEST.stat().st_size,
           "sha256": _sha256(DEST), "processing_status": "downloaded"}
    MANIFEST.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(f"manifest -> {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
