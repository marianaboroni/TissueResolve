#!/usr/bin/env python
"""
Import predictions from an external deconvolution tool run outside TissueResolve.

The benchmark cannot run every external tool automatically, but you can run them
yourself and import their predictions so they are compared as *executed
(imported)* — distinct from export-only tools.

Usage
-----
    python benchmarks/import_external_results.py \\
        --method RCTD --modality spatial \\
        --predictions path/to/rctd_predictions.tsv \\
        --out benchmarks/outputs/spatial/external/RCTD.tsv

Then:
    python benchmarks/run_all.py --use-existing-real-data --include-imported

The predictions file must be a TSV with obs (samples/spots) as rows and
cell types as columns (a header row + an index column).  Rows are renormalised
to sum to 1 on load.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402

from benchmarks.shared.io import OUTPUTS_DIR, write_tsv  # noqa: E402

# External tools whose results can be imported for a fair, executed comparison.
KNOWN_EXTERNAL = {
    "bulk": ["MuSiC", "BisqueRNA", "DWLS", "SCDC", "CIBERSORTx", "BayesPrism"],
    "spatial": ["RCTD", "cell2location", "stereoscope", "SPOTlight", "Tangram",
                "CARD", "DestVI"],
}


def main(argv=None) -> int:
    epilog = ("Supported external tools to import:\n"
              f"  bulk:    {', '.join(KNOWN_EXTERNAL['bulk'])}\n"
              f"  spatial: {', '.join(KNOWN_EXTERNAL['spatial'])}\n"
              "(any --method name is accepted; these are the recognised ones.)")
    ap = argparse.ArgumentParser(
        description=__doc__ + "\n" + epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", required=True, help="Method name, e.g. RCTD")
    ap.add_argument("--modality", required=True, choices=("bulk", "spatial"))
    ap.add_argument("--predictions", required=True, type=Path,
                    help="TSV of obs × cell_type predictions from the external tool")
    ap.add_argument("--spot-id-col", default=None,
                    help="Column holding spot/sample IDs (default: first/index column).")
    ap.add_argument("--out", type=Path, default=None,
                    help="Destination (default: benchmarks/outputs/<modality>/imported/<method>.tsv)")
    args = ap.parse_args(argv)

    if not args.predictions.exists():
        print(f"ERROR: predictions file not found: {args.predictions}", file=sys.stderr)
        return 1
    sep = "\t" if args.predictions.suffix.lower() in (".tsv", ".txt") else ","
    if args.spot_id_col:
        df = pd.read_csv(args.predictions, sep=sep, comment="#").set_index(args.spot_id_col)
    else:
        df = pd.read_csv(args.predictions, sep=sep, index_col=0, comment="#")
    if df.empty or df.shape[1] < 2:
        print("ERROR: predictions must be obs × cell_type with ≥2 cell-type columns.",
              file=sys.stderr)
        return 1
    if args.method not in KNOWN_EXTERNAL.get(args.modality, []):
        print(f"NOTE: {args.method!r} is not in the recognised {args.modality} list "
              f"({', '.join(KNOWN_EXTERNAL[args.modality])}); importing anyway.")
    out = args.out or (OUTPUTS_DIR / args.modality / "imported" / f"{args.method}.tsv")
    write_tsv(df, out, comment=[
        f"imported external predictions: {args.method} ({args.modality})",
        "executed_or_exported: executed_imported",
        f"source: {args.predictions}"])
    print(f"Imported {args.method} → {out} ({df.shape[0]} obs × {df.shape[1]} cell types).")
    print("Re-run: python benchmarks/run_all.py --use-existing-real-data --include-imported")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
