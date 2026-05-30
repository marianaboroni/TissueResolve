#!/usr/bin/env python
"""
01 — Build the TissueResolve reference from the single-cell breast-cancer atlas.

Loads the downloaded reference .h5ad, auto-detects the cell-type column, drops
cell types with too few cells, builds the reference with NB overdispersion
(so it serves both bulk and spatial), and writes the saved reference object
plus summary tables.

Usage
-----
    python scripts/01_prepare_reference.py
    python scripts/01_prepare_reference.py --cell-type-col cell_type_major
    python scripts/01_prepare_reference.py --min-cells 50
"""
from __future__ import annotations

import argparse
import sys

import _harness as H


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-h5ad", default=str(H.REFERENCE_H5AD))
    parser.add_argument("--cell-type-col", default=None,
                        help="Override the auto-detected cell-type column.")
    parser.add_argument("--min-cells", type=int, default=H.MIN_CELLS_PER_TYPE,
                        help="Drop cell types with fewer than this many cells.")
    args = parser.parse_args(argv)

    from pathlib import Path

    ref_path = Path(args.reference_h5ad)
    if not ref_path.exists():
        print(f"ERROR: reference not found: {ref_path}\n"
              "Run 00_download_data.py --run-real-data first (or place the "
              "file manually).", file=sys.stderr)
        return 1

    import anndata as ad

    print(f"Loading reference: {ref_path}")
    adata = ad.read_h5ad(ref_path)
    print(f"  {adata.n_obs} cells × {adata.n_vars} genes")

    prep = H.prepare_reference(
        adata, cell_type_col=args.cell_type_col, min_cells=args.min_cells,
        estimate_overdispersion=True,
    )
    print(f"  cell-type column: {prep.cell_type_col}  (source: {prep.counts_source})")
    print(f"  reference: {prep.reference.n_genes} genes × "
          f"{prep.reference.n_cell_types} cell types")

    H.ensure_dirs()
    prep.reference.save(H.SAVED_REFERENCE_DIR)
    H.write_tsv(prep.summary.set_index("metric"),
                H.OUT_REFERENCE_DIR / "reference_summary.tsv")
    H.write_tsv(prep.cell_type_counts.to_frame("n_cells"),
                H.OUT_REFERENCE_DIR / "cell_type_counts.tsv")
    (H.OUT_REFERENCE_DIR / "selected_annotation_column.txt").write_text(
        prep.cell_type_col + "\n", encoding="utf-8"
    )
    print(f"Saved reference object -> {H.SAVED_REFERENCE_DIR}")
    print(f"Saved summaries        -> {H.OUT_REFERENCE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
