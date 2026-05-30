#!/usr/bin/env python
"""
02 — Generate pseudobulk mixtures with known ground-truth proportions.

Builds >= 12 pseudobulk samples (easy / medium / hard mixing regimes) by
summing raw counts of cells sampled per target composition.  The recorded
ground truth is **mRNA proportions** (count-fraction per cell type), matching
what the bulk deconvolver estimates; cell-count proportions are also recorded
in the metadata.

Usage
-----
    python scripts/02_make_pseudobulk.py
    python scripts/02_make_pseudobulk.py --n-per-regime 5 --seed 1
"""
from __future__ import annotations

import argparse
import sys

import _harness as H


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-h5ad", default=str(H.REFERENCE_H5AD))
    parser.add_argument("--cell-type-col", default=None)
    parser.add_argument("--n-per-regime", type=int, default=4,
                        help="Mixtures per regime (>=4 → >=12 total).")
    parser.add_argument("--n-cells", type=int, default=200,
                        help="Cells sampled per pseudobulk mixture.")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    from pathlib import Path

    ref_path = Path(args.reference_h5ad)
    if not ref_path.exists():
        print(f"ERROR: reference .h5ad not found: {ref_path}\n"
              "Run 00_download_data.py --run-real-data first.", file=sys.stderr)
        return 1

    import anndata as ad

    print(f"Loading single-cell data: {ref_path}")
    adata = ad.read_h5ad(ref_path)
    col = H.detect_cell_type_col(adata.obs, args.cell_type_col)
    print(f"  cell-type column: {col}")

    # Harmonise gene identifiers the SAME way as the reference (script 01), so
    # pseudobulk gene names are symbols matching the saved reference, not the
    # numeric CELLxGENE soma_joinid var_names.
    adata, gene_info = H.harmonize_reference_genes(adata)
    print(f"  gene-id source  : {gene_info['gene_id_source']} "
          f"(duplicate strategy: {gene_info['duplicate_strategy']}, "
          f"{gene_info['n_duplicate_labels']} duplicate symbol(s))")
    print(f"  gene names head : {[str(g) for g in adata.var_names[:5]]}")

    counts_df, true_props_df, meta_df = H.generate_pseudobulk(
        adata, col, n_per_regime=args.n_per_regime,
        n_cells=args.n_cells, seed=args.seed,
    )
    print(f"  generated {counts_df.shape[1]} mixtures × {counts_df.shape[0]} genes")

    H.ensure_dirs()
    H.write_tsv(counts_df, H.PSEUDOBULK_COUNTS,
                comment=["pseudobulk raw counts (genes × samples)",
                         f"seed: {args.seed}", f"cell_type_col: {col}",
                         f"gene_id_source: {gene_info['gene_id_source']}"])
    H.write_tsv(true_props_df, H.PSEUDOBULK_TRUE_PROPS,
                comment=["GROUND TRUTH: mRNA proportions (count-fraction per "
                         "cell type), NOT cell fractions.",
                         "Rows sum to 1."])
    H.write_tsv(meta_df, H.PSEUDOBULK_META,
                comment=["Per-sample regime, counts, dominant type, and "
                         "cell-count proportions (cellprop_*)."])
    print(f"Saved pseudobulk -> {H.DERIVED_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
