#!/usr/bin/env python
"""
04 — Spatial validation on the real 10x Visium breast-cancer section.

Loads the Visium .h5ad (using ``layers['counts']`` when present), reuses the
same TissueResolve reference as the bulk validation, explicitly intersects
genes, builds the spatial graph from array coordinates, and runs the spatial
pipeline.  Visium has no per-spot ground-truth composition, so this validation
is qualitative: spatial coherence (Moran's I), QC, and surfaced warnings.

If the spatial workflow (Stage 4) is unavailable, the script exits with a
clear non-zero message rather than failing obscurely.

Usage
-----
    python scripts/04_run_spatial_validation.py
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

import _harness as H


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial-h5ad", default=str(H.SPATIAL_H5AD))
    parser.add_argument("--neighbourhood", action="store_true",
                        help="Also compute neighbourhood co-occurrence stats.")
    args = parser.parse_args(argv)

    from pathlib import Path

    # Stage-4 availability check (clear message if missing).
    try:
        import tissueresolve as tr  # noqa: F401
        from tissueresolve.results import ReferenceSignature  # noqa: F401
        import tissueresolve.spatial.pipeline  # noqa: F401
    except Exception as exc:
        print("ERROR: the spatial workflow (Stage 4) is required but could not "
              f"be imported: {exc}", file=sys.stderr)
        return 2

    sp_path = Path(args.spatial_h5ad)
    if not sp_path.exists():
        print(f"ERROR: Visium .h5ad not found: {sp_path}\n"
              "Run 00_download_data.py --run-real-data first.", file=sys.stderr)
        return 1
    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"ERROR: reference not found at {H.SAVED_REFERENCE_DIR}.  "
              "Run 01_prepare_reference.py first.", file=sys.stderr)
        return 1

    import anndata as ad
    import scipy.sparse as sp
    import pandas as pd
    from tissueresolve.results import ReferenceSignature
    import tissueresolve as tr

    print(f"Loading Visium: {sp_path}")
    adata = ad.read_h5ad(sp_path)
    print(f"  {adata.n_obs} spots × {adata.n_vars} genes")

    # counts: prefer layers['counts']
    mat, source = H.resolve_counts_matrix(adata)
    print(f"  counts source: {source}")
    Y = mat if sp.issparse(mat) else np.asarray(mat)

    # array coordinates (required by the hex-graph builder)
    obs = adata.obs
    if "array_row" in obs.columns and "array_col" in obs.columns:
        array_row = obs["array_row"].to_numpy()
        array_col = obs["array_col"].to_numpy()
    else:
        print("ERROR: array_row / array_col not found in obs.  The hex graph "
              "needs Visium array coordinates; obsm['spatial'] holds pixel "
              "coordinates only.  Re-export the Visium object with array "
              "coordinates (load_visium attaches them).", file=sys.stderr)
        return 1

    lib_sizes = (np.asarray(Y.sum(axis=1)).ravel()
                 if sp.issparse(Y) else Y.sum(axis=1)).astype("float32")
    gene_names = list(adata.var_names)
    spot_ids = list(adata.obs_names)

    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    overlap = H.gene_overlap(gene_names, ref.gene_names)
    print(f"Gene overlap (spatial vs reference): {overlap}")
    if overlap["n_shared"] == 0:
        print("ERROR: no shared genes between Visium and reference.", file=sys.stderr)
        return 1

    print("Running spatial deconvolution …")
    result = tr.deconv_spatial(
        Y, ref, array_row, array_col, lib_sizes, gene_names, spot_ids,
        run_neighbourhood=args.neighbourhood,
    )

    H.ensure_dirs()
    H.write_tsv(result.deconv.proportions,
                H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv",
                comment=["estimate_type: spot_rna_composition (NOT cell counts)",
                         f"lambda_spatial: {result.deconv.lambda_spatial}"])
    if getattr(result, "spot_qc", None) is not None:
        H.write_tsv(result.spot_qc, H.OUT_SPATIAL_DIR / "spatial_qc.tsv")
    if getattr(result, "morans_i", None) is not None:
        H.write_tsv(result.morans_i.to_frame("morans_i"),
                    H.OUT_SPATIAL_DIR / "morans_i.tsv")

    warnings = {
        "estimate_type": result.deconv.ESTIMATE_TYPE,
        "estimate_type_note": "Spot-level RNA-derived composition, not single-cell counts.",
        "converged": bool(result.deconv.converged),
        "lambda_spatial": float(result.deconv.lambda_spatial),
        "gene_overlap": overlap,
    }
    H.write_json(warnings, H.OUT_SPATIAL_DIR / "spatial_warnings.json")
    H.write_json(result.run_metadata, H.OUT_SPATIAL_DIR / "spatial_run_metadata.json")
    print(f"Saved spatial validation outputs -> {H.OUT_SPATIAL_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
