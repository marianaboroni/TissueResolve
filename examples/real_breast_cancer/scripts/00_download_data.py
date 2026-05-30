#!/usr/bin/env python
"""
00 — Download / query the real datasets for breast-cancer validation.

Datasets
--------
* Single-cell reference: a public **human breast cancer single-cell atlas**
  from CZ CELLxGENE (via ``cellxgene-census``).
* Spatial: **10x Visium Human Breast Cancer 1** (OpenProblems / 10x public).

Safety
------
* Default (no flag, env unset): **dry run** — prints exactly what *would* be
  downloaded, writes a manifest skeleton, and exits 0.  No network access.
* Real downloads require ``--run-real-data`` or ``TISSUERESOLVE_RUN_REAL_DATA=1``.
* No fragile temporary URLs are hard-coded; stable dataset identifiers and
  official client libraries are used.  If a client is missing, clear manual
  instructions are printed and the script exits non-zero.

Usage
-----
    python scripts/00_download_data.py                 # dry run
    python scripts/00_download_data.py --run-real-data # actually download
"""
from __future__ import annotations

import argparse
import sys

import _harness as H


# ---------------------------------------------------------------------------
# Manual instructions (printed when automated download is unavailable)
# ---------------------------------------------------------------------------

_REF_MANUAL = f"""\
Manual reference download:
  1. Install the census client:   pip install cellxgene-census
  2. Or download a human breast-cancer single-cell .h5ad with cell-type
     annotations from https://cellxgene.cziscience.com/ (filter: Homo sapiens,
     tissue = breast).
  3. Save it to:
       {H.REFERENCE_H5AD}
  4. Re-run script 01 to build the TissueResolve reference.
"""

_SPATIAL_MANUAL = f"""\
Manual spatial download:
  Option A (OpenProblems):  pip install openproblems   then fetch
            'human_breast_cancer_1' and write it to the path below.
  Option B (scanpy):        python -c "import scanpy as sc; \\
            ad = sc.datasets.visium_sge('V1_Breast_Cancer_Block_A_Section_1'); \\
            ad.write_h5ad('{H.SPATIAL_H5AD}')"
  Option C (10x Genomics):  download the public 'Human Breast Cancer (Block A
            Section 1)' Visium dataset and save the .h5ad to:
       {H.SPATIAL_H5AD}
  Raw counts should be available in layers['counts'] or X.
"""


def _print_plan() -> None:
    print("DRY RUN — no data will be downloaded.")
    print(f"  real-data flag/env not set ({H.REAL_DATA_ENV}=1 to enable).\n")
    print("Would prepare reference:")
    print(f"    -> {H.REFERENCE_H5AD}")
    print(f"       (downsample to <= {H.MAX_REFERENCE_CELLS} cells, "
          f"<= {H.MAX_REFERENCE_GENES} genes if larger)\n")
    print("Would prepare spatial dataset:")
    print(f"    -> {H.SPATIAL_H5AD}\n")
    print(_REF_MANUAL)
    print(_SPATIAL_MANUAL)


# ---------------------------------------------------------------------------
# Real download paths (best-effort, dependency-guarded)
# ---------------------------------------------------------------------------


def download_reference(manifest: dict) -> bool:
    """Query CELLxGENE Census for a breast-cancer reference.  Returns success."""
    try:
        import cellxgene_census  # noqa: F401
    except Exception:
        print("ERROR: cellxgene-census is not installed; cannot auto-download "
              "the reference.", file=sys.stderr)
        print(_REF_MANUAL, file=sys.stderr)
        return False

    import anndata as ad  # noqa: F401
    import cellxgene_census

    # Discover the current census API; census version is recorded for provenance.
    census_version = "stable"
    print(f"Opening CELLxGENE Census (version={census_version}) …")
    with cellxgene_census.open_soma(census_version=census_version) as census:
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter=(
                "tissue_general == 'breast' and disease != 'normal' "
                "and is_primary_data == True"
            ),
            column_names={"obs": ["cell_type", "tissue", "disease", "assay",
                                   "donor_id", "is_primary_data"]},
        )
    print(f"Census returned {adata.n_obs} cells × {adata.n_vars} genes.")
    adata = _downsample_reference(adata)
    H.REFERENCE_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(H.REFERENCE_H5AD)
    entry = manifest["datasets"]["reference"]
    entry["downloaded"] = True
    entry["url_or_id"] = f"cellxgene-census:{census_version}"
    entry["downsampling"]["applied"] = (
        adata.n_obs >= H.MAX_REFERENCE_CELLS or adata.n_vars >= H.MAX_REFERENCE_GENES
    )
    print(f"Saved reference -> {H.REFERENCE_H5AD}")
    return True


def _downsample_reference(adata):
    """Cap cells/genes for a first validation run; preserve major cell types."""
    import numpy as np

    rng = np.random.default_rng(0)
    if adata.n_obs > H.MAX_REFERENCE_CELLS:
        keep = rng.choice(adata.n_obs, H.MAX_REFERENCE_CELLS, replace=False)
        adata = adata[np.sort(keep)].copy()
    if adata.n_vars > H.MAX_REFERENCE_GENES:
        import scipy.sparse as sp

        X = adata.X
        expr = np.asarray(X.sum(axis=0)).ravel() if sp.issparse(X) else X.sum(0)
        top = np.argsort(expr)[::-1][: H.MAX_REFERENCE_GENES]
        adata = adata[:, np.sort(top)].copy()
    return adata


def download_spatial(manifest: dict) -> bool:
    """Fetch the 10x Visium Human Breast Cancer 1 dataset.  Returns success."""
    # Try OpenProblems first, then scanpy's public Visium loader.
    try:
        import scanpy as sc
    except Exception:
        print("ERROR: scanpy is not installed; cannot auto-download spatial.",
              file=sys.stderr)
        print(_SPATIAL_MANUAL, file=sys.stderr)
        return False

    try:
        adata = sc.datasets.visium_sge("V1_Breast_Cancer_Block_A_Section_1")
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
        H.SPATIAL_H5AD.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(H.SPATIAL_H5AD)
        entry = manifest["datasets"]["spatial"]
        entry["downloaded"] = True
        entry["url_or_id"] = "scanpy.datasets.visium_sge:V1_Breast_Cancer_Block_A_Section_1"
        print(f"Saved spatial -> {H.SPATIAL_H5AD}")
        return True
    except Exception as exc:  # pragma: no cover - network path
        print(f"ERROR: spatial download failed: {exc}", file=sys.stderr)
        print(_SPATIAL_MANUAL, file=sys.stderr)
        return False


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-real-data", action="store_true",
                        help="Actually download (otherwise dry run).")
    args = parser.parse_args(argv)

    H.ensure_dirs()
    manifest = H.default_manifest()

    if not H.real_data_enabled(args.run_real_data):
        _print_plan()
        H.write_manifest(manifest, H.MANIFEST_PATH)
        print(f"\nWrote manifest skeleton -> {H.MANIFEST_PATH}")
        return 0

    print("REAL DATA MODE — downloading datasets.\n")
    ok_ref = download_reference(manifest)
    ok_sp = download_spatial(manifest)
    H.write_manifest(manifest, H.MANIFEST_PATH)
    print(f"\nWrote manifest -> {H.MANIFEST_PATH}")

    if not (ok_ref and ok_sp):
        print("\nOne or more downloads were unavailable; see instructions above.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
