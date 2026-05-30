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

CELLxGENE Census compatibility
------------------------------
``cellxgene-census`` pins a ``tiledbsoma`` range.  Newer Census *stable*
releases are encoded with a SOMA object-encoding version that an older
``tiledbsoma`` cannot read (``Unsupported SOMA object encoding version``).
This script therefore does **not** hard-code ``"stable"``: it defaults to a
pinned Census release known to be readable by the pinned client, and falls
back across older pinned releases.  If none work, it stops with an actionable
message (installed versions + recommended command + manual fallback).

Usage
-----
    python scripts/00_download_data.py                       # dry run
    python scripts/00_download_data.py --run-real-data       # download
    python scripts/00_download_data.py --run-real-data \\
        --census-version 2023-07-25                          # pick a release
"""
from __future__ import annotations

import argparse
import os
import platform
import sys
import time
from datetime import datetime, timezone

import _download_utils as du
import _harness as H

# A pinned LTS Census release readable by older tiledbsoma builds, then
# progressively older fallbacks.  These are stable dated identifiers, not
# fragile URLs.  "stable" is tried last because newer stable releases may use
# a SOMA encoding the installed tiledbsoma cannot read.
DEFAULT_CENSUS_VERSION = "2024-07-01"
CENSUS_FALLBACK_VERSIONS = ["2023-12-15", "2023-07-25", "2023-05-15", "stable"]


# ---------------------------------------------------------------------------
# Manual instructions
# ---------------------------------------------------------------------------

_REF_MANUAL = f"""\
Manual reference download:
  1. Upgrade the client so it can read the current Census:
       python -m pip install -U cellxgene-census tiledbsoma
     or pick an older Census release with --census-version (e.g. 2023-07-25).
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


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _print_plan(census_version: str) -> None:
    print("DRY RUN — no data will be downloaded.")
    print(f"  real-data flag/env not set ({H.REAL_DATA_ENV}=1 to enable).\n")
    print("Would prepare reference (CELLxGENE Census):")
    print(f"    census_version : {census_version} "
          f"(fallbacks: {', '.join(CENSUS_FALLBACK_VERSIONS)})")
    print(f"    cellxgene-census: {du.get_package_version('cellxgene-census')}")
    print(f"    tiledbsoma      : {du.get_package_version('tiledbsoma')}")
    print(f"    -> {H.REFERENCE_H5AD}")
    print(f"       (downsample to <= {H.MAX_REFERENCE_CELLS} cells, "
          f"<= {H.MAX_REFERENCE_GENES} genes if larger)\n")
    print("Would prepare spatial dataset:")
    print(f"    -> {H.SPATIAL_H5AD}\n")
    print(_REF_MANUAL)
    print(_SPATIAL_MANUAL)


# ---------------------------------------------------------------------------
# Reference download (CELLxGENE Census, with version fallback + progress)
# ---------------------------------------------------------------------------

_CENSUS_OBS_FILTER = (
    "tissue_general == 'breast' and disease != 'normal' "
    "and is_primary_data == True"
)
_CENSUS_OBS_COLS = ["cell_type", "tissue", "disease", "assay",
                    "donor_id", "is_primary_data"]


def _soma_help(entry: dict, requested: str, exc: BaseException | None) -> str:
    return (
        "ERROR: could not read the CELLxGENE Census with the installed client.\n"
        f"  a. cellxgene-census installed : {entry['cellxgene_census_version']}\n"
        f"  b. tiledbsoma installed       : {entry['tiledbsoma_version']}\n"
        f"  c. census_version requested   : {requested} "
        f"(attempted: {entry.get('fallback_versions_attempted')})\n"
        f"  d. recommended fix:\n"
        "       python -m pip install -U cellxgene-census tiledbsoma\n"
        "     or choose an older release:  --census-version 2023-07-25\n"
        f"     last error: {exc}\n\n" + _REF_MANUAL
    )


def download_reference(
    manifest: dict,
    *,
    census_version: str = DEFAULT_CENSUS_VERSION,
    force: bool = False,
    clock=time.monotonic,
    now=_now_iso,
) -> bool:
    """Query CELLxGENE Census for a breast-cancer reference.  Returns success."""
    entry = manifest["datasets"]["reference"]
    entry["method"] = "cellxgene-census"
    entry["source"] = "CZ CELLxGENE Census (cellxgene-census)"
    entry["cellxgene_census_version"] = du.get_package_version("cellxgene-census")
    entry["tiledbsoma_version"] = du.get_package_version("tiledbsoma")
    entry["python_version"] = platform.python_version()
    entry["start_time"] = now()
    entry["fallback_versions_attempted"] = []

    # Reuse an already-downloaded reference unless --force.
    if H.REFERENCE_H5AD.exists() and not force:
        size = os.path.getsize(H.REFERENCE_H5AD)
        entry["downloaded"] = True
        entry["status"] = "already_exists"
        entry["method"] = "reuse_existing"
        entry["file_size_bytes"] = int(size)
        entry["end_time"] = now()
        print(f"Reference already exists ({du.format_bytes(size)}); skipping "
              f"Census query (use --force to re-download).\n  {H.REFERENCE_H5AD}")
        return True

    try:
        import cellxgene_census
    except Exception:
        print("ERROR: cellxgene-census is not installed; cannot auto-download "
              "the reference.  Install with: pip install -e \".[realdata]\"",
              file=sys.stderr)
        print(_REF_MANUAL, file=sys.stderr)
        entry["end_time"] = now()
        return False

    log = du.ProgressLogger("reference", clock=clock)
    log.stage(5, f"checking package versions "
                 f"(census={entry['cellxgene_census_version']}, "
                 f"tiledbsoma={entry['tiledbsoma_version']})")

    # Build the ordered candidate list (requested first, then fallbacks).
    candidates: list[str] = []
    for v in [census_version, *CENSUS_FALLBACK_VERSIONS]:
        if v not in candidates:
            candidates.append(v)

    adata = None
    used_version: str | None = None
    last_exc: BaseException | None = None
    for v in candidates:
        entry["fallback_versions_attempted"].append(v)
        try:
            log.stage(10, f"opening census (version={v})")
            with cellxgene_census.open_soma(census_version=v) as census:
                log.stage(20, "querying datasets table")
                log.stage(35, "selecting breast-cancer cells")
                log.stage(50, "retrieving AnnData subset")
                adata = cellxgene_census.get_anndata(
                    census,
                    organism="Homo sapiens",
                    obs_value_filter=_CENSUS_OBS_FILTER,
                    column_names={"obs": _CENSUS_OBS_COLS},
                )
            used_version = v
            break
        except Exception as exc:  # noqa: BLE001 - need to classify
            last_exc = exc
            if du.is_soma_encoding_error(exc):
                log.info(f"census version {v!r} incompatible (SOMA encoding); "
                         "trying next release…")
                continue
            # An unexpected error — stop with a clear message.
            print(_soma_help(entry, census_version, exc), file=sys.stderr)
            entry["end_time"] = now()
            return False

    if adata is None:
        print(_soma_help(entry, census_version, last_exc), file=sys.stderr)
        entry["end_time"] = now()
        return False

    log.info(f"loaded {adata.n_obs} cells × {adata.n_vars} genes "
             f"(census_version={used_version})")
    log.stage(70, "downsampling cells/genes if needed")
    n_before = (adata.n_obs, adata.n_vars)
    adata = _downsample_reference(adata)
    downsampled = (adata.n_obs, adata.n_vars) != n_before
    log.info(f"after downsampling: {adata.n_obs} cells × {adata.n_vars} genes")

    log.stage(85, "writing h5ad")
    H.REFERENCE_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(H.REFERENCE_H5AD)
    size = os.path.getsize(H.REFERENCE_H5AD)
    log.done(f"saved reference h5ad ({du.format_bytes(size)})")

    entry["downloaded"] = True
    entry["status"] = "downloaded"
    entry["census_version"] = used_version
    entry["url_or_id"] = f"cellxgene-census:{used_version}"
    entry["file_size_bytes"] = int(size)
    entry["content_length_available"] = False  # census API has no byte total
    entry["downsampling"]["applied"] = bool(downsampled)
    entry["end_time"] = now()
    entry["elapsed_seconds"] = round(log.elapsed(), 3)
    print(f"Saved reference -> {H.REFERENCE_H5AD} ({du.format_bytes(size)})")
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


# ---------------------------------------------------------------------------
# Spatial download (scanpy public Visium, Python-3.9-compatible, with progress)
# ---------------------------------------------------------------------------

_VISIUM_SAMPLE = "V1_Breast_Cancer_Block_A_Section_1"


def ensure_tarfile_data_filter(mod=None) -> bool:
    """Make ``tarfile.data_filter`` available on Python < 3.12 (PEP 706).

    scanpy's Visium loader extracts a tar archive and references
    ``tarfile.data_filter``, which only exists on Python 3.12+.  On older
    Pythons this raises ``module 'tarfile' has no attribute 'data_filter'``.
    We install a permissive shim (the 10x public dataset is trusted) so the
    extraction succeeds.

    Parameters
    ----------
    mod:
        Module to patch (defaults to the real ``tarfile``); injectable for tests.

    Returns
    -------
    bool
        Whether ``data_filter`` was **already** available before patching.
    """
    import tarfile

    mod = mod if mod is not None else tarfile
    available = hasattr(mod, "data_filter")
    if not available:
        def _permissive(member, dest_path):  # PEP 706 filter signature
            return member

        mod.data_filter = _permissive
        if not hasattr(mod, "fully_trusted_filter"):
            mod.fully_trusted_filter = _permissive
    return available


def _is_count_like(X) -> bool:
    """Heuristic: non-negative, integer-valued matrix (a sample of it)."""
    import numpy as np
    import scipy.sparse as sp

    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return False
    sample = np.asarray(data[:10000], dtype=float)
    return bool(np.all(sample >= 0) and np.allclose(sample, np.round(sample)))


def _validate_and_prepare_spatial(adata) -> dict:
    """Validate the Visium AnnData and ensure a counts layer when safe."""
    has_counts = "counts" in adata.layers
    has_spatial = "spatial" in adata.obsm
    copied = False
    if not has_counts and _is_count_like(adata.X):
        adata.layers["counts"] = adata.X.copy()
        has_counts = True
        copied = True
    return {
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "has_counts": has_counts,
        "has_spatial": has_spatial,
        "copied_counts_from_X": copied,
    }


def _default_visium_loader(sample: str):
    """Real loader: patch tarfile for 3.9, then call scanpy's Visium loader."""
    import scanpy as sc

    ensure_tarfile_data_filter()
    return sc.datasets.visium_sge(sample)


def download_spatial(
    manifest: dict,
    *,
    force: bool = False,
    loader=None,
    clock=time.monotonic,
    now=_now_iso,
) -> bool:
    """Fetch the 10x Visium Human Breast Cancer 1 dataset.  Returns success.

    *loader* is an injectable ``sample -> AnnData`` callable (defaults to the
    scanpy path) so the download can be tested fully offline.
    """
    entry = manifest["datasets"]["spatial"]
    entry["python_version"] = platform.python_version()
    entry["scanpy_version"] = du.get_package_version("scanpy")
    entry["tarfile_data_filter_available"] = ensure_tarfile_data_filter()
    entry["start_time"] = now()

    log = du.ProgressLogger("spatial", clock=clock)
    log.stage(10, "checking existing file")

    if H.SPATIAL_H5AD.exists() and not force:
        size = os.path.getsize(H.SPATIAL_H5AD)
        info = _read_existing_spatial_flags(H.SPATIAL_H5AD)
        entry.update(
            downloaded=True,
            spatial_status="already_exists",
            spatial_download_method="reuse_existing",
            spatial_file_size_bytes=int(size),
            spatial_has_counts_layer=info.get("has_counts"),
            spatial_has_spatial_obsm=info.get("has_spatial"),
            end_time=now(),
        )
        print(f"Spatial dataset already exists ({du.format_bytes(size)}); "
              f"skipping (use --force to re-download).\n  {H.SPATIAL_H5AD}")
        return True

    log.stage(20, "resolving download method")
    load = loader if loader is not None else _default_visium_loader

    try:
        log.stage(40, f"fetching Visium sample {_VISIUM_SAMPLE!r}")
        adata = load(_VISIUM_SAMPLE)
    except Exception as exc:
        print(f"ERROR: spatial download failed: {exc}", file=sys.stderr)
        print(_SPATIAL_MANUAL, file=sys.stderr)
        entry["spatial_status"] = "failed"
        entry["spatial_download_method"] = "scanpy.datasets.visium_sge"
        entry["end_time"] = now()
        return False

    log.stage(70, "validating AnnData")
    info = _validate_and_prepare_spatial(adata)
    log.info(f"{info['n_spots']} spots × {info['n_genes']} genes; "
             f"counts_layer={info['has_counts']} "
             f"(copied_from_X={info['copied_counts_from_X']}), "
             f"spatial_obsm={info['has_spatial']}")

    log.stage(85, "writing h5ad")
    H.SPATIAL_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(H.SPATIAL_H5AD)
    size = os.path.getsize(H.SPATIAL_H5AD)
    log.done(f"saved spatial h5ad ({du.format_bytes(size)})")

    entry.update(
        downloaded=True,
        spatial_status="downloaded",
        spatial_download_method="scanpy.datasets.visium_sge",
        url_or_id=f"scanpy.datasets.visium_sge:{_VISIUM_SAMPLE}",
        spatial_file_size_bytes=int(size),
        spatial_has_counts_layer=info["has_counts"],
        spatial_has_spatial_obsm=info["has_spatial"],
        content_length_available=False,
        end_time=now(),
        elapsed_seconds=round(log.elapsed(), 3),
    )
    print(f"Saved spatial -> {H.SPATIAL_H5AD} ({du.format_bytes(size)})")
    return True


def _read_existing_spatial_flags(path) -> dict:
    """Cheaply read an existing Visium h5ad to record counts/obsm flags."""
    try:
        import anndata as ad

        a = ad.read_h5ad(path)
        return {"has_counts": "counts" in a.layers, "has_spatial": "spatial" in a.obsm}
    except Exception:
        return {"has_counts": None, "has_spatial": None}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-real-data", action="store_true",
                        help="Actually download (otherwise dry run).")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if the target file already exists.")
    parser.add_argument("--census-version", default=DEFAULT_CENSUS_VERSION,
                        help="CELLxGENE Census release to try first "
                             f"(default {DEFAULT_CENSUS_VERSION}; falls back "
                             f"to {', '.join(CENSUS_FALLBACK_VERSIONS)}).")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    H.ensure_dirs()
    manifest = H.default_manifest()
    manifest["run_started"] = _now_iso()

    if not H.real_data_enabled(args.run_real_data):
        _print_plan(args.census_version)
        H.write_manifest(manifest, H.MANIFEST_PATH)
        print(f"\nWrote manifest skeleton -> {H.MANIFEST_PATH}")
        return 0

    print("REAL DATA MODE — downloading datasets.\n")
    t0 = time.monotonic()
    ok_ref = download_reference(
        manifest, census_version=args.census_version, force=args.force)
    ok_sp = download_spatial(manifest, force=args.force)
    manifest["run_finished"] = _now_iso()
    manifest["run_elapsed_seconds"] = round(time.monotonic() - t0, 3)
    H.write_manifest(manifest, H.MANIFEST_PATH)
    print(f"\nWrote manifest -> {H.MANIFEST_PATH}")

    if not (ok_ref and ok_sp):
        print("\nOne or more downloads were unavailable; see instructions above.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
