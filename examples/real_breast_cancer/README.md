# Real-data validation: breast cancer

Validation tooling for TissueResolve using **one** real human breast-cancer
single-cell reference for both:

1. **Bulk** validation via pseudobulk mixtures with known ground-truth
   proportions, and
2. **Spatial** validation on a real 10x Visium breast-cancer section.

This is **validation tooling only**. It does not modify or weaken the core
algorithms. The default TissueResolve test suite stays offline; nothing here
downloads data unless you explicitly opt in.

## Safety model

- **Dry run is the default.** Scripts never touch the network unless
  `--run-real-data` is passed or `TISSUERESOLVE_RUN_REAL_DATA=1` is set.
- All resolved dataset IDs, versions, filters, and downsampling are written to
  `data/download_manifest.json`.
- No fragile temporary URLs are hard-coded; official client libraries
  (`cellxgene-census`, `scanpy`) and documented dataset IDs are used. If a
  client is missing, the script prints exact manual instructions and exits
  non-zero.
- Real data files (`*.h5ad`, derived TSVs) are git-ignored under `data/`.

## Datasets

| Role | Dataset | Source | Saved to |
|---|---|---|---|
| Reference | Human breast cancer single-cell atlas | CZ CELLxGENE Census | `data/reference/breast_cancer_sc_reference.h5ad` |
| Spatial | 10x Visium Human Breast Cancer 1 (~4263 spots × 14906 genes) | OpenProblems / 10x public | `data/spatial/human_breast_cancer_1.h5ad` |

If the reference is large it is downsampled to ≤ 30,000 cells and ≤ 5,000
highly-expressed genes, preserving major cell types.

## Install extras

```bash
source ../../.venv/bin/activate
python -m pip install -e "../..[all]"
```

Real-data validation additionally requires `cellxgene-census` (a heavy,
opt-in dependency kept out of `[all]`). Install it via the `realdata` extra:

```bash
python -m pip install -e ".[realdata]"
# or, equivalently:
python -m pip install -U cellxgene-census
```

## Run order

```bash
cd examples/real_breast_cancer

# 0) Dry run — prints the plan, writes a manifest skeleton, no network:
python scripts/00_download_data.py

# 0') Real download (opt-in):
python scripts/00_download_data.py --run-real-data
#   or:  TISSUERESOLVE_RUN_REAL_DATA=1 python scripts/00_download_data.py
#   pick a specific Census release:
python scripts/00_download_data.py --run-real-data --census-version 2023-07-25
#   re-download even if the target files already exist:
python scripts/00_download_data.py --run-real-data --force

# 1) Build the TissueResolve reference (NB overdispersion for spatial too):
python scripts/01_prepare_reference.py            # --cell-type-col / --min-cells

# 2) Make pseudobulk mixtures (easy/medium/hard, >=12 samples):
python scripts/02_make_pseudobulk.py

# 3) Bulk validation (estimated vs true mRNA proportions):
python scripts/03_run_bulk_validation.py

# 4) Spatial validation (Moran's I, QC, warnings):
python scripts/04_run_spatial_validation.py

# 5) Summary report (works after a partial run):
python scripts/05_summarize_results.py
```

## Progress reporting

`00_download_data.py` prints progress as it works:

- **Reference (Census)** — stage-level progress, since the Census API exposes
  no byte total:
  ```
  [  5%] reference: checking package versions (census=1.15.0, tiledbsoma=1.11.4) (elapsed 00:00:00)
  [ 10%] reference: opening census (version=2024-07-01) (elapsed 00:00:01)
  [ 35%] reference: selecting breast-cancer cells (elapsed 00:00:09)
  [ 50%] reference: retrieving AnnData subset (elapsed 00:00:20)
  [ 70%] reference: downsampling cells/genes if needed (elapsed 00:01:05)
  [ 85%] reference: writing h5ad (elapsed 00:01:10)
  [100%] reference: saved reference h5ad (84.2 MB) (elapsed 00:01:12)
  ```
  Cell/gene counts are printed after loading and after downsampling, and the
  output file size after saving.
- **Direct HTTP downloads** — byte-level progress when `Content-Length` is
  known, else downloaded MB + elapsed:
  ```
  [42.3%] 317.5 MB / 750.0 MB elapsed 00:03:21
  [ --- ] 317.5 MB downloaded elapsed 00:03:21      # no Content-Length
  ```

## Re-running / skipping downloads

By default the downloader **reuses any file that already exists** and records
`status: already_exists` (reference) / `spatial_status: already_exists`
(spatial) in the manifest — so a partial run is resumable and the reference is
never re-fetched needlessly. Pass `--force` to re-download.

## Python 3.9 and the spatial `tarfile.data_filter` error

scanpy's Visium loader extracts a tar archive and references
`tarfile.data_filter`, which only exists on Python 3.12+ (PEP 706). On Python
3.9–3.11 this surfaces as:

```
module 'tarfile' has no attribute 'data_filter'
```

The downloader installs a permissive compatibility shim before calling the
loader (the 10x public dataset is trusted), so the spatial download works on
Python 3.9. Whether the attribute was natively available is recorded as
`tarfile_data_filter_available` in the manifest. The downloaded AnnData is
validated (spots, genes, `layers['counts']`, `obsm['spatial']`); if no counts
layer exists and `X` is count-like, `X` is copied into `layers['counts']` and
the decision is recorded.

## CELLxGENE Census compatibility

If you see:

```
ValueError: Unsupported SOMA object encoding version 1.1.0
```

it means the installed `tiledbsoma` (pinned by your `cellxgene-census`) cannot
read the **current `stable`** Census release. The downloader handles this:

- It does **not** default to `stable`. It defaults to a pinned LTS release
  (`2024-07-01`) and **falls back** across older releases
  (`2023-12-15`, `2023-07-25`, `2023-05-15`, then `stable`), skipping any that
  raise the SOMA-encoding error.
- Pick a release explicitly with `--census-version YYYY-MM-DD`.
- If **all** releases fail, the script stops with an actionable message listing
  the installed `cellxgene-census` / `tiledbsoma` versions, the requested
  version, and the recommended fix:
  ```
  python -m pip install -U cellxgene-census tiledbsoma
  ```
- Manual fallback: download a human breast-cancer single-cell `.h5ad` (with
  cell-type annotations) from https://cellxgene.cziscience.com/ and save it to
  `data/reference/breast_cancer_sc_reference.h5ad`, then run script 01.

The resolved `census_version`, both package versions, timings, file sizes, and
any fallback versions attempted are recorded in `data/download_manifest.json`.

## Expected outputs

```
outputs/reference/       reference_summary.tsv, cell_type_counts.tsv,
                         selected_annotation_column.txt, breast_cancer_reference/
data/derived/            pseudobulk_counts.tsv, pseudobulk_true_proportions.tsv,
                         pseudobulk_metadata.tsv
outputs/bulk/            bulk_estimated_proportions.tsv, bulk_validation_metrics.tsv,
                         bulk_per_celltype_metrics.tsv, bulk_qc.tsv, bulk_warnings.json
outputs/spatial/         spatial_spot_proportions.tsv, spatial_qc.tsv, morans_i.tsv,
                         spatial_warnings.json, spatial_run_metadata.json
outputs/validation_summary/  validation_summary.md, validation_summary.json
```

## What the numbers mean

- **Bulk ground truth is mRNA proportions** (count-fraction per cell type), not
  cell fractions — because the deconvolver estimates mRNA proportions. See
  `../../docs/output_interpretation.md`.
- **Spatial** has no per-spot ground truth; validation is qualitative
  (spatial coherence via Moran's I, QC, surfaced warnings).

## Tests

Offline-only tests live in `tests/` and use tiny synthetic AnnData; they never
download data. Run from the repo root:

```bash
python -m pytest examples/real_breast_cancer/tests -v
```
