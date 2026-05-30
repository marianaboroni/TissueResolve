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
# For real downloads of the reference:
python -m pip install cellxgene-census
```

## Run order

```bash
cd examples/real_breast_cancer

# 0) Dry run — prints the plan, writes a manifest skeleton, no network:
python scripts/00_download_data.py

# 0') Real download (opt-in):
python scripts/00_download_data.py --run-real-data
#   or:  TISSUERESOLVE_RUN_REAL_DATA=1 python scripts/00_download_data.py

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
