# Real-data TNBC benchmark (bulk TCGA + spatial Visium + single-cell atlas)

This documents a real-data run of TissueResolve on triple-negative breast
cancer (TNBC), using **one shared single-cell reference** for both a **bulk**
and a **spatial** workflow. All data is downloaded by the harness scripts and is
git-ignored; nothing here is committed.

> **Scientific framing — read first.** Real bulk tumours and real Visium have
> **no ground-truth cell proportions**, so this is a **concordance / plausibility
> / robustness / runtime** benchmark, **not** an accuracy benchmark. Accuracy
> (Pearson/RMSE vs known truth) is only meaningful on the **pseudobulk** harness
> (scripts 02–03), which has ground truth by construction.

## Datasets

| Role | Dataset | Source | Size used |
|---|---|---|---|
| Single-cell **reference** | Human breast cancer cell atlas (32 cell types, 5000 genes) | CZ CELLxGENE Census | 30k cells (capped) |
| **Bulk** TNBC | TCGA-BRCA, triple-negative (ER−/PR−/HER2− by IHC), STAR gene counts, primary tumour | NCI GDC API | 40 of 115 samples (this run) |
| **Spatial** | 10x Visium human breast cancer (counts in `layers['counts']`) | 10x / OpenProblems | 3798 spots × 36601 genes |

TNBC selection is reproducible from IHC receptor status in the TCGA-BRCA clinical
Biotab: **116 TNBC patients** of 1097; **115** have a primary-tumour STAR-Counts
RNA-seq file (~488 MB for the full cohort).

## Commands

```bash
export TISSUERESOLVE_RUN_REAL_DATA=1
cd examples/real_breast_cancer

# reference + spatial (already downloaded by 00_download_data.py)
python scripts/00_download_data.py --run-real-data
python scripts/01_prepare_reference.py            # builds outputs/reference/breast_cancer_reference

# NEW: real bulk TCGA TNBC (GDC). Omit --max-samples for all 115 samples (~488 MB)
python scripts/00b_download_tcga_tnbc_bulk.py --dry-run
python scripts/00b_download_tcga_tnbc_bulk.py --run-real-data --max-samples 40
python scripts/03b_run_bulk_tcga.py --run-real-data

# spatial deconvolution on the real Visium
python scripts/04_run_spatial_validation.py
```

## Results (this run — 40 TNBC bulk samples, 3798 Visium spots)

### Bulk TCGA TNBC — no-ground-truth method benchmark
`outputs/bulk_tcga/tcga_bulk_method_benchmark.tsv`

| Method | Runtime (s) | Mean reconstruction Pearson* | Concordance vs `auto` (Pearson / Spearman) |
|---|---|---|---|
| TissueResolve `auto` | 131 | 0.39 | 1.00 / 1.00 |
| TissueResolve `nnls` | 30 | 0.46 | 0.44 / 0.33 |
| TissueResolve `weighted_nnls` | 31 | 0.46 | 0.44 / 0.27 |

\* observed vs reference-reconstructed CPM (log1p), per sample, averaged — a
self-consistency check, **not** accuracy.

**Honest observations:**
- **Methods disagree substantially** on real bulk (cross-method Pearson ~0.27–0.44).
  Without ground truth there is nothing to adjudicate the disagreement — solver
  choice materially changes the answer on real data. This is the central caveat.
- Reconstruction concordance is **modest** (~0.39–0.46): a breast scRNA atlas
  only partially explains TCGA bulk expression (platform/protocol differences,
  tumour purity, batch). Plain/weighted NNLS reconstruct the observed profile
  slightly better than `auto` here, while `auto` (gene-masking CV) is ~4× slower.
- Mean predicted composition is **biologically plausible** for TNBC tissue:
  epithelial (mammary/luminal/basal) + endothelial + immune (NK/T, monocytes).

### Spatial Visium — no-ground-truth structure benchmark
`outputs/spatial/` (`spatial_spot_proportions.tsv`, `morans_i.tsv`, `report.html`)

- 3798 spots × 32 cell types; gene overlap 4934/5000 reference genes.
- Mean composition epithelial-dominant (mammary/luminal epithelial) with
  endothelial + immune — plausible breast-tumour architecture.
- **Moran's I ≈ 0.96–0.97** for epithelial populations → strong spatial
  structure (tumour epithelium forms coherent regions), exactly what a spatial
  model should recover. **No accuracy is claimed** (no spot-level ground truth).

## Reproducibility / provenance

- `data/download_manifest.json` carries a `bulk_tcga_tnbc` entry (selection
  criteria, file count, package versions, access date).
- `data/bulk_tcga_tnbc/tcga_tnbc_download.json` lists every GDC `file_id`.
- All inputs/outputs under `data/` and `outputs/` are git-ignored.

## Scale-up

For the full TNBC cohort, drop `--max-samples` (downloads all 115, ~488 MB).
The bulk benchmark scales linearly in samples; `auto` is the slow path
(gene-masking CV). For a ground-truth accuracy benchmark, use the pseudobulk
harness (`02_make_pseudobulk.py` → `03_run_bulk_validation.py`).
