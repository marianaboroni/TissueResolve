# External-Tool Benchmark Plan

Plan for a fair benchmark of TissueResolve against established deconvolution tools,
**after** the full-workflow validation (see `FULL_WORKFLOW_REPORT_VALIDATION.md`).
Nothing is committed; large inputs/outputs are git-ignored.

## 1. Datasets
* **Breast** single-cell reference (CELLxGENE), 30k cells, 41 donors, 41 fine → 8 families.
* **HLCA/lung** subset, 11.4k cells, 40 donors, 59 fine → 10 families.
* (Spatial real, no truth) breast Visium `human_breast_cancer_1.h5ad` if used.

## 2. Ground-truth availability
* **Bulk** (Category A): donor-held-out pseudobulk simulated from labelled cells → truth known.
* **Spatial synthetic** (Category B): pseudo-spots from labelled cells → truth known.
* **Real Visium** (Category C): **no ground truth** → only marker recovery / concordance /
  structure / stability / runtime. Accuracy is NOT reported (rule 7).

## 3. Reference construction
One shared reference per tissue, built from **reference donors only** (50% donor split,
seed-disjoint from query). Same reference, fine labels, and broad→fine hierarchy fed to
every tool where the tool's API allows it (rule 8). Min cells/type = 20.

## 4. Train / calibration / test splits
Donor- AND seed-disjoint: reference donors (build), query/test donors (pseudobulk).
Gene selection / weighting use reference donors only; the soft-gating confidence model is
fit on calibration seeds only; metrics come from held-out test seeds (rules 5, 6, 10).
Test truth is never shown to any tool.

## 5. Tool list (see `tool_registry.tsv` for the real installed state)
* Bulk: NNLS/WNNLS baselines, TissueResolve flat, TissueResolve hierarchical+soft,
  MuSiC, BisqueRNA, BayesPrism, (SCDC/DWLS if installed), CIBERSORTx (web/licensed).
* Spatial: NNLS-per-spot, TissueResolve spatial, RCTD/spacexr, cell2location, CARD,
  (SPOTlight/Tangram if installed).

## 6. Installation strategy
Use official packages / documented routes (CRAN/Bioconductor/GitHub for R via
`envs/Rlib`; pip/conda for cell2location via `envs/c2l_py39`). Record exact versions;
never silently skip — every unavailable/failed tool is logged with a reason
(`EXTERNAL_TOOL_INSTALLATION_NOTES.md`, `run_manifest.tsv`). See rule 9.

## 7. Harmonized inputs
Consistent gene identifiers (gene symbols; intersection reported in
`gene_overlap_summary.tsv`), harmonized cell-type labels and broad/fine levels
(`label_mapping.tsv`), the same donor-held-out test mixtures and query matrices,
normalization only as each tool requires (`normalization_summary.tsv`). Differences
documented (`input_manifest.tsv`). Large matrices not committed.

## 8. Harmonized outputs
All tools emit per-sample/per-spot proportions over the same cell-type set, reindexed and
0-filled to the shared label set before scoring. TissueResolve additionally emits
resolution-aware columns (unresolved mass, trusted resolution).

## 9. Metrics
* Broad: Pearson, Spearman, RMSE, MAE, JSD, Aitchison, dominant-family accuracy.
* Fine: abs RMSE, fine Pearson/Spearman, conditional within-family RMSE/Pearson,
  macro subtype error, rare sensitivity, rare precision, false-positive detection,
  pairwise spillover.
* Resolution-aware (TissueResolve only): unresolved precision/recall, false-resolution,
  false-abstention, correct-resolution-level rate, effective-N, entropy, richness,
  dominant fraction.

## 10. Runtime / memory capture
Wall-clock per tool/run + peak RSS via subprocess accounting, recorded in
`*_runtime.tsv` and `run_manifest.tsv`.

## 11. Failure handling
Each run is `executed | imported | exported_only | skipped | failed` with a reason
(unavailable / install failure / license / runtime failure / incompatible input /
excessive memory-runtime). Heavy tools (BayesPrism ~2180s prior; cell2location MPS/CPU
~160s+; RCTD/CARD) are bounded by a per-tool timeout; a timeout is recorded as
`failed:timeout`, never hidden.

## 12. Reporting strategy
Fresh internal bulk comparison (TR variants + NNLS/WNNLS baselines) on both tissues with
paired stats; external R/Python tools attempted with bounded timeouts and otherwise
documented from their real installed state + prior-run evidence. Consolidated in
`EXTERNAL_TOOL_BENCHMARK_REPORT.md`. Claims use careful language (outperformed on X /
matched / faster than / competitive but not superior); no "best overall" claim unless
supported across metrics and tissues (rule 10).

## Categories
* **A — Bulk with truth**: donor-held-out pseudobulk (this plan's primary, fully run).
* **B — Spatial synthetic with truth**: pseudo-spots; reuses validated spatial scripts.
* **C — Real Visium without truth**: marker recovery / concordance / structure only.
