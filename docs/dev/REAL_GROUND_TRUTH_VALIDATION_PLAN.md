# Real ground-truth validation plan (bulk↔flow, spatial↔Xenium/CosMx)

The single most important gate to move TissueResolve from "alpha with promise" to
"scientifically validated": measure accuracy on **real data with orthogonal ground
truth**, not pseudobulk/synthetic. This document defines candidate datasets, the local
data contract, metrics, and gates. A harness
(`benchmarks/dev/real_ground_truth_validation.py`) consumes local files and refuses to
run without them — **no data is downloaded automatically**; nothing here is fabricated.

> Author note: dataset identifiers below are from prior literature knowledge and
> **must be verified** at acquisition time (accessions/URLs change). The harness does
> not depend on any specific accession — it reads documented local files.

## What ground truth actually validates (honest scope)

- **Flow/CyTOF** on bulk gives proportions at the **broad immune-lineage** level
  (T, CD4/CD8, B, NK, monocyte, …), *not* fine collinear subtypes. So this validates
  **broad/family accuracy** (exactly the level our synthetic work found reliable) — and
  will *not* rescue collinear fine states (consistent with all prior findings).
- **Xenium/CosMx** segmentation gives per-cell type calls → per-spot proportions after
  co-registration with Visium — a real **per-spot** truth at the panel's cell-type
  resolution (typically 100s of genes → intermediate resolution).

## Candidate datasets (verify accessions before use)

### Bulk RNA-seq + flow/CyTOF proportions
- **Finotello et al. 2019 (quanTIseq)** — PBMC/tumour bulk RNA-seq with matched flow.
- **Newman et al. (CIBERSORT/CIBERSORTx)** — PBMC bulk with flow (whole-blood mixtures).
- **Racle et al. 2017 (EPIC)** — melanoma/PBMC bulk with flow (CD4/CD8/B/NK).
- **ImmPort SDY311 / SDY420 (AbbVie)** — bulk + CyTOF, widely used in deconvolution
  benchmarks (e.g. Monaco, Vallania).
- Any GEO series with **matched bulk RNA-seq + flow** for the same subjects.

Requirement: subjects with BOTH a bulk expression profile AND a flow/CyTOF proportion
vector, plus a single-cell reference whose labels can be mapped to the flow labels.

### Spatial (Visium) + Xenium/CosMx segmentation truth
- **10x Genomics matched Visium + Xenium human breast** (Janesick et al. 2023) —
  canonical matched dataset.
- **10x public Visium + Xenium** (other tissues as released).
- **NanoString CosMx** public datasets (e.g. lung) with cell segmentation + Visium of
  an adjacent/registered section.

Requirement: a Visium section AND a co-registered per-spot cell-type-proportion truth
derived from Xenium/CosMx segmentation (the registration/segmentation step is done
upstream; the harness consumes the resulting per-spot truth TSV).

## Local data contract (all gitignored)

```
examples/real_bulk_flow/data/
  bulk_counts.tsv         genes × samples (raw/normalised counts)
  flow_truth.tsv          samples × flow_labels (fractions, rows sum ~1)
  reference.h5ad          single-cell reference (or a prepared reference dir)
  label_map.tsv           columns: flow_label <TAB> reference_cell_type(s, '|'-sep)
  dataset_manifest.json   name, source, accession, access_date, versions, filters

examples/spatial_ground_truth/data/
  visium.h5ad             (or counts.tsv + coords.tsv)
  spot_truth.tsv          spot × cell_type (fractions from segmentation), rows sum ~1
  reference.h5ad
  label_map.tsv           segmentation_label <TAB> reference_cell_type(s)
  dataset_manifest.json
```

The harness prints these exact paths and an acquisition checklist if data is missing,
and writes a `dataset_manifest.json` template to fill.

## Metrics (truth-based, cross-subject / cross-spot)

**Bulk↔flow** (aggregate reference estimates up to the flow-label level via `label_map`):
- per-flow-label Pearson & **CCC** across subjects; overall RMSE / MAE / bias;
- compare **wNNLS (default)** vs **Poisson GLM (recommended experimental)**;
- report at the resolution the flow supports (broad lineages) — no fine-state claims.

**Spatial↔segmentation**:
- per-cell-type Pearson & CCC across spots; local RMSE; oversmoothing score vs truth;
- rare-niche behaviour where applicable; compare default / weak / edge_aware presets.

## Gates (configurable; validation = pass, not tuning)

- Bulk: per-major-lineage Pearson ≥ 0.6 (CCC ≥ 0.5), bias within ±0.05, on ≥2 datasets;
  Poisson GLM ≥ wNNLS on mean CCC (to justify promoting it to default).
- Spatial: per-cell-type Pearson ≥ 0.5 vs segmentation truth; local RMSE competitive
  with CARD/RCTD on the same section; no worse oversmoothing than reported.
- **Cross-cohort/protocol generalisation is mandatory** (≥2 datasets/platforms) — the
  conditional-estimator gate failed exactly on protocol shift, so single-dataset
  success is not sufficient.

## Non-negotiables

No auto-download; `--run-real-data` required; clear manual-acquisition instructions;
`dataset_manifest.json` recorded; outputs gitignored; **no fabricated results**; if a
tool/dataset is unavailable it is reported **deferred**, not invented; defaults
unchanged until gates pass on real data.

## Local data search (2026-07) — no usable orthogonal ground truth on-disk

A filesystem search for real ground-truth datasets found:
- Real **Visium breast** (`data/V1_Breast_Cancer_Block_A_Section_1/`,
  `examples/real_breast_cancer/data/spatial/human_breast_cancer_1.h5ad`) — **no**
  matched Xenium/CosMx segmentation truth.
- Real single-cell references only: breast atlas, HLCA lung, CRC (GSE200997). No
  ground-truth proportions.
- `~/Downloads/Planilhas/s3_cell_by_gene.csv` + `s3_mapped_cell_table.csv` — a real
  imaging-spatial dataset, but **BARISTAseq mouse cortex** (70 genes, annotated by
  cortical *layer*, not cell type): wrong species/tissue, no matching reference →
  **not usable**.
- `Xenium_hIO_v1_metadata.csv` — Xenium **panel metadata** (380 genes), not data.
- TCGA-TNBC bulk — real but **no** fine ground truth. All repo `*_truth.tsv` are
  **synthetic pseudobulk**.

**Conclusion: no orthogonal real ground truth is available locally.** Highest-value
acquisition (matches an existing reference, so the harness runs end-to-end):
**10x Genomics matched Visium + Xenium human breast (Janesick et al. 2023)** — pairs
with the existing breast single-cell reference. Fetch from the 10x portal, derive
`spot_truth.tsv` from the Xenium segmentation, place under
`examples/spatial_ground_truth/data/`. For bulk↔flow, additionally obtain a
PBMC/blood bulk+flow cohort **with a matching immune scRNA reference**.
