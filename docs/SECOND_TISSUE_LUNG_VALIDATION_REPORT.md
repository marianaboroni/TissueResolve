# Second-Tissue Validation — HLCA Lung (Final Report)

Independent second-tissue validation on the **Human Lung Cell Atlas (HLCA) core**
— the open + labeled atlas secured after the CRC routes failed (Pelka gated;
GSE200997 unlabeled). **Both promotion blockers are now addressed.** Nothing
committed; the 5.9 GB atlas is git-ignored.

## Dataset (open + labeled, verified from primary sources)
HLCA core (CZ CELLxGENE, CC BY 4.0): 584,944 cells × 27,402 genes, **raw counts**
(`adata.raw.X`, integer-verified), **107 donors**, manual consensus annotations
`ann_level_1..5` + `ann_finest_level`. Download was **truncated once (2.16 of
5.87 GB)** — caught by a size check (h5py `truncated file`), the downloader
hardened to verify Content-Length + resume, then completed cleanly. Manifest:
`examples/second_tissue_lung/data/download_manifest.json`.

**Validation atlas (subsampled, donor-disjoint capable):** 11,369 cells, **40
donors**, **broad = ann_level_2 (10 families)** → **fine = ann_finest_level (59
subtypes; 2 low-support flagged, not dropped)**, **9 families with ≥2 fine
subtypes** — structurally analogous to the breast hierarchy (8 families / ~30
subtypes). Hierarchy documented in `data/derived/hlca_hierarchy.tsv`.

## Part A — identifiability pattern GENERALIZES
Per within-family pair (donor-aware; deconvolution-relevant metrics):

| family | n pairs | median signature Pearson | median shared-lineage frac | median ratio-recovery MAE |
|---|---|---|---|---|
| Smooth muscle | 3 | 0.93 | 0.983 | 0.003 |
| Blood vessels | 10 | 0.90 | 0.974 | 0.004 |
| Fibroblast lineage | 15 | 0.88 | 0.970 | 0.004 |
| Myeloid | 78 | 0.76 | 0.940 | 0.003 |
| Lymphoid | 15 | 0.69 | 0.928 | 0.003 |
| Airway epithelium | 120 | 0.70 | 0.925 | 0.002 |
| Submucosal Gland | 6 | 0.23 | 0.871 | 0.001 |

**Median shared-lineage fraction = 0.940** (range 0.85–0.98) — matches breast
(0.86–0.99). Pairwise mixture recovery is again near-perfect (MAE ~0.003), i.e.
subtypes are collinear-but-recoverable in isolation; the difficulty is the
shared-lineage-dominated full-panel problem. **The breast diagnosis generalizes to
lung.**

## Part B — soft gating REPLICATES (8/8 prospective gates on lung)
Held-out-donor pseudobulk, cal seeds {0,1} / **test {2,3,4}** (donor- + seed-disjoint),
3 scenarios. Mode means (test):

| mode | fine Pearson | fine RMSE | broad RMSE | false-abstention | eff-N err | rare sens. |
|---|---|---|---|---|---|---|
| flat | 0.663 | 0.031 | 0.035 | — | 5.41 | 0.420 |
| ungated | 0.431 | 0.049 | 0.127 | 0.000 | 8.57 | 0.454 |
| **soft_gate** | 0.431 | **0.046** | 0.127 | 0.098 | 8.56 | 0.437 |
| hard_gate (current) | 0.216 | 0.050 | 0.127 | 0.634 | 15.66 | 0.118 |

Gates: vs ungated — broad/fine RMSE non-inferior (fine *better*), false-resolution
non-inferior, mass conserved ✅; vs hard — false-abstention 0.098 vs 0.634, rare
sensitivity 0.44 vs 0.12 (3.7×), eff-N error 8.6 vs 15.7, fine RMSE lower ✅.
**8/8 PASS** — the exact breast pattern (breast 9/9). **The 2-tissue replication
requirement for soft-gating promotion is now MET.**

## Part C — weaker-λ spatial finding GENERALIZES (modestly)
HLCA-derived synthetic spatial (5 seeds), weak λ=0.02 vs default λ=0.1:

| λ | fine RMSE | broad RMSE | oversmoothing | boundary F1 |
|---|---|---|---|---|
| default 0.1 | 0.074 | 0.144 | 1.64 | 0.271 |
| **weak 0.02** | **0.071** | **0.116** | **1.42** | **0.310** |

Weak λ is **≤ default on every metric** (broad RMSE and over-smoothing clearly
better; boundary F1 better; fine RMSE marginally better, Wilcoxon p=0.31 at n=5 —
not significant). Consistent with breast. **The weaker-λ default change is
supported on a second tissue** (at worst non-inferior; reduces over-smoothing on
both tissues). (One spatial seed hit the 200-iter cap at δΠ=6.7e-5 — minor.)

## Promotion status (both blockers cleared)
1. **Soft gating: PROMOTABLE.** Passed prospective gates on **two tissues**
   (breast 9/9, lung 8/8), donor-disjoint, mass-conserving. Recommend promoting it
   to the hierarchical default (re-specifying the original false-resolution gate as
   soft-vs-ungated, as recorded). **Not auto-applied** — awaiting explicit go-ahead
   (rule: no default change without sign-off).
2. **Weaker spatial λ (0.1→~0.02): SUPPORTED.** Generalizes to a second tissue
   (synthetic spatial). Recommend as a candidate default; ideally confirm on a
   real Visium with boundary ground truth before finalizing (modest fine-RMSE gain).

## Limitations
HLCA is healthy lung scRNA (no native Visium) — the spatial check uses synthetic
spatial; both tissues are synthetic-pseudobulk/synthetic-spatial for ground truth.
Lung subset capped at 40 donors / 200 cells per fine type for tractability. n=5
seeds for the spatial λ check (fine-RMSE improvement not individually significant).
Real-Visium-with-truth and image-based pseudo-spots remain the strongest future
confirmation.

## Deliverables (nothing committed; data git-ignored)
`examples/second_tissue_lung/` (README, `scripts/00_download_hlca.py`,
`01_audit_and_subsample.py`, `02_validate_generalization.py`,
`03_validate_spatial_lambda.py`, manifest, derived hierarchy/audit tables).
Outputs: `lung_identifiability_{pairs,family_summary}.tsv`,
`lung_softgating_{modes,gates}.tsv`, `lung_spatial_lambda.tsv`. Docs:
`PELKA_CRC_DATASET_VERIFICATION.md`, `SECOND_TISSUE_CRC_AUDIT.md`, this report.
