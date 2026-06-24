# Second-Tissue CRC — Metadata/Label Audit (Stage 3) — GATE FAIL

Audit of the committee-approved open CRC substitute **GSE200997** (downloaded
openly; manifest in `examples/second_tissue_crc/data/download_manifest.json`).
Per the rule "do not assume published annotations are suitable; audit first."

## What the data actually contains
- **Counts:** raw UMI matrix, 41,364-gene × 49,859-cell, integer — ✅ raw counts OK.
- **Donors:** `samples` column, **23 samples** (B_cac… normal ×7, T_… tumor ×16) — ✅ donor IDs present.
- **Annotation columns:** `samples`, `Condition` (Normal/Tumor), `Location`
  (Left/Right), `MSI_Status` (MSS/MSI-H), `bulk_prediction`+`prediction` (CMS
  molecular subtype). The cell-ID prefix is **B/T = benign/tumor**, not a cell type.

## Critical failure
**There are NO cell-type annotations (broad or fine).** The only per-cell labels
are tumor-level molecular metadata (CMS subtype, MSI, location) — not cell
identities. Therefore this dataset **cannot** support:
- the subtype-identifiability ceiling (needs within-family fine labels);
- the broad/fine reference build;
- the soft-gating broad+fine validation.

| Requirement | Verified | Decision |
|---|---|---|
| Raw counts | ✅ yes (integer UMI) | pass |
| Donor IDs | ✅ yes (23) | pass |
| **Broad cell-type labels** | ❌ **absent** | **FAIL** |
| **Fine cell-type labels** | ❌ **absent** | **FAIL** |
| Donor count adequacy | ⚠ 23 (modest) | marginal |

## Decision (rules: stop on critical failure; no fabrication; no silent self-annotation of fine types)
**GSE200997 fails the second-tissue gate at the label-audit stage.** The two CRC
routes are mutually exclusive in practice:
- **Pelka 2021 *Cell*** — has deep cTNI broad+fine labels, but **credential-gated**
  (SCP sign-in + EGA-controlled) → unacquirable in this no-credentials session.
- **GSE200997 (open)** — anonymously downloadable raw counts, but **no cell-type
  labels** → unusable for subtype validation.

De-novo fine-subtype annotation of GSE200997 would introduce **non-gold-standard
labels** and risks exactly the tissue-specific/circular overfitting the program
forbids; coarse self-annotation could give *broad* compartments only, not the
within-family fine resolution the validation targets. Not done.

## Options forward
1. **Credentialed Pelka acquisition** (you provide SCP/EGA access; I validate on the
   real labeled atlas). Restores the intended fine-label validation.
2. **Re-verify a different OPEN + LABELED atlas** — the audit lesson is to require
   *cell-type labels present AND open* up front. Re-check NSCLC (e.g. an open HLCA/
   NSCLC atlas on CELLxGENE with fine labels) or an open HGSOC atlas; verify both
   open raw counts AND broad+fine annotations before download.
3. **Defer the second tissue; proceed to Task D** (level-specific spatial smoothing)
   which needs no external data. Soft gating stays experimental (promotion blocked).

## Status
Acquisition scaffold (`examples/second_tissue_crc/`) + open download workflow are
built and reusable for any labeled dataset. **Validation NOT run** (no labels).
Manifest `processing_status` updated to `unsuitable_no_celltype_labels`. Nothing
committed; data git-ignored.
