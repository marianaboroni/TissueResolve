# Second-Tissue Validation Harness (CRC)

Acquisition + audit scaffold for independent second-tissue validation of
TissueResolve. Built for the committee-approved **open CRC substitute GSE200997**.

**Status: blocked at metadata audit.** GSE200997 downloads openly with raw counts +
donor IDs but has **no cell-type annotations** (only sample/condition/CMS metadata)
→ unusable for subtype-identifiability/soft-gating validation. See
`docs/SECOND_TISSUE_CRC_AUDIT.md`. The scaffold is **reusable** for any labeled +
open atlas (point the scripts at a new accession).

## Scripts
- `scripts/00_download_data.py` — open, resumable, checksummed download + manifest
  (`--dry-run`/`--force`; no credentials; data git-ignored). Works now for GSE200997.
- `scripts/01_validate_metadata.py` — Stage-3 audit (donors, labels, raw-count
  integrity, low-support flags, proposed fine→broad mapping). **Run; found no labels.**
- `scripts/02_prepare_reference.py` — build TissueResolve reference (to add once a
  labeled dataset is chosen; mirrors `examples/real_breast_cancer/scripts/01_*`).
- `scripts/03_create_splits.py` — donor-disjoint train/calibration/validation/test
  (to add; no donor crosses splits; no test-donor gene selection/tuning).

## Data
`data/` is git-ignored (rule 11). `data/download_manifest.json` records source,
accession, checksums, sizes, date, license, and processing status.

## Next step (committee decision)
1. credentialed Pelka acquisition (labeled, gated), or
2. re-verify an OPEN + LABELED atlas (require cell-type labels up front), or
3. defer second tissue → Task D (spatial smoothing).
