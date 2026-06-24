# Second-Tissue Selection & Validation Plan (Stage 6)

A real second tissue is **required before promoting soft gating or any new
signature strategy**. None is available in this environment (only the breast
atlas; PBMC68k is unusable). This plan specifies acquisition + validation; it does
not download data (rule 15) and changes nothing.

## Why PBMC68k is excluded
No donor labels (cannot do donor-disjoint splits, rule 4); not raw counts (scaled,
X max ≈ 28 — breaks count-level pseudobulk); 700 cells / 765 genes (too small);
single biological compartment (no epithelial/stromal/immune diversity). Using it
would violate the validation design.

## Required properties for the second tissue
raw UMI counts · ≥10 donors · broad + fine labels · ≥50 cells per fine subtype in
≥3 donors · multiple compartments (epithelial/stromal/immune) · independent tissue
or cancer type · gene symbols overlapping typical bulk/Visium · ideally both
pseudobulk and Visium validation potential.

## Candidate domains (in rough priority)
1. **Ovarian cancer** (HGSOC atlases; rich TME, CELLxGENE-available).
2. **Lung cancer** (NSCLC atlases; many donors, fine immune/epithelial labels).
3. **Colorectal cancer** (epithelial + stroma + immune; matched Visium exists).
4. Kidney / pancreas / adipose (well-annotated normal references).

## Acquisition protocol (when authorised)
- Pull from CELLxGENE Census (stable dataset IDs, donor_id + cell_type present).
- Record: dataset name, source URL/ID, version, access date, n donors, n cells,
  label columns (broad + fine), protocol(s), license, filters, downsampling.
- Build a `download_manifest.json` (per CLAUDE.md real-data rules).
- Harmonise to gene symbols; define a fine→broad hierarchy mapping analogous to
  the breast one.

## Splits (donor-disjoint; rules 4,10,11)
reference donors → signatures; calibration donors → confidence/threshold;
**test donors** → final evaluation only. Mirror the breast harness
(`synthetic_holdout`, `phase2b_stage1_validation`) with the new atlas.

## Validation to run on the second tissue
1. **Soft-gating replication** (the promotion gate): re-run the Stage-1 prospective
   gates (soft vs ungated/hard) on the second tissue. Soft gating is promoted to the
   hierarchical default **only if it replicates** (broad/fine RMSE non-inferior to
   ungated, false-abstention ≪ hard, mass conserved, rare sensitivity ↑).
2. **Identifiability ceiling** (Stage 1 probes) on the new families — confirm
   whether shared-lineage dominance generalises.
3. (If Stage 3 built) contrastive-signature gates on the second tissue.

## Stop/go
- Soft gating: **not promoted** until replication on ≥1 real second tissue passes.
- Any new signature method: experimental until it improves on ≥2 tissues.
- If replication fails: keep soft gating experimental; report the tissue-specificity.
