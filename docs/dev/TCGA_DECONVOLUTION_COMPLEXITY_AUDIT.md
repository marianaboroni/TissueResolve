# TCGA Deconvolution Complexity Audit

Purpose

This file summarises the targeted audit to determine whether the observed
low compositional complexity (4–5 detected populations per sample) in TCGA
bulk predictions is an artefact of the pipeline (solver, thresholds,
post-processing) or a plausible biological signal.

Procedure

1. Run the audit script in `benchmarks/bulk/tcga_sparsity_audit.py` to
   compute per-sample and per-method metrics from available prediction TSVs.
2. Inspect `benchmarks/outputs/tcga_prediction_complexity_audit_by_sample.tsv`
   to identify samples with extremely low effective numbers or high Gini.
3. For flagged samples, examine the raw prediction TSV (no postprocessing) to
   verify whether mass was zeroed or aggregated into `unresolved_` or `Other`.
4. Validate that the solver output saved by TissueResolve (if available) matches
   the prediction TSV used by the report generator (no silent post-processing).

Key checks

- Count of non-zero raw coefficients vs displayed counts in reports
- Presence of `unresolved_` columns and their mass
- Marker-panel size per family (are too few genes used?)
- Reconstruction R² per sample (low R² suggests mismatch, not sparsity)
- Condition number of the panel (collinearity can drive sparse NNLS)

Next steps after the audit

- If sparsity is due to post-processing (thresholding/top-N/Other): adjust
  report settings and rerun reporting only; do not change solver defaults.
- If sparsity arises from hierarchical unresolved gating: inspect
  `hierarchical` config thresholds and per-family separability metrics.
- If sparsity arises from solver (regularisation/normalisation): capture
  raw W_hat matrix before row L1-normalisation and compare.
# TCGA Deconvolution Complexity Audit (PART 17)

**Scope.** A dedicated analysis of the real TCGA-BRCA TNBC bulk predictions to
determine whether TissueResolve is *uniquely* sparse or whether the ~4–5
population observation is a generic property of constrained deconvolution on this
reference. Companion to `docs/BULK_PREDICTION_SPARSITY_AUDIT.md` (mechanism) — this
document is the per-sample, per-method TCGA evidence.

**Source data.** `benchmarks/outputs/tcga_prediction_complexity_audit.tsv`
(320 rows = 8 methods × 40 samples), generated deterministically by
`benchmarks/audit/tcga_sparsity_audit.py --run-real-data`. Cohort: 40 TCGA-TNBC
samples × 59 427 genes, 4 982 shared with the 32-population reference.

> **No ground truth.** Real TCGA tumours have no measured cell proportions. This
> audit reports *composition complexity* (how concentrated each prediction is),
> **never accuracy.** Accuracy requires the pseudobulk harness (Phase 2).

---

## 1. Per-method complexity across the 40-sample cohort

| Method | mean non-zero | median | eff-N (mean) | Gini | dominant frac | top-5 cum | unresolved mass |
|---|---|---|---|---|---|---|---|
| `flat_nnls`          | 2.8  | 3  | 2.1  | high | 0.64 | 1.00 | 0.00 |
| `flat_weighted_nnls` | 3.9  | 4  | 2.9  | high | 0.56 | 1.00 | 0.00 |
| `flat_marker_nnls`   | 11.4 | 12 | 7.0  | mid  | 0.32 | 0.83 | 0.00 |
| `flat_ridge_nnls`    | 22.9 | 23 | 17.5 | low  | 0.13 | 0.45 | 0.00 |
| `flat_auto` (→ridge) | 22.9 | 23 | 17.5 | low  | 0.13 | 0.45 | 0.00 |
| `hier_family` (8)    | 3.0  | 3  | 2.4  | —    | 0.58 | 1.00 | 0.27 |
| `hier_fine` (32)     | 2.8  | 3  | 1.8  | —    | 0.77 | 1.00 | 0.27 |
| `hier_combined` (36) | 4.0  | 4  | 2.7  | —    | 0.57 | 1.00 | 0.27 |

(Per-sample detail incl. n>1e-4/1e-3/5e-3/1e-2 and Shannon entropy in the TSV.)

## 2. Is TissueResolve *uniquely* sparse?

**Not as a tool — it is the chosen solver that is sparse.** Within TissueResolve
itself the spread is enormous on identical TCGA input: `nnls` → 2.1 effective
populations vs `ridge`/`auto` → 17.5. So "TissueResolve" is neither inherently
sparse nor inherently diffuse; the answer is entirely set by the solver and the
near-collinear reference (condition ≈ 486). The **shipped `auto` TCGA output is
the diffuse 23-population one**, not the sparse 4–5.

NNLS-driven sparsity on collinear references is a **known, general** behaviour of
constrained deconvolution (it is the reason MuSiC weights genes, DWLS reweights,
BayesPrism uses a full Bayesian prior, and CIBERSORTx uses ν-SVR). We therefore
*expect* unregularised per-gene NNLS to be the sparsest method and Bayesian /
regularised methods to be denser — but this must be **measured**, not assumed.

## 3. External-method comparison — STATUS: NOT YET EXECUTED

The brief asks to compare TissueResolve against **NNLS, MuSiC, BisqueRNA,
BayesPrism, CIBERSORTx**. Honest status (per the rule "Do not present skipped or
exported-only tools as executed"):

| Method | Status | Notes |
|---|---|---|
| NNLS (internal) | **Executed** | `flat_nnls`, this audit. |
| weighted/marker/ridge NNLS | **Executed** | internal baselines, this audit. |
| MuSiC | **Pending (Phase 2)** | R wrapper present (`benchmarks/bulk/methods/run_music.R`); not yet run in a controlled env. |
| BisqueRNA | **Pending (Phase 2)** | wrapper present (`run_bisque.R`); not run. |
| BayesPrism | **Pending (Phase 2)** | wrapper present (`run_bayesprism.R`); not run. |
| CIBERSORTx | **Export-only** | license-gated; `cibersortx_export.py` prepares inputs only — will be labelled *imported* if results are returned, never as *executed* here. |

The cross-method complexity comparison (does every method concentrate mass on
TCGA, or only NNLS?) will be completed in Phase 2 once these run in an isolated
environment with recorded versions. **Until then, no claim is made that
TissueResolve is more or less sparse than external tools.**

## 4. Where the report would (and would not) hide populations

Audited explicitly (see sparsity audit §3): on TCGA the report does **not** hide
populations via top-N display, abundance threshold, or "Other" aggregation — all
are default-off / display-only, and `auto`'s 23 populations pass through intact.
The only place mass is *moved out of fine subtypes* is the hierarchical
**unresolved-family gating** (27% of mass on TCGA), which is by design and is
reported as `unresolved_<family>` columns — not silently dropped.

## 5. Verdict

- The TCGA "4–5 populations" is **reproducible and explained**: it is the
  flat-`nnls`/`weighted_nnls` solution (≈3–4) and/or the hierarchical
  family+gating result (≈4) — a **solver/conditioning** phenomenon, not a
  threshold, "Other"-merge, or display artefact.
- It is **not yet established as biological.** The shipped `auto` output is in
  fact diffuse (≈23). Which level of concentration is correct is a Phase-2
  ground-truth question.
- **No default changed** at this stage (`CLAUDE.md` rule 10).

---

*Source: `benchmarks/outputs/tcga_prediction_complexity_audit.tsv`. Raw
per-method proportions: `benchmarks/outputs/raw/`. Generator:
`benchmarks/audit/tcga_sparsity_audit.py`.*
