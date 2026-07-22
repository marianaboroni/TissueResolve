# Stage 2B — Calibrated Identifiability Certificate (report)

## OBJECTIVE
Turn the identifiability certificate from an experimental diagnostic (Stage 2A) into a *calibrated*
output, along five axes: (1) nonnegative/simplex-aware identifiability; (2) swap-based recoverability
outcome; (3) donor-level uncertainty; (4) realistic depth/panel stress; (5) a third tissue.
Calibration only — no new solver, no cohort covariance, no adaptive/multipanel work; no default changed.

## AUDIT FINDINGS
- The 2A linear certificate is an unconstrained SVD bound: conservative (over-flags confounding the
  nonnegative solver resolves) with an abundance-dominated outcome (absolute RMSE) → false-merge ~1.0
  and gene-bootstrap uncertainty that under-covered 3–4×.
- `ReferenceSignature.donor_cv` (cross-donor CV, shape (G,K)) is the natural donor-uncertainty input,
  **but only populated when the reference is built with `donor_col`** — the shared harness builder
  omitted it, so 2A ran with `donor_cv=None` (donor uncertainty silently off). Fixed here by building
  the reference with `ReferenceConfig(donor_col=...)`.
- Third tissue: CRC (GSE200997) annotates cells by **CMS molecular subtype (CMS1–4), not cell type**
  (37% unlabelled) — unusable as a deconvolution reference; no other multi-donor cell-typed reference
  is available offline (pbmc68k_reduced is single-donor).

## IMPLEMENTATION
- `src/tissueresolve/reference/identifiability_calibration.py` —
  `calibrated_identifiability_certificate(ref, *, query_detectable_genes, library_size, n_sim, noise,
  donor_uncertainty, shift_scale, tau_cluster, ...)`. Geometry (cosine ≥ `tau_cluster`) proposes
  candidate confusable clusters; a **forward simulation** disposes: draw nonnegative compositions
  (Dirichlet + cluster-stress), generate query counts from a **donor-perturbed** generating profile
  (multiplicative lognormal with the reference CV, mean-corrected) at `library_size`, add Poisson/NB
  counting noise, and solve with the **production nonnegative solver**. Outputs per-type simulated
  RMSE/bias, donor-aware predicted SD (`pred_sd`), within-cluster **swap** RMSE, group-sum RMSE, and a
  nonneg-/swap-based recoverability class. `shift_scale` inflates the donor CV to match the observed
  train→held-out profile shift (calibrated on a donor-disjoint split). No new solver: it reuses the
  existing solver inside a simulation for analysis only.
- Benchmarks: `run_stage2b_calibration.py` (breast+lung, train/cal/test, realistic depth sweep) and
  `run_stage2b_third_tissue.py` (CRC audit + reproducible readiness path).

## TESTS ADDED
`tests/test_identifiability_calibration.py` (8): orthogonal→all RESOLVABLE/no clusters; identical
pair→GROUP_ONLY (high swap, recoverable group sum); **near-collinear resolved by nonnegativity**;
three-identical group; donor uncertainty widens `pred_sd`; depth improves recovery; determinism; uses
no truth. Bug fixed with a regression: large effective CV overflowed the Poisson `lam` → correct
lognormal σ = √log1p(CV²), median-corrected, with `mu` clipped to the library size.

## TEST RESULTS
8/8 new tests pass. Full default suite: **1519 passed, 1 skipped, 0 failed** (+8, no regressions); no
production algorithm modified — only a new experimental module + tests + benchmarks were added.

## FOCUS 1 — NONNEGATIVE / SIMPLEX-AWARE
Recoverability is measured by the nonnegative solver on simulated mixtures. Loosely-correlated large
clusters are resolved (breast 21-member cluster, cosine 0.988 → RESOLVED_BY_NONNEG, swap 0.073) while
tight clusters are flagged (breast macrophage 4-cluster, cosine 0.970 → GROUP_ONLY, swap 0.180). The
certificate no longer over-flags what nonnegativity resolves (2A false-merge ~1.0 → 0.0 here).

## FOCUS 2 — SWAP-BASED OUTCOME
Recoverability outcome = within-cluster conditional (swap) RMSE, not absolute RMSE. Predicted vs
empirical swap RMSE (held-out, donor-disjoint) is tightly correlated: **breast Pearson 0.991
(Spearman 0.80, 4 clusters); lung Pearson 0.995 (Spearman 1.0, 9 clusters)**. Examples (predicted /
empirical): lung AT2;AT2-proliferating 0.230/0.295; lung Club;Goblet;pre-TB 0.279/0.346; breast
macrophage 0.180/0.184. The certificate now predicts the *magnitude* of member confusion, not just its
presence.

## FOCUS 3 — DONOR-LEVEL UNCERTAINTY
Predicted SD propagates donor CV through the simulation; a single donor-shift inflation `shift_scale`
is fit on the **calibration** donors (grid, target nominal 90%) and validated on **test**. Coverage at
nominal 90%:

| condition | gene-bootstrap (2A-style) | predicted, uncalibrated s=1 | predicted, CAL-calibrated |
|---|---|---|---|
| breast full/shallow/panel | 0.15 / 0.19 / 0.22 | 0.84 / 0.83 / 0.93 | **0.91 / 0.91 / 0.92** (s=2,2,1) |
| lung full/shallow/panel | 0.18 / 0.21 / 0.32 | 0.78 / 0.78 / 0.91 | **0.89 / 0.88 / 0.92** (s=2,2,1) |

Donor-aware predicted intervals reach ~nominal coverage; the gene-subsampling bootstrap under-covers
3–5×. The CAL-fit inflation consistently landed at s≈2.

## FOCUS 4 — REALISTIC DEPTH / PANEL
Swept full-depth (~1e6), shallow (×0.1) and a 1000-gene panel, with the certificate's `library_size`
matched to each. Swap prediction and coverage are stable across depth (full≈shallow); the panel is the
main driver of which clusters appear. (Pseudobulk caps depth at ~1e6; real bulk 1e7+ is only easier.)

## FOCUS 5 — THIRD TISSUE
**NOT EXECUTED.** CRC = CMS molecular subtypes (verified programmatically), no other donor-annotated
cell-typed reference offline. `run_stage2b_third_tissue.py` prints an exact reproducible path (drop a
multi-donor, cell-typed `.h5ad` at `examples/third_tissue/data/reference.h5ad`; it then builds the
reference and runs the calibrated certificate) — the pipeline is tissue-agnostic and ready.

## BREAST RESULTS
swap Pearson 0.991 / Spearman 0.80 (4 clusters); swap-based false-merge 0.0, false-resolution 0.0,
Brier 0.006, ECE 0.037; calibrated coverage 0.91/0.91/0.92. Canonical confounders flagged: macrophage
subset; tight epithelial/endothelial pairs at the panel.

## LUNG RESULTS
swap Pearson 0.995 / Spearman 1.0 (9 clusters); swap-based false-merge 0.0, false-resolution 0.0,
Brier 0.0004, ECE 0.02; calibrated coverage 0.89/0.88/0.92. Flagged: AT2/AT2-proliferating,
Club/Goblet/pre-TB secretory, alveolar-macrophage subsets.

## CALIBRATION RESULTS
Under the swap-based outcome with donor uncertainty and a CAL-fit shift, the certificate is calibrated
on both tissues: swap prediction near-exact (r≈0.99), zero false-merge and zero false-resolution,
Brier ≤0.006, and ~nominal (0.88–0.92) uncertainty coverage. This closes every 2A gap
(false-merge 1.0→0.0; coverage 0.15–0.32 → 0.88–0.92).

## FALSE-MERGE / FALSE-RESOLUTION RATE
Swap-based: **0.0 / 0.0** on both tissues (vs 2A absolute-RMSE 0.97–1.0 / ~0.0). The swap outcome
removes the abundance artifact; the calibrated certificate's merge and resolve decisions match held-out
recoverability.

## NEGATIVE RESULTS / LIMITATIONS
- Only two tissues; third NOT EXECUTED (no data).
- Near-nominal coverage requires a CAL-fit `shift_scale` (≈2 here); a shipped default needs
  cross-tissue validation of that default.
- `donor_cv` is only available when the reference is built with donors (≥2); single-donor references
  fall back to counting-noise-only uncertainty (under-covers).
- Simulation uses one solver family (production Poisson) and Poisson/NB noise; clusters depend on the
  cosine threshold (pre-specified 0.90, not tuned on test).

## DECISIONS
**1. Is the certificate calibrated enough to become a core output? → B.**
It clears the calibration bar on both tissues (swap r≈0.99, zero false-merge/false-resolution,
~nominal coverage) — a decisive improvement over the 2A experimental diagnostic. Promote it to a
first-class, shipped output (donor-aware uncertainty + swap-based classes), reported alongside
estimates. Reserve full **default-on core** status (decision A) until (i) a third tissue confirms and
(ii) a validated default `shift_scale` (≈2) ships so it calibrates without a user calibration split.
Not C — the evidence is strong and reproducible. Per instruction, no default was changed this stage.

**2. Does cohort covariance justify future cohort-aware bulk modeling? → E.**
Unchanged from Stage 2A and out of 2B scope: with a known reference, cohort cross-sample covariance
adds no identifiability (null-space gain 0.0; beaten by the nonnegative per-sample solver). Do not
build cohort-aware bulk modeling as an identifiability mechanism.

## NEXT RECOMMENDED STAGE
Stage 2C — promotion to default core: validate a default `shift_scale` across tissues; add a third
tissue (drop-in path ready); wire the calibrated certificate as a reported core output of the bulk
pipeline (recoverability class + swap risk + donor-aware interval per type) without changing
deconvolution estimates; document the recommended_merge guidance for GROUP_ONLY clusters.
