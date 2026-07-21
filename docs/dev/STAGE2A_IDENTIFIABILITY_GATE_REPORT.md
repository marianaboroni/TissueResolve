# Stage 2A — Identifiability Gate Report (finalized)

## SCIENTIFIC QUESTION
Can TissueResolve predict, **before** deconvolution, which populations are recoverable and which
combinations are structurally indistinguishable — from the (reference, query) pair alone, no ground
truth? And: does modelling a **cohort's** cross-sample covariance recover information that
single-sample Poisson deconvolution cannot when the reference R is known? This is a *gate* deciding
whether the differentiating direction is **Pillar A** (identifiability-first + calibration) alone, or
Pillar A **+ Pillar B** (cohort / distribution-aware modelling). No core algorithm, no default changed.

---

## AUDIT FINDINGS
- The reference exposes `as_R_cpm()` → (K, G) means, `gene_names`, `cell_types`, `selected_genes`,
  plus `donor_cv` and `n_cells_per_type` (usable as reference-uncertainty inputs).
- No existing identifiability/effective-rank analysis; separability is only pairwise Bhattacharyya on
  mean profiles (`reference/separability.py`) — panel-dependent, not noise-aware, not θ-space.
- The production per-sample solver (`PoissonGLMSolver`) has no built-in uncertainty; a
  gene-subsampling bootstrap was added **in the benchmark only** (validation tooling).
- `split_donors` is 2-way; a donor-disjoint 3-way (train/calibration/test) split was derived from it.
- Third tissue (CRC) has only raw CSVs, no donor-annotated h5ad → **NOT EXECUTED**.

## IMPLEMENTATION (experimental, feature-flagged, wired into no default)
- `src/tissueresolve/reference/identifiability.py` — `identifiability_certificate(...)`: SVD of R
  restricted to query-detectable genes; condition number, stable rank, effective rank; noise-aware
  effective rank; near-null → confounded groups; per-type recoverability class; detection floor.
  Uses **no truth** (asserted by test).
- `src/tissueresolve/experimental/cohort_covariance.py` — 4 estimators (single / first_moment /
  cross_covariance [empirical Bayes] / oracle) for Experiment B. Not registered; changes nothing.
- Benchmarks: `benchmarks/signatures/run_stage2a_certificate.py` (Experiment A, 3-way split),
  `benchmarks/signatures/run_stage2a_cohort.py` (Experiment B; synthetic grid + real cohort).

## TESTS ADDED
16 synthetic regression tests. `tests/test_identifiability.py` (9): orthogonal→RESOLVABLE,
identical→UNRESOLVABLE/GROUP + correct merge, near-collinear depth-dependence, 3-type group,
query-gene-loss collapse, depth raises noise-adjusted rank, few-donors→NOT_TESTABLE, determinism,
no-truth-used. `tests/test_cohort_covariance.py` (7): null_direction finds confounded axis, no
information without variation, **no gain on null space even with variation**, gain only where
mathematically supported (shrinkage under noise), oracle=truth, bad-method rejected, determinism.

## TEST RESULTS
16/16 new tests pass. Full default suite after this stage's source additions: **1511 passed, 1
skipped, 0 failed**; no production source was modified during Stage 2A (only new experimental modules
+ benchmarks + tests were added), so the pipeline's behaviour is unchanged.

## IDENTIFIABILITY MODEL
Deconvolution solves y ≈ R θ. The recoverable information is the projection of θ onto the
well-conditioned, above-noise subspace of R (restricted to query-detectable genes); the rest is
structurally indistinguishable for any solver. Right singular vectors of R live in θ-space; near-null
directions are the confounded cell-type combinations.

## NOISE-AWARE EFFECTIVE RANK
Per singular direction v_i, a matched-filter SNR at the query library size L over baseline counts
b = L·p̄/1e6: snr = √Σ(s_i²/b), s_i = L·(R_q v_i)/1e6; detection floor = z/snr (minimal detectable
abundance shift). A direction is recoverable iff floor < 10%. Noise-adjusted rank = # recoverable
directions. This makes recoverability **depth dependent**: at ~1e6 pseudobulk counts breast keeps all
30 and lung all 49 directions above noise (condition numbers 222 / 147); non-resolvability appears
only when stressed (spot-like depth, or a 1000-gene panel).

## NULL-SPACE RESULTS
Confounded groups (connected components over shared near-null directions) are **biologically
canonical**: breast (spot-like) → an 8-member **lymphoid** group (CD4 helper, CD4 αβ, CD8 memory,
effector-memory CD8, lymphocyte, NK T, NK, Treg); lung (panel/spot) → {EC general capillary; EC
venous pulmonary}, {Club (non-nasal); Migratory DCs}, {AT0; Monocyte-derived Mph; Smooth muscle}, and
a large alveolar/macrophage group. Within-group **swap** RMSE (member/group-sum, abundance-normalised)
confirms the split among members is unrecoverable: breast 0.179 (vs 0.033 individual); lung small
groups 0.34–0.36; the diluted 24-member group 0.08.

## RECOVERABILITY CLASSES
RESOLVABLE (rec-fraction ≥ 0.80, not in a group) · WEAKLY_RESOLVABLE (≥ 0.50) · GROUP_ONLY (individual
axis near-null but group sum recoverable) · UNRESOLVABLE (type and group near-null) · NOT_TESTABLE
(too few query genes or reference donors). Per type: closest_confounders, recommended_merge,
recoverable_group, detection floor, reason.

## BREAST RESULTS (TEST donors; RESOLVABLE=108, GROUP_ONLY=8 across conditions)
Class predicts held-out error only weakly in **absolute** terms (held-out RMSE GROUP_ONLY 0.033 ≈
RESOLVABLE 0.036 — abundance-dominated, non-monotone), but the **direction** metrics behave: absent-type
false positives GROUP_ONLY 0.146 > RESOLVABLE 0.098; within-family conditional RMSE uniformly high
(~0.15–0.19, not class-discriminating — hierarchy families ≠ certificate groups); rare recall flat
(~0.62 both). RESOLVABLE reliable: false-resolution 0.019.

## LUNG RESULTS (TEST donors; RESOLVABLE=152, GROUP_ONLY=31, WEAKLY=7, UNRESOLVABLE=2)
Cleaner monotone ordering by class: held-out RMSE 0.024 < 0.025 < 0.040 < 0.060; **absent-type false
positives 0.051 < 0.064 < 0.139 < 0.389** (strong); spillover 0.011 < 0.013 < 0.022 < 0.050. Detection
recall higher for RESOLVABLE (0.569 vs 0.412 GROUP_ONLY); rare recall noisy/non-monotone. RESOLVABLE
reliable: false-resolution 0.000.

## CALIBRATION RESULTS
Class→P(recoverable) fit on **calibration** donors, evaluated on **test**: Brier 0.017 (breast) /
0.014 (lung); ECE 0.003 / 0.024. These look excellent but partly **trivially**: under the absolute-RMSE
recoverability outcome (RMSE < 0.10), almost every type is "recoverable" at natural abundances (even
confounded ones, which are individually low-abundance), so a near-1 probability is well-calibrated
against an almost-always-true label. The scientifically meaningful outcome is the within-group **swap**,
not absolute RMSE.

## FALSE-MERGE RATE
**1.00 (breast) / 0.97 (lung)** under the absolute-RMSE criterion — i.e., nearly all types flagged
UNRESOLVABLE/GROUP_ONLY are "recoverable" by absolute RMSE. This is a **measurement artifact** of the
abundance-dominated criterion: the abundance-normalised swap metric (Null-space results) shows the
within-group split is genuinely unrecoverable (0.18–0.36). The certificate over-flags because the
*nonnegative* solver recovers more than an unconstrained-linear analysis predicts.

## FALSE-RESOLUTION RATE
**0.019 (breast) / 0.000 (lung)** — when the certificate says RESOLVABLE it is essentially always
right. This is the safe direction and the certificate's strongest property.

## RESOLVED COVERAGE
Of the empirically recoverable test types, the fraction flagged RESOLVABLE/WEAKLY: **0.93 (breast) /
0.83 (lung)**. The certificate captures most recoverable types while (per false-resolution) rarely
promoting an unrecoverable one.

## UNCERTAINTY COVERAGE
Gene-subsampling bootstrap (nominal 90%) **severely under-covers**: empirical ~0.24–0.38, consistent
between calibration and test (breast overall cal 0.38 / test 0.26; lung 0.31 / 0.28). Class ordering
is monotone on lung (RESOLVABLE 0.37 > GROUP_ONLY 0.25 > WEAKLY 0.24 > UNRESOLVABLE 0.13) and flat on
breast. Interpretation: gene-resampling captures only a fraction of the true error (donor/biological
variation and model misspecification dominate); solver uncertainty is **not calibrated**, and the
certificate class predicts coverage only weakly. Uncertainty quantification is an open problem.

## COHORT COVARIANCE EXPERIMENT
Synthetic factor grid (collinearity × variation × depth × N): cross_covariance − single gain on R's
**null direction is exactly 0.0 in every confounded condition** (constant, clustered, or independently
varying θ; low/high noise). Overall gain is positive only under high noise + low true variance
(well-conditioned, σ=0.5 → +0.04–0.06 = pure variance reduction/shrinkage) and ≤0 when θ genuinely
varies. Real cohort (breast+lung, 1000-gene panel so confounding exists): cross_covariance ≈
first_moment ≈ single (breast ~0.048, lung ~0.038); paired cross_covariance − single within ±0.001,
win-rate → 0% by N=50 (no gain, no trend with N); confounded-pair swap RMSE identical (lung 0.061 =
0.061); the production **single-sample Poisson GLM beats every linear cohort method** (nonnegativity
> covariance). N capped at 50 by pool size; N=100 not reached.

## NEGATIVE RESULTS
(1) Cohort covariance adds **no** recoverable information when R is known. (2) Absolute held-out RMSE,
within-family conditional RMSE, and rare recall do **not** cleanly track the recoverability class —
they are abundance/depth-dominated, not identifiability-dominated, at these abundances. (3)
Gene-bootstrap uncertainty is badly **under-calibrated**. (4) At full depth both references are
well-conditioned (no confounded groups) — identifiability is not the full-depth bottleneck.

## KNOWN LIMITATIONS
Pseudobulk only (no real bulk with fine ground truth); two tissues (CRC not prepared); the certificate
is an unconstrained-linear **lower bound** on the nonnegative solver's recoverability (hence
conservative); the absolute-RMSE recoverability outcome is ill-posed (swap error is the right
criterion); uncertainty via gene subsampling is uncalibrated; cohort estimators are linear-Gaussian;
real cohort N ≤ 50.

## PARADIGM DECISION
**Certificate → B.** It is *sound where it matters* (false-resolution ≈ 0; predicts absent-type false
positives and within-group swaps; biologically canonical groups) but **not yet calibrated enough to be
a core, on-by-default output**: uncertainty is uncalibrated, the recoverability outcome needs the
swap-based definition, and behaviour is depth/panel dependent. Ship as an **experimental / opt-in
diagnostic** while the calibration refinement is done; do not make it core or change any default.
**Cohort covariance → E.** With a known reference it does not add identifiability (null-space gain 0.0
everywhere; no real-data gain; beaten by the nonnegative per-sample solver); any benefit is variance
reduction a simpler per-sample regulariser would capture. **Do not build Pillar B for identifiability.**
Net paradigm: pursue **Pillar A (identifiability-first + calibration)**; do not pursue Pillar B.

## NEXT RECOMMENDED STAGE
Certificate calibration refinement: (1) nonnegativity/simplex-aware identifiability so the certificate
stops over-flagging confounding the nonnegative solver actually resolves; (2) adopt within-group swap
RMSE (not absolute RMSE) as the recoverability outcome and re-fit the class→probability map on it; (3)
a calibrated uncertainty estimator (the gene bootstrap under-covers 3–4×) with donor-level resampling;
(4) an explicit realistic-bulk-depth sweep; (5) a third tissue once a donor-annotated reference exists.
Only after that re-audit whether the certificate qualifies as a core output (decision A).
