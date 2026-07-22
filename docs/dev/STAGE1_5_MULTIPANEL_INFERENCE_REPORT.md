# Stage 1.5 — multipanel inference & budget calibration report

## QUESTION
Does TissueResolve benefit from the multipanel architecture (broad / fine_global / sibling /
rare_confirmation), or only from the `fine_global` donor-aware panel? And are the earlier
panel comparisons confounded by budget/implementation bugs?

## IMPLEMENTATION AUDIT
`ReferenceSignatureOptimizer` audited. Found: (1) `optimize()` permanently set `self.dc=None`
on a no-donor input (state mutation); (2) gene budget applied as a flat `genes[:N]` on
type-grouped panels (starves later-listed cell types); (3) `rare_confirmation` was sibling-only
(not combined evidence); (4) panels labelled bulk_or_spatial despite only pseudobulk/Poisson
validation; (5) intended-use safeguard existed but the four-panel roles/manifest were incomplete.

## BUGS FOUND
- **§3.1 state mutation** — `self.dc` permanently disabled after a no-donor run (later
  donor-aware runs silently degraded).
- **§3.3 flat budget truncation** — `genes[:budget]` on a per-type-grouped panel dropped the
  last cell types' markers; this **confounded the Stage-1 lung equal-budget result**.

## CORRECTIONS
- state mutation → try/finally restore of `self.dc` (regression test: no-donor then donor).
- flat truncation → **stratified round-robin** allocation (min-per-type coverage + budget
  redistribution + per-type `budget_shortfall` records).
- `rare_confirmation` → **combined evidence** (sibling-supported ∩ fine-global marker set).
- honest labels → `validated_modality="bulk_experimental"`, `spatial_status="spatial_unvalidated"`.
- manifest → `budget_strategy` + `fine_global_shortfall` recorded.

## BENCHMARK IMPLEMENTATION
New permanent, versioned benchmark `benchmarks/signatures/` (replaces /tmp scripts):
`validation.py` (guards that ABORT), `run_multipanel_benchmark.py` (donor-disjoint, stratified
equal budget, Poisson, paired-by-mixture), `tests/` (offline), README, config. Outputs
(gitignored) in `benchmarks/results/signatures/`.

## DATA LEAKAGE CHECKS
Guards raise on: wrong-tissue hierarchy (regression-tested: breast hierarchy on lung types
aborts), any fine type unmapped, train/test donor overlap, selection-donor leakage, unrecorded
budget, mis-oriented matrices, truth/pred label mismatch. All wired into the runner.

## GENE BUDGET RESULTS (stratified, budget 300, 5 seeds, donor-disjoint)
| tissue | strategy | RMSE | cond_RMSE | rare_fpr | rare_recall |
|---|---|---|---|---|---|
| breast | current | 0.0365 | 0.210 | 0.280 | 0.430 |
| breast | fine_global | **0.0340** | **0.183** | **0.268** | **0.480** |
| lung | current | **0.0313** | 0.242 | 0.219 | 0.384 |
| lung | fine_global | 0.0340 | 0.242 | 0.222 | 0.350 |

The stratified fix improved lung fine_global from **0.043→0.034 RMSE** vs the flat-truncation
Stage-1 number — the truncation was a real confound. Still, at a **fixed 300** budget lung needs
its natural (~550) budget to beat current ⇒ **budget must scale with #fine types**.

## STRATEGY A–E RESULTS
- A (current) and B (fine_global): implemented + benchmarked (above).
- **C (union), D (hierarchical multi-panel inference), E (fine_global + independent sibling/rare
  confidence): NOT implemented this session.** The central question — do sibling/rare add
  predictive value over fine_global alone — is therefore **not yet answered**.

## PAIRED STATISTICAL RESULTS (B_fine_global − A_current; negative ⇒ B better)
- breast RMSE −0.0025, CI [−0.0048, −0.0005], 73% mixtures B better → **B wins**.
- breast cond_RMSE −0.027, CI [−0.043, −0.012], 77% → **B wins**.
- lung RMSE +0.0027, CI [−0.0005, +0.0054], 23% → **B worse** at 300.
- lung cond_RMSE −0.0002, CI [−0.023, +0.021], 53% → **tie**.

## PER-CELL-TYPE FAILURE MODES
`metrics_by_celltype.tsv` saved. Lung's fixed-budget deficit concentrates in fine-grained
epithelial families (many types starved at 300 genes) — the type-count/budget interaction.

## RARE-CELL EXPLORATORY RESULTS
rare_fpr/recall reported at threshold 0.01 (exploratory, not calibrated). fine_global improves
rare recall on breast (0.48 vs 0.43) at slightly lower FPR; on lung recall dips (0.35 vs 0.38)
at ~equal FPR — consistent with lung's budget starvation. Calibrated rare detection is Stage 2.

## RUNTIME AND MEMORY
Benchmark completed both tissues × 5 seeds in a few minutes (optimizer + Poisson). Optimizer
per reference: ~1–4 min (donor DE + minimal LODO). No memory issues.

## NEGATIVE RESULTS
- At a FIXED budget, fine_global is tissue-dependent (wins breast, not lung) — a fixed cap is the
  wrong policy.
- Multipanel *inference* (sibling/rare adding value) remains **undemonstrated** — not a positive.

## KNOWN LIMITATIONS
Pseudobulk only (no fine ground truth); exploratory rare thresholds; Strategy C/D/E unimplemented;
spatial unvalidated; real-bulk untested for these panels.

## PROMOTION DECISION
- **Corrections (state mutation, stratified budget, combined rare, honest labels): merged**
  (experimental module; no default changed; tests added; suite green).
- **fine_global: experimental_candidate_pending_budget_and_rare_validation at a PER-TYPE budget** (natural ~top_n×n_types), where
  prior + this evidence show a two-tissue win; **NOT at a fixed budget**.
- **Multipanel architecture: NOT declared validated** — sibling/rare have not demonstrated added
  value. Per §9 this awaits Strategy D or E.

## NEXT STAGE
Implement Strategy E first (sibling + rare_confirmation as independent **reliability/confidence
scores** layered on fine_global, without changing proportions) and test whether they reduce
rare-FPR-when-absent / flag confounded subtypes without hurting recall or breast/lung accuracy.
If they add no measurable value, record the negative and simplify. Only then Strategy D, then
Stage 2 (calibrated rare detection).

---

## UPDATE — Strategy E implemented and evaluated (the central question)

`benchmarks/signatures/run_evidence_value.py` (+ `src/tissueresolve/reference/signature_evidence.py`,
per-subtype sibling markers in the optimizer). Leakage-safe: sibling/rare evidence scores + the
fine_global abundance are collected per (mixture, subtype-with-a-panel); mixtures split into dev
(seeds 0–2) and test (3–4); logistic (FP detection among predicted-present) and linear (|error|)
models fit on DEV ONLY, evaluated on TEST. Pre-specified presence threshold 0.01.

| tissue | FP-AUC abundance | FP-AUC abund+evidence | ΔAUC | err-R² abund | err-R² +evid | ΔR² | corr(rare,|err|) | corr(sib,|err|) |
|---|---|---|---|---|---|---|---|---|
| breast | 0.766 | 0.707 | **−0.060** | 0.435 | 0.440 | +0.005 | 0.268 | 0.086 |
| lung | 0.494 | 0.602 | **+0.108** | 0.552 | 0.550 | −0.002 | 0.159 | 0.228 |

**Finding:** sibling/rare evidence adds **no reproducible cross-tissue value** over fine_global
abundance. It *helps* FP detection on lung (where abundance alone is ≈random, AUC 0.49→0.60) but
*hurts* on breast (abundance already strong, 0.77→0.71), and adds ~0 to error-R² either way. The
raw evidence↔error correlations are modest and tissue-swapped. Signals are opposite-signed across
the only two tissues and rest on a single dev/test split with small test FP-subsets.

## SCIENTIFIC CONCLUSION (§16): **E — evidence remains inconclusive** (leaning negative)
No reproducible benefit of the multipanel evidence was demonstrated. There is a real but
**tissue-dependent** FP-detection gain (lung only, where abundance is unreliable), not a general
one. This does not meet "multipanel adds reproducible value" (A); nor is fine_global proven
strictly sufficient (D), given lung's genuine gain. Strategy C (union) remains the prior negative;
**Strategy D (hierarchical multi-panel inference) was not implemented this session** — so the
architecture cannot be declared validated.

## ARCHITECTURE SIMPLIFICATION DECISION
Keep `sibling`/`rare_confirmation` **experimental, opt-in, NOT wired into the default pipeline**.
Do not add their complexity to inference on current evidence. `fine_global` remains the one panel
with (budget-scaled) deconvolution value — still `experimental_candidate_pending_budget_and_rare_validation`.
Before keeping the extra panels long-term, a larger evaluation (≥2 more tissues, cross-platform,
per-abundance-band, paired ΔAUC CIs, plus Strategy D) must show consistent value; otherwise drop them.

---

## FINALIZATION UPDATE — fallback, safeguards, budget curves, Strategies C/D

**Fallback chain (§3, implemented+tested):** donor→batch→pooled-cell→inadequate with provenance
(`requested/actual_method`, donor_aware, batch_aware, confidence_level, limitations) and states
PASS/PASS_WITH_RESTRICTIONS/EXPERIMENTAL_FINE/BROAD_ONLY/REFERENCE_INADEQUATE/SINGLE_SUBTYPE.
Non-donor groupings never claim donor-validated PASS; pooled-cell emits a CRITICAL warning.

**Safeguards (§5, integrated):** `SignaturePanel` + `validate_panel_use` (strict raises
`PanelRoleError`; permissive warns+records) + `deconvolve_with_signature_panel` (records
`panel_use_warnings`); `SignatureModel.as_panel`. Now CALLED at a real inference entry, not just
defined.

**Gene-budget curves (§7, run — settles the ~550 hypothesis with data):** `fine_global` wins
conditional within-family RMSE at nearly every budget on both tissues. It is gene-efficient on
breast (~110–200 genes ≈ current at 500–1000). On lung (~53 fine types) it needs **~400–550 genes
to beat current** (crossover: fine_global@553 RMSE 0.0299 vs current@500 0.0305) and is starved
below ~380 — so the fixed-300 deficit was budget starvation, and **budget must scale with #fine
types (~10/type)**, confirmed (not a fixed 550). Outputs: `budget_curve*.tsv`, `pareto_front.tsv`.

**Strategies B/C/D (§8-10, run):**
| tissue | B_fine_global | C_union | D_hierarchical |
|---|---|---|---|
| breast RMSE / cond / coverage | 0.036 / 0.190 / 1.00 | 0.040 / 0.220 / 1.00 | 0.047 / 0.195 / 0.13 |
| lung RMSE / cond / coverage | 0.030 / 0.238 / 1.00 | 0.030 / 0.230 / 1.00 | 0.034 / 0.251 / 0.12 |
D reuses the existing hierarchical inference (broad→fine→gate→unresolved, mass-conserving,
verified sum=1.0). **B is best on both tissues; C worse; D worse and over-abstains (~87% unresolved).**

## REVISED SCIENTIFIC CONCLUSION (§16): **D — fine_global alone is sufficient**
Across Strategy C (union worse), Strategy D (hierarchical multi-panel worse + over-abstains), and
Strategy E (sibling/rare evidence: no reproducible cross-tissue value), **no multipanel variant
beats the single flat `fine_global` panel** for deconvolution. The only positive signal (lung FP
detection in E) is non-reproducible (opposite on breast) and pertains to reliability flagging, not
deconvolution.

## ARCHITECTURE SIMPLIFICATION DECISION
**Simplify.** For the deconvolution path, `fine_global` (donor-aware fine one-vs-rest, per-type
budget) is the one panel that matters. Keep `broad` for broad-only references; keep
`sibling`/`rare_confirmation` **experimental, NOT wired into the inference pipeline** (their value
is undemonstrated) — candidates for removal if a larger multi-tissue/cross-platform/spatial study
still shows nothing. Do not build Strategy D's hierarchical-multipanel gate into defaults
(over-abstains, underperforms flat).
