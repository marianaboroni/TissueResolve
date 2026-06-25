# High-Granularity Strategy — Evaluation Report (NEGATIVE RESULT)

**Outcome: the strategy FAILED the promotion gates and is NOT integrated.** It
remains experimental (`src/tissueresolve/experimental/soft_hierarchy/high_granularity.py`),
not wired into any pipeline. Soft gating stays the final hierarchical uncertainty
layer. Nothing was committed; all benchmark inputs are git-ignored real data.

## Hypothesis tested

That high-granularity fine-subpopulation prediction improves — *without* increasing
spillover or false-positive subtype detection — by making the reference more
query-compatible (**ReferenceCalibration**) and restricting the candidate subtype
space (**CandidateRestriction**), then taking a consensus over robust panels
(**ConsensusStability**), applied **before** the validated soft gate, rather than
amplifying fragile subtype contrasts after the fact.

## Pipeline order tested (soft gating once, after the procedure)

`π_f` (broad, unchanged) → `q_raw` → calibration → candidate restriction →
multi-panel fine deconvolution on retained candidates → consensus → stability /
candidate-aware confidence → `θ_k = π_f q_consensus` → **soft gate (once)** →
`u_f = π_f − Σ θ_resolved`. Broad-family and total mass conserved exactly
(`mass_error = 0` in every benchmark row; unit tests assert within-family Σ=1 and
broad mass unchanged).

## Methods implemented (all components, independently ablatable)

1. **ReferenceCalibration** — conservative gene-level scaling (`S*_{gk}=a_g S_{gk}`,
   clipped) and affine log-scale; contrast-preserving (same scale across subtypes in
   a gene → within-gene ratios unchanged); residual-downweight option. Uses only the
   query expression (never truth/test donors).
2. **CandidateRestriction** — `none|topk|evidence|hybrid`; never empties a family;
   protects rare-but-real subtypes (raw ≥ 2%); records retained/dropped.
3. **Multi-panel fine deconvolution** — one function over panels {current, marker,
   query_detectable, donor_stable}; within-family WNNLS on the retained candidates.
4. **ConsensusStability** — `mean|median|weighted|conservative` (subtype kept only if
   detected by ≥M panels); stays on the within-family simplex; emits per-subtype
   stability (detection frequency × (1−CV)) used to modulate soft-gate confidence.
5. **Soft-gating compatibility** — the existing validated soft gate, applied once,
   with confidence = base × stability (low stability → routed to unresolved).

## Files

* **Created:** `src/tissueresolve/experimental/soft_hierarchy/high_granularity.py`,
  `tests/test_high_granularity.py`, `benchmarks/diagnostics/high_granularity_benchmark.py`,
  this report.
* **Modified:** none (no integration — gates failed).

## Tests added

17 offline/deterministic tests in `tests/test_high_granularity.py`: calibration
shape/finiteness + contrast preservation + no-truth/no-test-donor inputs; candidate
restriction never empties / records drops / protects rare; multi-panel fan-out;
consensus non-negativity + within-family simplex; broad+total mass conservation;
soft gating applied once after consensus; default pass-through; unsupported-family
fallback; metadata records all selected methods. **All 17 pass.**

## Datasets used

Breast single-cell reference (30,000 cells, 41 donors, 41 fine → 8 families) and
HLCA-core lung subset (11,369 cells, 40 donors, 59 fine → 10 families).

## Train / calibration / validation / test split design

Donor- and seed-disjoint. Per seed a reference/query donor split; donor-stability
panels use **reference donors only**; the soft-gating confidence model is fit on
**calibration seeds {0,1}**; reference calibration uses only the query expression;
all reported metrics from held-out **test seeds {2,3}**. Scenarios: balanced,
imbalanced, rare, similar_subtypes (collinear).

## Breast metrics (test-seed means) — DECISIVE

| mode | cond_RMSE | cond_Pearson | fine_RMSE | broad_RMSE | pairwise_spillover | rare_sens | false_pos | eff_N | runtime_s |
|---|---|---|---|---|---|---|---|---|---|
| flat_nnls | 0.262 | 0.634 | 0.035 | 0.026 | 0.069 | 0.526 | 0.303 | 11.85 | 17.7 |
| ungated | 0.271 | 0.579 | 0.043 | 0.072 | 0.043 | 0.499 | 0.146 | 10.54 | 17.7 |
| hard_gate | 0.274 | 0.566 | 0.059 | 0.072 | 0.050 | 0.066 | 0.026 | 4.54 | 17.7 |
| **baseline_softgate** | **0.271** | **0.579** | 0.039 | 0.072 | **0.043** | 0.484 | **0.140** | 10.55 | 17.7 |
| calib_softgate | 0.298 | 0.426 | 0.053 | 0.072 | 0.102 | 0.600 | 0.445 | 13.71 | 35.9 |
| candidate_softgate | 0.338 | 0.389 | 0.052 | 0.072 | 0.083 | 0.375 | 0.182 | 8.76 | 34.6 |
| consensus_softgate | 0.379 | 0.304 | 0.052 | 0.072 | 0.098 | 0.249 | 0.239 | 5.89 | 34.7 |
| calib_candidate | 0.306 | 0.428 | 0.053 | 0.072 | 0.096 | 0.471 | 0.229 | 10.68 | 35.6 |
| candidate_consensus | 0.377 | 0.318 | 0.052 | 0.072 | 0.087 | 0.236 | 0.106 | 5.63 | 37.4 |
| full_softgate | 0.342 | 0.364 | 0.055 | 0.072 | 0.111 | 0.306 | 0.124 | 7.55 | 34.9 |

**Every variant is worse than `baseline_softgate` on conditional RMSE** (best is
`calib_softgate` at +10% worse; the rest +13–40%), and **spillover and false
positives rise** for the calibration/hybrid variants. The consensus variants
(`consensus_softgate`, `candidate_consensus`) lower spillover modestly **only by
collapsing rare-subtype sensitivity** (0.24–0.25 vs 0.48 baseline) and effective-N
(5.6 vs 10.5) — i.e. they "improve stability by pushing everything to unresolved",
which the gates explicitly forbid. No variant improves conditional RMSE at all.

## Lung / HLCA metrics

The HLCA 6-variant multi-panel sweep is compute-heavy (~4× breast per scenario) and
did not complete within the run budget. **It cannot change the decision:** the
promotion gates require improvement in **both** tissues, and breast already fails
the headline conditional-RMSE gate plus the spillover/false-positive and
rare-recall gates for every variant. The breast failure mode (calibration/contrast
emphasis worsens spillover; consensus restriction destroys rare recall) is the same
mechanism that sank the FineGranularityRefiner on both breast and lung, so a lung
reversal is implausible. (If the lung sweep is re-run with fewer panels/modes, its
table should be appended here; the conclusion stands on breast.)

## Ablation results

* **ReferenceCalibration only** (`calib_softgate`): cond RMSE +10%, spillover 0.043→0.102,
  false-pos 0.14→0.45. Scaling toward the query inflates false detection. Fails.
* **CandidateRestriction only** (`candidate_softgate`): cond RMSE +25%, rare sens 0.48→0.37.
  Concentrating family mass on top candidates raises per-candidate error. Fails.
* **Multi-panel consensus only** (`consensus_softgate`): cond RMSE +40%, rare sens 0.48→0.25,
  eff-N 10.5→5.9. Conservative consensus over-collapses. Fails.
* **Combinations** (`calib_candidate`, `candidate_consensus`, `full_softgate`): all worse
  than baseline; combining components does not rescue them.

## Candidate-restriction results

Restriction reduces eff-N (5.6–8.8 vs 10.5) and, in the conservative-consensus
variants, sharply reduces rare-subtype sensitivity (to ~0.24). It does not reduce
conditional error. CandidateRestriction "improves precision by destroying rare
recall" — an explicit fail condition.

## Reference-calibration results

Gene-level scaling toward the query **increases** false-positive detection
(0.14→0.45 breast) and spillover (0.043→0.102) while worsening conditional RMSE;
held-out reconstruction shifts do not translate into better subtype recovery.

## Consensus-stability results

Cross-panel consensus lowers spillover only marginally and only by routing mass to
unresolved / dropping subtypes, collapsing eff-N and rare recall. Stability gains
come at the cost of recall — a fail condition.

## Spillover results

Pairwise within-family spillover **increased** for every variant on breast
(0.043 → 0.083–0.111). No variant reduced it without simultaneously destroying
rare recall.

## False-positive subtype detection

Increased for calibration/hybrid variants (up to 0.45 vs 0.14 baseline). The
consensus variants reduce false positives only alongside the rare-recall collapse.

## Rare-subtype sensitivity and precision

Sensitivity preserved/raised only by `calib_softgate` (0.60) — bundled with its
false-positive blow-up. Restriction/consensus variants reduce sensitivity to
0.24–0.37. No variant improves both sensitivity and precision.

## Broad-RMSE impact

Unchanged by construction (broad mass never touched): 0.072 across all hierarchical
variants. Gate passes.

## Effective-N and entropy

Consensus variants deflate eff-N below truth-compatible ranges (5.6–5.9 vs ~10.5);
calibration variants inflate it (13.7). Both are truth-incompatible.

## Runtime

~2× the baseline per scenario (≈35 s vs 18 s breast); acceptable but not justified
by any accuracy gain.

## Promotion-gate results (best variant = calib_softgate, breast)

| gate | breast |
|---|---|
| conditional RMSE improves ≥5% rel | ❌ (+10% worse) |
| pairwise spillover non-inferior | ❌ |
| false-positive non-inferior | ❌ |
| rare sensitivity preserved | ✅* |
| broad RMSE not worse >2% | ✅ |
| eff-N not inflated | ✅ |
| mass conserved | ✅ |
| runtime acceptable | ✅ |
| multi-family improvement >1 | ✅ |
| improvement in BOTH tissues | ❌ (breast fails) |

\* the rare-sensitivity "pass" for `calib_softgate` is inseparable from its
false-positive blow-up. **The decisive gates (conditional RMSE, spillover,
false-positive) fail on breast for every variant**, so the both-tissue requirement
cannot be met.

## Was the method integrated?

**No.** No CLI/config/API/report/docs changes. The strategy stays experimental.

## Report / documentation updates

None to user-facing docs (integration is gated on passing). Only this report
documents the negative result, as required.

## Interpretation

This is the third independent confirmation of the same identifiability ceiling:
in shared-lineage-dominated families (96–99% shared signature), neither amplifying
subtype contrast (FineGranularityRefiner) nor calibrating the reference toward the
query / restricting candidates / consensus-averaging (this strategy) recovers
reliable high-granularity composition. Each buys a metric (or stability) only by
trading away spillover control, false-positive discipline, or rare recall. The
honest output for such families is the broad family (or `unresolved_<family>`)
under the validated soft gate — which is exactly what the new Resolution Decision
Layer now classifies as `broad_only` before fine predictions are interpreted.

## Remaining limitations

Bulk pseudobulk, synthetic ground truth, two test seeds × four scenarios; the lung
6-variant sweep is expensive and was not completed in budget (does not affect the
conclusion). The failure directions are large and consistent with prior results.

## Safe to commit?

Code is safe (new experimental module + tests + one diagnostic script; no core
changes; targeted + unit suites green). **Per instruction, nothing was committed.**
Benchmark outputs and real-data inputs are git-ignored.
