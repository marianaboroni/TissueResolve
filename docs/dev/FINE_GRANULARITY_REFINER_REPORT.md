# FineGranularityRefiner — Evaluation Report (NEGATIVE RESULT)

**Outcome: the refiner FAILED the promotion gates on both tissues and is NOT
integrated.** It remains experimental (`src/tissueresolve/experimental/`), not
wired into any pipeline. Soft gating stays the final hierarchical uncertainty
layer. Nothing was committed; all benchmark inputs are git-ignored real data.

## Hypothesis tested

That a compact fine-level refinement layer — emphasising donor-stable
subtype-contrast genes (and optionally down-weighting shared-lineage genes and
calibrating within-family spillover) **before** the validated partial
confidence-weighted soft gate — improves *conditional within-family* subtype
prediction in families where subtypes are biologically distinguishable but highly
collinear in full-panel deconvolution.

## Pipeline order tested (soft gating applied exactly once, after the refiner)

`broad mass π_f` → `raw conditional q^raw_{k|f}` → **refiner** `q^refined_{k|f}` →
`θ_k = π_f q^refined_{k|f}` → confidence `c_k` → soft gate
`θ^resolved_k = θ_k c_k` → `u_f = π_f − Σ θ^resolved_k`. The refiner changes only
the within-family conditional; broad-family and total mass are conserved exactly
(verified: `mass_error = 0` in every benchmark row; unit tests assert
within-family Σ = 1 and per-family broad mass unchanged).

## Methods implemented (all four components, independently ablatable)

1. **Contrast-weighted WNNLS** — within-family re-solve with genes weighted by
   `w_{g,f} = (C·D·Q·P)/(1+S+V+R+T)` (subtype contrast, donor stability, query
   detectability, pairwise support / shared-lineage dominance, donor instability,
   redundancy, technical risk). No L1 sparsity; weights floored so no gene is
   hard-zeroed.
2. **Residual subtype signatures** — `S = B + D` decomposition with centered
   contrast (`Σ_k ω_k D_{k}=0`, exact); refine the conditional along the
   contrast-weighted residual directions. `D` used only for weighting/ranking;
   outputs remain non-negative simplex conditionals.
3. **Family-specific spillover calibration** — ridge map `q_corr = project(C_f
   q_pred)` learnt on donor-held-out calibration pseudobulk, projected back to the
   simplex; small/low-support families skipped.
4. **Conservative query-adaptive gene reweighting** — one reconstruction-residual
   down-weight + re-solve (implemented; off by default; not needed once 1–3 failed).

## Files

* **Created:** `src/tissueresolve/experimental/soft_hierarchy/fine_refiner.py`,
  `tests/test_fine_refiner.py`, `benchmarks/diagnostics/fine_refiner_benchmark.py`,
  this report.
* **Modified:** none (no integration — gates failed). The shipped breast hierarchy
  config is left untouched; 9 reference cell types absent from it are mapped into
  its existing broad families inside the benchmark only.

## Tests added

16 offline/deterministic tests in `tests/test_fine_refiner.py`: simplex projection;
residual-contrast centering; finite/non-negative contrast weights; missing-donor
handling; refined conditional non-negativity + within-family Σ=1 (both modes);
broad-family mass unchanged; refiner output is conditional (gate not applied
inside); canonical order conserves mass with soft gating applied once; default
mode is pass-through; insufficient-support fallback; calibration train/test
separation + simplex projection + skip-small-family + identity-when-none; metadata
records refinement status. **All 16 pass.**

## Datasets used

* **Breast** — real single-cell breast-cancer reference (CELLxGENE),
  30,000 cells × 5,000 genes, 41 donors, 41 fine types → 8 broad families.
* **Lung** — HLCA-core subset, 11,369 cells, 40 donors, 59 fine types → 10 broad
  families (the same open+labeled second tissue used for soft-gating validation).

## Train / calibration / test split design

Donor- **and** seed-disjoint. Per seed, `split_donors` partitions donors into
reference/query (50/50). Contrast weights and donor-stability residuals use
**reference donors only**; the confidence model and the spillover calibrator are
fit on **calibration seeds {0,1} only**; all reported metrics come from **held-out
test seeds {2,3}**. Gene weighting never sees test donors (rules 12–13).
Scenarios: balanced, imbalanced, rare, similar_subtypes (collinear).

## Metrics — breast (test-seed means)

| mode | cond_RMSE | cond_Pearson | fine_RMSE | broad_RMSE | pairwise_spillover | rare_sens | false_pos | eff_N | mass_err |
|---|---|---|---|---|---|---|---|---|---|
| flat_nnls | 0.262 | 0.634 | 0.035 | 0.026 | 0.069 | 0.526 | 0.303 | 11.85 | 0 |
| ungated | 0.271 | 0.579 | 0.043 | 0.072 | 0.043 | 0.499 | 0.146 | 10.54 | 0 |
| hard_gate | 0.274 | 0.566 | 0.059 | 0.072 | 0.050 | 0.066 | 0.026 | 4.54 | 0 |
| **baseline_softgate** | **0.271** | **0.579** | 0.039 | 0.072 | **0.043** | 0.484 | **0.140** | 10.55 | 0 |
| contrast_softgate | 0.336 | 0.377 | 0.057 | 0.072 | 0.089 | 0.395 | 0.451 | 9.23 | 0 |
| residual_softgate | 0.298 | 0.495 | 0.047 | 0.072 | 0.060 | 0.453 | 0.216 | 10.69 | 0 |
| spillover_softgate | 0.267 | 0.495 | 0.048 | 0.072 | 0.104 | 0.788 | 0.700 | 18.75 | 0 |

## Metrics — lung / HLCA (test-seed means)

| mode | cond_RMSE | cond_Pearson | fine_RMSE | broad_RMSE | pairwise_spillover | rare_sens | false_pos | eff_N | mass_err |
|---|---|---|---|---|---|---|---|---|---|
| flat_nnls | 0.289 | 0.409 | 0.028 | 0.031 | 0.100 | 0.486 | 0.225 | 17.23 | 0 |
| ungated | 0.259 | 0.512 | 0.043 | 0.104 | 0.078 | 0.473 | 0.158 | 13.57 | 0 |
| hard_gate | 0.267 | 0.490 | 0.043 | 0.104 | 0.074 | 0.095 | 0.012 | 5.41 | 0 |
| **baseline_softgate** | **0.259** | **0.512** | 0.041 | 0.104 | **0.078** | 0.465 | **0.148** | 13.57 | 0 |
| contrast_softgate | 0.307 | 0.292 | 0.044 | 0.104 | 0.120 | 0.431 | 0.248 | 15.73 | 0 |
| residual_softgate | 0.253 | 0.513 | 0.040 | 0.104 | 0.092 | 0.452 | 0.192 | 14.24 | 0 |
| spillover_softgate | 0.201 | 0.567 | 0.040 | 0.104 | 0.135 | 0.759 | 0.374 | 31.92 | 0 |

## Ablation results

* **Contrast-weighted WNNLS (Component 1): clearly harmful on both tissues** —
  conditional RMSE +23.9% (breast) / +18.7% (lung), Pearson down, pairwise
  spillover and false positives up. Up-weighting the small (1–4% of signature),
  donor-fragile subtype-contrast genes amplifies noise rather than signal in the
  full-panel problem — exactly the failure mode the prior diagnosis predicted.
* **Residual-contrast (Component 2): negative-to-neutral** — breast +10.0%, lung
  −2.2% (below the 5% bar); pairwise spillover increases on both.
* **Spillover calibration (Component 3): a metric artifact, not real improvement**
  — it is the only mode that lowers conditional RMSE (−1.6% breast / −22.3% lung),
  but it does so by **spreading mass across subtypes**: false-positive detection
  jumps to 0.70 (breast) / 0.37 (lung) vs ~0.14 baseline, pairwise spillover rises,
  and effective-N inflates to 18.8 (breast) / 31.9 (lung) vs ~13.6 truth-level. The
  ridge calibrator regresses conditional toward the average over collinear
  subtypes, which lowers RMSE while manufacturing exactly the **false fine-level
  precision** the gates forbid (gate 12).
* **Component 4** not required (1–3 failed).

## Spillover results

Pairwise within-family spillover **increased** for every refiner mode on both
tissues (breast baseline 0.043 → 0.060–0.104; lung 0.078 → 0.092–0.135). No mode
reduced spillover. Gate 2 fails universally.

## Rare-subtype results

Contrast/residual modes reduced rare sensitivity (e.g. contrast 0.40 vs 0.48
breast). Spillover-calibration "raised" rare sensitivity (0.79 breast / 0.76 lung)
only by detecting mass everywhere — inseparable from its false-positive blow-up,
so not a genuine gain.

## Broad-RMSE impact

Unchanged by construction (the refiner never touches broad-family mass): broad RMSE
identical across all hierarchical modes (0.072 breast / 0.104 lung). Gate 5 passes.

## Effective-N and entropy

eff-N collapses are avoided by contrast/residual modes, but spillover calibration
**inflates** eff-N well above the truth level (mass-spreading), the opposite
failure — degenerate over-resolution.

## Runtime

Acceptable. Per scenario the refiner adds <0.5 s over the base hierarchical solve
(breast ~8.5 s, lung ~69 s per scenario, dominated by the existing deconvolution).

## Promotion-gate results (best refiner = spillover_softgate, both tissues)

| gate | breast | lung |
|---|---|---|
| conditional RMSE improves ≥5% rel | ❌ (−1.6%) | ✅ (−22.3%) |
| pairwise spillover decreases | ❌ | ❌ |
| rare sensitivity non-inferior | ✅* | ✅* |
| false-positive not materially worse | ❌ | ❌ |
| broad RMSE not worse >2% | ✅ | ✅ |
| eff-N not collapsed | ✅ | ✅ |
| mass conserved | ✅ | ✅ |
| runtime acceptable | ✅ | ✅ |
| improvement in >1 family | ✅ | ✅ |

\* rare-sensitivity "pass" is driven by indiscriminate mass-spreading, not real
recovery. **The two decisive gates — pairwise spillover and false-positive
detection — fail on BOTH tissues, and the headline conditional-RMSE gate fails on
breast.** Gate 9 (improvement in both tissues) is therefore not met.

## Was the method integrated?

**No.** No CLI/config/API/report/docs were changed. The refiner stays experimental.

## Report / documentation updates

None to the QC-first report or user docs (integration is gated on passing, which
did not happen). Only this report documents the negative result, as required.

## Interpretation

This corroborates and sharpens the existing identifiability finding. Fine subtypes
in these families are **cell-classifiable** but their deconvolution-relevant
contrast is a tiny (1–4%), donor-fragile fraction of a shared-lineage-dominated
signature (96–99%). Any layer that amplifies that contrast (contrast weighting) or
regresses toward the family average (spillover calibration) buys lower conditional
RMSE only by trading away spillover control and false-positive discipline — i.e. it
manufactures false fine-level precision. The honest result for these families is
the broad family (or `unresolved_<family>`) under the validated soft gate, not a
forced sub-resolution. **Soft gating should remain the final hierarchical layer.**

## Remaining limitations

* Bulk pseudobulk on two tissues, two test seeds, four scenarios; a synthetic-truth
  ground truth (no native paired spatial truth). The direction of the failures is
  large and consistent, so additional seeds are unlikely to reverse it.
* Donor-stability residuals require ≥2 reference donors with ≥10 cells per member;
  very small families fall back (recorded), so donor stability is not exercised in
  the smallest families.
* The contrast-weight formula uses simple robust factor definitions; a more elaborate
  parameterisation was deliberately avoided (the negative signal is unambiguous and
  not a tuning artifact — the failures are in spillover/false-positives, which more
  aggressive weighting would worsen).

## Safe to commit?

The code is safe (new experimental module + tests + one diagnostic script; no core
changes; all tests green). **Per instruction, nothing was committed.** Benchmark
outputs and the real-data inputs are git-ignored and must not be committed.
