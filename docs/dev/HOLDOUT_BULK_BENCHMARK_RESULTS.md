# Held-Out-Donor Bulk Benchmark — First Ground-Truth Results (Phase 2, TissueResolve-only)

**What this is.** The first *accuracy* benchmark with known ground truth, run with
**no external tools**. Donors are split into disjoint reference (20) and query
(21) sets; the signature is built from reference donors only; count-level
pseudobulk is generated from the held-out query donors; predictions are scored
against the **mRNA-proportion truth** (what an RNA-based deconvolver estimates).
Runner: `benchmarks/bulk/run_holdout_bulk_benchmark.py`. Outputs:
`benchmarks/outputs/holdout_bulk/` (metrics TSVs, raw predictions, manifest, figures).

**Design.** 5 scenarios × 12 samples, 30 reference cell types, 600 cells/sample,
seed 0. Methods: flat `nnls` / `weighted_nnls` / `ridge_nnls` / `auto`, and
`hierarchical`. One donor split, one tissue — so these are **point estimates,
not yet CI-backed** (bootstrap + multi-split is the Phase-2 statistical layer).

---

## Headline numbers (mean over 5 scenarios)

| Method | fine Pearson | fine CCC | fine JSD↓ | broad Pearson | rare detection (true 0.7%) | complexity error↓ |
|---|---|---|---|---|---|---|
| `flat_nnls`          | **0.80** | **0.78** | 0.170 | 0.96 | 92 % | 3.3 (under-detects rich) |
| `flat_auto`          | 0.75 | 0.69 | 0.184 | **0.96** | 92 % | scenario-dependent |
| `flat_weighted_nnls` | 0.73 | 0.71 | 0.188 | 0.92 | 33 % | low |
| `flat_ridge_nnls`    | 0.66 | 0.53 | **0.161** | 0.88 | **100 %** | over-disperses concentrated |
| `hierarchical`       | 0.28 | 0.22 | 0.728 | 0.59–0.87 | **0 %** | **collapses to ~2.7** |

## What this says about the "4–5 populations" concern

1. **The true composition is diffuse (effective-N 10–25), not 4–5.** So *some*
   detected sparsity is genuine under-detection, not biology.
2. **No solver is universally best — the trade-off is real and measured:**
   - `nnls` (the sparse one) has the **best mean fine accuracy and rare-pop
     recall**, because it puts mass on the right dominant types — but it
     **under-detects** when the truth is genuinely rich (balanced: pred eff-N
     16.7 vs true 24.7).
   - `ridge` **best preserves complexity** on rich mixtures and has 100 % rare
     recall, but **over-disperses** concentrated mixtures (imbalanced: pred 24.8
     vs true 10.5) and has the **worst fine CCC** (spreads mass onto wrong types).
   - So sparsity-vs-diffuseness should be judged against the *actual* truth
     complexity; `auto`'s data-dependent choice is defensible.
3. **The clearest defect is `hierarchical` mode.** It collapses every scenario to
   ~2.7 effective populations (vs true 10–25), scores worst on fine **and broad**
   accuracy, and **detects the rare population in 0 % of samples** (it abstains it
   into unresolved T/NK mass). This corroborates the sparsity-audit finding that
   the within-family gate (separability ≥ 0.10, spillover ≤ 0.30) is **too
   aggressive** on near-collinear signatures — here it measurably *hurts accuracy*.

## Per-scenario complexity preservation (pred vs true effective-N)

See `bulk_complexity_metrics.tsv` / `figures/figA_richness_preservation`. Pattern:
ridge tracks the diffuse truth; nnls/auto under-detect rich mixtures by ~30 %;
hierarchical is flat at ~2–3 regardless of the truth.

## Honest caveats

- **One donor split, one tissue (breast), 12 samples/scenario.** No CIs yet;
  Phase 2 adds bootstrap + paired tests + ≥1 more tissue.
- **Hierarchical fine metrics strip `unresolved_*` columns**, which penalises
  honest abstention. But its **broad-level** accuracy (fair) is *also* lower and
  rare recall is 0 %, so the underperformance is real, not purely an artifact.
  A follow-up should also report `hierarchical_fair_metrics` (resolvable-only).
- **Scenario `balanced` (eff-N 24.7) is an extreme** near-uniform stress test;
  real tumours are closer to `imbalanced` (eff-N 10.5), where `nnls` wins and
  ridge over-disperses — i.e. for realistic compositions the *sparse* solver may
  be the more accurate one.
- All flat solvers **over-estimate the rare population** (bias +0.015–0.025): the
  detection floor is optimistic, the abundance is inflated.
- Truth = mRNA proportions (not cell fractions); cell-fraction truth also saved.

## Tentative implications (NOT yet acted on — CLAUDE.md rule 10)

- `auto` defaulting to `nnls`-like behaviour is reasonable for fine accuracy.
- **The hierarchical within-family gating thresholds are the strongest candidate
  for revision**, but only after: (a) external-method comparison (do MuSiC/
  BayesPrism also collapse?), (b) bootstrap CIs, (c) a fair resolvable-only
  re-score. This is the top Phase-2 question.

---

## Addendum — multi-split + bootstrap CIs + fair re-score (PART 14)

To harden the single-split point estimates, the benchmark was repeated across
**5 donor splits × 5 scenarios = 25 paired replicates** with bootstrap 95% CIs,
paired Wilcoxon tests, and a **fair hierarchical re-score** (fine accuracy on
only the families it chose to resolve). Runner:
`benchmarks/bulk/run_holdout_bulk_multisplit.py`; outputs:
`benchmarks/outputs/holdout_bulk_multisplit/` (`metrics_ci.tsv`,
`pairwise_wilcoxon.tsv`, `rank_stability.tsv`, figures `figD`/`figE`).

**Fine-level Pearson (mean [95% CI], n=25):**

| Method | fine Pearson | broad Pearson | complexity error↓ |
|---|---|---|---|
| `hierarchical_resolvable`* | 0.894 [0.857, 0.924] | 0.901 [0.864, 0.930] | 0.49 [0.34, 0.66] |
| `flat_nnls` | 0.747 [0.702, 0.791] | **0.959 [0.946, 0.970]** | **3.25 [2.20, 4.31]** |
| `flat_weighted_nnls` | 0.741 [0.692, 0.788] | 0.946 | 3.39 |
| `flat_auto` | 0.736 [0.692, 0.781] | 0.949 | 3.66 |
| `flat_ridge_nnls` | 0.678 [0.644, 0.714] | 0.887 | 7.86 [6.14, 9.60] |
| `hierarchical` (strict)* | 0.324 [0.251, 0.393] | (artifact, see below) | 12.54 |

*`hierarchical_resolvable` is scored only on the **resolvable subset** of
families (an easier panel) — **not directly comparable** to the flat full-panel
solvers. It quantifies *precision when the method commits*, not full-panel accuracy.

**What the statistics establish:**

1. **The earlier "hierarchical is worst" was largely an abstention-penalty
   artifact.** Strict 0.324 → fair 0.894: when judged only on the families it
   resolves, hierarchical mode is highly accurate. It is **high-precision,
   low-coverage** — right when it commits, but it abstains on most families
   (effective-N collapses to ~2.7). The trade-off, not a defect, but it does
   under-report diversity.
2. **Among the full-panel flat solvers, `nnls`/`auto`/`weighted` lead and `ridge`
   trails** on fine accuracy (nnls vs ridge: Wilcoxon p≈0.06, effect 0.43),
   broad accuracy, **and** complexity preservation. Crucially, averaged over
   realistic (non-uniform) scenarios, **`nnls` preserves complexity *better* than
   `ridge`** (error 3.25 vs 7.86) — ridge only wins on the extreme near-uniform
   `balanced` case. So the single-split "ridge best complexity" was driven by one
   extreme scenario; for realistic compositions the sparse solver is both more
   accurate and more complexity-faithful.
3. **Rank stability (fine Pearson):** among flat solvers `nnls` has the best mean
   rank (2.56) and is never beaten by `ridge` (4.92) across bootstrap resamples.

**Two honesty notes on the multi-split table:**

- `hierarchical` (strict) **broad** Pearson reads as ≈0 here because this runner
  strips `unresolved_*` columns *before* family aggregation, zeroing abstained
  families. The *correct* fair broad for hierarchical (unresolved mass mapped
  back to families) is 0.59–0.87 (see the single-split section). Treat the strict
  broad row as a known artifact, not a result.
- A fully apples-to-apples fair comparison would score the flat solvers on the
  *same* resolvable subset as `hierarchical_resolvable`. That refinement is the
  next iteration; today's takeaway is the strict-vs-fair gap and the flat-solver
  ranking.

**Updated implication.** The accuracy-optimal *flat* default is `nnls`/`auto`
(not ridge). Hierarchical mode is not inaccurate — it is conservative; the open
question is whether its **coverage** (how many families it resolves) is too low,
which is a threshold-tuning question for Phase 2, to be settled against external
methods + a same-subset fair comparison before any change (CLAUDE.md rule 10).

---

## Addendum 2 — external methods (MuSiC, BisqueRNA) on the SAME data (PART 3A / PART 17)

External bulk methods were executed in an **isolated R 4.1.2 library**
(`benchmarks/envs/r_lib`) on the **identical** held-out reference + per-scenario
bulk as TissueResolve (fair-input contract). Honest execution status
(`benchmarks/outputs/holdout_bulk/method_status.tsv`):

| Method | Status | Version | Note |
|---|---|---|---|
| MuSiC | **executed** | 1.0.0 | binary deps; DESCRIPTION TOAST pin (≥1.10.1) relaxed — core `music_prop` unaffected |
| BisqueRNA | **executed** | 1.0.5 | CRAN binary |
| BayesPrism | **failed** | — | C++ source compilation broken in this env (`fatal error: 'vector' file not found`); cannot build. Needs a working toolchain/conda env. |
| CIBERSORTx | not attempted | — | license-gated; export-only, never executed here |
| DWLS | not attempted | — | deferred |

> Environment note: this sandbox's C/C++ toolchain cannot compile R source
> packages (missing C++ SDK headers). MuSiC/Bisque succeeded only via **binary**
> installs; BayesPrism has no R-4.1 binary and failed to compile, exactly as in
> the prior `tool_installation_status.tsv`. Reported as failed, not faked.

**Combined fine-level accuracy (mean Pearson over 5 scenarios, same inputs):**

| Method | family | fine Pearson | fine CCC | fine JSD↓ |
|---|---|---|---|---|
| `flat_nnls` | TissueResolve | **0.796** | **0.777** | 0.170 |
| `flat_auto` | TissueResolve | 0.752 | 0.691 | 0.184 |
| `flat_weighted_nnls` | TissueResolve | 0.733 | 0.711 | 0.188 |
| **MuSiC** | external | 0.669 | 0.653 | **0.127** |
| `flat_ridge_nnls` | TissueResolve | 0.661 | 0.532 | 0.161 |
| **BisqueRNA** | external | 0.389 | 0.350 | 0.375 |
| `hierarchical` (strict) | TissueResolve | 0.282 | 0.224 | 0.728 |

**Two external-evidence conclusions:**

1. **TissueResolve is competitive — its flat solvers (nnls/auto/weighted) lead
   MuSiC and BisqueRNA on fine Pearson + CCC** on this held-out-donor breast
   benchmark. MuSiC has the best *distributional* shape (lowest JSD 0.127) but a
   lower correlation; Bisque is weak here. (Single dataset, one split — not a
   general claim; see caveats.)
2. **The "sparse composition" concern is NOT unique to TissueResolve — PART 17
   answered.** Effective-N vs the true 16.2: **MuSiC 17.3** (tracks truth),
   **BisqueRNA 9.4** (*more* concentrated than TissueResolve's nnls), TissueResolve
   nnls ~13–16, ridge over-disperses. So composition concentration is a
   **method/conditioning property shared across deconvolvers** — Bisque collapses
   more than TissueResolve — not a TissueResolve-specific defect.

Figure: `benchmarks/outputs/holdout_bulk/figures/figF_external_vs_tissueresolve`.
Tables: `external_fine_metrics.tsv`, `combined_fine_metrics.tsv`,
`combined_method_summary.tsv`. Runtimes (per scenario): MuSiC ~25–46s, Bisque
~13–36s (in `external_inputs/external_method_status.tsv`).

**Caveats (unchanged + new):** single tissue, one donor split, 12 samples/scenario
— no CIs across methods yet (the multi-split CI machinery exists and should next be
applied to MuSiC/Bisque too); MuSiC's relaxed TOAST pin affects only `music2` (not
used); external methods got 5000 genes (their own internal selection applies).

---

## Addendum 3 — external methods across 5 splits, with CIs + paired Wilcoxon (PART 14 + 3A)

MuSiC and BisqueRNA were executed on **all 25 (split × scenario) replicates** —
the identical held-out inputs as the TissueResolve multi-split run — enabling
bootstrap CIs and paired Wilcoxon. All 50 external runs **executed** (0 failed).
Scripts: `export_external_multisplit.py`, `run_external_multisplit.R`,
`score_external_multisplit.py`. Outputs: `benchmarks/outputs/holdout_bulk/multisplit_combined/`
(`metrics_long/metrics_ci.tsv`, `pairwise_wilcoxon_vs_nnls.tsv`, figG).

**Fine Pearson, mean [95% CI], n=25:**

| Method | fine Pearson | broad Pearson | complexity error↓ |
|---|---|---|---|
| `hierarchical_resolvable`* | 0.894 [0.857, 0.924] | — | — |
| `flat_nnls` | **0.747 [0.702, 0.791]** | **0.959 [0.946, 0.970]** | 3.25 [2.20, 4.31] |
| `flat_weighted_nnls` | 0.741 | — | — |
| `flat_auto` | 0.736 | — | — |
| `flat_ridge_nnls` | 0.678 | — | 7.86 [6.14, 9.60] |
| **MuSiC** | 0.614 [0.511, 0.709] | 0.818 [0.782, 0.857] | **1.94 [1.51, 2.43]** |
| **BisqueRNA** | 0.377 [0.313, 0.440] | 0.602 [0.558, 0.642] | 5.98 [4.18, 7.82] |
| `hierarchical` (strict) | 0.324 | — | — |

*resolved-subset score, not full-panel comparable.

**Paired Wilcoxon vs `flat_nnls` (fine Pearson, n=25 paired replicates):**

| Comparison | p | effect | verdict |
|---|---|---|---|
| nnls vs **MuSiC** | **0.0002** | +0.78 | nnls significantly better |
| nnls vs **BisqueRNA** | **<0.0001** | +1.00 | nnls better on *every* replicate |
| nnls vs flat_auto | 0.75 | +0.10 | tie |
| nnls vs flat_weighted | 1.00 | 0.00 | tie |
| nnls vs flat_ridge | 0.059 | +0.43 | nnls better (borderline) |

**Statistically-backed conclusions:**

1. **On this held-out-donor breast benchmark, TissueResolve (nnls/auto/weighted)
   significantly outperforms MuSiC and BisqueRNA on fine- AND broad-level
   accuracy** — nnls beats MuSiC (p=2e-4) and Bisque (p<1e-4, every replicate),
   with broad Pearson 0.96 vs MuSiC 0.82 vs Bisque 0.60. MuSiC's fine CI is wide
   ([0.51, 0.71]) — it is also less *stable* across donor splits.
2. **But MuSiC best preserves composition complexity** (effective-N error 1.94 vs
   nnls 3.25, ridge 7.86): MuSiC recovers the *right number* of populations even
   when its per-type correlation is lower. A fair, non-trivial nuance — "diversity
   recovery" and "per-type accuracy" are different axes, and different methods win
   each. BisqueRNA is weakest on both here.
3. **PART 17 reconfirmed with CIs:** composition concentration is method-specific
   and shared across deconvolvers (Bisque off by ~6 effective populations; ridge
   by ~8); it is not a TissueResolve-unique artifact.

**Caveats:** still a single tissue (breast), one atlas, 5 splits; external methods
received the same 5000-gene universe (their internal gene selection applies);
BayesPrism/CIBERSORTx remain unexecuted (toolchain/license). A second tissue is
the next generalization step.

---

*Source: `benchmarks/outputs/holdout_bulk/` (single split),
`benchmarks/outputs/holdout_bulk_multisplit/` (TissueResolve 5 splits + stats),
and `benchmarks/outputs/holdout_bulk/multisplit_combined/` (external 5 splits +
CIs/Wilcoxon). Generators/scorers:
`benchmarks/bulk/run_holdout_bulk_benchmark.py`,
`benchmarks/bulk/run_holdout_bulk_multisplit.py`; stats helpers:
`benchmarks/shared/stats.py`; figures: `benchmarks/bulk/plot_holdout_results.py`.
No core algorithm modified; nothing committed.*
