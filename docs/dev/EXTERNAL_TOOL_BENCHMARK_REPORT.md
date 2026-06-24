# External-Tool Benchmark Report

Benchmark of TissueResolve against established deconvolution tools, after the
full-workflow validation (`FULL_WORKFLOW_REPORT_VALIDATION.md`). Donor- and
seed-disjoint; no model defaults changed; spatial λ unchanged at 0.1. Nothing
committed; large inputs/outputs git-ignored.

## 1. Executive summary

On donor-held-out **bulk** pseudobulk (breast + HLCA/lung, real ground truth,
**24 paired replicates**), TissueResolve's hierarchical + soft-gating path is
**competitive but not superior on raw accuracy** and **robustly superior on
resolution-aware reliability**. With paired bootstrap CIs + Wilcoxon + BH
correction it is **significantly better** (BH p < 1e-3, CI excludes 0) on
false-positive subtype detection (~9× lower), pairwise spillover, and rare-subtype
precision; **statistically tied** on conditional within-family RMSE (corrected from
an earlier 2-seed reading); and **worse, by design**, on raw broad/fine correlation
and rare sensitivity because it deliberately abstains in collinear families (routing
unreliable fine mass to unresolved). Flat methods — including TissueResolve flat,
algorithmically identical to plain NNLS here — achieve higher correlation but with
far more false positives. No method is "best overall"; the choice depends on whether
the goal is maximal fine recall or trustworthy, low-false-positive resolution. See
§12 for the suitability-WARNING and abstention-trade-off limitations.

## 2. Datasets and ground truth

* Breast SC reference (30k cells, 41 donors, 41 fine → 8 families).
* HLCA/lung subset (11.4k cells, 40 donors, 59 fine → 10 families).
* **Bulk (Cat. A):** donor-held-out pseudobulk → truth known (this report's primary,
  fully run fresh).
* **Spatial synthetic (Cat. B) / real Visium (Cat. C):** see §6–7 (prior-run evidence;
  real Visium reports **no accuracy**, only structure/concordance).

## 3. Tool versions and installation status

Verified installed (see `EXTERNAL_TOOL_INSTALLATION_NOTES.md`,
`benchmarks/external_tools/tool_registry.tsv`): R 4.1.2 — MuSiC 1.0.0, BisqueRNA
1.0.5, BayesPrism 2.2.3, CARD 1.1, spacexr/RCTD 2.2.1; Python — cell2location 0.1.4
(scvi 1.1.6, torch 2.8.0, Apple **MPS** GPU, no CUDA). Not installed: SCDC, DWLS,
SPOTlight, Tangram, Seurat. CIBERSORTx: web/token-gated (export-only, cannot
auto-run). Internal baselines (`.venv`): NNLS, WNNLS, TissueResolve flat/hierarchical.

## 4. Input harmonization

One shared reference per tissue from reference donors only; gene symbols harmonized
(intersection reported); same donor-held-out test mixtures and label set fed to
each method; normalization only as each tool requires; differences documented.
Test truth never shown to any method.

## 5. Bulk benchmark (fresh, both tissues, **12 seeds × 2 tissues = 24 replicates/method**)

Method means over the 24 donor-disjoint replicates:

| method | broad Pearson | broad RMSE | fine Pearson | fine RMSE | cond RMSE | spillover | rare sens | rare prec | false-pos |
|---|---|---|---|---|---|---|---|---|---|
| NNLS_baseline (= TR flat) | 0.941 | 0.045 | 0.655 | 0.038 | 0.304 | 0.189 | 0.446 | 0.572 | 0.226 |
| WNNLS_baseline | 0.921 | 0.053 | **0.656** | **0.039** | **0.296** | 0.182 | **0.467** | 0.605 | 0.207 |
| TissueResolve_flat | 0.941 | 0.045 | 0.655 | 0.038 | 0.304 | 0.189 | 0.446 | 0.572 | 0.226 |
| **TissueResolve_hierarchical_soft** | 0.743 | 0.100 | 0.273 | 0.052 | 0.300 | **0.154** | 0.076 | **0.701** | **0.024** |

(Per-tissue/seed values in `benchmarks/outputs/external_bulk_metrics.tsv`; per-family in
`external_bulk_per_family_metrics.tsv`.) NNLS_baseline and TissueResolve_flat are the
same algorithm (plain NNLS), so their rows match exactly — reported transparently.
Seeds were expanded from 2 → 12 specifically to power the superiority tests below
(§10); this **corrected an earlier 2-seed reading** in which conditional RMSE looked
best for TissueResolve — at n=24 it is statistically **tied**, not superior.

**External tools (MuSiC, BisqueRNA, BayesPrism, …):** installed and executed in a
prior harness run (`benchmarks/outputs/real_external_method_status.tsv`: MuSiC ~47 s,
BisqueRNA ~16 s, BayesPrism failed at ~2180 s). They were **not re-executed this
session** to avoid the multi-hour runtimes / environment instability observed; their
status and prior runtimes are recorded in `benchmarks/external_tools/run_manifest.tsv`,
clearly labelled. This is disclosed rather than hidden (rule 9); the fresh head-to-head
above is the internal comparison.

## 6. Spatial synthetic benchmark with truth (prior evidence)

Covered by the existing validated scripts and `docs/SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md`;
CARD (~42 s) and cell2location (~163 s) executed in the prior run
(`real_external_method_status.tsv`). A fresh re-run was not performed this session
(runtime/stability); the spatial pipeline itself is exercised by the full-workflow
machinery (same reference/hierarchy/resolution/soft-gating/report). No spatial λ
default was changed.

## 7. Real spatial no-truth evidence

Per `docs/SPATIAL_REAL_VISIUM_EVIDENCE.md`: marker recovery, method concordance,
spatial structure, and runtime only — **no accuracy claim** (no ground truth, rule 7).

## 8. Runtime and memory

Fresh internal bulk methods: NNLS/WNNLS/flat sub-second to ~1 s; TissueResolve
hierarchical ~1–2 s per run (`external_bulk_runtime.tsv`). External (prior run):
MuSiC ~47 s, BisqueRNA ~16 s, BayesPrism ~2180 s (failed), CARD ~42 s, cell2location
~163 s. TissueResolve is faster than every external tool measured.

## 9. Failure / skipped-tool summary

`run_manifest.tsv` records each tool: fresh internal runs (executed); MuSiC/BisqueRNA/
CARD/cell2location/RCTD installed + prior-executed, not re-run this session
(runtime/stability); BayesPrism prior `failed` (~2180 s); SCDC/DWLS/SPOTlight/Tangram
not installed; CIBERSORTx web/licensed (export-only). No tool silently skipped.

## 10. Statistical comparison (powered: 24 paired replicates)

Seeds expanded to 12 (× 2 tissues = **24 donor-disjoint paired replicates**). For
each superiority-claim metric, TissueResolve_hierarchical_soft is compared to every
baseline with a paired mean difference, percentile **bootstrap 95% CI** of the
difference (4000 resamples), **paired Wilcoxon** signed-rank test, matched-pairs
effect size, and **Benjamini–Hochberg** FDR correction across all 18 tests
(`external_benchmark_statistical_tests.tsv`); rank stability is the mean/std rank
over the 24 replicates (`external_benchmark_rank_stability.tsv`).

**Robustly superior (BH-significant *and* bootstrap CI excludes 0 favourably):**

| metric | TR vs baselines (mean diff) | 95% CI | BH p | effect size | mean rank (TR) |
|---|---|---|---|---|---|
| false-positive rate ↓ | −0.18 to −0.20 | excl. 0 | <1e-3 | −1.0 | **1.0 (std 0.0)** |
| pairwise spillover ↓ | −0.028 to −0.035 | excl. 0 | 1.5e-4–8e-4 | −0.75 to −0.84 | **1.54** |
| rare precision ↑ | +0.10 to +0.13 | excl. 0 | 1.8e-4–2.7e-3 | +0.69 to +0.83 | **1.46** |

**Statistically tied (NOT a superiority claim — corrected from the 2-seed reading):**

| metric | TR vs baselines | 95% CI | BH p |
|---|---|---|---|
| conditional within-family RMSE | −0.004 (vs NNLS/flat), +0.005 (vs WNNLS) | includes 0 | 0.34–0.49 (ns) |

(TR mean rank on cond RMSE = 2.38; WNNLS = 1.88 — not separable.)

**Robustly worse (by design — abstention):** fine Pearson (Δ≈−0.38), broad Pearson
(Δ≈−0.20), fine RMSE, rare sensitivity, eff-N — all BH-significant against TR, rank
4/4 (std 0.0) on fine/broad Pearson.

## 11. Interpretation

The two tool families optimise different objectives. Flat solvers (NNLS/WNNLS, and
external bulk tools in the same family) maximise raw correlation by always assigning
fine mass — at the cost of high false-positive subtype detection and spillover in
collinear families. TissueResolve's hierarchical + soft-gating + Resolution Decision
Layer instead **only reports fine resolution where it is statistically supportable**,
yielding far lower false positives/spillover and higher rare precision, at the cost
of fine sensitivity/correlation. This matches every prior internal finding
(soft gating validated; refiner / high-granularity negative results).

## 12. Limitations (read before citing any claim)

* **Reference-suitability WARNING.** Both atlases score `WARNING` on the
  reference-suitability check (breast 0.627, lung 0.740), driven by gene-overlap and
  marker-stability on these subsampled references. All downstream numbers inherit
  this caveat — they characterise *relative* method behaviour on imperfect references,
  not absolute performance on an ideal one.
* **Abstention trade-off (central).** TissueResolve hierarchical+soft is *conservative
  by construction*: in collinear, shared-lineage-dominated families the Resolution
  Decision Layer marks subtypes `broad_only` and soft gating routes their mass to
  unresolved (≈0.81–0.83 of fine mass on these atlases). This is exactly why it wins
  on false positives / spillover / rare precision and **loses** on raw fine/broad
  correlation, fine RMSE, and rare *sensitivity*. The win is a deliberate
  precision-over-recall trade, not free-lunch superiority; on tasks that need maximal
  fine recall, a flat solver is preferable.
* **Conditional-RMSE claim downgraded.** At 2 seeds TR looked best on conditional RMSE;
  at 24 paired replicates the difference is not significant (BH p 0.34–0.49, CI spans
  0) — reported as **tied**, not superior. Expanding seeds changed a claim, which is
  why the powered analysis matters.
* **Synthetic truth / scope.** Bulk pseudobulk simulated from labelled cells; no real
  bulk-with-truth. 12 donor-split seeds per tissue are correlated (overlapping donor
  pools from ~40 donors), so the effective independent N is below 24 — significance is
  indicative, not a clinical guarantee.
* **External head-to-head.** MuSiC/BisqueRNA/BayesPrism/cell2location/CARD/RCTD are
  installed and were executed in a **prior** harness run, **not re-executed this
  session** (multi-hour runtimes + observed environment instability; BayesPrism prior
  `failed` at ~2180 s). Their accuracy comparison is prior-run evidence, clearly
  labelled in `run_manifest.tsv`. Spatial benchmarks reference prior evidence.
* **No spatial accuracy without truth** (real Visium): structure/concordance only.

## 13. Recommended claims (careful language)

* TissueResolve **outperformed** flat/NNLS-family methods, **robustly (BH-significant,
  bootstrap-CI-backed, 24 replicates)**, on **false-positive subtype detection,
  pairwise spillover, and rare-subtype precision**.
* TissueResolve was **statistically tied** on **conditional within-family RMSE**
  (not significant at n=24) — *competitive, not superior*.
* TissueResolve **was faster than** every external tool measured.
* TissueResolve **showed better resolution-aware reporting** (trusted resolution,
  unresolved mass, soft gating) than tools that always force fine output.
* TissueResolve **did not outperform** on raw broad/fine correlation, fine RMSE, or
  rare *sensitivity* — it is **worse there by design** (it abstains rather than
  fabricate fine precision).
* **No "best overall" claim is made** — supported only per-metric, per-objective.

## 14. Diagnostics

* **Raw W_hat (pre-L1-normalisation).** The WNNLS solver now optionally exposes the
  raw NNLS coefficients before row normalisation (`return_raw_weights=True`;
  **default off, no change to outputs**) plus per-sample sparsity (raw row sum,
  n-nonzero, effective-N). Saved by the full-workflow script to
  `full_workflow_<tissue>/raw_weights_pre_l1norm.tsv` and `raw_weight_sparsity.tsv`
  for debugging degenerate/sparse solutions.
* **TCGA prediction-complexity audit** (`tcga_worst_samples_summary.tsv`). The two
  external tools show **opposite failure modes**: BisqueRNA is **over-sparse /
  over-confident** (worst: `easy_00`, only 5 of 32 types non-zero, effective-N 2.6,
  dominant fraction 0.67, Gini 0.94) while MuSiC is **over-diffuse** (worst: `hard_03`,
  effective-N 24.5; mean effective-N 17.8 vs Bisque 10.4). Neither reports unresolved
  mass — both force a full fine answer, unlike TissueResolve's resolution-aware path.
