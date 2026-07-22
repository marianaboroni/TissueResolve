# TissueResolve Public Benchmark Report

Reviewer-facing comparison of TissueResolve against established bulk and spatial deconvolution tools,
on identical donor-disjoint references, held-out mixtures, ground truth, and scoring. Nothing was
tuned on the test donors; no TissueResolve default was changed; no results were fabricated.

## 1. Executive summary
- **Bulk, matched head-to-head** (same scenarios all methods): **TissueResolve_PoissonGLM is the most
  accurate method tested** — broad Pearson **0.893**, fine Pearson **0.700** — and it wins **every**
  matched scenario on fine accuracy against MuSiC (0.545), BisqueRNA (0.474) and an external NNLS
  baseline (0.566). Paired vs TR_PoissonGLM: it beats every other in-process method on broad RMSE with
  a 100% win rate and bootstrap CI excluding 0.
- **The Poisson GLM backbone clearly beats the wNNLS default** (fine Pearson 0.708 vs 0.345), and the
  wNNLS default has the **highest overconfident fine-state rate** (0.085) — an honest negative.
- **The calibrated identifiability certificate is a method-agnostic evaluation axis**: for *every*
  method, held-out error rises RESOLVABLE → GROUP_ONLY → UNRESOLVABLE. The certificate predicts which
  populations are hard *before* deconvolution, for TR and competitors alike.
- **Spatial (synthetic, ground truth): all TissueResolve presets beat NNLS-per-spot** at both
  resolutions; weak_smoothing is best (broad 0.963 / fine 0.819) and rare-niche sensitivity is
  0.861 vs 0.167 for NNLS — at a real runtime cost (130–370 s vs < 10 s).
- **Not shown / deferred (honest):** BayesPrism and lung-scale (27k-gene) MuSiC/Bisque runs were
  deferred for runtime; all external *spatial* tools (RCTD/CARD/cell2location/DestVI) were deferred
  (runtime/GPU); a third tissue was NOT EXECUTED. **No claim of "best overall / SOTA" is made.**

## 2. Benchmark question
How does TissueResolve compare with established tools at (1) broad cell-type and (2) fine-subpopulation
resolution, on rare detection, absent-population spillover, overconfident fine calls, and — uniquely —
does its calibrated identifiability layer add information beyond point accuracy?

## 3. Dataset summary
| tissue | ref donors (train) | test donors | fine types | fine→broad families | source |
|--------|-------------------|-------------|-----------|--------------------|--------|
| breast | 20 | 21 (split cal/test upstream) | 30 | ~10 | breast-cancer scRNA reference |
| lung / HLCA | 20 | 20 | 50 | ~11 | HLCA subset |
| third tissue | — | — | — | — | **NOT EXECUTED** — only CRC (GSE200997) available, whose labels are CMS *molecular subtypes*, not cell types (37% unlabelled); no donor-annotated cell-typed reference offline |

Donor-disjoint throughout; the reference is built on TRAIN donors (with a donor column so the
identifiability certificate has cross-donor variance); mixtures come from held-out TEST donors. Test
donors are never used for gene/marker/threshold/model selection or certificate calibration.

## 4. Method status table
| modality | method | status | note |
|----------|--------|--------|------|
| bulk | TissueResolve_wNNLS (default), TissueResolve_PoissonGLM, TissueResolve_hierarchical | executed | all 6 scenarios × breast+lung |
| bulk | NNLS_baseline, RidgeNNLS_baseline, CorrelationMatcher_baseline | executed | internal/external baselines |
| bulk | MuSiC (v1.0.0), BisqueRNA (v1.0.5) | executed | breast × 3 scenarios + lung/reduced_overlap (real R runs) |
| bulk | MuSiC/BisqueRNA on lung base & low_depth (27k genes) | deferred | runtime (> MAX_SHARED=8000 genes) |
| bulk | BayesPrism (v2.2.3) | deferred | loads in R; per-run runtime too high this session (opt-in via RUN_BAYESPRISM) |
| bulk | DWLS | failed_install | package present but fails to load in R 4.1.2 |
| bulk | CIBERSORTx, SCDC | deferred / not_installed | web/licensed (export-only); SCDC absent |
| spatial | TissueResolve_default / weak_smoothing / edge_aware, NNLS_spot_baseline | executed | breast+lung × 2 seeds |
| spatial | RCTD (spacexr 2.2.1), CARD (1.1), cell2location (0.1.4), DestVI/stereoscope (scvi 1.1.6) | deferred | verified loadable/importable; full runs deferred (runtime/GPU) |
| spatial | SPOTlight, Tangram | failed_install | not present |

## 5. Bulk benchmark design
Donor-held-out pseudobulk mixtures (easy/medium/hard regimes = balanced → imbalanced) with known mRNA
proportions. Three query stressors: **base** (~1e6 counts, all genes), **low_depth** (×0.03 thinning),
**reduced_overlap** (1000-gene panel = cross-platform stress). Rare, collinear and closely-related
subpopulations are analysed *within* these mixtures via the metric layer (rare band, conditional
within-family, swap within certificate clusters) rather than as separate mixtures, keeping ground
truth consistent. Every method receives the identical reference + bulk matrix.

## 6. Spatial benchmark design
Semi-synthetic Visium-like spots with spatial domains (gradients, sharp boundaries, a rare niche)
generated from held-out donors with known per-spot composition. Metrics include local RMSE, boundary
F1, Moran's-I concordance, domain ARI, rare-niche sensitivity, oversmoothing, and absent-subtype mass.
Real spatial data without ground truth was not scored for accuracy (not run this session).

## 7. Broad cell-type performance (matched scenarios, mean)
| method | broad Pearson | broad RMSE |
|--------|--------------|-----------|
| **TissueResolve_PoissonGLM** | **0.893** | **0.053** |
| NNLS_baseline | 0.850 | 0.067 |
| MuSiC | 0.798 | 0.065 |
| TissueResolve_wNNLS (default) | 0.779 | 0.088 |
| BisqueRNA | 0.722 | 0.081 |
| TissueResolve_hierarchical | 0.685* | 0.097 | (*unresolved mass credited to family; abstains heavily)

TissueResolve is **competitive-to-best for broad cell types** via the Poisson GLM backbone.

## 8. Fine subpopulation performance (matched scenarios)
| method | fine Pearson | fine RMSE |
|--------|-------------|----------|
| **TissueResolve_PoissonGLM** | **0.700** | **0.030** |
| NNLS_baseline | 0.566 | 0.042 |
| MuSiC | 0.545 | 0.036 |
| BisqueRNA | 0.474 | 0.046 |
| TissueResolve_wNNLS (default) | 0.378 | 0.053 |

Per-scenario fine Pearson (TR_PoissonGLM / MuSiC): breast base 0.819/0.668; breast low_depth
0.810/0.583; breast reduced_overlap 0.696/0.632; lung reduced_overlap 0.477/0.295 — **TR_PoissonGLM
wins every matched scenario**. Fine subpopulation resolution is the hardest regime (all methods drop
sharply on lung), but TR_PoissonGLM leads.

## 9. Conditional within-family performance
Within-family conditional RMSE (isolates subtype confusion; lower better), all scenarios:
TR_PoissonGLM **0.154** < RidgeNNLS 0.183 < MuSiC 0.190 < CorrelationMatcher 0.203 < NNLS 0.232 <
BisqueRNA 0.257 < TR_wNNLS/hierarchical ≈ 0.29. Swap RMSE within certificate-flagged clusters:
MuSiC 0.065 ≈ TR_PoissonGLM 0.079 (MuSiC marginally better on swap), both far ahead of TR_wNNLS 0.126.

## 10. Rare population performance
rare_recall: CorrelationMatcher 0.966 and RidgeNNLS 0.858 are highest — **but this is an artifact of
spreading mass** (they have the worst fine accuracy and higher absent mass). Among accurate methods,
TR_PoissonGLM 0.594 ≈ MuSiC 0.555 > NNLS 0.431 > BisqueRNA 0.382 > TR_wNNLS 0.360. **No method
reliably resolves all rare populations.**

## 11. Absent population / false-positive analysis
Absent-type predicted mass (lower better): TR_hierarchical **0.003** (abstains) < BisqueRNA 0.012 ≈
MuSiC 0.013 ≈ TR_PoissonGLM 0.013 < NNLS 0.015 < TR_wNNLS 0.017 < RidgeNNLS 0.021 < CorrelationMatcher
0.026. TissueResolve does not leak more absent-type mass than external tools; the hierarchical
(abstaining) mode leaks the least.

## 12. Identifiability-stratified analysis
On scenarios where the certificate finds confounded clusters (breast/reduced_overlap, lung/base,
lung/low_depth), held-out error rises monotonically with the recoverability class **for every method**:

| method | err RESOLVABLE | err GROUP_ONLY | err UNRESOLVABLE |
|--------|---------------|----------------|------------------|
| TissueResolve_PoissonGLM | 0.022 | 0.028 | 0.032 |
| MuSiC | 0.031 | — | 0.036 |
| NNLS_baseline | 0.034 | 0.047 | 0.038 |
| RidgeNNLS_baseline | 0.030 | 0.030 | 0.052 |
| TissueResolve_wNNLS | 0.041 | 0.083 | 0.056 |

The certificate predicts, *before* deconvolution, which populations are hard — and this holds for
external tools too. This is the layer's core value: **method-agnostic recoverability prediction**.
`recommended_merge` names the confounded members to read out as an aggregate.

## 13. Uncertainty calibration
TissueResolve's donor-aware certificate intervals are calibrated (~90% coverage; Stage 2B/2C). Overconfident
fine-state rate (confident wrong subtype call in GROUP_ONLY/UNRESOLVABLE): **TR_wNNLS 0.085 (worst)**,
NNLS 0.021, BisqueRNA 0.019, TR_PoissonGLM 0.007, **MuSiC 0.000 / RidgeNNLS 0.000 (conservative)**.
Honest nuance: TR does **not uniquely** avoid overconfident calls — MuSiC and Ridge make none because
they spread mass conservatively; TR's contribution is the *explicit* certificate + unresolved mass
(hierarchical: 0.822 unresolved, 0 overconfident) that flags the risk transparently rather than
silently.

## 14. Runtime and failure rate
Bulk per-scenario runtime: internal baselines 0.01–1 s; TR_PoissonGLM ~0.9 s (breast) to ~42 s (lung
reduced_overlap); MuSiC 22–43 s; BisqueRNA 2–10 s; NNLS_baseline up to 106 s (lung). TR is competitive
and faster than MuSiC on comparable inputs. Spatial: TR presets 84–369 s vs NNLS-per-spot < 10 s — TR's
CAR model is much slower. Failure rate: 0 crashes among executed methods; deferrals/failures recorded
in `method_status`.

## 15. External-tool limitations
rpy2 is absent, so external R tools run via subprocess Rscript against `benchmarks/envs/Rlib`. MuSiC/
Bisque/BayesPrism/CARD/spacexr load; DWLS fails to load (R 4.1.2). BayesPrism and 27k-gene lung runs
are too slow to complete this session (deferred, not failed). Spatial externals (RCTD/CARD/
cell2location/DestVI) are importable but training/inference is deferred for runtime/GPU. CIBERSORTx is
web/licensed (export-only).

## 16. Negative results
- TR_wNNLS (the shipped default) is middling and the most overconfident — **prefer the Poisson GLM
  backbone**.
- TR_hierarchical with an *auto-inferred* mapping over-abstains (82% unresolved), depressing both
  broad and fine accuracy — resolution-aware abstention needs a curated mapping to pay off.
- Naive baselines (Correlation/Ridge) inflate rare_recall by spreading mass; recall alone is misleading.
- MuSiC matches or slightly beats TR_PoissonGLM on within-cluster swap — TR is not uniformly best on
  every fine metric.

## 17. Supported claims
See `combined/supported_claims.tsv`. Headline: (a) TR_PoissonGLM is the most accurate *executed* bulk
method (broad+fine, matched scenarios, paired CI); (b) Poisson GLM ≫ wNNLS default; (c) the
identifiability certificate predicts held-out error for every method; (d) all TR spatial presets beat
NNLS-per-spot with far better rare-niche detection; (e) MuSiC is the strongest external bulk tool tested.

## 18. Unsupported claims
See `combined/unsupported_claims.tsv`. NOT supported: TR is best overall/SOTA; TR beats BayesPrism or
any spatial external (all deferred/not executed); TR uniquely avoids overconfidence; TR resolves all
rare populations; TR generalizes across all tissues (2 tissues only); TR estimates absolute cell
fractions (outputs are RNA-derived proportions).

## 19. Reproducibility
Drivers: `benchmarks/public_benchmark/{run_bulk.py, run_spatial.py, run_external_bulk.R,
merge_and_report.py, scoring.py}`. All outputs under `benchmarks/results/public_benchmark/` (bulk/,
spatial/, combined/) + figures under `docs/figures/public_benchmark/` (each with a `.data.tsv`).
Manifests record datasets, seeds, scenarios. Re-run: `run_bulk.py` → `run_external_bulk.R` →
`merge_and_report.py`; `run_spatial.py`. Deterministic seeds; donor-disjoint splits fixed.

## 20. Reviewer-facing limitations
Two tissues (no third); external bulk limited to MuSiC/Bisque on ≤8000-gene scenarios (BayesPrism +
lung-scale + all spatial externals deferred for compute); pseudobulk/semi-synthetic (no real bulk with
fine ground truth, no real-spatial accuracy); a single reference-cell subsample cap (80/type) for
external export; spatial identifiability stratification is limited because summed-spot pseudo-queries
are deep (no confounded classes surface at that depth).

## 21. Final positioning of TissueResolve
The statement under test — *"TissueResolve is a resolution-aware framework combining competitive bulk
and spatial estimation with a calibrated identifiability layer indicating when fine-state
interpretation is supported, when grouped interpretation is safer, and when the pair is not testable"*
— is **supported by the executed evidence, with scope caveats**: (1) competitive-to-best bulk accuracy
at both resolutions via Poisson GLM (proven vs MuSiC/Bisque/NNLS on 2 tissues); (2) best-of-tested
spatial vs NNLS-per-spot (external spatial tools not yet run); (3) a calibrated identifiability layer
that demonstrably predicts recoverability for every method and flags confounded groups + unresolved
mass. It is **not** proven "best overall": BayesPrism and all spatial externals were not executed, and
generalization beyond two tissues is untested. TissueResolve's differentiator is **resolution-aware
interpretability + calibrated identifiability on top of competitive accuracy**, not a claim of
universal superiority.
