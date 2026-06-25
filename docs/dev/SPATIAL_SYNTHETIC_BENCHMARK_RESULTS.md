# Synthetic Spatial Benchmark — Ground-Truth Results (Phase 3, TissueResolve-only)

**What this is.** The first **spatial** accuracy benchmark with known ground
truth, offline / no external tools. A 20×20 = **400-spot** grid with real spatial
structure (sharp vertical domain border + smooth top→bottom gradient + a rare-cell
niche) is realised **count-level** by sampling **held-out query-donor** cells per
spot; the reference is built from the disjoint reference donors. Predictions are
scored against the **per-spot mRNA-proportion truth** with per-spot accuracy
**and** spatial-fidelity metrics (PART 12). Runner:
`benchmarks/spatial/run_synthetic_spatial.py`; outputs:
`benchmarks/outputs/spatial_synthetic/`.

Methods: TissueResolve spatial **flat** and **hierarchical** (spatial-aware
NB-CAR) vs an **NNLS-per-spot** (spatial-naive) baseline.

---

## Results (400 spots, 30 cell types, seed 0)

| Method | fine Pearson | broad Pearson | dominant acc | local RMSE↓ | oversmoothing | domain ARI |
|---|---|---|---|---|---|---|
| **TissueResolve_spatial_flat** | **0.799** | **0.933** | **0.708** | **0.026** | 1.95 | 0.271 |
| NNLS_per_spot (naive) | 0.678 | 0.913 | 0.507 | 0.040 | 0.42 | 0.270 |
| TissueResolve_spatial_hierarchical | 0.112 | −0.275 | 0.020 | 0.057 | 2.75 | 0.181 |

Mean true Moran's I = 0.30 (predicted: NNLS 0.13, flat 0.58, hier 0.45).

## Conclusions

1. **The spatial-aware model adds real value.** TissueResolve_spatial_flat beats
   the spatial-naive NNLS-per-spot baseline on every accuracy axis — fine Pearson
   0.799 vs 0.678, dominant-type accuracy 0.71 vs 0.51, and **local
   (neighbourhood) RMSE 0.026 vs 0.040**. The NB-CAR spatial prior measurably
   improves per-spot and local accuracy over treating spots independently.
2. **Accuracy/smoothness tradeoff, surfaced honestly (PART 13 caveat).**
   TissueResolve **over-smooths ~2×** (predicted Moran's I 0.58 vs true 0.30;
   oversmoothing score 1.95), while NNLS-per-spot **under-smooths** (0.13, score
   0.42). The spatial prior buys accuracy at the cost of inflated spatial
   autocorrelation — exactly the effect the brief warns about ("high spatial
   smoothness may indicate oversmoothing"). It is a tunable trade (λ_spatial), not
   a free win; reported, not hidden.
3. **Hierarchical spatial underperforms badly here** (fine 0.11, broad −0.275).
   This mirrors — and exceeds — the bulk hierarchical issue. As in bulk, part is
   the abstention penalty + the strip-`unresolved`-before-broad-aggregation
   artifact (so the negative broad correlation is partly an evaluation artifact,
   not purely the method). **Flagged for investigation**: the spatial within-family
   gating + a fair resolved-subset re-score are needed before any conclusion.
4. **Domain recovery is modest for all methods** (ARI ≈ 0.27): predicted
   composition clusters into the true L/R/niche domains only weakly — a hard task
   on these compositionally-overlapping domains; no method is rescued by it.

## Honest caveats

- One scenario, one seed, one tissue — point estimates (the multi-split + CI
  machinery from bulk should next wrap the spatial runner).
- The `oversmoothing_score` is the **aggregate** mean-Moran ratio (a per-cell-type
  ratio is numerically unstable when a type's true Moran's I ≈ 0; the per-type
  ratio is *not* used). `morans_i_mae` / `morans_i_corr_true_pred` are the
  companion autocorrelation-preservation metrics.
- Synthetic spots are mini-pseudobulks of query-donor cells (donor-disjoint), not
  real Visium; segmentation/optical effects are not modelled.

## External spatial methods — executed (RCTD, cell2location); CARD failed

Attempted in isolated environments. Honest status
(`benchmarks/outputs/spatial_synthetic/external_spatial_status.tsv`):

| Method | Status | Version | Note |
|---|---|---|---|
| RCTD (spacexr) | **executed** | 2.2.1 | pure-R; installed (binary deps), full mode, 2078 s |
| cell2location | **executed** | 0.1.4 | pip into isolated py3.9 venv; scvi-tools 1.1.6, torch 2.8, **CPU**; 1561 s (ref 250 + spatial 6000 epochs, reduced for CPU) |
| CARD | **failed** | — | C++ `.so` fails to load (broken toolchain, as for BayesPrism) |
| SPOTlight | skipped | — | not attempted |
| Tangram | skipped | — | not attempted |

> Fix applied so cell2location would run: scvi-tools 1.1 removed the `use_gpu=`
> argument from `train()` (the pip resolver paired old cell2location 0.1.4 with a
> newer scvi-tools); the wrapper was updated to drop it (CPU is the default).

**Combined results (synthetic spatial, ground truth):**

| Method | family | fine Pearson | broad Pearson | dom. acc | local RMSE↓ | oversmoothing | domain ARI | runtime |
|---|---|---|---|---|---|---|---|---|
| **TissueResolve_spatial_flat** | TR | **0.799** | 0.933 | **0.708** | **0.026** | 1.95 | 0.271 | **12 s** |
| cell2location | external | 0.681 | 0.819 | 0.528 | 0.032 | **0.96** | 0.295 | 1561 s |
| NNLS_per_spot | TR | 0.678 | 0.913 | 0.507 | 0.040 | 0.42 | 0.270 | 6 s |
| RCTD | external | 0.538 | **0.971** | 0.428 | 0.041 | 0.87 | **0.313** | 2078 s |
| TissueResolve_spatial_hierarchical | TR | 0.112 | −0.275 | 0.020 | 0.057 | 2.75 | 0.181 | 13 s |

**Balanced, axis-by-axis conclusions (no single winner):**

1. **Fine-level accuracy: TissueResolve_spatial_flat leads** (0.799) — above
   cell2location (0.681), NNLS-per-spot (0.678), RCTD (0.538) — and also best on
   dominant-type accuracy and local RMSE.
2. **Broad-level + domain structure: RCTD leads** (broad Pearson 0.971, domain
   ARI 0.313). RCTD is strong at the family level and recovering tissue domains,
   weaker at fine resolution.
3. **Spatial calibration: the external methods win.** cell2location (oversmoothing
   0.96) and RCTD (0.87) match the true autocorrelation closely; **TissueResolve
   over-smooths ~2×** (1.95). An honest weakness — TissueResolve buys fine
   accuracy at the cost of inflated spatial autocorrelation.
4. **Runtime: TissueResolve is ~130–170× faster** (12 s vs RCTD 2078 s,
   cell2location 1561 s on CPU) — a large efficiency advantage at this scale.
5. **Hierarchical spatial remains the weak spot** (flagged above).

So the methods occupy different points on the accuracy/calibration/speed frontier:
TissueResolve = most accurate (fine) + fastest but over-smooths; RCTD = best broad
/ domain / calibration but slow + weaker fine; cell2location = balanced + best
calibrated. This is reported as a frontier, not a single ranking.

**Caveats:** one scenario / seed / tissue; CPU-only cell2location with **reduced
epochs** (6000 vs default 30000 — may understate its accuracy); RCTD/c2l output
absolute abundances renormalised to proportions; CARD/SPOTlight/Tangram not
executed. Multi-seed CIs are the next hardening step.

Figure: `figures/figJ_spatial_external`. Tables: `combined_spatial_metrics.tsv`,
`external_spatial_status.tsv`. Wrappers: `methods/run_rctd_synthetic.R`,
`methods/run_cell2location_synthetic.py`; export `export_spatial_external.py`;
scorer `score_external_synthetic.py`.

---

## Addendum — multi-seed CIs + paired Wilcoxon (PART 14)

The single run above used one 400-spot layout. To test robustness, the benchmark
was repeated across **5 seeds** (each = a new donor split + spatial layout) on
**144-spot grids** (smaller, to bound the ~15-min-per-run external cost — all 5
seeds × {RCTD, cell2location} executed, 0 failures). External runtime: RCTD ~868 s,
cell2location ~812 s/seed (CPU); TissueResolve ~6 s. Outputs:
`benchmarks/outputs/spatial_multiseed/` (`metrics_ci.tsv`, `pairwise_wilcoxon.tsv`,
figK). Replicate unit = seed; **n = 5 (small — CIs are wide and Wilcoxon cannot
reach p<0.05 below the 0.0625 floor at n=5)**.

| Method | fine Pearson [95% CI] | broad Pearson | oversmoothing (1=ideal) | domain ARI | runtime |
|---|---|---|---|---|---|
| RCTD | **0.794 [0.656, 0.873]** | **0.979 [0.975, 0.982]** | 0.78 | **0.311** | 868 s |
| TissueResolve_spatial_flat | 0.770 [0.715, 0.819] | 0.884 | 1.89 [1.71, 2.08] | 0.217 | **6 s** |
| NNLS_per_spot | 0.718 [0.665, 0.770] | 0.924 | 0.37 | 0.253 | 4 s |
| cell2location | 0.650 [0.603, 0.712] | 0.785 | **0.99 [0.89, 1.11]** | 0.278 | 812 s |
| TissueResolve_spatial_hierarchical | 0.066 [−0.01, 0.13] | −0.338 | 36 (degenerate) | 0.137 | 7 s |

Paired Wilcoxon vs `TissueResolve_spatial_flat` (fine Pearson, n=5): vs RCTD
**p=0.625** (effect −0.33 — *no significant difference*); vs cell2location
**p=0.0625, effect +1.00** (TissueResolve higher on *every* seed); vs NNLS
p=0.125 (+0.87); vs hierarchical p=0.0625 (+1.00).

**Robust, CI-backed conclusions (these supersede the single-run numbers):**

1. **The single-run "TissueResolve ≫ RCTD on fine accuracy" did NOT replicate.**
   Across 5 seeds, RCTD (0.794) and TissueResolve_flat (0.770) are **statistically
   tied** on fine Pearson (Wilcoxon p=0.625). RCTD has the **higher mean but much
   wider CI** (0.66–0.87 — high seed-to-seed variance); TissueResolve is **lower
   mean but more stable** (0.72–0.82). This is the multi-seed methodology doing its
   job: the earlier single-layout advantage was not robust. (Grids also differ —
   144 vs 400 spots — so absolute values are not directly comparable to the single
   run; the multi-seed estimate is the more reliable one.)
2. **RCTD is the broad-level and domain-structure leader** (broad 0.979, domain
   ARI 0.311), with tight CIs — its strength is family-level / domain recovery.
3. **cell2location is the best spatially-calibrated** (oversmoothing CI 0.89–1.11
   brackets the ideal 1.0) but trails on fine accuracy — TissueResolve_flat beats
   it on *every* seed.
4. **TissueResolve_flat: competitive + stable + ~130–150× faster**, but
   **over-smooths ~2×** robustly (CI 1.71–2.08) — a consistent, real weakness.
5. **Hierarchical spatial is reproducibly broken** (negative broad on every seed;
   degenerate autocorrelation) — a confirmed defect to investigate, not noise.

**Takeaway.** No single method dominates spatially: RCTD wins broad/domain and ties
on fine (but slow + high-variance); TissueResolve ties on fine, is the most
stable and far the fastest, but over-smooths; cell2location is best-calibrated but
less accurate. n=5 is small — more seeds would tighten these CIs (each external
seed is ~15 min CPU, so this is the cost/robustness trade-off made explicit).

Figure: `benchmarks/outputs/spatial_multiseed/figures/figK_spatial_multiseed_ci`.

---

*Source: `benchmarks/outputs/spatial_synthetic/`. Generator:
`benchmarks/shared/synthetic_spatial.py`; spatial metrics:
`benchmarks/shared/spatial_metrics.py`; runner:
`benchmarks/spatial/run_synthetic_spatial.py`. No core algorithm modified; nothing committed.*
