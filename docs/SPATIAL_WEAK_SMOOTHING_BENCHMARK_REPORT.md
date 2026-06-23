# Spatial Weak-Smoothing Benchmark Report

Experimental evaluation of a reduced spatial-smoothing preset for TissueResolve.
Run 2026-06-23. No default changed; no Redeconve-inspired regularization; no new
solver.

Reproduce:

```bash
PYTHONPATH=src:. python benchmarks/spatial/run_weak_smoothing_grid.py --run-real-data
```

## 1. Motivation

The gold-truth spatial benchmark
(`docs/GOLD_TRUTH_SPATIAL_EXTERNAL_BENCHMARK_REPORT.md`) showed TissueResolve's
flat spatial mode leads on fine accuracy and local RMSE but **oversmooths** more
than CARD/RCTD/cell2location (oversmoothing score ≈ 1.95 vs 0.6–0.96). The NB-CAR
spatial prior at the default λ=0.1 spreads mass across neighbouring spots. This
report tests whether a **lower** λ reduces oversmoothing while preserving
TissueResolve's strengths.

## 2. Current default behavior (Task 1 findings)

- Default smoothing strength: `SpatialSolverConfig.lambda_spatial = 0.1`
  (`src/tissueresolve/config.py:174`). `0` disables spatial smoothing (pure
  NB-MAP per spot).
- Already settable from CLI (`tissueresolve spatial --lambda-spatial`) and config
  YAML; `deconv_spatial(..., config=cfg)` honours it.
- λ is recorded in `SpatialDeconvResult` metadata and figure captions.
- Oversmoothing / local RMSE / boundary / Moran's-I metrics live in
  `benchmarks/shared/spatial_metrics.py`.
- **Change site:** a new opt-in preset layer
  `src/tissueresolve/experimental/spatial_presets.py` + a CLI
  `--spatial-preset` option. The default lambda and default behaviour are
  untouched (the `default` preset is a no-op).

## 3. Lambda grid

`[0.0, 0.01, 0.02, 0.03, 0.05, 0.1]` where 0.1 is the package default, 0.0 a
no-smoothing control, and 0.02 the candidate `weak_smoothing` preset.

## 4. Datasets / scenarios

- **One dataset:** breast (`real_breast_cancer` single-cell reference, donor-
  disjoint split). **HLCA/lung spatial-like synthetic was not built** → results
  are single-dataset (see §15).
- **Three scenario seeds** (0,1,2): each is a distinct realisation of the
  structured scenario embedding **sharp boundary** (L/R domains), **top→bottom
  gradient**, **rare niche** (top-right block enriched for the rare type),
  **mixed spots**, and **high-collinearity fine subtypes** (within-family).
  400 spots (20×20), 40 cells/spot. Reference fixed; only the spatial scenario
  varies across seeds (isolates spatial-structure variation).
- Ground truth: per-spot RNA-derived (mRNA) proportions.

## 5. Metrics

Accuracy (broad/fine Pearson & RMSE, conditional within-family RMSE, local
RMSE), spatial structure (oversmoothing, boundary F1, edge blurring, Moran's-I
preservation, domain recovery ARI, rare-niche sensitivity/precision/FPR),
composition (entropy, effective-N, sparsity, dominant fraction, absent-subtype
mass, pairwise spillover), and runtime/warnings. Outputs:
`benchmarks/outputs/spatial_weak_smoothing/{lambda_grid_metrics,per_scenario_metrics,per_family_metrics,runtime_metrics}.tsv`.

## 6. Results by lambda (mean over 3 seeds)

| λ | preset | fine Pear | broad Pear | local RMSE | cond RMSE | oversmooth | boundary F1 | rare sens | FP-subtype | runtime s |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.00 | no_smoothing | 0.765 | 0.916 | 0.038 | 0.246 | 0.36 | 0.317 | 0.370 | 0.0011 | 1.1 |
| 0.01 | grid | 0.844 | 0.953 | 0.027 | 0.167 | 0.93 | 0.354 | 0.481 | 0.0013 | 12.7 |
| **0.02** | **weak_smoothing** | **0.859** | **0.952** | **0.025** | **0.149** | **1.32** | **0.352** | **0.593** | **0.0018** | 8.5 |
| 0.03 | grid | 0.860 | 0.948 | 0.024 | 0.143 | 1.49 | 0.340 | 0.630 | 0.0019 | 6.6 |
| 0.05 | grid | 0.859 | 0.942 | 0.024 | 0.138 | 1.65 | 0.307 | 0.630 | 0.0019 | 6.1 |
| 0.10 | default | 0.854 | 0.934 | 0.024 | 0.133 | 1.81 | 0.290 | 0.593 | 0.0019 | 4.5 |

## 7. Boundary and oversmoothing results

- **Oversmoothing decreases monotonically with λ** and is highly consistent
  across seeds. weak_smoothing (0.02) vs default (0.1): **1.32 vs 1.81 → 27%
  reduction** (gate: ≥20% ✓).
- **Boundary F1 improves** with weak smoothing: 0.290 → 0.352; edge blurring is
  lower (0.005 → 0.003). Sharper domain boundaries.
- λ=0.0 (no smoothing) over-corrects: oversmoothing 0.36 but fine Pearson
  collapses to 0.765 and boundary F1 drops — pure NB-MAP is too noisy.

## 8. Rare niche results

- Rare-niche sensitivity is **noisy** (the niche is a 9-spot block): per seed it
  varies 0.22–0.67. Averaged over 3 seeds, weak_smoothing **equals** the default
  (0.593 = 0.593) and is in fact **more stable** (no single-seed dip to 0.44 that
  the default shows at seed 1). Rare-niche precision is preserved (0.042 vs
  0.036) and false-positive rate not increased.
- λ=0.0 hurts rare-niche recall (0.370). Some smoothing helps the rare niche.

## 9. Fine/broad accuracy tradeoff

- **Fine Pearson is preserved/slightly improved** (0.854 → 0.859); fine RMSE
  unchanged (0.025).
- **Broad Pearson improves** (0.934 → 0.952, well within the 2% gate).
- **Conditional within-family RMSE is slightly worse** with less smoothing
  (0.133 → 0.149): smoothing helps stabilise collinear within-family subtypes.
  This is the clearest accuracy cost of weak smoothing.
- Local RMSE is essentially unchanged (+0.0001, within noise).
- Effective-N: truth 17.04; default predicts 17.97 (slightly inflated),
  weak_smoothing 16.92 (closer to truth, not inflated).

## 10. Runtime

weak_smoothing ~8.5 s vs default ~4.5 s per 400-spot scenario (lower λ → more
solver iterations to converge). Both are practical; no failures, warnings are the
usual separability notices.

## 11. Recommended preset

Expose **`weak_smoothing` (λ=0.02) as an experimental, opt-in preset** for
analyses that prioritise **boundary sharpness and reduced oversmoothing**
(e.g. delineating domain edges or rare niches), accepting a small increase in
within-family conditional RMSE. λ=0.03 is a near-equivalent alternative (slightly
higher rare-niche sensitivity but only ~17% oversmoothing reduction, below the
20% gate).

## 12. Whether default should remain unchanged

**Yes — the default (λ=0.1) is unchanged.** Although weak_smoothing passes 10 of
11 promotion gates, the evidence is from a **single dataset** (gate 11 not met).
Per the rules, single-dataset evidence does not justify changing the default.

## 13. Supported claims

- On the breast synthetic spatial benchmark (3 seeds), the experimental
  `weak_smoothing` preset (λ=0.02) **reduces oversmoothing by ~27% and improves
  boundary F1 while preserving fine Pearson, local RMSE, broad Pearson, and
  rare-niche detection**, without increasing spillover or false-positive subtype
  detection, and brings effective-N closer to truth.
- It may be useful for boundary / rare-niche-focused spatial analyses.

## 14. Unsupported claims

- **NOT** claimed: weak smoothing is universally better — within-family
  conditional RMSE worsens slightly (0.133 → 0.149).
- **NOT** claimed: it generalises beyond breast — only one dataset tested.
- **NOT** claimed: the default should change.
- **NOT** claimed: any real-Visium effect (synthetic ground truth only).

## 15. Remaining limitations

- **Single dataset** (breast); HLCA/lung spatial-like synthetic was not built →
  promotion gate 11 (multi-dataset replication) is unmet, so the preset stays
  experimental and the default is unchanged.
- Three seeds of one structured scenario; distinct scenario *types* are embedded
  in one generator rather than run as separate matrices.
- Rare-niche metrics are coarse (9-spot niche); seed-level variance is high.
- Conditional within-family RMSE degrades slightly under weak smoothing — the
  collinearity-stabilisation that a future Redeconve-inspired state regulariser
  targets is NOT addressed here (out of scope by instruction).

## Promotion-gate summary (weak_smoothing 0.02 vs default 0.1; mean over 3 seeds, breast)

| # | gate | result |
|---|---|---|
| 1 | oversmoothing −≥20% | ✅ −27% |
| 2 | fine Pearson preserved/improved | ✅ 0.854→0.859 |
| 3 | local RMSE preserved/improved | ✅ +0.0001 (noise) |
| 4 | broad Pearson within 2% | ✅ improved 0.934→0.952 |
| 5 | boundary preservation improved/preserved | ✅ F1 0.290→0.352 |
| 6 | rare-niche detection improved/preserved | ✅ 0.593=0.593 (more stable) |
| 7 | FP subtype detection not increased | ✅ 0.0019→0.0018 |
| 8 | pairwise spillover not increased | ✅ 0.0020→0.0015 |
| 9 | effective-N not inflated beyond truth | ✅ 16.92 vs truth 17.04 |
| 10 | stable in >1 scenario | ✅ 3 seeds |
| 11 | replicate in >1 dataset | ❌ breast only |

**Verdict:** 10/11 gates pass → **keep `weak_smoothing` as an experimental
opt-in preset; do NOT change the default** (single-dataset evidence; minor
within-family conditional cost).
