# Spatial Weak-Smoothing Benchmark Report

Experimental evaluation of a reduced spatial-smoothing preset for TissueResolve,
on **two datasets** (breast + HLCA/lung) across **5 seeds**, with an external-tool
(CARD/RCTD/cell2location) reference. No default changed; no Redeconve-inspired
regularization; no new solver.

Reproduce:

```bash
PYTHONPATH=src:. python benchmarks/spatial/run_weak_smoothing_grid.py \
    --run-real-data --datasets breast lung --seeds 0 1 2 3 4
# external (lung):
PYTHONPATH=src:. python benchmarks/spatial/export_lung_spatial_external.py --run-real-data
Rscript benchmarks/spatial/methods/run_card_synthetic.R \
    benchmarks/outputs/spatial_synthetic_lung/external_inputs
```

## 1. Motivation

The gold-truth spatial benchmark showed TissueResolve's flat spatial mode leads
on fine accuracy / local RMSE but **oversmooths** more than CARD/RCTD/cell2location
(oversmoothing score ≈ 1.8–1.95 vs 0.6–0.97). The NB-CAR prior at the default
λ=0.1 spreads mass across neighbouring spots. This report tests whether a lower λ
reduces oversmoothing while preserving TissueResolve's strengths — and whether the
effect **replicates across datasets** before any strong recommendation.

## 2. Current default behavior (Task 1)

- Default: `SpatialSolverConfig.lambda_spatial = 0.1` (`src/tissueresolve/config.py`).
  `0` disables spatial smoothing.
- Settable via CLI (`--lambda-spatial`) and config; `deconv_spatial(config=cfg)`
  honours it. λ recorded in result metadata + captions.
- **Change site:** opt-in preset layer
  `src/tissueresolve/experimental/spatial_presets.py` + CLI `--spatial-preset`.
  Default behaviour untouched (`default` preset is a no-op).

## 3. Lambda grid

`[0.0, 0.01, 0.02, 0.03, 0.05, 0.1]`; 0.1 = package default, 0.02 = candidate
`weak_smoothing` preset, 0.0 = no-smoothing control.

## 4. Datasets / scenarios

- **Breast** (`real_breast_cancer` reference, ~25 ref types) and **HLCA/lung**
  (`hlca_subset`, 49 ref types after the min-cell filter, broad mapping from
  `cell_type_broad`). Both donor-disjoint.
- The spatial generator was generalised (`domain_families`) so the structured
  scenario (sharp boundary + gradient + rare niche + mixed/collinear spots) maps
  to each tissue's broad families. Lung rare niche = `NK cells`; breast =
  `regulatory T cell`.
- **5 scenario seeds** per dataset (each a distinct realisation), 400 spots
  (20×20), 40 cells/spot. Ground truth = per-spot RNA-derived proportions.

## 5. Metrics

Accuracy, spatial structure, composition, runtime (same panel as before).
Outputs in `benchmarks/outputs/spatial_weak_smoothing/`
(`lambda_grid_metrics.tsv` now has a `dataset` column;
`per_scenario_metrics.tsv`, `per_family_metrics.tsv`, `runtime_metrics.tsv`).

## 6. Results by lambda (mean over 5 seeds)

**Breast**

| λ | fine Pear | broad Pear | local RMSE | oversmooth | boundary F1 | rare sens | FP-subtype |
|---|---|---|---|---|---|---|---|
| 0.00 | 0.773 | 0.920 | 0.038 | 0.36 | 0.327 | 0.356 | 0.001 |
| 0.02 (weak) | 0.865 | 0.953 | 0.024 | 1.31 | 0.353 | 0.667 | 0.002 |
| 0.10 (default) | 0.859 | 0.935 | 0.024 | 1.80 | 0.260 | 0.711 | 0.002 |

**Lung**

| λ | fine Pear | broad Pear | local RMSE | oversmooth | boundary F1 | rare sens | FP-subtype |
|---|---|---|---|---|---|---|---|
| 0.00 | 0.486 | 0.868 | 0.051 | 0.52 | 0.326 | 0.622 | 0.054 |
| 0.02 (weak) | 0.591 | 0.923 | 0.039 | 1.52 | 0.283 | 0.822 | 0.087 |
| 0.10 (default) | 0.571 | 0.899 | 0.039 | 1.85 | 0.270 | 0.756 | 0.116 |

## 7. Boundary and oversmoothing results

- Oversmoothing decreases monotonically with λ on both datasets. weak vs default:
  **breast −27%** (1.80→1.31, passes the ≥20% gate); **lung −18%** (1.85→1.52,
  **fails** the ≥20% gate at λ=0.02). The reduction is dataset-dependent in
  magnitude.
- Boundary F1: improved on breast (0.260→0.353) and marginally on lung
  (0.270→0.283).
- Even at weak λ, TissueResolve oversmooths more than CARD (see §16).

## 8. Rare niche results — DATASET-DEPENDENT (key finding)

- **Breast: weak smoothing HURTS rare-niche sensitivity** (0.711 → 0.667).
- **Lung: weak smoothing HELPS rare-niche sensitivity** (0.756 → 0.822).

The rare-niche effect points in **opposite directions** across the two datasets.
This single inconsistency is decisive against a strong/universal recommendation.

## 9. Fine/broad accuracy tradeoff

- Fine Pearson improves on **both** (breast 0.859→0.865; lung 0.571→0.591).
- Broad Pearson improves on **both** (breast 0.935→0.953; lung 0.899→0.923).
- Local RMSE preserved on both (Δ ≤ 0.0002).
- Conditional within-family RMSE **slightly worsens** on both (breast
  0.133→0.148; lung 0.276→0.283) — smoothing stabilises collinear within-family
  subtypes.

## 10. Composition / effective-N

- Effective-N: breast truth 17.0, weak 16.7 (not inflated, ✓); lung truth 12.2,
  weak 16.7, default 19.3 — **lung weak is closer to truth than default but still
  inflated beyond it** (gate 9 fails on lung).
- False-positive subtype detection / spillover **decreases** with weak smoothing
  on both (notably lung 0.116→0.087).

## 11. Recommended preset

`weak_smoothing` (λ=0.02) stays **experimental and opt-in**. It is a reasonable
choice for **breast-like** boundary/oversmoothing-focused analyses, but its
benefits are **not uniform across tissues** (see §8, §7). λ=0.03 is similar.

## 12. Whether default should remain unchanged

**Yes — default (λ=0.1) unchanged.** Across two datasets the promotion gates do
**not** pass uniformly (§ gate summary), so there is no basis to change the
default.

## External-tool comparison (Task 4)

Same seed-0 scenarios. RCTD (not currently installable) and cell2location
(~26 min/run) were **not** re-run for the new lung scenario; breast external
numbers are from the gold-truth spatial benchmark; CARD was run on both.

**Breast (seed 0):** TR_flat fine_r 0.799 / oversmooth 1.95; CARD 0.674 / 0.59;
RCTD 0.538 / 0.87; cell2location 0.681 / 0.96.

**Lung (seed 0):**

| method | fine Pear | broad Pear | oversmooth | boundary F1 | FP-subtype |
|---|---|---|---|---|---|
| TR default (0.1) | 0.594 | 0.928 | 1.88 | 0.351 | 0.097 |
| TR weak (0.02) | 0.614 | 0.951 | 1.52 | 0.286 | 0.072 |
| CARD | 0.518 | 0.935 | 0.97 | 0.328 | 0.055 |

- TissueResolve **beats CARD on fine Pearson** on both datasets.
- **CARD oversmooths much less** than TissueResolve on both; weak_smoothing
  **partially closes** the gap (lung 1.88→1.52, breast 1.95→1.31) but never
  reaches CARD's level (~0.6–0.97).
- CARD has the lowest false-positive spillover; weak_smoothing reduces
  TissueResolve's spillover toward it.

## 13. Supported claims

- Across breast + lung (5 seeds), `weak_smoothing` (λ=0.02) **improves fine and
  broad Pearson, reduces oversmoothing, and reduces false-positive spillover,
  while preserving local RMSE** — a consistent, modest improvement on those axes.
- It moves TissueResolve toward CARD's lower-oversmoothing regime without losing
  TissueResolve's fine-accuracy lead over CARD.

## 14. Unsupported claims

- **NOT** claimed: weak smoothing is universally better — the **rare-niche effect
  reverses between datasets** (helps lung, hurts breast).
- **NOT** claimed: it should be the default or a strong recommendation.
- **NOT** claimed: it reaches CARD-level oversmoothing.
- **NOT** claimed: within-family conditional accuracy improves (it slightly
  worsens on both).
- **NOT** claimed: any real-Visium effect (synthetic ground truth only).

## 15. Remaining limitations

- Two datasets, 5 seeds each, one structured scenario type per dataset.
- Rare-niche metric is coarse (9-spot niche) with high seed variance.
- Lung effective-N remains inflated even at weak λ.
- External comparison on lung is CARD-only (RCTD not installed; cell2location
  runtime-infeasible for a new scenario in-session); breast external is seed-0.
- Within-family conditional RMSE degrades slightly under weak smoothing — the
  collinearity-stabilisation a future Redeconve-inspired regulariser targets is
  out of scope here.

## Promotion-gate summary (weak 0.02 vs default 0.1; mean over 5 seeds)

| # | gate | breast | lung |
|---|---|---|---|
| 1 | oversmoothing −≥20% | ✅ −27% | ❌ −18% |
| 2 | fine Pearson preserved/improved | ✅ | ✅ |
| 3 | local RMSE preserved/improved | ✅ | ✅ |
| 4 | broad Pearson within 2% | ✅ improved | ✅ improved |
| 5 | boundary preservation improved/preserved | ✅ | ✅ |
| 6 | rare-niche detection improved/preserved | ❌ 0.711→0.667 | ✅ 0.756→0.822 |
| 7 | FP subtype not increased | ✅ | ✅ |
| 8 | pairwise spillover not increased | ✅ | ✅ |
| 9 | effective-N not inflated beyond truth | ✅ 16.7≤17.0 | ❌ 16.7>12.2 |
| 10 | stable in >1 scenario | ✅ 5 seeds | ✅ 5 seeds |
| 11 | replicate in >1 dataset | partial — see below | partial |
| 13 | improvement in BOTH datasets | mixed (accuracy ✅; rare niche ✗ on breast) | mixed |

**Verdict — does weak_smoothing become a strong recommendation?**

**No.** The gates do **not** pass uniformly across datasets: breast fails the
rare-niche gate, lung fails the oversmoothing-≥20% and effective-N gates, and the
rare-niche effect **reverses direction** between tissues. weak_smoothing delivers
consistent but modest gains (fine/broad Pearson, oversmoothing reduction, lower
spillover) and is retained as an **experimental, opt-in preset** for
boundary/oversmoothing-focused analyses — **but it is not promoted to a strong
recommendation and the default λ=0.1 is unchanged.**
