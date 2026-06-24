# TissueResolve — Performance Benchmark Report

Publishable performance summary for the TissueResolve alpha. All numbers come from
the repository's gold-truth, donor-disjoint benchmarks; deferred tools are reported
as deferred, never fabricated. Estimates are **RNA-derived (mRNA) proportions**.

## 1. Overview

TissueResolve is benchmarked for bulk and spatial deconvolution against available
external tools on synthetic-but-real-derived gold-truth (pseudobulk / pseudo-spots
built from held-out donor cells, so the truth is known). The headline results:
the opt-in **Poisson GLM bulk solver** is the strongest validated improvement; the
default spatial solver is locally accurate but oversmooths relative to CARD; and
state-regularization experiments are negative (kept as diagnostics).

## 2. Benchmark design

Donor-disjoint splits (reference donors vs held-out query donors). Bulk: pseudobulk
mixtures with known mRNA-proportion truth. Spatial: pseudo-spots realized from
held-out cells with a structured scenario (sharp boundary + gradient + rare niche +
mixed spots). Each method receives the **identical** reference and query. Truth is
the per-sample/per-spot mRNA-derived proportion — the same estimate type all scored
methods produce.

## 3. Datasets

- **Breast** — real breast-cancer single-cell atlas, 30 fine types, donor-disjoint.
- **HLCA / lung** — HLCA subset, 49 fine types (many collinear), donor-disjoint.

Bulk external comparison: breast (5 scenarios, seed 0) and lung (5 seeds × 3 scenarios).
Spatial external comparison: breast synthetic scenario (seed 0; CARD this session,
RCTD/cell2location from prior exact-scenario runs) + multi-seed internal metrics.

## 4. Methods compared

Bulk: TissueResolve `wNNLS` (default), Poisson GLM, NB GLM; external NNLS control;
MuSiC; BisqueRNA. **Deferred:** BayesPrism, DWLS, SCDC (install failed — deps/
compilation). CARD is spatial and is **excluded** from bulk.

Spatial: TissueResolve default (flat NB-CAR), `weak_smoothing`, `edge_aware_smoothing`;
CARD; RCTD; cell2location; NNLS-per-spot control. **Deferred:** SPOTlight, Tangram,
DestVI (not installed).

## 5. Metrics

Bulk: fine/broad Pearson & RMSE, conditional within-family RMSE, rare recall/
precision/FPR, false-positive subtype rate, absent-subtype mass, effective-N, runtime.
Spatial: fine/broad Pearson, local RMSE, oversmoothing score, boundary F1, domain ARI,
rare-niche sensitivity/precision/FPR, effective-N, runtime.

Source tables: `docs/figures/bulk_performance_summary.tsv`,
`docs/figures/spatial_performance_summary.tsv`,
`docs/figures/benchmark_method_status.tsv`, `docs/figures/negative_results_summary.tsv`.

## 6. Bulk benchmark results

Identical donor-disjoint inputs; means over scenarios. (effective-N truth: breast 16.2,
lung 16.6.)

### Breast
| method | fine r | broad r | cond-RMSE↓ | rare rec | rare prec | rare FPR↓ | effN |
|---|---|---|---|---|---|---|---|
| **TR Poisson GLM** (rec.) | **0.814** | **0.980** | **0.156** | 0.86 | 0.59 | 0.55 | **16.1** |
| TR NB GLM | 0.814 | 0.980 | 0.156 | 0.86 | 0.59 | 0.55 | 16.1 |
| TR wNNLS (default) | 0.671 | 0.786 | 0.233 | 0.37 | 0.69 | 0.17 | 11.9 |
| external NNLS | 0.796 | 0.964 | 0.219 | 0.59 | 0.57 | 0.41 | 13.4 |
| MuSiC | 0.669 | 0.854 | 0.176 | 0.92 | 0.54 | 0.62 | 17.3 |
| BisqueRNA | 0.389 | 0.610 | 0.280 | 0.28 | 0.73 | 0.07 | 9.4 |

### Lung (49 types, 5 seeds)
| method | fine r | broad r | cond-RMSE↓ | rare rec | rare prec | rare FPR↓ | effN | runtime |
|---|---|---|---|---|---|---|---|---|
| **TR Poisson GLM** (rec.) | **0.758** | **0.957** | **0.198** | 0.60 | **1.00** | **0.00** | 19.1 | ~7 s |
| TR NB GLM | 0.760 | 0.957 | 0.196 | 0.70 | 1.00 | 0.00 | 19.1 | ~8 s |
| TR wNNLS (default) | 0.566 | 0.860 | 0.297 | 0.00 | – | 0.03 | 11.6 | ~8 s |
| external NNLS | 0.611 | 0.919 | 0.297 | 0.40 | 0.07 | 0.21 | 13.0 | ~5 s |
| MuSiC | 0.619 | 0.851 | 0.237 | 0.93 | 0.20 | 0.31 | 22.9 | ~150–300 s |
| BisqueRNA | 0.378 | 0.657 | 0.300 | 0.60 | 0.06 | 0.31 | 15.2 | ~100–330 s |

**Interpretation.** The Poisson GLM leads on fine/broad accuracy, conditional
within-family RMSE, and effective-N calibration on **both** tissues, beating wNNLS
(within-tissue conditional RMSE −39% breast / −33% lung) and the external tools, at
20–40× lower runtime than MuSiC/Bisque. NB ≈ Poisson (the gain is the count
likelihood, not overdispersion). On lung the rare type (NK) is identifiable and the
GLM detects it at precision 1.0 / FPR 0; on breast the rare type (regulatory T cell)
is collinear (see §8). (Figures: `bulk_fine_pearson.png`, `bulk_conditional_rmse.png`.)

## 7. Spatial benchmark results

Synthetic gold-truth (breast, seed 0; multi-seed internal). No real-Visium accuracy
is claimed (no ground truth exists for real Visium).

| method | fine r | broad r | local RMSE↓ | oversmoothing (1=truth) | domain ARI | runtime |
|---|---|---|---|---|---|---|
| **TR spatial (default/flat)** | 0.770 | 0.884 | **0.026** | 1.95 | 0.271 | 10.2 s |
| CARD | 0.674 | 0.859 | 0.033 | **0.59** | 0.253 | 34.6 s |
| RCTD | 0.538 | **0.971** | 0.041 | 0.87 | **0.313** | ~35 min |
| cell2location | 0.681 | 0.819 | 0.032 | 0.96 | 0.295 | ~26 min |
| NNLS per-spot (control) | 0.718 | 0.924 | 0.041 | 0.37 | 0.253 | 5.3 s |

Experimental smoothing presets (synthetic breast/lung): `weak_smoothing` reduces
oversmoothing (breast 1.80→1.31, lung 1.85→1.52) but does **not** reach CARD's level;
`edge_aware_smoothing` improves boundary F1 but fails the strict ≥20% oversmoothing
gate; `combined_weak_edge_smoothing` is redundant with `weak_smoothing`.

**Interpretation.** TissueResolve's default spatial solver has the **best local RMSE**
and competitive fine accuracy, but **oversmooths the most** (1.95 vs CARD 0.59) — the
clearest spatial weakness. (Figure: `spatial_oversmoothing_vs_accuracy.png`.)

## 8. Rare-state detection

The bulk count GLM raises rare *recall*; whether precision/FPR can be controlled
depends on **identifiability**. The opt-in rare-detection layer
(`experimental/rare_detection.py`) makes this explicit via average-precision (AP) vs
base rate:

- **NK cell (lung; distinct markers):** detectable — AP ≈ 2× base rate; GLM achieves
  precision 1.0 / FPR 0.
- **regulatory T cell (breast; collinear):** not identifiable — AP 0.27 < base rate
  0.36; the gate can only trade recall for FPR (FPR 0.44→0.17 at recall 0.15→0.08).

(Figure: `rare_detection_precision_recall.png`.) The layer reduces FPR where the rare
type is identifiable and honestly flags when it is not — it does not over-call.

## 9. Adaptive resolution / identifiability diagnostics

Every run emits (reporting only; estimates unchanged) `adaptive_resolution.tsv` and,
when the raw reference is available, `reference_uncertainty.tsv` /
`state_reliability.tsv` / `family_reliability.tsv`. Families are classified as
`resolved_fine`, `partially_resolved_group`, `broad_only`, `unresolved_family`, or
`diagnostic_only` from separability + within-family conditioning + marker support +
reference reliability. This distinguishes identifiable fine states from those that
should be interpreted as grouped/broad/unresolved.

## 10. Negative-result experiments (kept deliberately)

| experiment | outcome | finding |
|---|---|---|
| NB overdispersion vs Poisson | no improvement | NB≈Poisson even with donor-level dispersion |
| combined_weak_edge_smoothing | redundant | ≈ weak_smoothing; edge term inert at low λ |
| post-fit state regularization | negative | did not improve conditional within-family RMSE |
| in-solver state regularization | negative | competition over-concentrates, laplacian over-smooths; worse than default |
| state-similarity regularized solver | negative | cannot recover non-identifiable collinear fine states |

**These are central, not incidental:** *regularization cannot recover non-identifiable
fine states when reference signatures are collinear*, and *simple smoothing cannot
fully solve spatial oversmoothing or fine-state collinearity*. TissueResolve's response
is honest abstention (unresolved mass, adaptive resolution), not overconfident splits.

## 11. Runtime

Bulk: TissueResolve solvers ~5–8 s/run vs MuSiC/Bisque ~100–330 s (20–40× faster).
Spatial: TissueResolve flat 10.2 s, CARD 34.6 s (practical); RCTD ~35 min,
cell2location ~26 min (run once on the fixed scenario).

## 12. Supported claims

- The opt-in Poisson GLM bulk solver is the strongest validated algorithmic
  improvement in TissueResolve to date and is **recommended as an experimental bulk
  solver**, while `wNNLS` remains the default.
- TissueResolve's default spatial solver has the best local RMSE and competitive fine
  accuracy on the synthetic scenario.
- Adaptive-resolution and reference-uncertainty diagnostics separate identifiable from
  non-identifiable fine states.

## 13. Unsupported claims

- ✗ best / state-of-the-art / validated on all tissues / solves all subtypes.
- ✗ NB > Poisson (they tie). ✗ cell fractions by default (RNA-derived proportions).
- ✗ real-Visium accuracy (no ground truth). ✗ closing the spatial oversmoothing gap to
  CARD. ✗ superiority over deferred tools (BayesPrism/DWLS/SCDC not run).

## 14. Limitations

Synthetic-but-real-derived gold truth (no real-Visium truth); breast + lung only;
bulk external comparison on MuSiC/Bisque only (others deferred); spatial external is a
single-scenario probe + multi-seed internal; spatial hierarchical multiseed harness is
degenerate and excluded (flat/default reported); outputs are mRNA proportions.

## 15. Reproducibility

Benchmarks: `benchmarks/diagnostics/nb_bulk_solver_benchmark.py`,
`benchmarks/diagnostics/external_bulk_comparison.py`,
`benchmarks/bulk/export_external_bulk_generic.py` + `run_external_bulk_generic.R` +
`benchmarks/diagnostics/score_external_bulk_generic.py` (lung),
`benchmarks/diagnostics/rare_detection_calibration.py`, and the spatial scripts under
`benchmarks/spatial/`. All require `--run-real-data` and local data; outputs are
gitignored. Method status: `docs/figures/benchmark_method_status.tsv`.
