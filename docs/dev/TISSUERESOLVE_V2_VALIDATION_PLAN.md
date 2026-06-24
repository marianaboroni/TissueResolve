# TissueResolve v2 — Validation Plan

Defines the gates a v2 experimental method must pass **before** any default
change (rules 3, 4, 10, 11). Margins are fixed here, before the final test set is
touched.

## Frozen baseline
`benchmarks/baselines/current_baseline/` (Stage 0). All v2 comparisons are paired
against these tables. The held-out **test** donors/seeds are reserved: never used
for gene selection, threshold tuning, or confidence calibration (rules 10, 11).

## Data splits (nested, donor-aware)
- **train donors** → build signatures, select genes, learn modules;
- **validation donors / held-out synthetic** → tune thresholds, calibrate
  confidence, select λ;
- **test donors + ≥2 tissues** → final evaluation only.
Bulk replicate unit = split×scenario (n=25); spatial = seed (n≥5, more if a gate
is borderline).

## Evaluation axes (reported separately — never a single flattened Pearson)
broad accuracy · fine accuracy (absolute) · within-family conditional accuracy ·
rare-population recall/precision · pairwise spillover · unresolved precision/recall ·
false-resolution / false-abstention rate · composition complexity (eff-N, entropy,
richness) · spatial boundary preservation / oversmoothing / local RMSE · runtime ·
memory.

## Statistics (every comparison)
paired bootstrap 95% CIs · paired Wilcoxon · effect sizes · multiple-testing
(Benjamini–Hochberg) · bootstrap rank stability · per-dataset + per-cell-type
breakdowns. (Reuses `benchmarks/shared/stats.py`.)

## Promotion gates (ALL must hold simultaneously; non-inferiority margins)
1. broad RMSE: no worse than baseline by >2% (relative).
2. fine RMSE: ≥5% relative improvement.
3. within-family conditional accuracy: improved.
4. pairwise spillover: ≥10% relative reduction.
5. rare-population recall: improved, with false-positive rate increase ≤ a
   predefined limit (set: ≤+5 percentage points).
6. false-resolution rate: reduced.
7. unresolved precision/recall: maintained or improved.
8. predicted richness & entropy: closer to truth than baseline (|eff-N error| ↓).
9. spatial over-smoothing: reduced (oversmoothing ratio closer to 1.0).
10. spatial boundary score / local RMSE: no worse than baseline.
11. mass conservation: rows sum to 1 within 1e-6.
12. seed stability: rank-1 retained across ≥80% of bootstrap resamples.
13. improvement reproduced in **≥2 tissues**.
14. runtime ≤ 3× baseline (unless a documented accuracy gain justifies more);
    memory profiled and acceptable.
15. ties or beats baseline on every **safety-critical** metric
    (broad RMSE, mass conservation, false-resolution, unresolved precision).

A method failing any gate stays **experimental** (rule 17). If a simpler variant
ties a complex one, the simpler wins (rule 16).

## Mandatory comparison methods
- Bulk: TissueResolve baseline, NNLS, MuSiC, BisqueRNA, BayesPrism (when a working
  toolchain is available), + one more executable method.
- Spatial: TissueResolve baseline, RCTD, cell2location, CARD, SPOTlight (toolchain
  permitting), + one more. Executed-only are ranked; failed/skipped reported, never
  ranked (already enforced by the harness).

## Per-component validation (rule 15)
For soft_hierarchy / rna_content / distribution_alignment / nbcar_vi:
mathematical spec + identifiability · unit tests · synthetic recovery tests (can
the estimator recover known parameters?) · profiling · ablation (each component's
marginal effect, measured before combination) · docs · feature-status label.

## Deliverables at promotion review
`docs/TISSUERESOLVE_V2_FINAL_REPORT.md`, `benchmarks/outputs/v2_ablation_metrics.tsv`,
`v2_method_comparison.tsv`, `v2_runtime_memory.tsv`, `v2_promotion_gates.tsv`,
`v2_failure_analysis.tsv`. The target claim may be used **only** if gates pass on
≥2 tissues with the stated margins.
