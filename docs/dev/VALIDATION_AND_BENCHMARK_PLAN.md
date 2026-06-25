# Validation and Benchmark Plan (summary)

High-level phases

Phase 1 — Pre-benchmark audit

- TCGA sparsity audit (script + metrics) — implemented under
  `benchmarks/bulk/tcga_sparsity_audit.py` and `docs/BULK_PREDICTION_SPARSITY_AUDIT.md`.
- Unit checks: solver raw coefficients, hierarchical assembly, unresolved mass.

Phase 2 — Synthetic bulk benchmarks

- Generate pseudobulk mixtures from held-out donors with varying rarity,
  batch effects, and sequencing depths. Use `benchmarks/synthetic` helpers.
- Compute broad/fine metrics, rare detection curves, spillover matrices.

Phase 3 — External method execution

- Run MuSiC, Bisque, BayesPrism, CIBERSORTx (when available), DWLS, NNLS
  baselines using `benchmarks/bulk/methods` adapters.
- Record versions, runtimes, memory, success/failure.

Phase 4 — Spatial synthetic & external methods

- Simulate Visium-like spatial mixtures and run TissueResolve spatial
  plus external wrappers (cell2location, RCTD, CARD, SPOTlight).

Phase 5 — Report generation and statistical analysis

- Aggregate per-method metrics, bootstrap confidence intervals, paired tests,
  rank stability plots, and the final `benchmarks/outputs/benchmark_report.html`.

Reproducibility

- Use deterministic seeds in all synthetic generators.
- Save raw outputs and metadata for every run in `benchmarks/outputs`.

Notes

- This plan intentionally avoids changing solver algorithms before the
  diagnostic audit completes. Config defaults may be suggested but not
  pushed without clear evidence and tests.
# TissueResolve Validation & Benchmark Plan

Phased execution plan for the bulk + spatial validation/benchmark study. This is
the controlling document; it maps the 21-part brief onto executable phases,
records what **already exists** in the repo (so we don't rebuild it), and marks
what is **missing**. Heavy external-tool installs are kept out of the default
test environment (`CLAUDE.md` environment + testing rules).

**Bulk and spatial are benchmarked and ranked separately throughout.** No mixed
ranking. No accuracy claims on real data without ground truth.

---

## Status legend
✅ exists & usable · 🟡 partial / needs extension · ❌ missing · ⏳ planned

## Existing infrastructure (audited 2026-06-08)

| Capability | Status | Location |
|---|---|---|
| Bulk pipeline + solvers (nnls/weighted/marker/ridge/auto/ensemble) | ✅ | `src/tissueresolve/bulk/`, `solver/` |
| Hierarchical bulk + within-family resolvability gating | ✅ | `bulk/hierarchical.py`, `reference/hierarchy.py` |
| Spatial pipeline (NB-CAR, graph, QC, Moran's I, auto-params) | ✅ | `src/tissueresolve/spatial/` |
| Accuracy metrics: Pearson/Spearman/RMSE/MAE/bias/calibration/dominant | ✅ | `benchmarks/shared/metrics.py:39–88` |
| Hierarchical fair metrics: family-level, fine-resolvable, unresolved precision/recall | ✅ | `metrics.py:169–246` |
| No-GT spatial metrics: concordance, entropy, near-zero, marker recovery, JSD | ✅ | `metrics.py:92–123`, `spatial_multimetric_ranking.py` |
| Method registry + status (executed/failed/skipped/exported/imported) | ✅ | `shared/method_registry.py`, `base.py`, `imported.py` |
| Composite scorecards (bulk + spatial), separated | ✅ | `composite_score.py`, `spatial_multimetric_ranking.py` |
| Benchmark HTML report + figures | 🟡 | `shared/report.py`, `plotting/*benchmark*` |
| Synthetic toy bulk + toy spatial w/ ground truth | 🟡 | `shared/synthetic.py` |
| External wrappers (MuSiC/Bisque/BayesPrism/cell2loc/RCTD/CARD/SPOTlight/Tangram) | 🟡 | `benchmarks/{bulk,spatial}/methods/`, not yet executed |
| **TCGA complexity audit (PART 1/17)** | ✅ | `benchmarks/audit/tcga_sparsity_audit.py`, this study |

## Gaps to build

❌ JSD / CCC / Aitchison / total-variation for **bulk**; ❌ per-sample Gini/effective-N
in the standard runners (now exist in the audit, to be promoted to `metrics.py`);
❌ rare-population & cross-donor & spillover-intensive synthetic scenarios;
❌ image-based pseudo-spot aggregation (Xenium/MERFISH); ❌ bootstrap CIs +
paired Wilcoxon + FDR + rank-stability; ❌ data-driven robustness sweeps;
❌ decision-oriented report redesign (PART 16 sections + bubble heatmap, rank
stability, detection curves, richness/entropy scatter, Pareto).

---

## Phases

### Phase 1 — Audit + metrics + synthetic bulk  *(complete)*
1. ✅ **PART 1** bulk sparsity audit → `docs/BULK_PREDICTION_SPARSITY_AUDIT.md`.
2. ✅ **PART 17** TCGA complexity audit → `docs/TCGA_DECONVOLUTION_COMPLEXITY_AUDIT.md` + `benchmarks/outputs/tcga_prediction_complexity_audit.tsv`.
3. ✅ Metrics in `benchmarks/shared/metrics.py`: JSD, total-variation, Aitchison
   (multiplicative zero replacement), Lin's CCC, `compositional_metrics`, plus
   effective-N / Gini / `composition_complexity` panel (**PART 5/6/10**).
   Tests: `benchmarks/tests/test_compositional_metrics.py` (14).
4. ✅ The 8 PART-1 diagnostic figures →
   `benchmarks/outputs/figures/part1_audit/fig1..8` (PNG+SVG+`.data.tsv`+caption),
   via `benchmarks/audit/plot_sparsity_audit.py`.
5. ✅ Held-out-donor synthetic bulk generator
   `benchmarks/shared/synthetic_holdout.py`: count-level, donor-disjoint
   (reference vs query), records **both** cell-fraction and mRNA-proportion
   truth; scenarios balanced/imbalanced/rare(0.001–0.05)/similar-subtypes/
   missing/extra + cross-platform split (**PART 2A.1, 8, 9**). Validated on the
   real 41-donor h5ad. Tests: `benchmarks/tests/test_synthetic_holdout.py` (10).
   *Remaining for Phase 2: wire the generator into a runner that builds the
   donor-disjoint reference, runs the methods, and applies the metrics.*

### Phase 2 — External bulk methods + broad/fine/rare/spillover benchmark
0. ✅ **TissueResolve-only ground-truth runner** (no installs):
   `benchmarks/bulk/run_holdout_bulk_benchmark.py` — donor-disjoint, count-level,
   5 scenarios; fine/broad/complexity/rare metrics →
   `benchmarks/outputs/holdout_bulk/*.tsv` + figures; results in
   `docs/HOLDOUT_BULK_BENCHMARK_RESULTS.md`. Key finding: hierarchical mode
   collapses to ~2.7 effective populations and is worst on accuracy + rare
   detection; no flat solver dominates (nnls best fine accuracy, ridge best
   complexity preservation). Top Phase-2 question: hierarchy gating thresholds.
6a. ✅ **External bulk methods executed** on the SAME held-out inputs (isolated
   R 4.1.2 lib `benchmarks/envs/r_lib`, binary installs): **MuSiC v1.0.0** and
   **BisqueRNA v1.0.5** executed; **BayesPrism failed** (C++ source compilation
   broken in this env — no R-4.1 binary). Scripts: `export_external_inputs.py`,
   `run_external_holdout.R`, `score_external_holdout.py`. Results:
   `benchmarks/outputs/holdout_bulk/combined_*.tsv`, `method_status.tsv`, figF;
   write-up in `docs/HOLDOUT_BULK_BENCHMARK_RESULTS.md` Addendum 2. Findings:
   TissueResolve flat nnls/auto > MuSiC > ridge > BisqueRNA on fine Pearson;
   PART 17 answered — Bisque is *more* concentrated than TissueResolve,
   so sparsity is not TissueResolve-specific.
6b. ✅ **Multi-split CIs + paired Wilcoxon for MuSiC/BisqueRNA** (all 25
   split×scenario replicates, 0 failures): `export_external_multisplit.py`,
   `run_external_multisplit.R`, `score_external_multisplit.py` →
   `benchmarks/outputs/holdout_bulk/multisplit_combined/` + Addendum 3. Result:
   nnls significantly beats MuSiC (p=2e-4) and BisqueRNA (p<1e-4, every
   replicate) on fine + broad accuracy; MuSiC best preserves complexity
   (eff-N error 1.94). Remaining: BayesPrism/CIBERSORTx in a working
   toolchain/conda env; a second tissue.
6. ⏳ (remaining external) Execute in **isolated envs** (R: MuSiC, Bisque,
   BayesPrism; CIBERSORTx export→import only). Record version, env, runtime,
   peak memory, GPU, **execution status**. Failures reported honestly; skipped /
   exported-only methods **excluded from rankings** (**PART 3A, 19**).
7. ⏳ Broad-level (**PART 5**) and fine-level conditional+absolute (**PART 6**)
   metrics; rare-population detection curves (**PART 8**); spillover matrix +
   AUROC (**PART 9**); hierarchy/unresolved metrics (**PART 7**); complexity
   preservation score (**PART 10**); robustness sweeps (**PART 11**).
8. 🟡 Statistical layer: bootstrap CIs, paired Wilcoxon, effect size,
   rank-stability — **implemented** in `benchmarks/shared/stats.py` (tests:
   `test_stats.py`) and applied via `benchmarks/bulk/run_holdout_bulk_multisplit.py`
   (5 splits × 5 scenarios = 25 paired replicates) →
   `benchmarks/outputs/holdout_bulk_multisplit/` + addendum in
   `docs/HOLDOUT_BULK_BENCHMARK_RESULTS.md`. Remaining: FDR correction across the
   pairwise tests; same-subset fair comparison vs `hierarchical_resolvable`.
   Key result: strict→fair re-score flips hierarchical (0.32→0.89 fine Pearson),
   so its apparent failure was largely an abstention penalty; among flat solvers
   nnls/auto > ridge on accuracy AND average complexity preservation.

### Phase 3 — Synthetic spatial + external spatial methods
0. ✅ **Synthetic spatial benchmark (ground truth)** — `benchmarks/shared/synthetic_spatial.py`
   (sharp border + gradient + rare niche, donor-disjoint, count-level),
   `benchmarks/shared/spatial_metrics.py` (Moran's I preservation, oversmoothing,
   local RMSE, domain ARI; tests `test_spatial_synthetic_and_metrics.py`),
   runner `benchmarks/spatial/run_synthetic_spatial.py` →
   `benchmarks/outputs/spatial_synthetic/`; write-up
   `docs/SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md`. Findings: TissueResolve spatial
   (flat) beats NNLS-per-spot on accuracy (fine 0.80 vs 0.68, local RMSE 0.026 vs
   0.040) but over-smooths ~2×; hierarchical spatial underperforms (flagged).
9. ⏳ Synthetic spatial generator: domains, borders, gradients, mixed/rare niches,
   missing pop, variable library, noise; known broad + fine + domain truth
   (**PART 2B.1**). Optional image-based pseudo-spot aggregation (**PART 2B.2**).
10a. ✅ **External spatial methods executed** on the synthetic scenario:
   **RCTD (spacexr 2.2.1)** and **cell2location 0.1.4** (isolated py3.9 venv, CPU)
   ran; **CARD failed** (C++ `.so` load — toolchain); SPOTlight/Tangram not
   attempted. Scripts: `export_spatial_external.py`, `methods/run_rctd_synthetic.R`,
   `methods/run_cell2location_synthetic.py`, `score_external_synthetic.py` →
   `benchmarks/outputs/spatial_synthetic/combined_spatial_metrics.tsv`,
   `external_spatial_status.tsv`, figJ; Addendum in
   `docs/SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md`. Findings (frontier, no single
   winner): TissueResolve_flat best fine accuracy + ~130× faster but over-smooths
   ~2×; RCTD best broad/domain + well-calibrated but slow + weak fine;
   cell2location balanced + best-calibrated. Remaining: multi-seed CIs;
   CARD/SPOTlight/Tangram in a working toolchain.
10b. ✅ **Multi-seed spatial CIs + paired Wilcoxon** (5 seeds × 144-spot grids;
   all RCTD + cell2location runs executed): `run_spatial_multiseed.py`,
   `methods/run_rctd_multiseed.R`, `methods/run_cell2location_multiseed.py`,
   `score_spatial_multiseed.py` → `benchmarks/outputs/spatial_multiseed/`
   (`metrics_ci.tsv`, `pairwise_wilcoxon.tsv`, figK) + Addendum in
   `docs/SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md`. **Key correction:** the
   single-run "TissueResolve ≫ RCTD on fine" did NOT replicate — across 5 seeds
   RCTD (0.794) and TissueResolve_flat (0.770) are statistically tied (p=0.625),
   RCTD higher-mean/high-variance, TR more stable + ~130× faster but over-smooths
   ~2×; cell2location best-calibrated; hierarchical spatial reproducibly broken.
   n=5 small (Wilcoxon floor 0.0625). Remaining: more seeds; CARD/SPOTlight.
10. ⏳ (remaining external) Execute spatial methods in isolated envs (cell2location, RCTD,
    CARD, SPOTlight, Tangram) + NNLS-per-spot baseline + TissueResolve spatial
    modes. Same status/version/runtime/memory recording (**PART 3B**).
11. ⏳ Spatial accuracy + fidelity metrics on synthetic truth (**PART 12**).

### Phase 4 — Real spatial (no GT) + report redesign + rankings
0. ✅ **Real Visium no-GT evidence** — `benchmarks/spatial/run_real_visium_evidence.py`
   on the breast Visium section (500/3798 spots): marker-recovery proxy, Moran's I/
   Geary's C, entropy, near-zero, reconstruction Pearson, cross-method concordance
   (3 internal methods) → `benchmarks/outputs/spatial_real_visium/` + figL; write-up
   `docs/SPATIAL_REAL_VISIUM_EVIDENCE.md`. Explicitly NON-accuracy. Findings:
   TissueResolve_flat best marker-recovery proxy; NNLS extremely sparse (near-zero
   0.91) on real data; hierarchical barely concords (flagged).
0b. ✅ **Real-Visium cross-method concordance with established tools** — RCTD
   (1279 s) + cell2location (1309 s) executed on the same 500 spots
   (`methods/run_rctd_real.R`, `methods/run_cell2location_real.py`,
   `score_real_visium_external.py`) → `spatial_real_metrics_all.tsv`,
   `method_concordance_*_all.tsv`, figM. RCTD first failed on duplicate Visium
   gene symbols → fixed (dedup) + re-ran. Findings: TissueResolve_flat sits inside
   the established-tool concordance cluster (↔cell2location 0.54, ↔RCTD 0.45, vs
   RCTD↔cell2location 0.51) and ~ties cell2location on marker recovery; the
   **hierarchical spatial path is the clear outlier** (concords 0.05–0.15) —
   triangulated with its synthetic-GT failure as the component to fix.
12. ⏳ Real Visium evidence: marker recovery, concordance, Moran's I/Geary's C,
    smoothness, runtime — **labelled non-accuracy** (**PART 13**).
13. ⏳ Eight separate rankings (**PART 15**); decision-oriented report (**PART 16**)
    with all required sections + visualisations; `benchmark_report.html`.
14. ⏳ Tests + reproducibility (**PART 18**); deterministic seeds; every figure
    exports PNG/SVG(/PDF)/`.data.tsv`/caption/metric-definitions.

---

## Fair-input contract (PART 4) — enforced for every method
One canonical benchmark input layer + per-method adapters. Every method receives
the **same** reference cells, genes, label hierarchy, normalization intent, and
donor split. Recorded per method: cells, genes, hierarchy, broad/fine labels,
normalization, transform, library type, gene-id type, duplicate handling, min
cells/label, filtering, marker selection, batch correction, parameters. No
silent cell-type removal, relabeling, or giving TissueResolve richer labels than
competitors. Hierarchy-unaware methods: evaluate fine directly, aggregate to
broad for broad-level eval.

## Hard rules carried from the brief
- No superiority claims without evidence; no ranking of skipped/failed/exported.
- Bulk vs spatial never co-ranked.
- No accuracy metrics on real (no-GT) spatial.
- No algorithm changes before the audit (**done**); no tuning on final test sets.
- Preserve all raw method outputs + metadata. Do not commit generated outputs
  unless explicitly instructed.

## Required output inventory (PART 20) — tracking
`docs/VALIDATION_AND_BENCHMARK_PLAN.md` ✅ · `docs/BULK_PREDICTION_SPARSITY_AUDIT.md` ✅ ·
`docs/TCGA_DECONVOLUTION_COMPLEXITY_AUDIT.md` ✅ ·
`docs/VALIDATION_AND_BENCHMARK_FINAL_REPORT.md` ⏳ ·
`benchmarks/outputs/tcga_prediction_complexity_audit.tsv` ✅ ·
`bulk_{broad,fine,rare_population,complexity,spillover,robustness}_metrics.tsv` ⏳ ·
`spatial_{synthetic,real}_metrics.tsv` ⏳ · `runtime_memory_metrics.tsv` ⏳ ·
`method_status.tsv` 🟡 (exists per-run) · `method_rankings.tsv` ⏳ ·
`benchmarks/outputs/benchmark_report.html` 🟡 (exists, needs PART-16 redesign).
