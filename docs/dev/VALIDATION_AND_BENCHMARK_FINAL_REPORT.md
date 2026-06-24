# TissueResolve — Validation & Benchmark Final Report (PART 16)

Decision-oriented synthesis of the bulk + spatial validation. **Bulk and spatial are reported and ranked separately.** Real-data spatial is evidence only (not accuracy). Only executed methods are ranked. Full HTML with figures: `benchmarks/outputs/benchmark_report.html`.

## Executive summary

- The TCGA **"4–5 populations" is a solver effect** (unregularised NNLS on a near-collinear reference), not post-processing/threshold/display; established methods vary similarly (BisqueRNA is *more* concentrated). See `BULK_PREDICTION_SPARSITY_AUDIT.md`.

- **Bulk (25 held-out-donor replicates):** TissueResolve `nnls`/`auto` significantly beat MuSiC (p=2e-4) and BisqueRNA (p<1e-4) on fine + broad accuracy, and are >100× faster.

- **Spatial synthetic (5 seeds):** RCTD ≈ TissueResolve_flat on fine accuracy (tie, p=0.62); RCTD leads broad/domain; cell2location best-calibrated; TissueResolve most stable + ~130× faster but over-smooths ~2×.

- **Spatial real (no GT):** TissueResolve_flat sits inside the established-tool concordance cluster; ties cell2location on marker-recovery proxy.

- **One clear defect:** the hierarchical path under-performs on bulk + spatial ground truth and is the concordance outlier on real data — fix, don't trust.


## Bulk performance — mean [95% CI], 25 split×scenario replicates

| method | fine Pearson | broad Pearson | complexity err↓ |
| --- | --- | --- | --- |
| flat_nnls | 0.747 [0.702,0.791] | 0.959 [0.946,0.970] | 3.247 [2.199,4.312] |
| flat_auto | 0.736 [0.692,0.781] | 0.949 [0.932,0.963] | 3.661 [2.503,4.975] |
| flat_weighted_nnls | 0.741 [0.692,0.788] | 0.946 [0.927,0.962] | 3.392 [2.349,4.444] |
| flat_ridge_nnls | 0.678 [0.644,0.714] | 0.887 [0.859,0.915] | 7.864 [6.143,9.597] |
| MuSiC | 0.614 [0.511,0.709] | 0.818 [0.782,0.857] | 1.944 [1.514,2.431] |
| BisqueRNA | 0.377 [0.313,0.440] | 0.602 [0.558,0.642] | 5.985 [4.177,7.819] |
| hierarchical | 0.324 [0.251,0.393] | -0.026 [-0.084,0.039] | 12.545 [10.416,14.696] |
| hierarchical_resolvable | 0.894 [0.857,0.924] | 0.901 [0.864,0.930] | 0.486 [0.339,0.657] |


## Spatial synthetic — mean [95% CI], 5 seeds

| method | fine Pearson | broad Pearson | oversmoothing | domain ARI |
| --- | --- | --- | --- | --- |
| RCTD | 0.794 [0.656,0.873] | 0.979 [0.975,0.982] | 0.779 [0.734,0.832] | 0.311 [0.294,0.334] |
| TissueResolve_spatial_flat | 0.770 [0.715,0.819] | 0.884 [0.872,0.896] | 1.894 [1.712,2.082] | 0.217 [0.167,0.267] |
| NNLS_per_spot | 0.718 [0.665,0.770] | 0.924 [0.915,0.934] | 0.374 [0.302,0.447] | 0.253 [0.209,0.295] |
| cell2location | 0.650 [0.603,0.712] | 0.785 [0.759,0.808] | 0.986 [0.891,1.106] | 0.278 [0.265,0.292] |
| TissueResolve_spatial_hierarchical | 0.066 [-0.010,0.134] | -0.338 [-0.438,-0.261] | 36.497 [2.675,103.635] | 0.137 [0.025,0.302] |


## Best method by scenario (bulk & spatial separate)

| Scenario | Best / leading method | Supporting metric | Caveat |
| --- | --- | --- | --- |
| Bulk — broad cell types | flat_nnls (TissueResolve) | broad Pearson 0.96 [0.95,0.97] | ties auto/weighted; > MuSiC 0.82, Bisque 0.60 |
| Bulk — fine subpopulations | flat_nnls / auto (TissueResolve) | fine Pearson 0.75 [0.70,0.79] | sig > MuSiC (p=2e-4) & Bisque (p<1e-4); tie vs weighted/auto |
| Bulk — rare population (0.7%) | ridge / nnls (TissueResolve) | detection 92–100% | all flat solvers over-estimate; hierarchical misses it (0%) |
| Bulk — complexity preservation | flat_nnls (TissueResolve) | eff-N err 3.25 | best on realistic mixtures; ridge over-disperses; MuSiC also good (1.94 on single split) |
| Bulk — fastest | flat solvers (TissueResolve) | <2 s/run | MuSiC ~25–46s, Bisque ~13–36s |
| Spatial — fine accuracy | RCTD ≈ TissueResolve_flat | 0.79 vs 0.77 (tie, p=0.62) | RCTD higher-mean/high-variance; TR more stable + ~130× faster |
| Spatial — broad / domain | RCTD | broad 0.98, domain ARI 0.31 | RCTD leads structure |
| Spatial — calibration | cell2location | oversmoothing ~1.0 | TR over-smooths ~2× |
| Spatial — fastest | TissueResolve_spatial | ~6 s | RCTD/cell2location ~13–35 min CPU |
| Spatial real (no GT) — marker proxy | cell2location ≈ TissueResolve_flat | 0.29 vs 0.29 | concordance, not accuracy; hierarchical is the outlier |


## Not performed (honest gaps)

- Dedicated spillover matrix / discrimination AUROC (Part 9) — partial coverage via within-family separability audit + `similar_subtypes` scenario.

- Parametric robustness sweeps (Part 11: depth, gene dropout, cell count) — partial coverage via cross-donor/missing/extra/rare scenarios.

- Peak-memory instrumentation; cross-tissue generalisation; CARD/BayesPrism/CIBERSORTx/SPOTlight/Tangram (toolchain/license).


## Limitations & readiness

Single tissue; small spatial n (5 seeds); cell2location CPU-only with reduced epochs; real Visium has no ground truth. **Strong enough for an alpha release and a methods-grade preprint of the breast-cancer benchmark**; broader claims (general superiority) need ≥1 more tissue, spillover + robustness analyses, and more seeds. See `VALIDATION_AND_BENCHMARK_PLAN.md` for the phase tracker.


*Source tables: `benchmarks/outputs/`. Companion docs: `BULK_PREDICTION_SPARSITY_AUDIT.md`, `TCGA_DECONVOLUTION_COMPLEXITY_AUDIT.md`, `HOLDOUT_BULK_BENCHMARK_RESULTS.md`, `SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md`, `SPATIAL_REAL_VISIUM_EVIDENCE.md`. No core algorithm modified; nothing committed.*
