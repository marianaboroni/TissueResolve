# TissueResolve feature status

Honest classification of each feature. **stable** = validated + default-safe;
**experimental** = usable behind an explicit flag, not validated across datasets;
**partial** = implemented but not fully wired/measured; **planned** = designed,
not implemented; **not implemented** = absent (and must not be claimed).

| Feature | Status | Notes |
|---|---|---|
| Bulk deconvolution (wNNLS / pipeline) | **stable** | default path; broad test coverage |
| Spatial deconvolution (NB-CAR) | **stable** | default path |
| Solver `auto` (gene-masking CV) | **stable** | matches NNLS on real bulk |
| Hierarchical broad→fine | **stable** | default `resolution-mode auto/hierarchical`; unresolved mass at family level |
| Resolution Decision Layer (trusted resolution per family) | **stable** | `src/tissueresolve/resolution.py`; classifies each family `broad_only`/`selected_fine`/`full_fine` from full-panel reliability + signature/query evidence (NOT cell-level AUROC) before fine predictions are interpreted; recorded in run metadata + QC table + report (table shown before fine predictions); backward-compatible (numeric estimates unchanged; soft gating still final) |
| Modality-aware pipeline metadata | **stable** | bulk vs spatial recorded explicitly (prediction unit sample/spot, gene-weighting mode, spatial smoothing used, H&E availability, fine trusted vs diagnostic); `compute_modality_aware_gene_weights` shared compatibility wrapper |
| Within-family soft gating (partial confidence-weighted) | **stable / default** | `--hierarchical-gating soft` (DEFAULT); mass-conserving; validated on breast + lung benchmarks |
| Within-family hard gating (binary threshold) | **legacy** | `--hierarchical-gating hard`; kept for reproducibility; over-abstains in collinear families; not default |
| Within-family ungated (diagnostic) | **diagnostic** | `--hierarchical-gating ungated`; no abstention; not calibrated |
| FineGranularityRefiner (contrast-weighted / residual / spillover-calibrated fine refinement) | **experimental — negative result; NOT integrated** | `experimental/soft_hierarchy/fine_refiner.py`; **failed promotion gates on breast + lung** — lowered conditional RMSE in some settings only by increasing spillover and false-positive subtype detection; not used by default, not wired into any pipeline; soft gating remains the final hierarchical layer. See `docs/FINE_GRANULARITY_REFINER_REPORT.md` |
| Granular signatures (broad / cell-type / state panels) | **experimental** | `reference/granular_signatures.py`, tested; built from per-cell AnnData at build time |
| State-aware hierarchy (broad→cell type→state) | **experimental** | `bulk/state_aware_hierarchical.py`; behind `--state-aware` / `deconv_bulk(state_aware=True)`; **not default**; falls back to 2-level when no state labels; **not yet validated on real data** |
| Spatial weak-smoothing preset (`--spatial-preset weak_smoothing`, λ=0.02) | **experimental — opt-in; default unchanged** | `experimental/spatial_presets.py`. In synthetic breast + lung benchmarks it reduced oversmoothing and improved broad/fine correlation while preserving local RMSE, but rare-niche behavior and effective-N calibration did not improve consistently across tissues, so it is not promoted and the default `lambda_spatial=0.1` is unchanged. See `docs/SPATIAL_WEAK_SMOOTHING_BENCHMARK_REPORT.md` |
| Reference suitability score | **stable** | PASS/CAUTION/WARNING/FAIL with worst-component override |
| Report generation (QC-first, embedded figures) | **stable** | unified report; figure manifest |
| Spatial benchmark (concordance/structure, no-truth) | **partial** | bulk/spatial split done; multi-metric ranking module added; marker-recovery/synthetic-truth/dashboard not yet wired into the run |
| External-tool benchmark (MuSiC/CARD/cell2location/…) | **partial** | several tools executed; others skipped/failed honestly |
| Reference adaptation (residual-driven gene-weight update) | **not implemented** | deferred; do not enable/claim |
| Cell-type-specific expression reconstruction | **not implemented** | deferred; do not enable/claim |
| Hyperparameter tuning | **not implemented** | `tuning/` is a quarantined stub that raises `NotImplementedError`; deferred to v0.2 |
| Full BayesPrism-like Bayesian model (Gibbs) | **not implemented** | deferred |

## Required wordings (use verbatim in README/CLI/report)

- **State-aware:** "Experimental: state-aware deconvolution runs behind
  `--state-aware`. It is not part of the default v0.1 workflow and has not been
  validated across real datasets."
- **State-aware report integration:** "State-aware outputs are experimental; full
  standard-report integration is limited."
- **Expression reconstruction:** "Cell-type-specific expression reconstruction is
  planned/deferred and not implemented in v0.1." Use "cell-type-level
  deconvolution estimates" / "RNA-derived composition estimates" instead.
- **External benchmark:** "Only executed or imported tools are ranked.
  Exported-only or skipped tools are not benchmarked."
- **Spatial benchmark:** "The real Visium benchmark does not measure accuracy
  without ground truth. It evaluates concordance, spatial structure, marker
  agreement, runtime, and output completeness."
- **Estimate type:** "Bulk estimates are RNA-derived proportions, not absolute
  cell fractions. Spatial estimates are spot-level RNA-derived composition, not
  single-cell labels."
- **Publication status:** "TissueResolve is currently alpha / early-access
  research software. It is not yet a fully publication-ready method until
  external benchmarks, multi-dataset validation, and final API stabilization are
  complete."

Deferred features (**reference adaptation**, **expression reconstruction**,
**hyperparameter tuning**, **full Bayesian model**) are intentionally absent and
must not be presented as available in the README, CLI help, or report.
