# Changelog

All notable changes to TissueResolve are documented here. This project is
**alpha / early-access research software** (pre-release v0.1); APIs and defaults
may change.

## [0.1.0] — alpha / early-access (unreleased)

First public alpha. Resolution-aware deconvolution for bulk RNA-seq and 10x
Visium spatial transcriptomics from a shared single-cell/-nucleus reference.

### Core (default-safe)
- Reference-based **bulk** and **spatial** deconvolution; outputs are
  RNA-derived (mRNA) composition, never reported as absolute cell counts.
- **Hierarchical broad→fine** mode by default (`--resolution-mode auto`) with
  `unresolved_<family>` mass where subtypes are not separable.
- **Partial confidence-weighted soft gating** (`--hierarchical-gating soft`),
  validated on two tissues (breast + lung benchmarks).
- Reference suitability, separability and spillover diagnostics; QC-first
  publication HTML report with source data saved for every figure.
- Honest internal + optional external-tool benchmarks (bulk and spatial kept
  separate); paired bootstrap CIs for the gold-truth bulk benchmark.

### Experimental (opt-in, behind explicit flags; not default)
- **Spatial weak-smoothing preset** (`--spatial-preset weak_smoothing`,
  λ_spatial=0.02). In synthetic breast + lung benchmarks it reduced
  oversmoothing and improved broad/fine correlation while preserving local RMSE,
  but rare-niche behaviour and effective-N calibration did not improve
  consistently across tissues — so it remains opt-in and the **default
  λ_spatial=0.1 is unchanged**. See
  `docs/SPATIAL_WEAK_SMOOTHING_BENCHMARK_REPORT.md`.
- State-aware three-level hierarchy (`--state-aware`), granular signatures,
  external-tool benchmark runners, spatial multi-metric ranking.

### Fixed
- **Preset/bootstrap bug:** the `quick` and `standard` presets declare
  `bootstrap: False` but the config default (`n_bootstrap=200`) was never
  cleared, so bootstrap CIs ran against the preset's intent. Presets that
  disable bootstrap now set `n_bootstrap=0`. Regression test added
  (`tests/cli/test_presets.py`).

### Not implemented / deferred (do not assume available)
- Reference adaptation, cell-type-specific expression reconstruction,
  hyperparameter tuning, and a full Bayesian (BayesPrism-like) model.

### Known limitations
- Not yet a publication-ready method: external benchmarks and multi-dataset
  validation are ongoing; the API is not stabilized.
- Spatial accuracy is demonstrated on **synthetic** ground truth only; real
  Visium has no ground truth (concordance/structure are reported, not accuracy).
