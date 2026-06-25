# TissueResolve — Alpha Release Notes

Status: **alpha** (`pyproject.toml`: version 0.1.0, Development Status :: 3 - Alpha).
This document records the public-release audit and the method taxonomy used to keep
the public surface honest. No defaults were changed for the alpha; no new methods
were added.

## Phase 0 — Repository audit

- **README:** present (422 lines), already states RNA-derived proportions and the
  hierarchical/soft-gating defaults. Rewritten for alpha (recommended vs experimental
  separation, Poisson-GLM recommendation, limitations, alpha status).
- **docs/:** ~70 internal reports/plans + audits (development history). These remain
  for provenance but are **not** the public entry points; the public-facing set is:
  `README.md`, `docs/PERFORMANCE_BENCHMARK_REPORT.md`, `docs/FEATURE_STATUS.md`,
  `docs/POISSON_GLM_RECOMMENDED_SOLVER.md`, `docs/tutorials/*`, this file,
  `docs/ALPHA_RELEASE_CHECKLIST.md`, `docs/FUTURE_WORK_DISTRIBUTION_AWARE_REFERENCE.md`.
- **CLI:** `tissueresolve run` exposes `--bulk-solver` (wNNLS default + Poisson/NB GLM
  experimental) and `--spatial-preset` (default + experimental smoothing/state presets).
  Experimental presets are retained but labelled experimental/benchmark-only in docs.
- **`.gitignore`:** already excludes `benchmarks/outputs/`, `*.h5ad`, `*.h5`,
  `benchmarks/envs/r_lib/`, `c2l_*`, example outputs. Verified: all generated
  benchmark artifacts, raw predictions, matrices, logs, envs, `.h5ad`/`.RDS` are
  ignored. Only small public summary TSVs/PNGs under `docs/figures/` are intended for
  version control.
- **tests:** full suite passes (see Phase 8); experimental modules are covered.
- **package metadata:** `pyproject.toml` correct (Python ≥3.9, alpha classifier).

### Large outputs that must remain gitignored
`benchmarks/outputs/**` (bulk/spatial metrics, external preds, exports),
`benchmarks/envs/r_lib`, `c2l_*`, all `*.h5ad`/`*.h5`/`*.RDS`, logs. Confirmed ignored.

## Phase 1 — Public method taxonomy

**Recommended / public-facing**

| Modality | Feature | Status |
|---|---|---|
| Bulk | `wNNLS` solver | **default** |
| Bulk | Poisson GLM solver (`--bulk-solver poisson_glm_experimental`) | **recommended experimental** |
| Bulk | bootstrap uncertainty | stable |
| Bulk/Spatial | RNA-derived proportions; `unresolved_<family>` mass; hierarchical broad→fine + soft gating | stable/default |
| Bulk/Spatial | adaptive resolution report; reference uncertainty diagnostics | stable (reporting) |
| Spatial | default NB-CAR solver | **default** |
| Spatial | `weak_smoothing`, `edge_aware_smoothing` | experimental |

**Diagnostic / interpretability-only (not estimation modes)**

state-similarity graph; adaptive grouping; state-ambiguity report; rare-detection
calibration. These inform interpretation; they do not change the recommended estimates.

**Benchmark-only / NOT recommended** (retained for reproducibility, demoted in docs)

`combined_weak_edge_smoothing`; post-fit state regularization; in-solver state
regularization; sparsity-state regularization; `state_regularized_solver_*` variants.
All tested and **negative** for improving conditional within-family recovery.

**Deferred external tools (attempted, not available here)**

BayesPrism, DWLS, SCDC (bulk) — install failed (dependency/compilation). SPOTlight,
Tangram, DestVI (spatial) — not installed. Reported as deferred, never fabricated.

## Alpha positioning

TissueResolve is a transparent, resolution-aware deconvolution framework for bulk and
spatial transcriptomics that reports **RNA-derived proportions** by default (not
absolute cell fractions), quantifies reference separability, exposes unresolved mass,
and avoids overconfident fine-state calls when the reference is not identifiable. The
strongest validated algorithmic improvement is the **opt-in Poisson GLM bulk solver**;
the default bulk solver remains `wNNLS` pending broader validation.

## Negative results are kept on purpose

The negative results (regularization cannot recover non-identifiable collinear fine
states; simple smoothing cannot fully close spatial oversmoothing) are a **feature** of
the project: they document where the tool correctly declines to promise impossible
fine resolution. They are preserved in `docs/PERFORMANCE_BENCHMARK_REPORT.md` §10 and
`docs/figures/negative_results_summary.tsv`.
