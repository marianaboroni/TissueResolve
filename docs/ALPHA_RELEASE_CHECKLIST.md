# TissueResolve — Alpha Release Checklist

Status of each item for the alpha (0.1.0). No GitHub release/tag is created here; no
commits are made.

| Item | Status | Notes |
|---|---|---|
| Tests passing | ✅ | Full suite green (see `docs/PERFORMANCE_BENCHMARK_REPORT.md` + Phase 8 of release prep) |
| README updated | ✅ | Recommended vs experimental separated; Poisson-GLM recommendation; limitations; alpha status |
| Tutorials | ✅ | `docs/tutorials/bulk_poisson_glm_quickstart.md`, `adaptive_resolution_and_reference_uncertainty.md`, `spatial_deconvolution_quickstart.md` |
| Performance/benchmark report | ✅ | `docs/PERFORMANCE_BENCHMARK_REPORT.md` (real numbers; deferred tools marked) |
| Figures/tables | ✅ | `docs/figures/*.tsv` (+ `*.png` where plotting available) |
| Outputs gitignored | ✅ | `benchmarks/outputs/**`, `*.h5ad`/`*.h5`/`*.RDS`, `r_lib`, `c2l_*`, logs — verified |
| No large files staged | ✅ | only small `docs/figures/*` (≤39 KB) are public artifacts |
| Package builds | ⚠️ verify | run `python -m build` if build tooling present (see Phase 8 result in final response) |
| CLI help reviewed | ✅ | `--bulk-solver` + `--spatial-preset` documented; experimental flags labelled |
| Default behaviour unchanged | ✅ | bulk default `wNNLS`; spatial default λ=0.1; no defaults changed |
| Recommended experimental features labelled | ✅ | Poisson GLM = recommended experimental; smoothing/state presets = experimental/benchmark-only |
| Limitations documented | ✅ | report §14 + README limitations + negative results §10 |
| Negative results preserved | ✅ | `docs/figures/negative_results_summary.tsv` + report §10 (kept on purpose) |
| Future-work roadmap | ✅ | `docs/FUTURE_WORK_DISTRIBUTION_AWARE_REFERENCE.md` (P2a/b, design only) |

## Pre-publication review notes

- Public entry points: `README.md`, `docs/PERFORMANCE_BENCHMARK_REPORT.md`,
  `docs/FEATURE_STATUS.md`, `docs/POISSON_GLM_RECOMMENDED_SOLVER.md`,
  `docs/tutorials/*`, `docs/ALPHA_RELEASE_NOTES.md`, this file.
- The ~70 internal `docs/*.md` reports are development provenance, not user docs.
- No claim of: best / state-of-the-art / all-tissue validation / cell fractions by
  default / real-Visium accuracy / superiority over deferred tools.
- Deferred tools (BayesPrism, DWLS, SCDC, SPOTlight, Tangram, DestVI) are reported as
  deferred with reasons, never fabricated.

## Not done (by instruction)

- No commit, no tag, no GitHub release.
- No new methods. No default changes.
