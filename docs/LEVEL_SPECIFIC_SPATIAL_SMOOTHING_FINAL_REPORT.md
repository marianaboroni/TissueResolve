# Level-Specific Spatial Smoothing — Final Report (Task D)

**Status: EXPERIMENTAL; NOT promoted.** Isolated module
`src/tissueresolve/experimental/spatial_smoothing/` (opt-in; default unchanged; no
NB-CAR VI; no PyTorch/JAX). Nothing committed.

## Hypothesis
Broad families tolerate stronger spatial smoothing than fine subtypes, so
**λ_broad > λ_fine** should reduce fine over-smoothing/edge blurring without
hurting broad accuracy.

## Implementation
`fit_level_specific_spatial()` composes the EXISTING NB-CAR solver twice — family
reference at λ_broad, fine reference at λ_fine — then combines mass-consistently
(`combined_k = broad_{f(k)} · within-family-share`). 5 unit tests (mass
conservation, separate λ, both-zero, boundary/edge-blurring, rare-niche).

## Benchmark — 10 seeds × 6 modes (structured synthetic: sharp border + gradient + rare niche)
`benchmarks/diagnostics/phase_taskD_level_spatial.py`.

| mode (λ_broad, λ_fine) | fine RMSE | broad RMSE | edge blur | boundary F1 | oversmoothing | rare-niche | local RMSE |
|---|---|---|---|---|---|---|---|
| none (0, 0) | 0.051 | 0.087 | 0.002 | 0.411 | 0.36 | 0.456 | 0.051 |
| fine_only (0, 0.1) | 0.041 | 0.087 | 0.000 | 0.411 | 1.03 | 0.756 | 0.041 |
| broad_only (0.1, 0) | 0.047 | 0.077 | 0.003 | 0.307 | 0.47 | 0.456 | 0.046 |
| **weak_equal (0.02, 0.02)** | **0.037** | **0.072** | 0.001 | 0.360 | 1.50 | **0.800** | **0.036** |
| default_equal (0.1, 0.1) | 0.038 | 0.077 | 0.002 | 0.307 | 1.80 | 0.744 | 0.037 |
| level_specific (0.1, 0.02) | 0.038 | 0.077 | 0.002 | 0.307 | 1.47 | 0.756 | 0.037 |

## Promotion gates (level_specific vs current default) — 5/6, critical gate FAILS

| gate | level_specific | default | pass |
|---|---|---|---|
| **fine RMSE improves vs default** | 0.0380 | 0.0376 | ❌ (Wilcoxon p=0.002, *worse* direction) |
| broad RMSE preserved (≤+2%) | 0.0774 | 0.0774 | ✅ |
| edge blurring reduced | 0.0019 | 0.0020 | ✅ |
| boundary F1 preserved | 0.307 | 0.307 | ✅ |
| rare-niche not reduced | 0.756 | 0.744 | ✅ |
| oversmoothing closer to 1.0 | 1.47 | 1.80 | ✅ |

## Verdict (rules 16, 17)
**The level-specific hypothesis is NOT supported.** Decoupling (strong broad / weak
fine) gives essentially the same fine RMSE as the default (0.038 vs 0.038) and is
**dominated by the simpler `weak_equal` (λ=0.02 everywhere)** on fine RMSE, broad
RMSE, local RMSE, and rare-niche sensitivity. Level-specific does reduce
over-smoothing (1.47 vs 1.80) but so does — more — weak_equal, without the extra
two-solve complexity. → **level-specific smoothing is not promoted and not
recommended for further development.**

## The actionable finding (confirms Phase-1)
The real lever is the **simpler change: lower the global λ from the default 0.1 to
~0.02.** `weak_equal` is the best or tied-best mode on every accuracy/structure
metric here and on the Phase-1 global-λ sweep (λ=0.02 fine Pearson 0.793 vs 0.775).
The effect is real but **modest** (fine RMSE 0.037 vs 0.038; the larger gains are
broad RMSE 0.072 vs 0.077, rare-niche 0.800 vs 0.744, over-smoothing 1.50 vs 1.80).

## Recommendation
- **Do not adopt level-specific smoothing** (no benefit over a single weaker λ).
- **Candidate default change:** lower `lambda_spatial` 0.1 → ~0.02. But per "do not
  change stable defaults before validation" + no validated second tissue, **hold**:
  validate the weaker λ on a second tissue (and, ideally, image-based pseudo-spots
  with sharp ground-truth boundaries) before changing the default.
- Both spatial smoothing changes remain **experimental**.

## Limitations
Single tissue; synthetic spatial only (no real-Visium boundary ground truth);
boundary F1 uses dominant-broad-family domain agreement (a proxy); the composed
two-solve level-specific path differs slightly from the production hierarchical
spatial path (used here to isolate the λ effect). Outputs:
`level_specific_spatial_smoothing.tsv`, `spatial_boundary_metrics.tsv`,
`taskD_promotion_gates.tsv`, `figures/taskD/figT_level_specific`.
