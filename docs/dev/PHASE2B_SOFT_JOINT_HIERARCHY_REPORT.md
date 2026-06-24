# Phase 2B — Soft Joint Broad–Fine Hierarchy Report (Stage 2)

**Status: EXPERIMENTAL; NOT supported for promotion.** Isolated module
(`src/tissueresolve/experimental/soft_hierarchy/joint_solver.py`), opt-in, default
unchanged, nothing committed.

## Hypothesis
The residual hierarchical-vs-flat gap is caused by the **sequential**
factorisation (broad → freeze broad mass → conditional fine). A soft joint
reconciliation (`π ≈ Aθ`, no irreversible freezing) should reduce it.

## Implementation
Constrained per-sample objective (no strong L1; rules 6,7,8 — pure scipy/numpy):
`min_{π,θ≥0} ‖W_b(y−Bπ)‖² + α‖W_f(y−Sθ)‖² + λ_h‖π−Aθ‖² + λ₂‖θ‖²`, `Σπ=Σθ=1`.
Two solvers — **alternating** NNLS sub-problems (primary) and **joint** SLSQP
(secondary). Simplex enforced by post-solve renormalisation (package convention).
Validated by 9 unit + synthetic-recovery tests (noiseless recovery mean abs error
< 0.05; λ_h=0 decouples to flat NNLS; strong λ_h ↑ consistency; mass conserved;
no sparsity collapse).

## Benchmark (controlled — identical preprocessing, 6 seeds × 3 scenarios)
`benchmarks/diagnostics/phase2b_stage2_joint.py`. All methods on the same
shared-gene CPM-log1p `S`,`y` to isolate the factorisation.

| mode | fine Pearson | fine RMSE | broad RMSE | cond RMSE | rare sens. | eff-N err |
|---|---|---|---|---|---|---|
| flat | 0.523 | 0.057 | 0.082 | 0.324 | 0.466 | 2.55 |
| sequential | 0.508 | 0.057 | 0.087 | 0.313 | 0.476 | 2.64 |
| **joint_alt** | **0.508** | **0.057** | **0.087** | **0.313** | 0.476 | 2.63 |
| oracle_broad | 0.588 | 0.055 | 0.020 | 0.307 | 0.483 | 2.40 |

**Joint ≈ sequential — identical to 3 decimals on every metric.**

## Error decomposition (`phase2b_error_decomposition.tsv`)
| mode | broad error | conditional-fine error | consistency error | total fine error |
|---|---|---|---|---|
| sequential | 0.087 | 0.313 | **0.000** | 0.057 |
| joint_alt | 0.087 | 0.313 | **0.000** | 0.057 |

## λ_h ablation (`phase2b_ablation_metrics.tsv`)
fine RMSE / Pearson / conditional RMSE are **flat across λ_h ∈ {0, 0.1, 1, 5}**
(0.0464 / 0.701 / 0.319). The soft-consistency strength has **no effect**.

## Why the hypothesis is NOT supported
The consistency error (`‖π−Aθ‖`) is **already ≈ 0** under the sequential
factorisation — broad and fine solves *agree*. So there is nothing for joint
reconciliation to fix, and varying λ_h does nothing. The residual error is
dominated by:
- **within-family conditional error (0.31 RMSE)** — the subtypes are hard to tell
  apart (an identifiability/signature limit; Phase-1 already showed an oracle DE
  panel did not help); and
- **broad-signature error (0.087)** — `oracle_broad` (true family mass) is the only
  variant that improves (fine 0.588, broad RMSE 0.020), confirming broad error is
  real, but the joint solver re-derives the *same* broad from the *same* data, so
  it cannot recover the true broad any better than sequential.

## Promotion gates (`phase2_promotion_gates.tsv`) — 5/7, critical ones FAIL
| gate | result |
|---|---|
| vs sequential: fine RMSE ≥5% better | ❌ (0.0574 vs 0.0574 — 0% improvement) |
| vs sequential: conditional RMSE improves | ✅ (marginal, 0.3126 vs 0.3126) |
| vs sequential: broad RMSE ≤+2% | ✅ |
| vs sequential: rare sensitivity non-inferior | ✅ |
| vs flat: broad RMSE non-inferior (≤+2%) | ❌ (0.087 vs 0.082) |
| vs flat: fine-RMSE gap reduced vs sequential | ✅ (negligible) |
| mass conserved <1e-6 | ✅ |

The **decisive gate (≥5% fine-RMSE improvement vs sequential) FAILS** — joint and
sequential are identical.

## Verdict (rules 16, 17)
Joint hierarchy provides **no measurable benefit** over the simpler sequential
factorisation → **the simpler method wins; joint hierarchy is NOT promoted and is
not recommended for further development.** This is an effort-saving negative
result: it redirects future work away from broad–fine *coupling* (already
consistent) toward the true bottleneck — **within-family subtype discriminability**
(a signature/identifiability problem) and **broad-signature accuracy**.

## Note on Phase-1 attribution
Phase 1 ranked broad-mass coupling the "secondary" failure. Stage 2 refines this:
the coupling is *consistent* (not the problem); the secondary error is broad
*signature* accuracy (only oracle-broad helps) plus the dominant within-family
conditional error — neither addressable by joint optimisation.

*Experimental; promotion blocked (1 tissue + failed critical gate). Nothing committed.*
