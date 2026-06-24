# TissueResolve v2 — Phase 2 Consolidated Final Report

Covers Stage 0 (verify + plan), Stage 1 (independent soft-gating validation),
Stage 2 (soft joint hierarchy). **Stage 3 (level-specific spatial smoothing) is
DEFERRED** (see end). Everything experimental; default unchanged; nothing committed.

## Stage 0 — baselines verified
`git diff src/tissueresolve/ d15a3da..HEAD` empty (core unchanged since freeze).
Phase 2A opt-in; production default = original stable path. Recorded the
originally-failed Phase-2A false-resolution gate and the revised *prospective*
criterion (soft non-inferior to **ungated**) declared before Stage-1 results
(rule 5). Plan: `PHASE2B_IMPLEMENTATION_AND_VALIDATION_PLAN.md`.

**Critical data limitation:** no adequate **second tissue** (only breast atlas;
PBMC68k has no donors, not raw counts, 700 cells). The **2-tissue replication
requirement is UNMET** → **no method can be promoted this phase**, regardless of
metrics.

## Stage 1 — independent soft-gating validation: **PASS (9/9 prospective gates)**
10 test seeds × 6 scenarios (donor- + seed-disjoint), breast. Soft gating
replicates the Phase-2A pattern at scale:

| mode | fine Pearson | fine RMSE | false-abstention | eff-N error | rare sens. |
|---|---|---|---|---|---|
| flat | 0.715 | 0.038 | — | 3.21 | 0.51 |
| ungated | 0.616 | 0.054 | 0.000 | 4.61 | 0.50 |
| **soft_gate** | 0.615 | **0.048** | 0.150 | 4.63 | 0.49 |
| hard_gate (current) | 0.324 | 0.065 | 0.778 | 12.43 | 0.10 |

vs ungated: broad/fine RMSE non-inferior (fine *better*), false-resolution
non-inferior, mass conserved. vs hard: false-abstention 0.150 vs 0.778, rare
0.49 vs 0.10 (5×), eff-N error 4.6 vs 12.4. **Soft gating is a validated, robust,
mass-conserving fix for the hard-gating collapse** — but stays experimental
(2-tissue replication unmet).

## Stage 2 — soft joint hierarchy: **NOT supported (negative result)**
Hypothesis (gap = sequential factorisation) **rejected**. Controlled comparison
(identical preprocessing, 6 seeds × 3 scenarios):
joint_alt **≡ sequential** (fine Pearson 0.508 = 0.508; RMSE 0.057 = 0.057).
Error decomposition: **consistency error ≈ 0** for both → broad and fine already
agree; nothing to reconcile. λ_h ablation: metrics flat across {0,0.1,1,5}.
Residual error = within-family conditional (0.31) + broad-signature (0.087);
only `oracle_broad` (true family mass) improves (fine 0.588). Critical promotion
gate (≥5% fine-RMSE improvement vs sequential) **FAILS** (0% improvement).
→ simpler sequential wins (rules 16,17); joint hierarchy not promoted, not
recommended for further development.

## Consolidated conclusions
1. **Soft gating (Phase 2A): the validated win.** It fixes the primary Phase-1
   failure (hard gating) — recovers accuracy, richness, rare detection; conserves
   mass; passes 9/9 prospective gates on 10 seeds × 6 scenarios. Blocked from
   promotion only by the missing second tissue.
2. **Joint hierarchy (Phase 2B): a clean negative.** The broad–fine *coupling* is
   already consistent; the real bottleneck is **within-family subtype
   discriminability** (identifiability/signature limit — an oracle DE panel did
   not help in Phase 1) and **broad-signature accuracy**. Future effort should
   target these, not joint optimisation.
3. Both remain **experimental**; the production default is unchanged.

## Recommended next actions
- **Obtain a genuine second tissue** (donor-labelled, raw-count single-cell atlas)
  to satisfy the replication requirement; if soft gating replicates, promote it as
  the hierarchical default (re-specify the false-resolution gate vs ungated).
- **Re-target the residual gap** at within-family discriminability (contrastive
  signatures / evidence modules) and broad-signature accuracy — NOT joint
  optimisation.
- **Execute Stage 3** (level-specific spatial smoothing) — deferred below.

## Stage 3 — level-specific spatial smoothing: DEFERRED (not executed this turn)
Independent of Stage 2. Phase-1 already established the direction (global λ sweep:
λ=0.02 beats the default λ=0.1 on accuracy *and* over-smoothing). A true
*level-specific* (λ_broad, λ_fine) experiment requires new experimental code in
the spatial hierarchical path (two λ at the two resolution levels) plus a
10-seed × 6×6-grid sweep — a full additional work unit. Design is specified in
`PHASE2B_IMPLEMENTATION_AND_VALIDATION_PLAN.md` §"Stage 3". **Not executed here to
avoid a rushed/fragile implementation;** `level_specific_spatial_smoothing.tsv`
and `LEVEL_SPECIFIC_SPATIAL_SMOOTHING_REPORT.md` are pending.

## Deliverables produced this turn
Docs: `PHASE2B_IMPLEMENTATION_AND_VALIDATION_PLAN.md`,
`PHASE2A_INDEPENDENT_VALIDATION_REPORT.md`, `PHASE2B_SOFT_JOINT_HIERARCHY_REPORT.md`,
this report. Code: `experimental/soft_hierarchy/joint_solver.py` (+ exports).
Tests: `tests/test_joint_hierarchy.py` (9) — full new Phase-2 suite 19 pass;
`benchmarks/tests/` 119 pass. Tables: `phase2a_independent_validation.tsv`,
`phase2a_tissue_summary.tsv`, `phase2a_resolution_metrics.tsv`,
`phase2a_stage1_gates.tsv`, `phase2b_joint_hierarchy_metrics.tsv`,
`phase2b_error_decomposition.tsv`, `phase2b_ablation_metrics.tsv`,
`phase2_promotion_gates.tsv`. Baseline preserved; nothing committed.
**Pending:** `level_specific_spatial_smoothing.tsv`,
`LEVEL_SPECIFIC_SPATIAL_SMOOTHING_REPORT.md` (Stage 3).
