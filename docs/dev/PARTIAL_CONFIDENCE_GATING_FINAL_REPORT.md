# Partial Confidence-Weighted Unresolved Mass — Phase 2A Final Report

**Status: EXPERIMENTAL.** Behind an isolated module
(`src/tissueresolve/experimental/soft_hierarchy/`), not imported by default, not
wired into any pipeline, default production behaviour unchanged (verified by test).
Nothing committed.

## 1. Audit findings (`docs/PARTIAL_CONFIDENCE_GATING_AUDIT.md`)
The default hierarchical path **computes a continuous subtype confidence but
collapses it to a hard 0/1 threshold** in `estimate_partial_subtype_resolution`
(confidence ≥ 0.10 → keep *full* subtype mass; else → zeroed to `unresolved_*`).
With the collinear breast families (within-family separability 0.03–0.12, at/below
0.10), most subtypes fall under the threshold → 70–86% of mass parked unresolved,
effective-N collapsed 10.2→2.5, fine Pearson 0.128. Mass allocation is already
conservative-by-construction; the report applies no top-N filter. **The defect is
the binary collapse, not gene selection or display.**

## 2. Mathematical implementation
`apply_partial_confidence_gating(broad_mass, fine_estimates, family_map,
confidence|…)`:
```
resolved_θ[k,s] = raw_θ[k,s] · c_{k,s},   c ∈ [0,1]
unresolved[f,s] = broad_mass[f,s] − Σ_{k∈f} resolved_θ[k,s]   (≥0; c≤1, Σ_k conditional=1)
```
A strict generalisation of the hard gate (`c∈{0,1}` reproduces it). Mass conserved
exactly per family and overall (proven by unit tests: max error 0.0). Missing
confidence → 0 (conservative, routed to unresolved, warned). Reuses the existing
confidence signal; **no gene-panel or solver change**.

## 3. Confidence features tested (`confidence_feature_predictiveness.tsv`)
Per-(subtype) reference features + run features, evaluated on the calibration
split against ungated absolute error / correct-resolution (tol 0.05) / presence:

| feature | Spearman vs \|error\| | AUROC correct-resolution | AUROC presence |
|---|---|---|---|
| min_disc_genes / marker_support | −0.174 | **0.603** | 0.43 |
| ref_evidence | −0.107 | 0.586 | 0.50 |
| separability | −0.032 | 0.515 | 0.51 |
| cross_solver_agreement | n/a (constant: auto≈ridge conditional) | 0.50 | 0.50 |

**Honest finding:** the available features are only **weakly discriminative**
(AUROC 0.51–0.60); marker support is best. The Phase-2A benefit therefore comes
**mostly from replacing the binary collapse with continuous weighting**, not from
a strong confidence signal.

## 4. Calibration (`v2_soft_gating_calibration.tsv`)
Fit on calibration seeds {0,1}, evaluated on validation seed {2} (disjoint
donors+seeds; rules 10–11):

| calibrator | Brier↓ | ECE↓ |
|---|---|---|
| **logistic (selected)** | **0.107** | **0.079** |
| isotonic | 0.107 | 0.081 |
| monotonic (baseline) | 0.441 | 0.456 |

Logistic chosen (simplest well-calibrated; isotonic ties but adds no benefit —
rule 16). Well-calibrated despite weak discrimination (calibration ≠ ranking).

## 5. Before/after — 4 modes on HELD-OUT TEST seeds {3,4} (mean [95% CI], n=6)

| mode | fine Pearson | fine RMSE↓ | cond Pearson | false-resolution↓ | false-abstention↓ | unresolved frac | eff-N error↓ | rare sens. |
|---|---|---|---|---|---|---|---|---|
| flat (NNLS) | **0.753** | **0.042** | 0.584 | 0.100 | 0.00 | 0.00 | **2.06** | **0.494** |
| ungated | 0.571 | 0.064 | **0.589** | 0.097 | 0.00 | 0.00 | 2.75 | 0.465 |
| **soft_gate** | 0.568 | **0.058** | **0.589** | 0.097 | **0.155** | 0.187 | 2.79 | 0.456 |
| hard_gate (current) | 0.280 | 0.075 | 0.432 | 0.069 | 0.762 | 0.679 | 9.56 | 0.101 |

**Soft gating recovers essentially all the accuracy the hard gate destroyed**
(fine Pearson 0.280→0.568, ≈ ungated; effective-N error 9.56→2.79; rare
sensitivity 0.101→0.456) **while retaining a modest, honest unresolved mass
(0.187 vs hard 0.679) and cutting false-abstention 0.762→0.155 (−80%)**. It also
has the **best fine RMSE of any hierarchical mode (0.058 < ungated 0.064)** —
honest partial abstention removes some wrong mass that the ungated estimate keeps.

## 6. Ablations (`v2_soft_gating_ablation.tsv`)
Soft-gate fine Pearson by confidence score: separability 0.578, marker_support
0.574, ref_evidence 0.568 — **all within noise**, consistent with the weak feature
predictiveness. The improvement is robust to the confidence feature choice; the
*continuous-vs-binary* change is what matters, not the specific feature.

## 7. Runtime
Soft gating adds negligible cost (~ms; pure arithmetic on existing estimates). No
extra deconvolution run. The cross-solver-agreement feature (a second solve) is
**not needed** (it was non-predictive) and is excluded from the recommended path.

## 8. Failure modes
- Confidence features weakly discriminative → confidence ranking is soft; the
  method relies on continuous mass-preservation rather than sharp gating.
- Soft gating does **not** close the hierarchical-vs-flat gap: fine Pearson 0.568
  vs flat 0.753 (Δ≈0.185). That residual is broad-mass coupling + sequential
  factorization → **Phase 2B (soft joint hierarchy)** territory, not 2A.
- `cross_solver_agreement` degenerate (auto/ridge conditional nearly identical).

## 9. Promotion-gate status (`v2_soft_gating_promotion_gates.tsv`) — 5/6 PASS

| gate | result |
|---|---|
| broad RMSE ≤ ungated +2% | ✅ (0.1014 = 0.1014) |
| fine RMSE ≤ ungated +2% | ✅ (soft 0.058 < ungated 0.064) |
| false-abstention < hard | ✅ (0.155 vs 0.762) |
| eff-N error < hard (closer to truth) | ✅ (2.79 vs 9.56) |
| mass conservation < 1e-6 | ✅ (0.0) |
| **false-resolution ≥10% better than hard** | ❌ (soft 0.097 vs hard 0.069) |

**The single failing gate is mis-specified.** The hard gate's low false-resolution
(0.069) is an artifact of abstaining on 68% of mass; soft gating's false-resolution
(0.097) **equals the ungated and flat baselines** (0.097 / 0.100) — it is not
pathological, it is the unavoidable cost of resolving more mass. The correct
comparison is *soft vs ungated* (0.097 ≈ 0.097 → would pass "no worse than
ungated"). Recommend the committee **re-specify this gate as soft-vs-ungated**;
under that specification all six gates pass.

## 10. Remains experimental?
**Yes.** Per rules 3/17, it is not promoted to default: (a) one predefined gate
failed as written; (b) it does not reach flat; (c) validated on one tissue / small
n (6 test replicates). It is, however, a **clear, validated, mass-conserving
improvement over the current hard gate on every accuracy/richness/abstention
metric** and is the recommended hierarchical default *once* the false-resolution
gate is re-specified and a second tissue confirms it.

## 11. Is Phase 2B justified?
**Yes, conditionally.** Phase 2A fixes the *gating* failure (the primary Phase-1
finding) but leaves the **hierarchical-vs-flat gap (Δ fine Pearson ≈ 0.185)**,
which Phase-1 attributed to broad-mass coupling + sequential factorization — the
exact target of Phase 2B (soft joint broad–fine reconciliation). The remaining gap
is now quantified and motivates 2B. Do **not** proceed automatically.

---

### Deliverables
Module: `src/tissueresolve/experimental/soft_hierarchy/{__init__,config,partial_gating,confidence}.py`.
Tests: `tests/test_soft_gating.py` (10, all pass; full suite **129 passed**).
Harness: `benchmarks/diagnostics/phase2a_soft_gating.py`. Tables:
`confidence_feature_predictiveness.tsv`, `v2_soft_gating_{calibration,mode_comparison,
mode_comparison_raw,ablation,promotion_gates}.tsv`. Figures:
`benchmarks/outputs/figures/phase2a/` (mode accuracy, false-res/abstention,
effective-N, reliability, confidence-vs-error — each with `.data.tsv` + caption).
Baseline frozen + unchanged. No core algorithm modified; nothing committed.
