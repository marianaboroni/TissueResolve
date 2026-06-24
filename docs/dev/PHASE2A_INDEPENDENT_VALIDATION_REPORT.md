# Phase 2A — Independent Validation Report (Stage 1)

Expanded, independent validation of partial confidence-weighted unresolved mass.
Harness: `benchmarks/diagnostics/phase2b_stage1_validation.py`. Soft gating stays
**experimental**. Nothing committed; default unchanged.

## Design (donor- AND seed-disjoint; rules 3,4)
- Calibration seeds {0,1} (confidence logistic fit), validation {2}, **TEST seeds
  {3..12} = 10 independent donor splits**.
- 6 scenarios: balanced, imbalanced, similar_subtypes, rare, missing_population,
  reduced_gene_overlap → **60 test replicates**.
- Modes: flat / hard_gate / ungated / soft_gate. Bootstrap 95% CIs over replicates.
- Reference cells restricted to mapped cell types (some donor subsets surface
  unmapped atlas types e.g. `mast cell`) — recorded.

## Data limitation (critical)
**No adequate second tissue is available.** The only non-breast dataset (scanpy
PBMC68k) has no donors, is not raw counts, and is 700 cells — unusable for
donor-disjoint count-level validation. The **2-tissue replication requirement is
UNMET**, so soft gating **cannot be promoted** regardless of metrics (rules 2,12).
This validation strengthens evidence *within one tissue* only.

## Results — mode means over 10 seeds × 6 scenarios (breast)

| mode | fine Pearson | fine RMSE | broad RMSE | false-resolution | false-abstention | eff-N error | rare sens. |
|---|---|---|---|---|---|---|---|
| flat | 0.715 | 0.038 | 0.034 | 0.094 | — | 3.21 | 0.511 |
| ungated | 0.616 | 0.054 | 0.091 | 0.089 | 0.000 | 4.61 | 0.499 |
| **soft_gate** | 0.615 | **0.048** | 0.091 | 0.091 | 0.150 | 4.63 | 0.486 |
| hard_gate (current) | 0.324 | 0.065 | 0.091 | 0.060 | 0.778 | 12.43 | 0.099 |

(Per-mode CIs in `phase2a_tissue_summary.tsv`; resolution metrics in
`phase2a_resolution_metrics.tsv`; raw per-replicate in
`phase2a_independent_validation.tsv`.)

## Prospective gates (declared in the plan, evaluated here) — **9/9 PASS**

vs **ungated**: broad RMSE ≤+2% ✅; fine RMSE ≤+2% ✅ (soft *better*: 0.048 vs
0.054); false-resolution ≤+10% relative ✅ (0.091 vs 0.089 — revised criterion,
non-inferior to ungated); eff-N error not worse ✅; mass conserved <1e-6 ✅.
vs **hard**: false-abstention materially lower ✅ (0.150 vs 0.778); rare
sensitivity materially higher ✅ (0.486 vs 0.099, ~5×); eff-N error substantially
lower ✅ (4.63 vs 12.43); fine RMSE lower ✅ (0.048 vs 0.065).

## Replication of the Phase-2A finding
The Phase-2A pattern reproduces at 10× the seeds and 2× the scenarios: soft gating
recovers the hard-gate collapse (fine Pearson 0.324→0.615 ≈ ungated), restores
rare detection and effective-N, cuts false-abstention ~80%, has the best
hierarchical fine RMSE, and conserves mass exactly. **Stable across all 10 seeds**
(CIs in the summary table exclude the null gap vs hard gating).

## Verdict
- **Prospective gates: PASS (9/9), single tissue.**
- **2-tissue replication: UNMET** → soft gating **remains experimental; NOT
  promoted**.
- Stop/go: prospective gates (the operational Stage-2 trigger) pass → proceed to
  **implement + experimentally validate Stage 2 (joint hierarchy)**, which also
  stays experimental and promotion-blocked. The original failed Phase-2A
  false-resolution gate is preserved (not rewritten); the revised soft-vs-ungated
  criterion was declared *before* these results (rule 5).

## Limitations
Single tissue; the `reduced_gene_overlap` scenario uses a 40% random gene subset
on an imbalanced base; confidence features remain weakly discriminative (Phase 2A)
so the gain is from continuous-vs-binary weighting; soft still trails flat on fine
Pearson (0.615 vs 0.715) — the residual broad-coupling/factorization gap that
Stage 2 targets.
