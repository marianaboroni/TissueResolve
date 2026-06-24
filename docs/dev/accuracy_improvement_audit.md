# Accuracy improvement audit

Read-only diagnosis of why TissueResolve currently trails `NNLS_baseline` on
fine-level pseudobulk accuracy, with **measured** numbers, before any code
changes. All numbers are from the existing real breast-cancer pseudobulk
(12 samples, 5000 shared genes) against known mRNA-proportion ground truth.

## 1. Why does NNLS_baseline outperform TissueResolve on fine-level accuracy?

Controlled experiment (plain per-sample NNLS unless noted):

| configuration | n_genes | fine Pearson |
|---|---|---|
| **NNLS, all shared genes** | 5000 | **0.768** |
| NNLS on TissueResolve's selected panel (no weights) | 500 | 0.705 |
| **TissueResolve_flat** (panel **+** wNNLS weights) | 500 | 0.675 |
| NNLS top-2000 by specificity | 2000 | 0.620 |
| NNLS top-500 by specificity | 500 | 0.403 |
| NNLS top-100 by specificity | 100 | 0.108 |

**Decomposition of the gap (0.768 → 0.675):**

- **Gene-panel restriction (~0.06 loss):** restricting NNLS to TissueResolve's
  500-gene panel drops 0.768 → 0.705. On clean pseudobulk, informative signal
  is spread across thousands of genes; a linear model benefits from using them.
  (TissueResolve's panel is far better than a naive top-500 — 0.705 vs 0.403 —
  so the selection is *good*, just smaller than optimal for clean data.)
- **Weighting (~0.03 loss):** wNNLS weighting on the panel drops 0.705 → 0.675.
  Protocol/specificity/stability weighting trades a little clean-data accuracy
  for robustness to protocol mismatch and spillover — beneficial on real noisy
  bulk, mildly harmful on clean pseudobulk.

**Conclusion:** the gap is a *design trade-off*, not a bug. TissueResolve's
marker selection + weighting are tuned for robustness; a full-gene plain NNLS is
near-optimal for clean linear mixtures. The fix is to let TissueResolve *adopt*
the strong backbone when cross-validation says it is best — i.e. `solver=auto`
selected by gene-masking CV — rather than competing against it externally.

## 2. Is TissueResolve using the same genes as NNLS_baseline?

No. `NNLS_baseline` uses **all 5000 shared genes**; TissueResolve selects a
**500-gene** protocol-aware marker panel (`GeneConfig.n_genes=500`). This is the
dominant driver of the fine-level gap.

## 3. Is TissueResolve applying weighting/filtering that reduces fine-level accuracy?

Yes — composite gene weights (specificity × stability × protocol × concordance)
feed a weighted NNLS. Measured cost on clean pseudobulk ≈ 0.03 Pearson. It is
protective on real noisy/protocol-mismatched data, so it should remain available
(and default for robustness), but `auto` should be able to drop it when CV shows
plain NNLS reconstructs better.

## 4. Is hierarchical mode all-or-nothing for unresolved families?

**Yes (current limitation).** `evaluate_within_family_resolvability` →
`decide_unresolved_families` marks an entire family resolvable or not; if not,
*all* of its mass goes to `unresolved_<family>` and *every* subtype is zeroed.
A family that is partially resolvable (one clear subtype, the rest ambiguous)
cannot assign confident mass to the clear subtype. Part 3 (partial resolution)
addresses this.

## 5. Are family-level metrics computed correctly?

Yes. `metrics.family_level_metrics` aggregates both truth and predictions to
families (mapping fine subtypes and `unresolved_<family>` columns to families)
and compares — mass-preserving. Real-data family-level Pearson: hierarchical
0.749 vs flat 0.661 vs NNLS 0.844.

## 6. Are fine-level metrics computed only on resolvable subtypes (for hierarchical)?

`metrics.fine_metrics_resolvable_only` exists and is reported, but the headline
`accuracy_metrics` still includes the zeroed abstained columns, which is why the
hierarchical "fine Pearson" looks very low (0.062). The benchmark already
reports the resolvable-only and family-level views; reporting should foreground
them for hierarchical methods.

## 7. Does auto mode use hierarchical broad→fine when broad/fine labels exist?

Yes (resolution-mode `auto`). Separately, there is currently **no solver `auto`**
— the solver backbone (panel + weighted NNLS) is fixed. Part 2 adds `solver=auto`.

## 8. Are broad/fine labels used consistently in bulk and spatial?

Yes — both call `build_cell_type_hierarchy` + `run_hierarchical_*`, which share
`aggregate_reference_by_family` / `assemble_hierarchical_estimates`.

## 9. Does the benchmark compare old vs improved TissueResolve?

Not yet. The benchmark compares flat vs hierarchical vs baselines, but has no
"improved" (solver-auto / partial-resolution / ensemble) TissueResolve modes.
Part 10 adds a before/after accuracy-improvement report.

## 10. Which parameters most likely affect accuracy?

In rough order of measured impact on clean pseudobulk:

1. **`GeneConfig.n_genes`** (panel size) — largest effect; 500 → all genes
   recovers most of the gap.
2. **Gene weighting** (`weight_*`) — ~0.03; disabling helps clean data, hurts
   robustness.
3. **Solver backbone** (NNLS vs ridge vs weighted) — interacts with conditioning.
4. **Hierarchical all-or-nothing threshold** — affects how much fine mass is
   abstained (family-level vs fine-level trade-off).
5. **`lambda_spatial`** (spatial only) — smoothing vs sharpness.

## Proposed changes (implemented in subsequent parts)

- **Solver backbone selection (`solver=auto`)** driven by **gene-masking CV**:
  evaluate plain NNLS (all genes), weighted NNLS (panel), marker NNLS, ridge,
  and an ensemble; pick the best by a multi-objective score (reconstruction +
  conditioning + spillover penalty), not reconstruction alone. This lets
  TissueResolve match/beat the NNLS baseline while keeping the resolution-aware
  layer on top.
- **Partial hierarchical resolution**: assign confident subtype mass to
  resolvable subtypes within a family and keep only the ambiguous remainder as
  `unresolved_<family>`.
- **Within-family pairwise markers** for subtype estimation.
- **Reference suitability score**, **ensemble**, and **cell-type expression
  reconstruction** as additional, clearly-labelled modules.

These preserve the reference-based, resolution-aware design: TissueResolve
estimates the **most reliable level of resolution** supported by the reference
and data, and now also adopts the strongest backbone CV selects.
