# Stage 1 design — ReferenceSignatureOptimizer (multilevel signature)

**Scope of this session's vertical slice.** The donor-aware pseudobulk substrate, broad and
within-family candidate marker scoring, deterministic selection, a donor-disjoint benchmark,
and comparison to the current marker path already exist and are validated
(`reference/gene_selection.py`, `benchmarks/dev/gene_selection_benchmark.py`; donor_de beats
`GeneSelector` on breast+lung+cross-platform). The missing Stage-1 §3.2 pieces — a first-class
multilevel optimizer, the **rare-confirmation** level, **minimal-within-tolerance** selection,
the **manifest**, and the **fallback status** vocabulary — are what this slice adds.

## Public abstraction

```python
from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
opt = ReferenceSignatureOptimizer(celltype_col, donor_col, broad_col_or_mapping,
                                  batch_col=None, tolerance=0.05, seed=0)
model = opt.optimize(adata)          # -> SignatureModel
model.save(out_dir)                  # manifest TSVs + JSON
```

`opt.optimize` composes the existing primitives — it does **not** reimplement DE/selection:
- **Broad signature**: `select_donor_aware_genes` with cells relabelled to broad families
  (one-vs-rest across families). The validated `donor_de` panel.
- **Within-family fine signatures**: per broad family with ≥2 fine members,
  `donor_aware_de(mode="sibling")` among siblings.
- **Rare-confirmation signature**: per fine subtype, a small **specificity-first** panel
  (rank by within-family specificity × donor support, NOT raw LFC); size-capped.

## Minimal-within-tolerance

For the broad panel (and optionally each family), rank genes by score and evaluate a
**donor-held-out** deconvolution metric (NNLS on donor pseudobulks; leakage-safe: never the
final benchmark truth, never selection on eval donors) at increasing sizes; choose the
smallest size whose primary metric is within `tolerance` (default 5%, configurable) of the
best size. Records the full size→metric trace.

## Fallback status (reference quality)

From donor/cell counts + within-family separability, per family and overall:
`PASS` / `PASS_WITH_RESTRICTIONS` / `BROAD_ONLY` / `EXPERIMENTAL_FINE` / `REFERENCE_INADEQUATE`.
No donor column → status capped at `EXPERIMENTAL_FINE` with a prominent warning (donor
stability unvalidated). A family with no separable siblings → `BROAD_ONLY` (fine panel still
emitted but flagged diagnostic). Poor everywhere → `REFERENCE_INADEQUATE` but the best broad
signature is still produced.

## Outputs (manifest)

`broad_signature.tsv`, `fine_signature.tsv`, `rare_confirmation_signature.tsv`,
`signature_gene_scores.tsv`, `signature_selection_trace.tsv`, `signature_manifest.json`
(status, sizes, tolerance, donor/cell counts, provenance, seed). Each gene row carries the
attributes available (family/subtype, log2FC, support fraction, single-donor-penalised,
specificity, selected level, reason).

## Vertical-slice result (honest — NEGATIVE for the combined panel)

Donor-disjoint Poisson deconvolution on breast (3 seeds), optimizer panels vs the current
`GeneSelector` markers:

| gene_set | n_genes | RMSE | Pearson | cond_RMSE | rare_recall |
|---|---|---|---|---|---|
| current_markers | 500 | 0.036 | 0.630 | 0.212 | 0.477 |
| optimizer broad (family-level DE) | 194 | 0.042 | 0.571 | 0.248 | 0.373 |
| optimizer union (broad+fine+rare) | 529 | 0.039 | 0.617 | 0.220 | 0.396 |

**The optimizer's panels do NOT beat current markers on fine whole-composition deconvolution.**
The broad-family signature is the wrong panel for a fine task (no within-family markers); the
naive union is comparable-but-slightly-worse. This does **not** contradict the earlier win: the
validated `donor_de` panel (already integrated via `build_reference(gene_selection="donor_de")`)
is **fine-level one-vs-rest** DE (each fine type vs *all* others), which is a different — and for
this task, better — set than `broad-family DE ∪ within-family sibling DE`. Sibling markers don't
separate across families and family markers don't separate siblings, so their union is not
equivalent to fine-vs-all-rest.

**Consequence / correction:** the optimizer is a useful scaffold for the *broad-only* signature,
the *rare-confirmation* panel (Stage 2), and the reference-quality status/manifest — but its
combined panel is **not promoted** as a fine deconvolution replacement. The recommended
fine-deconvolution gene selection remains `donor_de` (fine one-vs-rest). A future revision should
expose a distinct fine-discriminative panel (= `donor_de`) alongside the broad/sibling levels,
rather than relying on the union.

## Non-goals / guarantees

- **No default change.** The optimizer produces panels + an optional `ref.selected_genes`;
  the production pipeline default (GeneSelector + wNNLS) is untouched.
- Deterministic (seeded); offline-testable on toy AnnData.
- Promotion to default requires the Stage-1 §3.8 criteria (improve/preserve broad, fine
  conditional RMSE, spillover, rare PR, absent-FPR on donor-disjoint) — donor_de already
  meets the broad/fine/spillover part; rare-PR/absent-FPR is Stage 2.
