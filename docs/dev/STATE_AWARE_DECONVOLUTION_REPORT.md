# State-aware hierarchical deconvolution — implementation report

Status report for the BayesPrism-inspired state-aware task. Honest about
implemented-and-tested vs documented next increment. No prediction values
changed; nothing committed.

## Audit (Part 1)
`docs/STATE_AWARE_DECONVOLUTION_AUDIT.md`. Confirmed TissueResolve is **two-level
only** (fine→broad), conflating cell type and cell state; gene panels are global
(within-family panels exist but not state-level); the solver aggregates
fine→family only; unresolved mass is single-level (`unresolved_<broad_family>`);
there is **no expression reconstruction**. The within-family layer (`43c9eda`)
and the resolvability/conditional/unresolved arithmetic in `hierarchy.py` are the
building blocks to extend to a third level.

## Implemented this increment (tested, green)
**Part 2 — three-level hierarchy support**:
`src/tissueresolve/reference/three_level_hierarchy.py` +
`tests/shared/test_three_level_hierarchy.py` (8 tests pass).
- `ThreeLevelHierarchy` (state→cell_type→broad) with upward aggregation helpers
  (`state_to_broad`, `states_of_celltype`, `celltypes_of_broad`, `frame`).
- `build_three_level_hierarchy` (structural validation — a state's cell type
  must map to a broad family, else `ValueError`; nothing silently inferred).
- `build_three_level_from_two_level` — **backward-compatible** lift of a
  `{fine→broad}` mapping (cell_type = fine, `has_states=False`).
- `build_three_level_from_obs(obs, broad_col, cell_type_col, state_col=None)` —
  derive from annotation columns; rejects inconsistent broad-per-celltype.
- `validate_three_level_hierarchy` — per-state table flagging low-support states
  (`< min_cells_per_state`) and single-donor cell types; **emits warnings, drops
  nothing**.
- `summarize_three_level_hierarchy` + `write_three_level_outputs` →
  `outputs/hierarchy/hierarchy_three_level.tsv`,
  `hierarchy_validation.tsv`, `hierarchy_summary.tsv`.

This is the foundation the state-aware solver and state-level gene panels build
on. It maps directly onto the task's Part 2 contract (three columns, validation,
the three output tables, two-level fallback).

## Answers to the final-report questions (current truth)
- **Three-level hierarchy supported?** **Yes (data model + validation + outputs),
  this increment.** Not yet consumed by the solver.
- **State-aware gene panels generated?** **No yet** — Part 3 next (the
  within-family machinery generalises from "within family" to "within cell
  type", but is not yet written/run).
- **State-aware solver uses them?** **No yet** — Part 4 next.
- **Effect on broad / cell-type / state accuracy and unresolved mass?** **Not
  measured** — the solver is unchanged this increment, so **no accuracy or
  unresolved-mass change is claimed** (no false precision).
- **Expression reconstruction implemented?** **No yet** — Part 6 next.

## Documented next increment (sequenced, each lands tested)
1. **Part 3 — state-level gene panels:** `build_broad_gene_panel`,
   `build_celltype_within_family_gene_panels`,
   `build_state_within_celltype_gene_panels` (generalise
   `within_family_markers` to within-cell-type), granularity weights, pairwise
   state marker support; `outputs/signatures/*`.
2. **Part 4 — state-aware solver:** broad → cell-type-within-family →
   state-within-celltype, with unresolved mass at the correct level
   (`unresolved_<broad>` / `unresolved_<cell_type>`); gate states on marker
   support / separability / spillover / overlap / stability.
3. **Part 5 — conservative reference adaptation:** residual-driven gene-weight
   updates (NOT full reference rewrite), CLI
   `--reference-adaptation none|gene_weight_update|experimental_reference_update`,
   before/after metrics.
4. **Part 6 — inferred cell-type expression reconstruction** + masked-gene QC
   (marker genes by default to bound output size).
5. **Parts 7–8 — old-vs-new state-aware benchmark + report section.**

## Also pending from prior tasks (uncommitted)
`benchmarks/shared/spatial_multimetric_ranking.py` (+12 tests) and the
granularity/spatial audit docs from the previous task; the within-family +
benchmark-split commits (`43c9eda`, `b657011`) are local but **the push was
blocked by the environment's permission policy**.

## Safe to commit?
The three-level hierarchy module + tests and the spatial ranking module + tests
are self-contained and green — safe to commit. The state-aware solver itself is
**not** implemented, so the granularity issue is **not** yet fixed on real data;
that requires the next increment.
