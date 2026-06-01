# Within-family marker selection — implementation report

Adds a two-level signature strategy: keep the existing **global** selection for
broad families, and add **within-family** HVG/DE/pairwise selection for the
fine-level step. Additive and backward-compatible (defaults unchanged).

## Audit result (Part 1)
`docs/WITHIN_FAMILY_MARKER_SELECTION_AUDIT.md`. Key findings: gene selection was
**global only**; HVGs were **never** computed; pairwise markers were
diagnostics-only and never used by the solver; the fine-level step reused the
**global** panel; within-family resolvability was judged on the **global** gene
set (the root cause of 6/8 families being unresolved); raw per-cell AnnData is
available **only at build time**, so within-family selection must happen there
and the panels must be threaded into the solver.

## Implemented (complete + tested, green)

### Parts 2–3 — `src/tissueresolve/reference/within_family_markers.py`
- `select_within_family_hvgs` — HVGs from one family's cells (within-family
  log-normalised variance ranking; detection + min-cells gating).
- `select_within_family_de_genes` — one-vs-rest subtype markers within a family.
- `select_pairwise_discriminative_genes_within_family` — pair-vs-pair
  discriminative genes with optional cross-sample stability.
- `build_family_specific_gene_panels` — per-family panels = HVG ∪ pairwise,
  minus mito/ribo/stress and minus genes undetected in the query; full Part-3
  scoring (variance, specificity (tau), pairwise-discriminatory, ref/query
  detection, sample-stability, stress penalty) → `final_within_family_weight`
  normalised 0–1 per family; returns panels + scores + excluded + marker-support
  tables; `.save()` writes all Part-3 TSVs.
- `within_family_resolution_summary` — global vs within-family separability /
  resolvability per family (Part 5), flags `separability_improved` /
  `verdict_changed`. **Improvement is reported, never assumed.**

### Part 4 — solver integration (additive, backward-compatible)
- `hierarchy.evaluate_within_family_resolvability` and
  `compute_within_family_subtype_confidence` gained an optional
  `family_gene_panels` arg; when given, separability/confidence are computed on
  each family's panel (restricting the reference via `subset_genes`). New
  `n_panel_genes` column records the panel size used.
- `hierarchy.assemble_hierarchical_estimates` and
  `bulk.hierarchical.run_hierarchical_bulk` accept and thread
  `family_gene_panels`; metadata records `within_family_panels`.
- Default (`None`) reproduces the original global-gene behaviour **exactly**
  (verified by `test_default_behaviour_unchanged_without_panels`).

### Part 8 — tests (all green)
`tests/shared/test_within_family_markers.py` (6) — per-family HVGs differ and
recover within-family subtype markers (not global broad markers); pairwise/DE
tables have family+pair columns; panels differ per family; mito/ribo/stress and
query-undetected genes excluded with reasons; weights ∈ [0,1]; `.save()` writes
the TSVs. `tests/shared/test_hierarchical.py` (+5) — family panel actually
changes separability (flat panel → unresolvable, discriminating panel →
resolvable); default unchanged without panels; `run_hierarchical_bulk` threads
panels and preserves mass; `within_family_resolution_summary` schema.
**Result: 37 passed (6 + 31 existing hierarchical + 5 new integration), 0
failures.** Existing global marker selection and all prior tests untouched.

## Deferred to the report/benchmark restructuring task (next)
The remaining integration layer is intentionally folded into the immediately
following "report summary + benchmark split" task to avoid wiring that would be
restructured anyway:
- **Part 6 (report subsection + before/after figures):** building family panels
  from the real-data AnnData at report time and rendering a "Within-family
  marker selection" subsection with before/after separability.
- **Part 7 (old-vs-new hierarchical benchmark):** the benchmark section is being
  split into separate bulk/spatial sections in the next task; the
  global-vs-within-family hierarchical comparison will be added there.
- **Part 9 (real-data regen):** run after the report wiring lands.

## Output files (Part 3) written by `FamilyPanels.save()`
`within_family_marker_scores.tsv`, `within_family_selected_genes.tsv`,
`within_family_pairwise_markers.tsv`, `within_family_gene_weights.tsv`,
`within_family_excluded_genes.tsv`, `marker_support_by_family.tsv`,
`marker_support_by_pair.tsv`; plus Part-5 `global_*` /
`within_family_pairwise_separability` / `within_family_resolution_summary`
frames from `within_family_resolution_summary`.

## Effect (to be quantified on real data in the regen step)
The mechanism is verified on synthetic data (panel restriction flips
resolvability). On the real breast-cancer reference, whether any of the 6
currently-unresolved families (T/NK, B/Plasma, Myeloid, Epithelial, Endothelial,
Mural) become resolvable will be **measured and reported** by
`within_family_resolution_summary` — not asserted in advance.

## Constraints honored
No core-algorithm change (NNLS/CAR untouched); no prediction-value change;
global selection kept; defaults unchanged; warnings not hidden; no fine-subtype
accuracy overclaim; nothing committed.
