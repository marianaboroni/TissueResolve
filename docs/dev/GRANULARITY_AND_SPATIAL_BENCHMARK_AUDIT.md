# Granularity & spatial-benchmark audit

Read-only audit (Part 1), written before any new code change. Reflects the state
**after** commits `43c9eda` (within-family markers) and `b657011`
(bulk/spatial benchmark split).

## Preflight
- `git status`: clean except this session's new docs; two prior-session audit
  docs intentionally untracked.
- Tests: the affected subset (report + plotting + hierarchical + within-family +
  examples + external benchmark) is green for all changed code; full suite was
  928 passed at the branch baseline.

## Answers

1. **Broad vs fine signatures separate?** *Partially.* Fine-level within-family
   panels now exist (`reference/within_family_markers.py`,
   `build_family_specific_gene_panels`). Broad estimation still uses the
   family-aggregated reference with the **global** panel — there is **no
   dedicated broad-signature panel builder** yet (`build_broad_signature_panel`
   does not exist).
2. **HVGs recalculated within each family?** **Yes** —
   `select_within_family_hvgs` (committed `43c9eda`).
3. **DE / pairwise markers within each family?** **Yes** —
   `select_within_family_de_genes`,
   `select_pairwise_discriminative_genes_within_family`.
4. **Within-family genes used by the hierarchical solver?** **Capable but not
   yet wired end-to-end.** `evaluate_within_family_resolvability`,
   `compute_within_family_subtype_confidence`, `assemble_hierarchical_estimates`
   and `run_hierarchical_bulk` accept an optional `family_gene_panels`; when
   passed, resolvability/confidence/fine-split use family genes. **But the report
   pipeline (`07_generate_reports.py`) does not yet BUILD panels from the
   reference AnnData and pass them**, so the real-data run still uses global
   genes. Spatial (`spatial/hierarchical.py`) does **not** yet accept panels.
5. **Only diagnostics?** No longer only diagnostics (the solver path consumes
   panels), but in the *current report run* the panels are not supplied, so in
   practice the fine split is still global until wired.
6. **Unrelated-family genes excluded in the fine step?** **Yes when panels are
   passed** (`subset_genes(panel)` restricts to family genes). Not enforced when
   panels are absent.
7. **Gene weights separate for broad vs fine?** Fine: yes
   (`within_family_gene_weights`, 0–1 per family). Broad: no dedicated
   broad-level weights (uses the global selector).
8. **Families with enough marker support?** Measured by
   `within_family_resolution_summary` / `marker_support_by_family`
   (`marker_support_status` PASS/CAUTION/WARNING/FAIL), but **not yet computed on
   the real reference** in the report run.
9. **Families unresolved due to low separability?** Currently 6/8 (T/NK,
   B/Plasma, Myeloid, Epithelial, Endothelial, Mural) — judged on global genes;
   re-evaluation with within-family panels is implemented but **not yet run on
   real data**.
10. **Benchmark compares global vs within-family marker strategy?** **No** — no
    `granularity_signature_benchmark` exists yet.
11. **Spatial benchmark distinguishes accuracy / concordance / marker recovery /
    structure / runtime?** *Partially.* `b657011` split bulk vs spatial and added
    spatial status / structure (Moran's I) / runtime / completeness with explicit
    no-ground-truth captions. **Missing:** marker-recovery metric,
    method-concordance heatmap (needs ≥2 executed spatial methods), a synthetic
    spatial accuracy benchmark, and a transparent multi-metric ranking.
12. **Report visually separates bulk vs spatial benchmark?** **Yes** — sections 9
    (bulk, accuracy) / 10 (spatial, concordance/structure) / 11 (scorecard),
    committed in `b657011`.

## Gap list → what this task must add
- **Granular signatures (Parts 2–3):** `build_broad_signature_panel`,
  `write_granular_signature_outputs` (the `outputs/signatures/*` tables +
  `family_solver_gene_panels.tsv`), and **wire panels end-to-end** into the
  real-data report build (build panels from AnnData → pass to
  `run_hierarchical_bulk`), plus spatial solver acceptance.
- **Before/after benchmark (Part 4):** global vs multi-granularity metrics.
- **Report subsection (Part 5):** "Multi-granularity signature construction".
- **Spatial multi-metric system (Parts 6–9):** marker recovery, JSD,
  concordance, structure, runtime/completeness; `spatial_multimetric_ranking.py`
  (real no-truth ranking, synthetic accuracy ranking, overall scorecard);
  bubble-heatmap dashboard; restructured spatial benchmark section; optional
  lightweight synthetic spatial benchmark.

## Constraints
No package rewrite; no reference-free methods; no prediction change before audit;
no hidden warnings; no spatial-accuracy overclaim; no mixed bulk/spatial ranking;
skipped/exported tools not ranked; nothing committed without instruction.
