# State-aware deconvolution audit

Read-only audit (Part 1), written before any code change. Reflects the codebase
after commits `43c9eda` (within-family markers) and `b657011` (benchmark split),
plus the uncommitted spatial ranking module.

## Answers

1. **More than two annotation levels?** **No.** The hierarchy is strictly
   two-level: a `{fine_cell_type → broad_family}` mapping
   (`reference/hierarchy.py:build_cell_type_hierarchy`,
   `validate_hierarchy`). `HierarchicalConfig` exposes only
   `broad_cell_type_col` and `fine_cell_type_col` (`config.py:317-318`);
   `io/validation.py` detects/validates only those two
   (`detect_broad_cell_type_col`, `detect_fine_cell_type_col`,
   `validate_hierarchical_annotations`). There is **no `state_col`**.

2. **Broad family / cell type / cell state treated separately?** **No.** Only
   "broad family" and "fine" exist; "fine" conflates cell type and cell state.

3. **Genes selected separately per level?** **Partially, two levels at most.**
   Global selection (`markers.py`) is single-panel; within-family HVG/DE/pairwise
   panels now exist (`within_family_markers.py`, committed `43c9eda`) and the
   solver can consume them via `family_gene_panels`, but (a) they are not yet
   built/passed in the real-data run, and (b) there is **no third (state-level)
   panel** and **no dedicated broad-only panel builder**.

4. **Flat fine estimation, or aggregate states → types → families?**
   **Two-level only.** `bulk/hierarchical.py:run_hierarchical_bulk` runs a
   family-level deconvolution and a flat fine-level deconvolution, then computes
   conditional within-family proportions and aggregates fine→family. There is
   **no state→cell-type→broad** three-level aggregation.

5. **Unresolved mass at the correct level?** **One level only** —
   `unresolved_<broad_family>` when within-family subtypes are not separable.
   There is no `unresolved_<cell_type>` (state-within-type) level.

6. **Reconstruct cell-type-specific expression?** **No.** The reference stores
   aggregated profiles (`R_cpm`, `R_log`); nothing reconstructs per-sample /
   per-spot cell-type expression or residuals (no `reconstruct_*expression`,
   no residual/QC outputs).

7. **Missing vs a BayesPrism-inspired state-aware strategy:**
   - a third **cell-state** annotation level (state → cell type → broad family);
   - **granularity-specific gene panels** at all three levels, **wired into the
     solver** (broad panel; cell-type-within-family panels; state-within-cell-type
     panels);
   - **state-aware hierarchical allocation** with unresolved mass assigned at the
     correct level (broad vs cell type);
   - a **conservative reference-adaptation** layer (residual-driven gene-weight
     updates, not full Bayesian reference updating);
   - **inferred cell-type-specific expression reconstruction** + reconstruction QC
     (masked-gene score);
   - an **old-vs-new state-aware benchmark** and report section.

## Existing building blocks to reuse (not rewrite)
- `within_family_markers.py`: HVG/DE/pairwise selection + scoring + panels —
  generalise from "within family" to "within cell type" for the state level.
- `hierarchy.py`: resolvability gating, conditional proportions, partial
  resolution, unresolved-mass arithmetic — extend to a third level.
- `bulk/hierarchical.py` `family_gene_panels` plumbing — extend to per-cell-type
  state panels.
- `spatial_multimetric_ranking.py` (uncommitted) — for benchmark scoring.

## Constraints
No package rewrite; no full Bayesian Gibbs sampling; keep NNLS/auto solvers;
never hide unresolved mass; no fine-state accuracy overclaim; no prediction
change before audit; nothing committed.

## Planned increments (honest sequencing)
1. **Part 2 — three-level hierarchy support** (this increment): validated
   `state → cell_type → broad` model with backward-compatible 2-level fallback,
   plus the `outputs/hierarchy/*` tables. Self-contained + tested.
2. Parts 3–4 — state-level gene panels + state-aware solver (extends the
   committed within-family layer to a third level).
3. Parts 5–6 — conservative reference adaptation + expression reconstruction.
4. Parts 7–8 — old-vs-new benchmark + report section.
Each lands tested; no real-data accuracy improvement is claimed until measured.
