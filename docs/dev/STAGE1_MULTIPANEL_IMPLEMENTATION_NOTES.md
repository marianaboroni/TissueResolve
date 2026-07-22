# Stage 1 multi-panel — pre-implementation inspection notes

Answers to the §1 questions, verified against the code (not the README).

1. **How `current_markers` are generated.** `reference/markers.py::GeneSelector.select(ref)`
   on the aggregated **mean-profile** `ReferenceSignature` (not per-cell): a composite
   log-additive score (specificity × stability × protocol-safety × bulk-concordance), two-phase
   (per-type + pairwise fill). Budget = `GeneConfig.n_genes` (**default 500**, total-panel cap).
   Stability uses `ref.donor_cv` but selection is on mean profiles, **not donor pseudobulk DE**.
2. **How `donor_de` does fine one-vs-rest.** `reference/gene_selection.py::select_donor_aware_genes`
   = union over fine types of `donor_aware_de(mode="one_vs_rest", top_n_per_type=15)`, i.e.
   **each fine type vs ALL other fine types**. This is the "fine-global" contrast.
3. **Donor pseudobulks vs pooled cells.** `donor_de` builds **donor × cell-type pseudobulks**
   (`donor_pseudobulk`) and runs a Welch test on donor pseudobulks — **donor-aware, not pooled**.
   Cells are never treated as independent replicates. ✓
4. **Selection strictly inside training donors?** The optimizer/`donor_de` select on whatever
   AnnData is passed. The benchmark passes **train (reference) donors only**; evaluation uses
   held-out **query donors**. So donor-disjoint and leakage-free **as invoked**. (Risk: a caller
   could pass all cells — documented; the optimizer warns when no donor column.)
5. **Gene-budget enforcement.** `GeneSelector` caps at `n_genes` (500). `donor_de` budget is the
   **union size** (top_n_per_type × #types, deduped ≈ 296 breast). The optimizer's broad/sibling
   use their own top_n. **Budgets are NOT harmonised across strategies.**
6. **Bulk vs spatial ranking.** Both bulk (`BulkPipeline`) and spatial (`SpatialPipeline`) select
   via the same `GeneSelector(config=self.cfg.genes)`. `donor_de` / `ref.selected_genes` is a
   bulk-path feature currently. **No bulk/spatial ranking difference in the default path.**
7. **Evaluation-donor leakage.** None in the benchmark: train (ref) donors and query donors are
   disjoint; selection sees only train donors.
8. **Equal gene budget?** **NO.** The Stage-1 table compared 500 (current) vs 194 (optimizer
   broad) vs 529 (union) vs ~296 (donor_de). **A fair comparison must equalise the budget** —
   part of the optimizer-broad deficit is simply fewer genes. This stage's benchmark controls it.

## Consequences for this stage

- **Reuse, do not duplicate:** `donor_aware_de` / `select_donor_aware_genes` already implement
  donor-aware fine one-vs-rest correctly. The **fine-global panel wraps them** (no new DE).
- **Equal-budget benchmark** is mandatory to fairly re-test the panels.
- The panels have **distinct statistical contrasts** and must not be interchanged:
  broad = family-vs-family; fine-global = fine-vs-all-other-fine; sibling = fine-vs-siblings;
  rare-confirmation = specificity-first per subtype.

## Equal-budget benchmark result (donor-disjoint, Poisson, 3 seeds)

Fair test (all panels capped to the SAME 300 genes) that Stage 1 lacked:

| tissue | strategy | RMSE | Pearson | cond_RMSE | rare_fpr | rare_recall |
|---|---|---|---|---|---|---|
| breast | current | 0.037 | 0.612 | 0.220 | 0.281 | 0.441 |
| breast | **fine_global** | **0.036** | **0.651** | **0.190** | **0.240** | **0.499** |
| lung | **current** | **0.032** | **0.562** | 0.258 | 0.179 | 0.371 |
| lung | fine_global | 0.043 | 0.370 | 0.276 | 0.197 | 0.378 |

**Finding (honest, mixed):** at a FIXED 300-gene budget, `fine_global` wins **every** metric on
breast but **loses** on lung. Cause: lung has ~2× the fine types, so a fixed cap starves the
per-type panel. At its **natural (~550) budget** `donor_de`/`fine_global` beat current on lung
too (prior `gene_selection_benchmark`: 0.031/0.558 vs 0.032/0.528). ⇒ **the gene budget must
scale with the number of fine types**, not be a fixed constant. `broad_only` and `union` are
never the best fine panel (confirming the multi-panel role separation).

**Promotion decision:** the **multi-panel architecture is validated** (fine_global is the correct
primary fine panel; breast equal-budget clean win; distinct-role panels + misuse safeguards).
But `fine_global` is **NOT promoted to default**: at a fixed budget it is tissue-dependent, its
optimal budget scales with #types, and its real-bulk behaviour is untested. Recommended next
config to benchmark cleanly: `fine_global` with a **per-type budget** (top_n_per_type × n_types),
breast+lung+cross-platform, 5 seeds — the setting under which prior evidence already shows a
two-tissue win.

## Per-type-budget confirmation (the corrected promotion benchmark)

`fine_global` at its NATURAL per-type budget = `select_donor_aware_genes(top_n_per_type=15)` =
the `donor_de_ovr` panel already benchmarked (5 seeds, breast+lung, cross-platform,
`benchmarks/results/gene_selection/`). Result vs current markers (500):

- breast/complete: RMSE 0.037→0.037, Pearson 0.630→**0.646**, cond-RMSE 0.214→**0.196**,
  rare_fpr 0.258→**0.224**, rare_recall 0.474→**0.491** (win on all).
- breast/cross_platform: RMSE 0.050→**0.046**, Pearson 0.560→**0.617**, cond-RMSE 0.283→**0.276** (win).
- lung/complete: RMSE 0.032→**0.031**, Pearson 0.528→**0.558**, cond-RMSE 0.243→**0.228** (win);
  **rare_fpr 0.177→0.214 (WORSE)**, rare_recall 0.401→0.392 (~).

**Promotion decision (revised from the fixed-budget stage):** at the correct per-type budget,
`fine_global` beats current markers on **RMSE / Pearson / conditional within-family RMSE across
both tissues + cross-platform (5 seeds)** — it is a **`experimental_candidate_pending_budget_and_rare_validation` for fine
deconvolution**. It is **not flipped to default yet** because of (1) a modest **lung rare-FPR
increase** (0.214 vs 0.177) — §3.8 forbids promotion that materially worsens rare FPR — and
(2) untested **real-bulk** behaviour. Both are addressable: the **rare_confirmation panel**
(Stage 2) is the designed mechanism to control rare FPR, and `fine_global` is already usable
opt-in via `build_reference(gene_selection="donor_de")` (identical panel). Recommended next:
Stage 2 (calibrated rare-cell detection with rare_confirmation) to clear the rare-FPR caveat,
then flip fine_global to default for the fine path.
