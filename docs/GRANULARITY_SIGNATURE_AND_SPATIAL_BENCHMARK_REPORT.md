# Granularity signature & spatial-benchmark — implementation report

Status report for the "signature granularity + spatial benchmark" task. Honest
about what is implemented-and-tested vs what is a documented next increment.

## 1. Audit result
`docs/GRANULARITY_AND_SPATIAL_BENCHMARK_AUDIT.md` (Part 1, read-only). Headline:
the **within-family signature layer already exists and is wired into the solver
API** (committed `43c9eda`), but is **not yet built from the real AnnData and
passed in the report run**; there is **no broad-specific panel builder**; and the
spatial benchmark (split in `b657011`) still lacks marker-recovery, concordance,
a synthetic accuracy scenario, and a transparent multi-metric ranking.

## 2–4. Within-family genes recalculated / used by solver?
**Yes, recalculated; yes, the solver consumes them when supplied** (commit
`43c9eda`): `reference/within_family_markers.py` recomputes HVGs + one-vs-rest DE
+ pairwise discriminative genes **inside each broad family**, excludes
mito/ribo/stress + query-undetected genes, scores and weights them 0–1, and
`hierarchy.evaluate_within_family_resolvability` /
`compute_within_family_subtype_confidence` / `assemble_hierarchical_estimates` /
`bulk.run_hierarchical_bulk` accept an optional `family_gene_panels` that
restricts the fine step to family genes (verified: a flat panel flips a
resolvable family to unresolvable; default `None` reproduces prior behaviour).
**Not yet done:** building the panels from the real reference AnnData and passing
them in `07_generate_reports.py` (so the *real-data* run still uses global genes
until wired), the spatial solver path, and a dedicated
`build_broad_signature_panel`.

## 10–11. Spatial benchmark — implemented this increment
**`benchmarks/shared/spatial_multimetric_ranking.py`** (Parts 6–7), 12 tests
green (`benchmarks/tests/test_spatial_multimetric_ranking.py`):
- **Primitives:** `jensen_shannon_divergence` (handles zeros, symmetric, ∈[0,1]),
  `marker_recovery_score` (predicted abundance vs marker-gene expression — a
  proxy, not validation), `normalize_metric_direction` (min–max + invert for
  error/runtime).
- **Scenario-separated scoring (never pooled):**
  `compute_real_spatial_score` (marker recovery / concordance / structure /
  stability / runtime-completeness / interpretability — **no accuracy**),
  `compute_synthetic_spatial_accuracy_score` (inverse error from RMSE/JSD/MAE,
  correlation from Pearson/Spearman/CCC, dominant accuracy, unresolved-aware,
  runtime), `compute_spatial_overall_scorecard` (labelled scorecard, not truth).
- **Ranking:** `rank_spatial_methods` excludes skipped/exported/failed, flags
  imported, and refuses to claim a ranking with <2 rankable methods;
  `explain_spatial_ranking` emits transparent Markdown (states "not accuracy" for
  the real scenario). Default weights match the task spec.

## What this delivers vs the prior split
`b657011` already split the report into 9 Bulk / 10 Spatial / 11 Scorecard with
explicit no-ground-truth captions. This increment adds the **measurement +
ranking foundation** the spatial section needs (marker recovery, JSD,
concordance-ready scoring, scenario-separated transparent ranking).

## Effects on accuracy / unresolved mass
**Not yet measured on real data** — the granular panels are not yet passed in the
real-data report run, so no real-data fine/family Pearson or unresolved-mass
delta is claimed here. The mechanism is verified on synthetic data only. **No
improvement is claimed where it has not been measured.**

## Files created (this increment)
- `benchmarks/shared/spatial_multimetric_ranking.py`
- `benchmarks/tests/test_spatial_multimetric_ranking.py`
- `docs/GRANULARITY_AND_SPATIAL_BENCHMARK_AUDIT.md`
- `docs/GRANULARITY_SIGNATURE_AND_SPATIAL_BENCHMARK_REPORT.md`

## Files relied on (already committed)
- `src/tissueresolve/reference/within_family_markers.py`,
  `hierarchy.py`, `bulk/hierarchical.py` (`43c9eda`)
- `src/tissueresolve/plotting/{bulk,spatial}_benchmark_plots.py`,
  `07_generate_reports.py` (`b657011`)

## Tests added / result
12 new ranking tests — **all pass**. (Within-family + benchmark-split tests from
the two prior commits remain green.)

## Remaining limitations / documented next increment
1. **Granular signatures (Parts 2–3, 5):** `build_broad_signature_panel`,
   `write_granular_signature_outputs` (the `outputs/signatures/*.tsv` +
   `family_solver_gene_panels.tsv`), **wire panels end-to-end** into the
   real-data report + spatial solver, and the "Multi-granularity signature
   construction" report subsection. (The algorithm exists; this is plumbing +
   outputs + figures.)
2. **Before/after benchmark (Part 4):** global vs within-family strategy
   metrics (`granularity_signature_benchmark.tsv`).
3. **Synthetic spatial benchmark (Part 6B):** lightweight simulator (coordinates,
   regions, known proportions, sharp-border + out-of-reference scenarios) →
   `spatial_synthetic_metrics.tsv`.
4. **Bubble-heatmap dashboard (Part 8)** + **spatial report restructuring
   (Part 9)** consuming the ranking module's outputs.
5. **Real-data run (Part 11)** once the wiring lands.

## Safe to commit?
The new ranking module + tests are self-contained and green; safe to commit.
The remaining parts are not yet implemented and should land as the next
increment before claiming the granularity issue is fully fixed on real data.
