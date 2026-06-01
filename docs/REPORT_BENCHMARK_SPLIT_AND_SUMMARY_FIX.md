# Report benchmark split & main-summary-figure fix

## What was wrong
1. **Composite "main summary" figures** (`bulk_main_summary_figure`,
   `spatial_main_summary_figure`) were crowded multi-panel diagnostic dumps shown
   **interleaved among the primary figures** (alphabetical order put them 3rd /
   8th), duplicating clearer dedicated plots.
2. **Benchmark was a single section** mixing bulk and spatial method status,
   runtime and the composite scorecard. Although the accuracy leaderboard already
   excluded spatial, there was **no separate spatial benchmark** with
   concordance/structure metrics, and the presentation read as one shared
   comparison. (Audit: `docs/REPORT_SUMMARY_AND_BENCHMARK_AUDIT.md`.)

## Figures removed / replaced / moved
- `bulk_main_summary_figure` → **moved** to a collapsible "Technical composite
  diagnostics (bulk)" appendix; **excluded** from the primary bulk grid.
- `spatial_main_summary_figure` and `spatial_spot_pie_charts` → **moved** to a
  collapsible "Technical / exploratory figures (spatial)" appendix; **excluded**
  from the primary spatial grid.
- Primary bulk section now leads with composition barplot/heatmap/QC; primary
  spatial section leads with H&E alignment / dominant-family / abundance maps.

## New report section order
1 Executive decision summary · 2 Reference quality · 3 Signature quality &
hierarchy · 4 Input data quality · 5 Bulk deconvolution · 6 Spatial
deconvolution · 7 Hierarchical broad→fine · 8 Resolution/separability/spillover ·
**9 Bulk benchmark: accuracy against pseudobulk ground truth** · **10 Spatial
benchmark: concordance and spatial structure** · **11 Benchmark scorecards &
method status** · 12 Warnings & limitations · 13 Methods · 14 Output files &
source data. (Sidebar nav auto-updates from the section list.)

## Bulk benchmark figures (`plotting/bulk_benchmark_plots.py`)
`bulk_method_status_summary`, `bulk_fine_accuracy_leaderboard`,
`bulk_family_accuracy_leaderboard`, `bulk_rmse_mae_comparison`,
`bulk_runtime_comparison`. Bulk-only; accuracy metrics justified by pseudobulk
ground truth. Section states: "Bulk pseudobulk mixtures have known ground truth,
so accuracy metrics are valid here."

## Spatial benchmark figures (`plotting/spatial_benchmark_plots.py`)
`spatial_method_status_summary`, `spatial_concordance_heatmap` (needs ≥2 spatial
methods; otherwise `missing_data`), `spatial_structure_metrics_summary` (Moran's
I / entropy / dominant / near-zero fraction), `spatial_runtime_comparison`,
`spatial_output_completeness_summary`. **No Pearson/RMSE-vs-truth.** Section
states: "Real Visium data do not have spot-level ground truth in this benchmark;
therefore these metrics evaluate concordance and spatial structure, not
accuracy." Every spatial caption carries the no-ground-truth caveat.

## Benchmark scorecard (section 11)
`composite_scorecard` (+ overall method status), explicitly labelled a weighted
multi-criteria **scorecard**, not an objective accuracy ranking; bulk and spatial
scored within their own modality. Placed **after** the measured bulk/spatial
sections.

## Manifest / source data
Manifest now distinguishes `section ∈ {bulk_benchmark, spatial_benchmark,
benchmark_scorecard}`; figures + `.data.tsv` written under
`benchmark/{bulk,spatial,scorecard}/figures/`. Bulk and spatial benchmark
figures never share a section (test-enforced).

## Tests added
- `tests/plotting/test_split_benchmark_plots.py` (10) — bulk plots exclude
  spatial; spatial plots carry no-ground-truth caveat and expose no
  accuracy/RMSE/Pearson function; concordance needs ≥2 methods.
- `examples/real_breast_cancer/tests/test_benchmark_split.py` (8) —
  `_SUMMARY_FIGS` constant; `_figure_cards` exclude/only relocation; bulk
  captions mention ground truth; spatial captions disclaim accuracy; generated
  report has separate sections + collapsible appendix; manifest distinguishes
  benchmark sections (real-report checks skip when outputs absent).

## pytest result
Targeted runs green (split benchmark plots + integration: 10 + 8 passed). Full
affected subset (report + plotting + hierarchical + within-family + examples +
external benchmark) run as the gate — see final summary.

## Report path
`examples/real_breast_cancer/outputs/report.html` (regenerated; 14 sections,
benchmark split confirmed, main summaries in collapsible appendices).

## Remaining limitations
- `spatial_concordance_heatmap` requires ≥2 executed/imported spatial methods
  with per-spot predictions; recorded as `missing_data` otherwise.
- `bulk_family_accuracy_leaderboard` / `bulk_rmse_mae_comparison` depend on the
  benchmark metric TSVs exposing family-accuracy / rmse-mae columns.
- A **synthetic** spatial benchmark (known simulated truth) would enable a
  clearly-separate "Spatial synthetic benchmark: accuracy" subsection; not yet
  provided.
- The composite scorecard's non-accuracy dimensions remain qualitative priors.
