# Report figure implementation report

Implementation of the missing informative figures and their integration into the
existing QC-first report (the structural redesign was **not** changed).

## Figures implemented (17 new figure generators)

**Reference QC** (`plotting/reference_qc_plots.py`, report section 2):
- `reference_broad_family_composition` — cells per broad family (sorted, family palette)
- `reference_fine_subpopulation_support` — fine-label counts grouped by family + min-cells line
- `reference_celltype_imbalance` — ranked counts + Gini imbalance score
- `gene_overlap_by_modality` — reference/query/shared genes per modality
- `reference_suitability_components` — traffic-light component bar

**Signature quality & hierarchy** (`plotting/signature_qc_plots.py`, section 3):
- `signature_matrix_heatmap` — top-variance genes × cell types (grouped by family, row z-scored, top-N only)
- `top_confusable_pairs` — lollipop of the 20 least-separable pairs (severity-coloured)
- `within_vs_between_family_separability` — box+points, within vs between family
- `hierarchy_map` — broad→fine treemap sized by reference cells, ⚠ unresolved flags
- `marker_support_by_family` — mean discriminating markers per within-family pair + threshold

**Resolution / spillover / unresolved** (`plotting/resolution_plots.py`, section 8):
- `separability_distribution` — histogram + high-risk threshold
- `unresolved_mass_by_family` — mean unresolved mass per family
- `trusted_resolution_summary` — per-family recommended interpretation level (broad→fine)

**Benchmark** (`plotting/benchmark_report_plots.py`, section 9):
- `benchmark_method_status_summary` — methods per status
- `bulk_accuracy_leaderboard` — measured accuracy, executed/imported bulk only (spatial excluded)
- `runtime_comparison` — runtime for executed/imported
- `composite_scorecard` — weighted scorecard, explicitly labelled "NOT objective accuracy", per modality

Every figure: writes `<name>.data.tsv` source data, uses the deterministic
hierarchical palette (families distinct, fine subtones, grey for Other/
unresolved), and is embedded in the report (PNG `<img>` or `<iframe>`).

## Figures deferred (documented, recorded as missing_data when applicable)

Already-existing bulk/spatial result plots were **not duplicated** (composition
barplot/heatmap, dominant-type map, abundance maps, mean composition, Moran's I,
H&E overlays, separability/spillover heatmaps). Genuinely-missing result extras
deferred to a follow-up pass (need per-spot coordinate plumbing or are
harness-dependent):
- bulk: `solver_selection_cv`, `bulk_prediction_entropy`, `bulk_unresolved_mass`
- spatial: `spatial_entropy_map`, `spatial_dominant_fraction_map`,
  `spatial_unresolved_mass_map`
- benchmark: `spatial_concordance_dashboard` (needs ≥2 executed spatial methods)

When their inputs are absent the report records a `missing_data` manifest row
with a reason rather than pretending the figure exists.

## Files created
- `src/tissueresolve/plotting/reference_qc_plots.py`
- `src/tissueresolve/plotting/signature_qc_plots.py`
- `src/tissueresolve/plotting/resolution_plots.py`
- `src/tissueresolve/plotting/benchmark_report_plots.py`
- `tests/plotting/test_reference_qc_plots.py`
- `tests/plotting/test_signature_qc_plots.py`
- `tests/plotting/test_resolution_plots.py`
- `tests/plotting/test_benchmark_report_plots.py`
- `examples/real_breast_cancer/tests/test_report_figures.py`
- `docs/REPORT_PLOT_IMPLEMENTATION_PLAN.md`, `docs/REPORT_FIGURE_IMPLEMENTATION_REPORT.md`

## Files modified
- `src/tissueresolve/plotting/palette.py` — `color_dicts_from_map` helper
- `src/tissueresolve/report/figures.py` — manifest gains `caption`,
  `reason_if_missing`; status vocab `generated/missing_data/skipped`
- `examples/real_breast_cancer/scripts/07_generate_reports.py` —
  `generate_diagnostic_figures`, `generate_benchmark_figures`, helpers, section
  `_figure_cards` wiring (signature/resolution/benchmark), 17 specific captions,
  manifest missing-data records + outputs metrics
- `tests/plotting/test_palette.py`, `tests/report/test_report_components.py` — new tests

## Tests added
- 4 plotting test modules (reference/signature/resolution/benchmark) — generate
  from toy data, save `.data.tsv`, missing input handled, hierarchical colours,
  top-genes cap.
- palette `color_dicts_from_map`; manifest `caption`/`reason_if_missing` +
  missing-data record; examples figure-integration (caption specificity, figure
  embedding, cross-module palette reuse).

## pytest result
- Targeted/affected subset (plotting + report + benchmark + examples report
  tests): **PASS** (see final summary).
- Full suite (~35 min incl. real-data examples) to be run as the final gate.

## Quantitative
- Embedded figures in the regenerated report: **39** (33 PNG + 6 interactive iframe)
- Empty figure bodies: **0**
- Generic captions: **0**
- Figure-manifest rows: **39** (39 generated, 0 missing_data on this dataset)
- Source-data tables: one `.data.tsv` per figure (39+)
- Hierarchical palette used: **yes** (all new figures)
- QC-first: **yes** (reference/signature/input precede results — unchanged)

## Remaining limitations
- Deferred result-map figures listed above.
- `signature_matrix_heatmap` shows top-variance genes (not curated markers).
- Benchmark figures depend on `benchmarks/outputs/*.tsv` being present.
- Spatial figures still come from the existing generators; new spatial maps
  (entropy / dominant-fraction / unresolved) are deferred.

## Suggested next commit message
```
Add reference/signature/resolution/benchmark QC figures to the report

Implement 17 informative plot generators (reference_qc_plots,
signature_qc_plots, resolution_plots, benchmark_report_plots), wire them into
the QC-first report sections, extend the figure manifest with caption/
reason_if_missing and a generated/missing_data status vocabulary, and add
offline tests. All figures save source data, use the deterministic hierarchical
palette, and are embedded. No core-algorithm or prediction changes.
```
