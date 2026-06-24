# Report plot implementation plan

Plan for implementing the missing informative figures and wiring them into the
existing **QC-first** report (do not redo the structural redesign). Written
before any code change (Part 1).

## Conventions to reuse (discovered)

- Plotly publication layer: `tissueresolve.plotting.export.export_figure(fig,
  out_dir, name, *, data=..., caption=..., data_comment=[...])` writes
  `<name>.html`, static PNG/SVG/PDF (when kaleido present) **and**
  `<name>.data.tsv` in the same dir. The report's `_figure_cards` globs
  `figures/*.html` and embeds each + links its sibling `.data.tsv`. **New figures
  must write into the section `figures/` dir** so the existing report consumes
  them automatically (no new report system).
- Layout/colors: `style.plotly_layout(title, subtitle=, height=)`,
  `style.color_sequence(n)`, `style.ESTIMATE_SUBTITLE`.
- Hierarchical palette: `palette.build_hierarchical_color_map(fine, mapping)` →
  DataFrame (`broad_cell_type, fine_cell_type, color_hex, color_role, …`);
  `OTHER_COLOR` (light grey), `UNRESOLVED_COLOR` (dark grey). A small helper
  `palette.color_dicts_from_map(df)` will be added to return `(fine→hex,
  broad→hex)`.
- Captions: `tissueresolve.plotting.captions` (add new caption builders).
- Reference object `ReferenceSignature`: `gene_names`, `cell_types`,
  `n_cells_per_type`, `R_cpm`/`R_log` (G×K signature matrices).

## Figures already generated (do NOT duplicate)

| section | existing figures |
|---|---|
| bulk | bulk_composition_clustered_barplot, bulk_composition_heatmap, bulk_qc_summary, bulk_uncertainty_plot, bulk_separability_heatmap, bulk_spillover_heatmap, bulk_main_summary_figure, spillover_network |
| spatial | spatial_mean_composition_barplot, spatial_dominant_cell_type_map, spatial_abundance_maps, spatial_morans_i_barplot, spatial_spot_pie_charts (exploratory), spatial_separability/spillover_heatmap, he_spots_check, he_dominant_cell_type, he_abundance_* |
| reference | **NONE** (empty — biggest gap) |
| signature/hierarchy | none (tables only) |
| resolution | none (tables only) |
| benchmark | none in report (external harness writes its own HTML) |

## Implemented in THIS pass (data available, high value)

New module `plotting/reference_qc_plots.py` (Part 3), source data alongside in
`reference/figures/`:
| figure_id | input table | source data |
|---|---|---|
| reference_broad_family_composition | `reference/cell_type_counts.tsv` + hierarchy mapping | `reference_broad_family_composition.data.tsv` |
| reference_fine_subpopulation_support | cell_type_counts + mapping (+ min-cells line) | `..._support.data.tsv` |
| reference_celltype_imbalance | cell_type_counts (ranked + Lorenz/Gini) | `..._imbalance.data.tsv` |
| gene_overlap_by_modality | bulk + spatial gene_overlap (from `*_warnings.json`) | `gene_overlap_by_modality.data.tsv` |
| reference_suitability_components | `reference_suitability_components.tsv` | `..._components.data.tsv` |

New module `plotting/signature_qc_plots.py` (Part 4), into `reference/figures/`:
| figure_id | input | source data |
|---|---|---|
| signature_matrix_heatmap | `ref.R_log` top-variance genes × cell types (grouped by family) | `signature_matrix_heatmap.data.tsv` |
| top_confusable_pairs | `resolution/pairwise_separability.tsv` (top 20) | `top_confusable_pairs.data.tsv` |
| within_vs_between_family_separability | pairwise_separability + mapping | `..._separability.data.tsv` |
| hierarchy_map | cell_type_counts + mapping + unresolved flags (nested barplot) | `hierarchy_map.data.tsv` |
| marker_support_by_family | pairwise_separability `n_discriminating_genes` by family | `marker_support_by_family.data.tsv` |

New module `plotting/resolution_plots.py` (Part 8), into `resolution/figures/`:
| figure_id | input | source data |
|---|---|---|
| separability_distribution | pairwise_separability (histogram + threshold) | `separability_distribution.data.tsv` |
| unresolved_mass_by_family | `hierarchical/bulk_unresolved_family_mass.tsv` | `unresolved_mass_by_family.data.tsv` |
| trusted_resolution_summary | mapping + separability + unresolved (recommended level) | `trusted_resolution_summary.data.tsv` |

New module `plotting/benchmark_report_plots.py` (Part 9), into `benchmark/figures/`
(only if `benchmarks/outputs/*.tsv` exist, else manifest `missing_data`):
| figure_id | input | source data |
|---|---|---|
| benchmark_method_status_summary | `real_external_method_status.tsv` | `..._status_summary.data.tsv` |
| bulk_accuracy_leaderboard | benchmark bulk metrics (executed/imported only) | `bulk_accuracy_leaderboard.data.tsv` |
| runtime_comparison | status `runtime_seconds` | `runtime_comparison.data.tsv` |
| composite_scorecard | `composite_scores.tsv` (labeled scorecard, per modality) | `composite_scorecard.data.tsv` |

## Manifest (Part 2)

Extend `report/figures.py` `FigureRecord`/manifest with `caption` and
`reason_if_missing`; status vocabulary `generated | missing_data | skipped`.
Keep existing columns (backward compatible). Report links it in the Outputs
section (already does). A new generator records `missing_data` rows when a
figure's inputs are absent — no crash.

## Tests to add

- `tests/plotting/test_reference_qc_plots.py`, `test_signature_qc_plots.py`,
  `test_resolution_plots.py`, `test_benchmark_report_plots.py`: each plot from
  toy data writes html + `.data.tsv`; missing input → handled (no crash);
  hierarchical colors used; top-genes cap on heatmap.
- `tests/plotting/test_palette.py`: add `color_dicts_from_map` test.
- `tests/report/`: manifest has new columns; missing figures carry
  `reason_if_missing`.
- Integration ordering already covered by `test_report_redesign_ordering.py`.

## Deferred this pass (documented, manifest `missing_data`/`skipped`)

- **Parts 6–7 bulk/spatial result plots** that essentially already exist
  (composition, heatmap, dominant map, abundance maps, mean composition,
  Moran's I). Genuinely-missing extras — `solver_selection_cv`,
  `bulk_prediction_entropy`, `bulk_unresolved_mass` (bulk),
  `spatial_entropy_map`, `spatial_dominant_fraction_map`,
  `spatial_unresolved_mass_map` — are **deferred** to a follow-up pass (they need
  per-spot coordinate plumbing already used by `histology.py`); they will be
  recorded in the manifest as `missing_data` with a reason so the report does not
  pretend they exist.
- `spatial_concordance_dashboard` (Part 9.5): only meaningful with multiple
  executed spatial methods; recorded as `missing_data` otherwise.

## Non-negotiables honored

No core-algorithm changes; no prediction changes; RNA-derived labelling kept; no
spatial-accuracy claims (no ground truth); deterministic hierarchical palette;
every figure saves source data; nothing committed; outputs not committed.
