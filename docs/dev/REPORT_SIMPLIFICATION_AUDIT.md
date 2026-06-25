# Report simplification audit (Part 3)

Read-only audit of the current generated report
(`examples/real_breast_cancer/outputs/report.html`, 14 sections) classifying
every section and figure for a simpler main report. Classes:

- **MAIN** — essential to assess quality or interpret final predictions
- **SECONDARY** — useful; collapse by default
- **APPENDIX** — technical; move to a technical appendix
- **SOURCE ONLY** — keep as `.data.tsv`, not shown as a figure
- **REMOVE** — redundant/confusing in the main body

## Section-level verdicts

| # | Section | Verdict | Rationale |
|---|---------|---------|-----------|
| 1 | Executive decision summary | **MAIN** | one-page interpretability verdict (status cards + checklist) |
| 2 | Reference quality | **MAIN** | can the reference support deconvolution? |
| 3 | Signature quality & hierarchy | **MAIN** (trim) | keep separability summary + recommended level; full heatmaps → appendix |
| 4 | Input data quality | **MAIN** | does the query match the reference? |
| 5 | Bulk deconvolution | **MAIN** (trim) | keep composition/heatmap/QC; heatmaps+network → appendix |
| 6 | Spatial deconvolution | **MAIN** (trim heavily) | keep H&E/dominant/top-abundance/mean/Moran; rest → appendix |
| 7 | Hierarchical broad→fine | **MERGE→SECONDARY** | fold into §3/§5–6; standalone section is redundant |
| 8 | Resolution, separability & spillover | **MAIN** (trim) | keep trusted-resolution + unresolved-mass; full pairwise table → appendix |
| 9 | Bulk benchmark | **MAIN** | accuracy vs pseudobulk truth (bulk only) |
| 10 | Spatial benchmark | **MAIN** | concordance/structure, no ground truth (spatial only) |
| 11 | Benchmark scorecards & status | **SECONDARY** | scorecard after measured metrics; collapse |
| 12 | Warnings & limitations | **MAIN** (summarize) | visible but summarized |
| 13 | Methods | **APPENDIX** | reproducibility detail |
| 14 | Output files & source data | **APPENDIX** | links to manifest/source data |

Target main body ≈ 9 logical sections (Part 4 structure); the rest collapse or
move to a technical appendix.

## Figure-level verdicts

### Reference QC (§2) — currently 5
| figure | verdict |
|---|---|
| reference_suitability_components | MAIN |
| reference_broad_family_composition | MAIN |
| reference_fine_subpopulation_support | MAIN |
| reference_celltype_imbalance | SECONDARY (collapse) |
| gene_overlap_by_modality | MAIN (belongs in Input §4) |

### Signature/hierarchy (§3) — currently 5
| figure | verdict |
|---|---|
| top_confusable_pairs | MAIN |
| within_vs_between_family_separability | MAIN |
| signature_matrix_heatmap | APPENDIX (large) |
| hierarchy_map | SECONDARY (collapse) |
| marker_support_by_family | SECONDARY (collapse) |

### Bulk (§5) — currently 8
| figure | verdict |
|---|---|
| bulk_composition_clustered_barplot | MAIN |
| bulk_composition_heatmap | MAIN |
| bulk_qc_summary | MAIN |
| bulk_uncertainty_plot | SECONDARY (collapse) |
| bulk_separability_heatmap | APPENDIX |
| bulk_spillover_heatmap | APPENDIX |
| spillover_network (bulk) | APPENDIX |
| bulk_main_summary_figure | REMOVE from main (already in collapsible; keep appendix-only) |

### Spatial (§6) — currently 14
| figure | verdict |
|---|---|
| he_spots_check | MAIN (alignment first) |
| he_dominant_cell_type / spatial_dominant_cell_type_map | MAIN (one dominant map) |
| spatial_abundance_maps (top 3–4 broad) | MAIN |
| spatial_mean_composition_barplot | MAIN |
| spatial_morans_i_barplot | MAIN |
| he_abundance_<each cell type> (3 here) | APPENDIX (keep only top 3–4 broad in main; per-fine-type → appendix) |
| spatial_separability_heatmap | APPENDIX |
| spatial_spillover_heatmap | APPENDIX |
| spillover_network (spatial) | APPENDIX |
| spatial_spot_pie_charts | APPENDIX (exploratory; never primary) |
| spatial_main_summary_figure | REMOVE from main (appendix-only) |

### Resolution (§8) — currently 3
| figure | verdict |
|---|---|
| trusted_resolution_summary | MAIN |
| unresolved_mass_by_family | MAIN |
| separability_distribution | SECONDARY (collapse) |

### Benchmark (§9–11)
Bulk: `bulk_fine_accuracy_leaderboard` + `bulk_rmse_mae_comparison` +
`bulk_runtime_comparison` MAIN; `bulk_method_status_summary` SECONDARY;
`bulk_family_accuracy_leaderboard` SECONDARY. Spatial:
`spatial_structure_metrics_summary` + `spatial_method_status_summary` MAIN;
`spatial_runtime_comparison`/`completeness`/`concordance` SECONDARY. Scorecard
(`composite_scorecard`, `benchmark_method_status_summary`) SECONDARY, after
measured metrics, clearly labelled.

## Per-question answers (applied across sections)
- **Needed in main / helps assess quality / helps interpret predictions?** MAIN
  figures above answer yes; APPENDIX/SOURCE figures answer no.
- **Redundant with another figure?** `*_main_summary_figure` (duplicates
  dedicated plots) → REMOVE-from-main; duplicate bulk/spatial separability &
  spillover heatmaps → APPENDIX (one shared diagnostic suffices).
- **Too technical for main?** full separability/spillover heatmaps,
  signature_matrix_heatmap, spillover_network, spot pies → APPENDIX.
- **Replace with simpler figure?** the §3 separability story is better told by
  `top_confusable_pairs` + `trusted_resolution_summary` than by full heatmaps.
- **Source-data-only?** large pairwise-separability and full benchmark tables →
  SOURCE ONLY (collapsible link, not rendered in full).

## Summary of changes the implementation should make
1. Main body = the Part-4 nine-section logic; merge §7 into §3/§5–6; move §13–14
   to a technical appendix.
2. Per results section, keep only the MAIN figures inline; route APPENDIX
   figures (full heatmaps, spillover networks, signature heatmap, spot pies,
   per-fine-type H&E maps, `*_main_summary_figure`) to a collapsible "Technical
   appendix" (or a separate `technical_appendix.html`).
3. Collapse SECONDARY figures by default; keep large tables as SOURCE ONLY
   (collapsible + `.data.tsv` link).
4. Keep bulk/spatial benchmarks separate (already done); scorecard last.
5. Every retained main figure keeps a specific caption + source-data link; no
   empty cards; no spatial-accuracy claim without ground truth; fine-level
   results labelled cautious when separability is low.
