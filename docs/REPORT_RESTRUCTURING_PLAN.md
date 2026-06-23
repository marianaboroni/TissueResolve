# Report Restructuring Plan (QC-first, decision-oriented)

Status: **plan only — no code changed in this document** (Part 1 deliverable).

Central rule applied throughout: *a figure stays in the main report only if it
helps the user decide whether a prediction is reliable (**reliability test**) or
helps interpret the final bulk/spatial prediction (**interpretation test**).*
Everything else → technical appendix, source-data only, or removed.

---

## 0. Current state (audited)

Source audited: the harness unified report
(`examples/real_breast_cancer/outputs/report.html`, **112 KB**),
`technical_appendix.html` (27 KB), and the ~45 figure HTML files under
`outputs/{reference,signature,resolution,bulk,spatial,benchmark}/figures/`.
Builder: `examples/real_breast_cancer/scripts/07_generate_reports.py::generate_unified_report`
(renders via the canonical `report/orchestration.py` → unified shell).

**Current sections (14):** 1 Executive decision summary · 2 Reference quality ·
3 Signature quality & hierarchy · 4 Input data quality · 5 Bulk deconvolution ·
6 Spatial deconvolution · 7 Hierarchical broad→fine · 8 Resolution, separability
& spillover · 9 Bulk benchmark · 10 Spatial benchmark · 11 Benchmark scorecards
& method status · 12 Warnings · 13 Methods · 14 Output files & source data.

**Problems:** 14 sections (target 10); predictions (5–6) appear **before** the
trusted-resolution/QC decision (7–8) — not QC-first; a whole benchmark-scorecard
section (11) of bookkeeping figures; methods + outputs split into two trailing
sections; ~45 figures available with several redundant/exploratory ones eligible
for the main grid; report is large and table/heatmap-heavy.

---

## 1. Per-figure classification

Columns: **Rel?** passes reliability test · **Int?** passes interpretation test
· **Clr?** visually clear · **Pub?** publication-quality · **Dup?** redundant
with another figure · **Decision**. Q1–Q6 (assess reference / input / resolution
/ bulk / spatial / benchmark) noted in *serves*.

### Reference & input
| Figure | Rel? | Int? | Clr? | Pub? | Dup? | Decision | serves |
|---|---|---|---|---|---|---|---|
| `reference_suitability_components` | ✅ | – | ✅ | ✅ | no | **MAIN** | Q1 reference good enough (PASS/CAUTION/WARN per component) |
| `reference_broad_family_composition` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** | Q1 balance; frames what families exist |
| `gene_overlap_by_modality` | ✅ | – | ✅ | ✅ | no | **MAIN** | Q3 query↔reference overlap (bulk+spatial) |
| `reference_celltype_imbalance` | ✅ | – | ✅ | ✅ | yes (suitability) | APPENDIX | imbalance detail already flagged by suitability |
| `reference_fine_subpopulation_support` | ✅ | – | ✅ | ✅ | yes (marker_support) | APPENDIX | per-fine support detail |

### Signature & resolution (the "trusted resolution" decision)
| Figure | Rel? | Int? | Clr? | Pub? | Dup? | Decision | serves |
|---|---|---|---|---|---|---|---|
| `trusted_resolution_summary` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** | Q4 broad / selected-fine / fine trust (the decision figure) |
| `top_confusable_pairs` | ✅ | – | ✅ | ✅ | no | **MAIN** | Q2 can signatures separate expected types |
| `within_vs_between_family_separability` | ✅ | – | ✅ | ✅ | no | **MAIN** | Q4 are fine types separable within family |
| `unresolved_mass_by_family` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** | Q4 how much mass stays at family level |
| `separability_distribution` | ✅ | – | ✅ | ✅ | yes (confusable/within-between) | APPENDIX | overview duplicated by the two MAIN separability views |
| `hierarchy_map` | – | ✅ | ✅ | ✅ | partial | APPENDIX | structural; decision conveyed by trusted_resolution_summary |
| `marker_support_by_family` | ✅ | – | ✅ | ✅ | yes | APPENDIX | marker counts; overlaps separability |
| `signature_matrix_heatmap` | – | – | partial | ✅ | – | APPENDIX | large gene×type heatmap (already appendix) |

### Bulk predictions & QC
| Figure | Rel? | Int? | Clr? | Pub? | Dup? | Decision | serves |
|---|---|---|---|---|---|---|---|
| `bulk_composition_clustered_barplot` | – | ✅ | ✅ | ✅ | no | **MAIN** | Q5 final bulk composition (top types + Other) |
| `bulk_qc_summary` | ✅ | – | ✅ | ✅ | no | **MAIN** | Q5 per-sample reconstruction QC |
| `bulk_uncertainty_plot` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** (or card if not computed) | Q5 bootstrap CIs / low-confidence flags |
| `bulk_main_summary_figure` | – | ✅ | partial | ✅ | yes (multi-panel of the above) | APPENDIX | duplicates the cleaner single-idea figures |
| `bulk_composition_heatmap` | – | ✅ | ✅ | ✅ | yes (clustered barplot) | APPENDIX | denser duplicate of predictions |
| `bulk_separability_heatmap` | – | – | partial | ✅ | yes | APPENDIX | full matrix (already appendix) |
| `bulk_spillover_heatmap` | – | – | partial | ✅ | yes | APPENDIX | full matrix (already appendix) |
| `spillover_network` (bulk) | – | – | partial | ✅ | yes | APPENDIX | exploratory network (already appendix) |

### Spatial predictions & structure
| Figure | Rel? | Int? | Clr? | Pub? | Dup? | Decision | serves |
|---|---|---|---|---|---|---|---|
| `spatial_mean_composition_barplot` | – | ✅ | ✅ | ✅ | no | **MAIN** | Q6 tissue-wide composition |
| `he_dominant_cell_type` | – | ✅ | ✅ | ✅ | no | **MAIN** | Q6 dominant type on tissue image |
| `spatial_abundance_maps` | – | ✅ | ✅ | ✅ | no | **MAIN** | Q6 spatial distribution of key types |
| `spatial_morans_i_barplot` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** | Q6 which populations are spatially structured |
| `spatial_dominant_cell_type_map` | – | ✅ | ✅ | ✅ | yes (he_dominant) | APPENDIX | coordinate-only fallback of he_dominant |
| `spatial_main_summary_figure` | – | ✅ | partial | ✅ | yes | APPENDIX | multi-panel duplicate (already appendix) |
| `he_abundance_<type>` ×3 | – | ✅ | ✅ | ✅ | yes (abundance_maps) | APPENDIX | per-type exploratory duplicates |
| `he_spots_check` | ✅ | – | partial | – | – | APPENDIX | overlay-alignment QC (debug) |
| `spatial_separability_heatmap` | – | – | partial | ✅ | yes | APPENDIX | full matrix (already appendix) |
| `spatial_spillover_heatmap` | – | – | partial | ✅ | yes | APPENDIX | full matrix (already appendix) |
| `spatial_spot_pie_charts` | – | – | partial | ✅ | – | APPENDIX | exploratory (already appendix) |
| `spillover_network` (spatial) | – | – | partial | ✅ | yes | APPENDIX | exploratory (already appendix) |

### Benchmark (bulk has ground truth; spatial does not)
| Figure | Rel? | Int? | Clr? | Pub? | Dup? | Decision | serves |
|---|---|---|---|---|---|---|---|
| `benchmark/bulk/bulk_fine_accuracy_leaderboard` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** | Q7 bulk accuracy vs pseudobulk truth |
| `benchmark/spatial/spatial_structure_metrics_summary` | ✅ | ✅ | ✅ | ✅ | no | **MAIN** | Q7 spatial concordance/structure (no accuracy) |
| `benchmark/bulk/bulk_rmse_mae_comparison` | ✅ | – | ✅ | ✅ | yes (leaderboard) | APPENDIX | second accuracy view |
| `benchmark/bulk/bulk_runtime_comparison` | – | – | ✅ | ✅ | – | APPENDIX | runtime (secondary) |
| `benchmark/bulk/bulk_method_status_summary` | – | – | ✅ | – | yes | APPENDIX | executed/skipped bookkeeping |
| `benchmark/spatial/spatial_runtime_comparison` | – | – | ✅ | ✅ | – | APPENDIX | runtime |
| `benchmark/spatial/spatial_method_status_summary` | – | – | ✅ | – | yes | APPENDIX | bookkeeping |
| `benchmark/spatial/spatial_output_completeness_summary` | – | – | ✅ | – | – | APPENDIX | completeness bookkeeping |
| `benchmark/figures/composite_scorecard` | – | – | ✅ | ✅ | yes (scorecard dir) | APPENDIX | opinionated composite, not objective accuracy |
| `benchmark/scorecard/composite_scorecard` | – | – | ✅ | ✅ | yes (dup) | APPENDIX | duplicate |
| `benchmark/figures/bulk_accuracy_leaderboard` | ✅ | ✅ | ✅ | ✅ | **yes (dup of bulk/ leaderboard)** | **REMOVE** | redundant top-level duplicate |
| `benchmark/figures/runtime_comparison` | – | – | ✅ | ✅ | yes | APPENDIX | duplicate runtime |
| `benchmark/figures/benchmark_method_status_summary` | – | – | ✅ | – | yes | APPENDIX | duplicate bookkeeping |
| `benchmark/scorecard/benchmark_method_status_summary` | – | – | ✅ | – | yes | APPENDIX | duplicate bookkeeping |

**MAIN total: 16 figures** (3 reference/input + 4 signature/resolution + 3 bulk +
4 spatial + 1 bulk-bench + 1 spatial-bench) — within the 15–22 target.
**REMOVE: 1** (`benchmark/figures/bulk_accuracy_leaderboard`, exact duplicate).
Everything else → APPENDIX (or already there). No new figures are introduced.

Full raw prediction tables, separability/spillover matrices, composite-score
tables, method-status tables, installation logs and debug metadata → APPENDIX /
`source_data/` only (never inline in the main report).

---

## 2. New main report structure (10 sections, QC-first)

| # | Section | From current | Main figures |
|---|---|---|---|
| 1 | Executive decision summary | 1 | (decision status grid + checklist; no plots) |
| 2 | Reference and signature quality | merge 2 + 3 | reference_suitability_components, reference_broad_family_composition, top_confusable_pairs, within_vs_between_family_separability |
| 3 | Input compatibility | 4 | gene_overlap_by_modality |
| 4 | Trusted resolution and uncertainty | merge 7 + 8 (**moved before predictions**) | trusted_resolution_summary, unresolved_mass_by_family |
| 5 | Final bulk predictions | 5 | bulk_composition_clustered_barplot, bulk_qc_summary, bulk_uncertainty_plot |
| 6 | Final spatial predictions | 6 | spatial_mean_composition_barplot, he_dominant_cell_type, spatial_abundance_maps, spatial_morans_i_barplot |
| 7 | Bulk benchmark summary | 9 | bulk_fine_accuracy_leaderboard |
| 8 | Spatial benchmark summary | 10 | spatial_structure_metrics_summary |
| 9 | Warnings and limitations | 12 | (severity warning box; no plots) |
| 10 | Methods and source data | merge 13 + 14 | (methods text + links; no plots) |

Key moves vs. current:
- **QC-first ordering enforced:** sections 2–4 (reference, signature, input,
  trusted resolution) all precede predictions (5–6); predictions precede
  benchmarks (7–8); warnings appear in §1 (decision summary) **and** §9.
- Merge current 2+3 → §2; merge 7 (hierarchical) + 8 (resolution/spillover) → §4
  and **move it above predictions**; merge 13+14 → §10.
- **Drop current §11 (Benchmark scorecards & method status) from the main report**
  → technical appendix (composite score is an opinionated scorecard, not
  objective accuracy).
- Hierarchical broad→fine is presented *inside* §4 as the trusted-resolution view
  (broad family composition + unresolved mass), not as its own section.

### Technical appendix (`technical_appendix.html`) contains
All APPENDIX figures above: full separability/spillover heatmaps (bulk+spatial),
spillover networks, signature_matrix_heatmap, spatial_spot_pie_charts, the
multi-panel `*_main_summary_figure`s, per-type `he_abundance_*`, hierarchy_map,
separability_distribution, marker_support_by_family, reference imbalance/support
detail, all benchmark runtime/method-status/output-completeness/composite-score
figures, full benchmark tables, full raw prediction tables, composite-score
tables, debug metadata and installation logs. The appendix links to/from the
main report.

---

## 3. Size target & adherence

- **10 sections** (was 14). ✅
- **~16 main figures** (target 15–22). ✅
- No full matrices / large raw tables / debug outputs inline → appendix or
  `source_data/`. ✅
- "When in doubt → appendix" applied to: hierarchy_map, separability_distribution,
  marker_support_by_family, bulk_composition_heatmap, spatial_dominant_cell_type_map.

---

## 4. Implementation principle (how to build it — Part 2/3/4)

Redesign by **subtraction and reorganisation, not addition**. No new modelling,
no new figures. Concretely, in `07_generate_reports.py::generate_unified_report`
(and matching captions/anchors):

1. **Re-route figures** by expanding the appendix sets to match §1 above:
   - `_BULK_APPENDIX_FIGS` += `bulk_composition_heatmap`.
   - `_SPATIAL_APPENDIX_FIGS` += `spatial_dominant_cell_type_map`,
     `he_abundance_*`, `he_spots_check`.
   - new `_SIGNATURE_APPENDIX_FIGS` += `hierarchy_map`, `marker_support_by_family`.
   - new `_RESOLUTION_APPENDIX_FIGS` = `separability_distribution`.
   - new `_REFERENCE_APPENDIX_FIGS` = `reference_celltype_imbalance`,
     `reference_fine_subpopulation_support`.
   - benchmark: keep only `bulk_fine_accuracy_leaderboard` (bulk) and
     `spatial_structure_metrics_summary` (spatial) in MAIN; route the rest
     (rmse/mae, runtime, method-status, output-completeness, composite-scorecard)
     to appendix; **drop** the duplicate top-level `bulk_accuracy_leaderboard`.
2. **Merge sections** 2+3 → "Reference and signature quality"; 7+8 → "Trusted
   resolution and uncertainty"; 13+14 → "Methods and source data".
3. **Reorder** so "Trusted resolution and uncertainty" is emitted **before** the
   bulk/spatial prediction sections.
4. **Improve captions/hierarchy**: each main figure gets a one-line "what to
   check" caption framed as a reliability or interpretation question; tables in
   main sections become compact summaries with full versions collapsed/appendix.
5. Keep the same `source_data/` + figure-manifest behaviour (every figure still
   writes its `.data.tsv`); only *placement* changes.

Constraints honoured: no algorithm/value changes; no new benchmark methods; no
reference adaptation or expression reconstruction; warnings never hidden (shown
in §1 and §9); spatial benchmark labelled concordance/structure (no accuracy);
fine-subtype resolution gated by §4; nothing committed.

### Tests to add/keep green (when implemented)
- main report has ≤10 sections in the exact QC-first order (§2 ordering test);
- "Trusted resolution" precedes bulk/spatial predictions;
- each appendix figure id is **absent** from `report.html` and **present** in
  `technical_appendix.html`;
- no `bulk_accuracy_leaderboard` duplicate in either file;
- main report figure count within 15–22;
- warnings appear in both the executive summary and the warnings section;
- spatial benchmark section contains the no-ground-truth caveat.
