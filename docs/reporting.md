# Reporting

TissueResolve can generate HTML reports from a results directory.

## v0.1 architecture: concise main report + technical appendix

The unified report is split into two files in the output directory:

- **`report.html`** — the concise, QC-first, user-facing report. One question
  per section, one main figure per idea, broad-level before fine-level, tables
  collapsed, warnings summarized. Links to the technical appendix and figure
  manifest.
- **`technical_appendix.html`** — everything intentionally kept out of the main
  report: full separability/spillover heatmaps, spillover networks, spot pie
  charts, the large signature-matrix heatmap, composite/technical summary
  figures, and the complete source-data table listing. Links back to
  `report.html`.

Both files are produced through the unified single-page shell
(`report/unified.py` + `report/components.py`); the appendix is a second set of
sections rendered to a second file. Figure routing lives in
`examples/real_breast_cancer/scripts/07_generate_reports.py`
(`_BULK_APPENDIX_FIGS`, `_SPATIAL_APPENDIX_FIGS`, `_SIGNATURE_APPENDIX_FIGS`,
`_SUMMARY_FIGS`).

**Canonical code path (single path):** `report/orchestration.py` is the one
entry point (`generate_report(modality, source, out)`); it builds ordered
`(title, body)` sections from a results directory **or** an in-memory pipeline
result and renders them through `report/unified.py` + `report/components.py`
(+ `style`, `glossary`, `figures`). Content/data come from `report/sections.py`
(results-dir) and `report/html.py`'s `*_result_sections` (in-memory), with prose
from `report/interpretation.py` and loading from `report/assets.py`.

`report/html.py`'s public `generate_report`/`generate_bulk_report`/
`generate_spatial_report` are **deprecated shims** that delegate to
`orchestration`; the separate per-modality page renderer has been removed (the
CLI/API and `tissueresolve run` all emit the unified report). See
[`docs/REPORT_PATH_CONSOLIDATION_PLAN.md`](REPORT_PATH_CONSOLIDATION_PLAN.md).

## Generate a report

```bash
tissueresolve report --modality bulk --results-dir results/bulk
```

```bash
tissueresolve report --modality spatial --results-dir results/spatial
```

## Combined bulk + spatial report

One `tissueresolve run` processes one modality, written to its own directory.
To present bulk and spatial together (they share a reference), merge two
existing run directories:

```bash
tissueresolve combine-report \
  --bulk-dir results/bulk --spatial-dir results/spatial \
  --out results/combined
```

`combine-report` reads the two run directories (no re-running of deconvolution)
and writes `results/combined/report.html`, `methods.txt`, `warnings.json` and
`run_metadata.json`. Sections: (1) executive summary, (2) shared reference
summary (when metadata is available), (3) bulk QC and predictions, (4) spatial
QC and predictions, (5) bulk benchmark summary, (6) spatial benchmark summary,
(7) warnings and limitations (bulk and spatial warnings, attributed), (8)
methods and source-data links back to each run directory. Bulk and spatial — and
their benchmark summaries — are kept **separate**; the combined `warnings.json`
includes both runs' warnings tagged by modality. The report states explicitly
that it summarises two separate runs sharing a reference, **not** a single joint
bulk+spatial model. It validates that `--bulk-dir` is a bulk run and
`--spatial-dir` is a spatial run (when their metadata records a modality).

## What the report includes

- Executive summary cards
- Key findings and interpretation
- Main publication figure
- QC diagnostics and warnings
- Detailed tables in collapsible sections
- File list for reproducibility

For hierarchical runs, the report's "Hierarchical resolution-aware
deconvolution" section opens with a **Trusted resolution by family** table (from
the Resolution Decision Layer) shown **before** any fine predictions: each family
is `full_fine`, `selected_fine`, or `broad_only`, decided from full-panel
deconvolution reliability and signature/query evidence (not cell-level AUROC).
`broad_only` families are flagged so their fine split reads as **diagnostic only**,
and the block restates that values are RNA-derived proportions, not cell fractions.
The section then includes a **Within-family gating mode** block that names
the active gating mode and states that **soft** gating is the default (partial
confidence-weighted, validated on breast + lung), **hard** is legacy (binary
threshold; over-abstains in collinear families), and **ungated** is
diagnostic-only. It also reminds readers that a high cell-classification AUROC is
not proof of deconvolution reliability and that values are RNA-derived
proportions, not cell fractions, and embeds the standard gating methods note.

The same block carries a **spillover & false-positive caution** (collinear-family
conditional estimates risk leaking mass onto the wrong subtype and detecting absent
subtypes — fine predictions are gated by trusted resolution) and a note that the
experimental **FineGranularityRefiner** is a documented negative result (it
increased spillover and false positives on breast + lung) and is **not applied**.
The report never presents fine refinement as a promoted/active feature.

## Figures and source data

Report figures are Plotly-based. When `kaleido` is installed, the package
also writes static exports in PDF, SVG, and PNG formats.

Every figure includes a `.data.tsv` file containing the underlying data used
for the plot.

## Run metadata and analysis plan

A report built from a results directory reads `run_metadata.json` if present.
This file contains resolved parameters and provenance.

If the run created `analysis_plan.json`, the report also surfaces the
analysis plan and the selected preset.

## Output structure

A typical report directory contains:

- `tables/` — data tables and QC outputs
- `figures/` — Plotly figures and static exports
- `run_metadata.json` — run provenance
- `analysis_plan.json` — selected mode and preset
- `warnings.json` — warnings and issues
- `methods.txt` — methods text for the report
