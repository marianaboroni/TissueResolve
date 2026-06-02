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

## What the report includes

- Executive summary cards
- Key findings and interpretation
- Main publication figure
- QC diagnostics and warnings
- Detailed tables in collapsible sections
- File list for reproducibility

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
