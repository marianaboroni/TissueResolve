# Reporting

TissueResolve can generate HTML reports from a results directory.

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
