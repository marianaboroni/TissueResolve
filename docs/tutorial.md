# TissueResolve Tutorial

## 1. Overview

TissueResolve is a unified deconvolution toolkit for:

- bulk RNA-seq deconvolution using protocol-aware weighted NNLS,
- 10x Visium spatial deconvolution using a negative-binomial CAR model.

It uses a shared single-cell reference abstraction, reports QC and warnings,
provides separability and spillover diagnostics, and exports publication-ready
reports and figures.

## 2. Installation

```bash
git clone https://github.com/marianaboroni/TissueResolve.git
dcd TissueResolve
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[all]"
```

Run the test suite:

```bash
python -m pytest -q
```

## 3. Input formats

### Single-cell reference

TissueResolve accepts a single-cell or single-nucleus reference as a `.h5ad`
file. The file should contain:

- `obs` with a cell-type label column such as `cell_type`
- `var` with gene names or identifiers
- raw counts or normalized expression values

The package also accepts a saved `ReferenceSignature` directory.

### Bulk counts table

A bulk query can be a tab-delimited counts table with genes in rows and
samples in columns. The first column is used as gene names.

### Visium input

Spatial data can be passed as either:

- a Visium `.h5ad` file containing spatial coordinates, image paths, and
  count data, or
- a Space Ranger-style output folder containing `filtered_feature_bc_matrix`
  and a `spatial/` directory.

## 4. Running a quick test

A fast smoke test uses the top-level `run` command with toy inputs:

```bash
tissueresolve run --reference reference.h5ad --query bulk_counts.tsv \
  --out results/bulk --mode bulk --preset quick
```

A spatial quick test is similar:

```bash
tissueresolve run --reference reference.h5ad --query visium.h5ad \
  --out results/spatial --mode spatial --preset standard
```

## 5. Building a reference

TissueResolve builds a reference from a single-cell `.h5ad` automatically
when passed with `--reference`. If you already have a saved reference
directory, pass that path instead.

There is no separate CLI command for reference construction in this version.
Use the `.h5ad` reference directly or build it in Python with the API.

## 6. Bulk analysis step-by-step

1. Prepare a single-cell reference `.h5ad` with a valid cell-type column.
2. Prepare a bulk counts table with genes as row names.
3. Run the top-level command:

```bash
tissueresolve run --reference reference.h5ad --query bulk_counts.tsv \
  --out results/bulk --mode bulk --preset standard
```

4. Generate an HTML report:

```bash
tissueresolve report --modality bulk --results-dir results/bulk
```

### What this produces

- `results/bulk/tables/`
- `results/bulk/figures/`
- `results/bulk/report.html`
- `results/bulk/run_metadata.json`
- `results/bulk/analysis_plan.json`

## 7. Spatial analysis step-by-step

1. Prepare the same single-cell reference.
2. Provide a Visium input file or folder.
3. Run spatial deconvolution:

```bash
tissueresolve run --reference reference.h5ad --query visium.h5ad \
  --out results/spatial --mode spatial --preset publication
```

Alternatively:

```bash
tissueresolve spatial run --visium visium.h5ad \
  --reference reference.h5ad --output results/spatial
```

4. Generate a report:

```bash
tissueresolve report --modality spatial --results-dir results/spatial
```

## 8. Generating HTML reports

The root report command reads a results directory and writes a self-contained
HTML report.

```bash
tissueresolve report --modality bulk --results-dir results/bulk
```

The report includes:

- executive summary
- main publication figure
- interpretation text
- QC and warnings
- collapsible detailed outputs

## 9. Understanding output folders

A typical run writes:

- `tables/` — TSV summary tables, QC metrics, and diagnostic reports
- `figures/` — Plotly figures and static exports
- `report.html` — publication-style report
- `run_metadata.json` — resolved parameters and provenance
- `analysis_plan.json` — the planned modality and preset
- `warnings.json` — warnings raised during the run
- `methods.txt` — methods text for reports

## 10. Understanding figures

Figures may include:

- Bulk clustered composition barplots
- Bulk QC summary plots
- Spatial mean composition barplots
- Spatial dominant cell-type maps
- Separability and spillover summaries

![Bulk clustered composition](assets/bulk_clustered_composition.png)
![Bulk QC summary](assets/bulk_qc_summary.png)
![Spatial mean composition](assets/spatial_mean_composition.png)
![Spatial dominant cell type](assets/spatial_dominant_cell_type.png)
![Separability and spillover](assets/separability_spillover.png)

## 11. Understanding QC metrics

TissueResolve reports:

- reconstruction quality for bulk samples
- gene overlap between query and reference
- protocol compatibility risk
- spatial Moran's I and neighborhood statistics

Low overlap, poor separability, or non-convergence are surfaced as warnings.

## 12. Understanding separability and spillover

- Separability evaluates whether two cell types can be distinguished reliably.
- Spillover quantifies confusion between cell types from expression similarity.
- If cell types are confusable, the system produces family-level recommendations
  instead of silently merging or hiding uncertainty.

## 13. Tuning parameters

### General settings

- `--preset` — controls defaults for plotting, bootstrap, and tuning.
- `--dry-run` — write the analysis plan without running pipelines.

### Spatial settings

- `--config` — path to a YAML configuration file.
- `--cell-type-col` — reference cell-type column name (default `cell_type`).
- `--marker-genes` — optional marker gene list.
- `--lambda-spatial` — override CAR smoothing strength.
- `--max-iter` — override solver iterations.
- `--random-state` — RNG seed.
- `--min-counts` — minimum UMI per spot.
- `--min-genes` — minimum genes detected per spot.
- `--neighbourhood / --no-neighbourhood` — compute spot co-occurrence stats.
- `--genome` — reference genome assembly (`hg38` or `mm10`).

## 14. Common problems and solutions

- `low gene overlap` — check gene identifiers in bulk/spatial and reference.
- `missing cell-type column` — set `--cell-type-col` or rename the obs column.
- `non-converged spatial solver` — inspect `run_metadata.json` and try
  `--lambda-spatial` or `--max-iter`.
- `bootstrap not computed` — enable bootstrap in the preset or configuration.

## 15. Example: real breast cancer validation

Run the example scripts under `examples/real_breast_cancer/scripts/`:

```bash
python examples/real_breast_cancer/scripts/00_download_data.py
python examples/real_breast_cancer/scripts/01_prepare_reference.py
python examples/real_breast_cancer/scripts/02_make_pseudobulk.py
python examples/real_breast_cancer/scripts/03_run_bulk_validation.py
python examples/real_breast_cancer/scripts/04_run_spatial_validation.py
python examples/real_breast_cancer/scripts/05_summarize_results.py
python examples/real_breast_cancer/scripts/07_generate_reports.py
```

This harness is designed to keep real-data validation separate from the default
package workflows.
