# TissueResolve

**TissueResolve:** spillover-aware deconvolution for bulk RNA-seq and 10x Visium
spatial transcriptomics.

TissueResolve estimates RNA-derived cell-type and cell-state composition using a
shared single-cell reference layer for both bulk and spatial workflows.
It combines protocol-aware gene selection, robust QC, separability diagnostics,
spillover-aware interpretation, and publication-ready HTML reports.

## What TissueResolve does

- Estimates RNA-derived composition from bulk RNA-seq and 10x Visium data.
- Builds a reference from single-cell `.h5ad` files or saved reference
  directories.
- Supports bulk deconvolution with protocol-aware weighting, bootstrap
  uncertainty, and compatibility checks.
- Supports spatial deconvolution on Visium data using a graph-aware NB-CAR
  model.
- Reports gene overlap, separability, spillover, warnings, and methods text.
- Generates publication-style reports with figures and source data.

## Main differentiators

- Unified bulk + spatial workflow with the same reference abstraction.
- Protocol-aware reference and query compatibility checks.
- Explicit bulk mRNA-proportion output warnings.
- Spatial graph-aware Visium modeling and neighborhood statistics.
- Family-aware handling of non-separable cell types.
- Plotly HTML reports with saved figure source data.
- Real-data breast cancer validation harness kept separate from the default
  offline test suite.

## Installation

```bash
git clone https://github.com/marianaboroni/TissueResolve.git
cd TissueResolve
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[all]"
python -m pytest -q
```

For finer dependency control:

```bash
python -m pip install -e ".[spatial,report,realdata]"
```

## Quickstart

### Bulk

```bash
tissueresolve run --reference reference.h5ad --query bulk_counts.tsv \
  --out results/bulk --mode bulk --preset standard
```

### Spatial

```bash
tissueresolve run --reference reference.h5ad --query visium.h5ad \
  --out results/spatial --mode spatial --preset publication
```

Or use the spatial subcommand:

```bash
tissueresolve spatial run --visium visium.h5ad \
  --reference reference.h5ad --output results/spatial
```

### Generate reports

```bash
tissueresolve report --modality bulk --results-dir results/bulk

tissueresolve report --modality spatial --results-dir results/spatial
```

## Supported CLI commands

- `tissueresolve run`
- `tissueresolve spatial run`
- `tissueresolve report`
- `tissueresolve bulk report`
- `tissueresolve spatial report`

> Note: `tissueresolve bulk run` is present in the CLI tree but not yet
> implemented; use `tissueresolve run --mode bulk` instead.

## Outputs

Runs produce:

- `tables/` — result tables and QC outputs
- `figures/` — Plotly figures and static exports
- `report.html` — publication-style report
- `run_metadata.json` — run provenance
- `analysis_plan.json` — selected mode and preset
- `warnings.json` — warnings and issues
- `methods.txt` — methods text for reports

## Reports and figures

The reporting layer embeds:

- executive summary cards
- key findings and interpretation
- a main publication figure
- detailed collapsible outputs
- warnings and QC notes

Figure outputs include Plotly HTML and, when `kaleido` is installed,
PDF/SVG/PNG static exports. Every figure writes a `.data.tsv` source data file.

## Interpretation

- Bulk outputs are RNA-derived **mRNA proportions**, not absolute cell counts.
- Spatial outputs are spot-level RNA-derived composition estimates, not direct
  single-cell counts.
- Non-separable cell-type pairs may be unreliable; family-level interpretation
  is provided when appropriate.
- Bootstrap uncertainty is available in publication and diagnostic modes.

## Real-data breast cancer validation

A dedicated validation harness lives in
`examples/real_breast_cancer/` and is designed to run only when explicitly
requested.

Scripts include data download, reference building, pseudobulk generation,
bulk validation, spatial validation, and report generation.

## Documentation

See the documentation pages in `docs/`:

- `docs/tutorial.md`
- `docs/quickstart.md`
- `docs/input_formats.md`
- `docs/advanced_parameters.md`
- `docs/reporting.md`
- `docs/output_interpretation.md`

## Status

Research software / pre-release.

## License

License metadata is declared as BSD-3-Clause in `pyproject.toml`.
A top-level `LICENSE` file is not currently present in this repository.
