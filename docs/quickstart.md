# Quickstart

## Installation

```bash
git clone https://github.com/marianaboroni/TissueResolve.git
cd TissueResolve
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[all]"
```

## Run a bulk analysis

```bash
tissueresolve run --reference reference.h5ad --query bulk_counts.tsv \
  --out results/bulk --mode bulk --preset standard
```

## Run a spatial analysis

```bash
tissueresolve run --reference reference.h5ad --query visium.h5ad \
  --out results/spatial --mode spatial --preset publication
```

Or use the spatial subcommand:

```bash
tissueresolve spatial run --visium visium.h5ad \
  --reference reference.h5ad --output results/spatial
```

## Generate reports

```bash
tissueresolve report --modality bulk --results-dir results/bulk

tissueresolve report --modality spatial --results-dir results/spatial
```

## Verify installation

```bash
python -m pytest -q
```

## Notes

- Bulk outputs are RNA-derived mRNA proportions, not absolute cell counts.
- Spatial outputs are spot-level RNA-derived composition estimates.
- `tissueresolve bulk run` is present in the CLI tree but is not yet implemented;
  use `tissueresolve run --mode bulk` instead.
