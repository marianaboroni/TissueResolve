# Installing external benchmark tools

External tools are **optional** and live in **separate environments** so they
cannot break the main TissueResolve `.venv`.

## R tools (MuSiC, BayesPrism, BisqueRNA, CARD, SPOTlight)

```bash
conda env create -f benchmarks/envs/benchmark_r.yml
conda activate tissueresolve-bench-r
bash benchmarks/envs/install_external_tools.sh   # installs the R packages, logs each
```

## Python spatial tools (cell2location, optional tangram/stereoscope)

```bash
conda env create -f benchmarks/envs/benchmark_python_cell2location.yml
conda activate tissueresolve-bench-c2l
```

## Status

`benchmarks/run_real_external_benchmark.py --install-tools` runs the installer
and writes `benchmarks/outputs/tool_installation_status.tsv`. Tools that fail to
install are recorded with their error and an install hint — they are **not**
counted as benchmarked. If you run a tool elsewhere, import its predictions with
`benchmarks/import_external_results.py` (status = imported).

> Installing Bioconductor/GitHub R packages and cell2location requires network
> access and can take a long time; cell2location benefits from a GPU. None of
> this is required for the internal benchmark, which always runs offline.
