# External-Tool Installation Notes (this environment)

Real, verified installation state (macOS / Darwin, Python 3.9, R 4.1.2). Versions
probed directly. No tool is silently skipped — every status has a reason.

## Installed and usable

| tool | modality | version | env | how |
|---|---|---|---|---|
| MuSiC | bulk | 1.0.0 | `benchmarks/envs/Rlib` | `R_LIBS=benchmarks/envs/Rlib Rscript` |
| BisqueRNA | bulk | 1.0.5 | `benchmarks/envs/Rlib` | Rscript bridge |
| BayesPrism | bulk | 2.2.3 | `benchmarks/envs/Rlib` | Rscript bridge (**very slow**: prior run ~2180 s, treated as failed:timeout) |
| CARD | spatial | 1.1 | `benchmarks/envs/Rlib` | Rscript bridge |
| spacexr (RCTD) | spatial | 2.2.1 | `benchmarks/envs/Rlib` | Rscript bridge (Seurat absent; core RCTD usable) |
| cell2location | spatial | 0.1.4 | `benchmarks/envs/c2l_py39` | scvi 1.1.6, torch 2.8.0, **MPS (Metal) GPU** available; no CUDA. ~160 s+/run |
| Rectangle (rectanglepy) | bulk | 1.5.0 | `benchmarks/envs/rectangle_py311` | **Python 3.11** (Homebrew python@3.11 3.11.15); project venv is 3.9 so Rectangle runs **out-of-process** via `benchmarks/dev/run_rectangle.py`, same pattern as cell2location. Needs single-cell reference (not a mean signature) → benchmarked in a dedicated dev script, not the signature-based product harness. Emits a scalar `Unknown` residual column. Install: `bash benchmarks/envs/install_rectangle.sh`. |

Internal baselines (always available, `.venv`): NNLS, WNNLS, ridge-NNLS, marker-NNLS,
TissueResolve flat, TissueResolve hierarchical + soft gating, NNLS-per-spot.

## Not installed / not runnable here (with reason)

| tool | modality | reason |
|---|---|---|
| SCDC | bulk | not installed (no CRAN/GitHub build in this env) |
| DWLS | bulk | not installed (prior run: skipped) |
| CIBERSORTx | bulk | web/token-gated academic license; cannot be auto-run. Inputs can be *exported* for manual upload only |
| SPOTlight | spatial | not installed (Bioconductor build absent) |
| Tangram | spatial | not installed (absent from the c2l env) |
| Seurat | (dep) | not installed; RCTD/spacexr core works without it but Seurat-based wrappers do not |

## Practical constraints observed

* **Runtime/stability**: BayesPrism and the spatial GPU/CPU tools (cell2location, RCTD,
  CARD) are minutes-to-tens-of-minutes per run; on this machine they have repeatedly
  pushed total runs into the multi-hour range and contributed to environment
  instability. They are attempted with **bounded per-tool timeouts**; a timeout is
  recorded as `failed:timeout`, not hidden.
* **GPU**: only Apple MPS is available (no CUDA); cell2location runs but slowly.
* **Reproducibility**: R tools are invoked via `Rscript` subprocesses with
  `R_LIBS=benchmarks/envs/Rlib`, exact versions captured in `tool_registry.tsv`.

## Install routes (documented, not re-run)

* R: `BiocManager::install(...)` / `remotes::install_github(...)` per
  `benchmarks/envs/install_external_tools.sh` (already used to populate `envs/Rlib`).
* Python: cell2location via `benchmarks/envs/benchmark_python_cell2location.yml`
  (already materialised in `envs/c2l_py39`).
