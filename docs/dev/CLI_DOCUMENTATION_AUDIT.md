# CLI Documentation Audit

Audit date: 2026-06-02. Scope: every CLI command documented in `README.md` and
`docs/` compared against the actual CLI (`src/tissueresolve/cli.py`) and the
benchmark scripts under `benchmarks/`. Captured `--help` for each command.

## Captured commands (all exist and respond to `--help`)

```
tissueresolve --help          → groups: bulk, info, report, run, spatial      (exit 0)
tissueresolve run --help      → top-level run                                  (exit 0)
tissueresolve bulk --help     → benchmark, check-compatibility, report, run    (exit 0)
tissueresolve spatial --help  → benchmark, info, report, run                   (exit 0)
python3 benchmarks/run_all.py --help                  (exit 0)
python3 benchmarks/import_external_results.py --help  (exit 0)
```

## Per-command audit table

| Documented command | Exists | Runs `--help` | Flags match docs | Output path documented | Status | Recommended fix |
|---|---|---|---|---|---|---|
| `tissueresolve run` | yes | yes | yes | partial | **mismatch** | Output claim overstated — see Critical check 2; correct README/tutorial "Outputs". |
| `tissueresolve run --mode bulk` | yes | yes | yes | partial | OK | none (mode choice `auto/bulk/spatial` present). |
| `tissueresolve run --mode spatial` | yes | yes | yes | partial | OK | none. |
| `tissueresolve spatial run` | yes | yes | yes | yes | OK | none (`--visium/--reference/--output/...` match `advanced_parameters.md`). |
| `tissueresolve spatial benchmark` | yes | yes | n/a (not in user docs) | yes | OK | synthetic Visium benchmark; not referenced by README workflow. |
| `tissueresolve spatial info` | yes | yes | n/a | n/a | OK | none. |
| `tissueresolve report --modality bulk/spatial` | yes | yes | yes | yes | OK | none (`--results-dir`, `--out`). |
| `tissueresolve bulk report` | yes | yes | yes | yes | OK | none. |
| `tissueresolve spatial report` | yes | yes | yes | yes | OK | none. |
| `tissueresolve bulk run` | yes (stub) | yes | n/a | n/a | **deprecated/stub** | Exits 2 "not yet implemented". README already steers users to `run --mode bulk`. Keep note; **tutorial example using it is wrong** (see below). |
| `tissueresolve bulk check-compatibility` | yes (stub) | yes | n/a | n/a | **stub** | Exits 2 "not yet implemented". Not advertised in README. OK as long as not documented as working. |
| `tissueresolve bulk benchmark` | yes (stub) | yes | n/a | n/a | **stub** | Exits 2. Benchmarking docs use `benchmarks/bulk/run_bulk_benchmark.py` (separate script), not this subcommand. OK. |
| `tissueresolve info` | yes | yes | n/a | n/a | OK | none. |
| `python benchmarks/run_all.py` | yes | yes | **no** | yes | **mismatch** | `--include-imported` documented (README:268) but NOT a `run_all.py` flag → command errors. Fixed by adding/forwarding the flag. |
| `python benchmarks/import_external_results.py` | yes | yes | yes | yes | OK | help text references the broken `run_all.py --include-imported` (fixed together). |
| `python benchmarks/bulk/run_bulk_benchmark.py` | yes | yes | yes | yes | OK | `--dry-run/--toy/--use-existing-real-data/--no-external/--include-imported`. |
| `python benchmarks/spatial/run_spatial_benchmark.py` | yes | yes | yes | yes | OK | same + `--max-spots/--seed/--spatial-auto`. |
| `python benchmarks/run_real_external_benchmark.py` | yes | yes | yes | yes | OK | `--dry-run/--prepare-inputs/--install-tools/--run-all/...`. |

## Critical checks

1. **Does `tissueresolve run --mode bulk` exist?** — **Yes.** `--mode` is
   `click.Choice(["auto","bulk","spatial"])`. Verified end-to-end (a synthetic
   bulk run completes and writes results).

2. **Does the documented bulk workflow match the current CLI?** — **Partly.**
   The *command* (`tissueresolve run --mode bulk ...`) is correct. The
   *documented outputs* are overstated: README "Outputs" and tutorial §6/§9 list
   `tables/`, `figures/`, `report.html`, `warnings.json`, `methods.txt` as
   products of a run. A real `tissueresolve run` writes only
   `analysis_plan.json`, `run_metadata.json`, `deconv/` and `qc/`. `report.html`
   is produced by the *separate* `tissueresolve report` step; the full
   `tables/`+`figures/`+`methods.txt`+`warnings.json` bundle is produced by the
   validation-harness script `07_generate_reports.py`, not by `run`. **Fix:**
   correct the documentation (Part 5/8) rather than expand `run`.

3. **Does the documented spatial workflow match the current CLI?** — **Yes** for
   the commands. `tissueresolve run --mode spatial` and `tissueresolve spatial
   run --visium ... --reference ... --output ...` both exist with the documented
   flags and run end-to-end on a synthetic Visium fixture. Same output-claim
   caveat as bulk applies to the `run` path.

4. **Are `--state-aware`, `--solver auto`, `--resolution-mode`, and benchmark
   flags documented accurately?**
   - `--state-aware` — present, **labelled EXPERIMENTAL** in CLI help, README,
     and `FEATURE_STATUS.md`; requires `--resolution-mode hierarchical`. ✅
   - `--solver auto` (+ `nnls/weighted_nnls/marker_nnls/ridge_nnls/ensemble_nnls/pipeline`)
     — present and matches `advanced_parameters.md`. ✅
   - `--resolution-mode auto/hierarchical/flat/none/suggest` — present and
     matches README/tutorial/advanced_parameters. ✅
   - Benchmark flags — **one mismatch**: `run_all.py --include-imported` is
     documented but the flag lived only on the per-modality scripts. Fixed.

5. **Are experimental flags clearly labelled experimental?** — **Yes.**
   `--state-aware` help begins "EXPERIMENTAL:" and states it is not part of the
   default v0.1 workflow and not validated across real datasets. Covered by
   `tests/test_v0_1_stabilization.py::test_cli_help_marks_state_aware_experimental`.

6. **Are any README examples obsolete?**
   - README: no obsolete *commands* (the `bulk run` note correctly steers to
     `run --mode bulk`). The `run_all.py --include-imported` example was broken
     (fixed).
   - `docs/tutorial.md` §"Bulk, flat" used `tissueresolve bulk run --bulk ...
     --cell-type-col ... --resolution-mode none --out ...` — **obsolete/broken**:
     `bulk run` is an unimplemented stub and accepts none of those flags. Fixed
     to the working `tissueresolve run --mode bulk ... --resolution-mode flat`.

## Mismatches found

| # | Severity | Location | Problem | Fix applied |
|---|---|---|---|---|
| 1 | HIGH | `docs/tutorial.md` "Bulk, flat" example | `tissueresolve bulk run` is a stub; flags don't exist → command fails | Rewrote to `tissueresolve run --mode bulk ... --resolution-mode flat`. |
| 2 | HIGH | `README.md:268`, `benchmarks/import_external_results.py:17,82`, `docs/benchmarking.md` | `run_all.py --include-imported` errors (`unrecognized arguments`) | Added `--include-imported` to `run_all.py` and forwarded it to both sub-benchmarks. |
| 3 | MED | `README.md` Outputs, `docs/tutorial.md` §6/§9, `docs/reporting.md` Output structure | `tissueresolve run` does not emit `tables/`, `figures/`, `report.html`, `warnings.json`, `methods.txt` | Corrected the documentation to describe actual `run` outputs and the separate report step. |
| 4 | LOW | `docs/output_interpretation.md` "run with `--n-bootstrap > 0`" | No `--n-bootstrap` CLI flag exists on `run` (bootstrap is preset-driven) | Documentation note; bootstrap is controlled by `--preset publication/diagnostic`. |
</content>
</invoke>
