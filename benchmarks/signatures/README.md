# Multipanel signature benchmark (Stage 1.5)

Permanent, leakage-guarded, donor-disjoint benchmark answering: **does TissueResolve
benefit from the multipanel architecture, or only from the fine_global donor-aware panel?**

- `validation.py` — guards that ABORT on hierarchy/dataset mismatch (incl. the wrong-tissue
  regression), donor overlap, selection-donor leakage, unrecorded budget, mis-orientation.
- `run_multipanel_benchmark.py` — real-data-gated runner (`--run-real-data`); equal STRATIFIED
  gene budget across strategies; Poisson solver; paired-by-mixture stats. Outputs (gitignored)
  under `benchmarks/results/signatures/`.
- `tests/` — offline guard tests (run in the default suite).

Run: `PYTHONPATH=src:. python benchmarks/signatures/run_multipanel_benchmark.py --run-real-data`

Strategies A (current) and B (fine_global) are implemented; C/D/E are the next increment.
