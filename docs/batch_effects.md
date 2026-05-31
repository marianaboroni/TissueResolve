# Batch-effect diagnostics

Reference single-cell/single-nucleus data often carry batch/donor/library
structure that can confound cell-type signatures. TissueResolve's philosophy is
**diagnose, don't blindly correct**: aggressive batch correction can remove real
cell-type signal, so the default surfaces the problem and scores marker
stability rather than altering the expression matrix.

## What is diagnosed

`benchmarks/shared/batch_effects.py` provides:

- `detect_batch_columns` — auto-detect `batch`/`donor`/`sample`/`library`/… columns.
- `compute_celltype_batch_confounding` — per-cell-type dominant-batch fraction,
  number of batches, and a `single_batch_only` / `confounded` flag.
- `compute_batch_mixing_score` — Shannon evenness of batch sizes.
- `compute_marker_batch_stability` — per-cell-type cross-donor marker CV
  (lower = more stable); reports when donor CV is unavailable rather than
  assuming stability.
- `batch_aware_reference_summary` — a compact summary for the report.

## Options (benchmark / reference building)

```
--batch-col auto|COLUMN|none
--donor-col auto|COLUMN|none
--library-type-col auto|COLUMN|none
--batch-aware-markers / --no-batch-aware-markers
--batch-correction none|diagnostic|marker_stability|combat|harmony|scanvi_export
```

Defaults: `--batch-correction diagnostic`, and batch-aware marker scoring is on
when donor/batch/library metadata exist. Raw expression is not modified unless a
correction mode is explicitly requested.

The benchmark compares batch-aware vs batch-naive references and writes
`batch_diagnostics.tsv`, `batch_marker_stability.tsv`,
`batch_aware_vs_naive_metrics.tsv`, and `batch_effect_summary.md`.
