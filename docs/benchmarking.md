# Benchmarking

TissueResolve ships an **optional** benchmark harness under `benchmarks/`. It
compares TissueResolve against internal baselines and (when installed) external
deconvolution tools, and evaluates robustness to normalization, protocol,
library type, batch effects, and reference composition.

The harness is offline-first: it runs with only the core dependencies using
internal baselines, and **skips any external tool that is not installed**,
recording the reason and an install hint. One missing tool never fails the run.

## Commands

```bash
# show the plan (which methods would run / be skipped) — no computation
python benchmarks/bulk/run_bulk_benchmark.py --dry-run
python benchmarks/spatial/run_spatial_benchmark.py --dry-run

# toy synthetic data with known ground truth (fast, fully offline)
python benchmarks/bulk/run_bulk_benchmark.py --toy
python benchmarks/spatial/run_spatial_benchmark.py --toy

# existing real breast-cancer harness outputs (no download)
python benchmarks/bulk/run_bulk_benchmark.py --use-existing-real-data
python benchmarks/spatial/run_spatial_benchmark.py --use-existing-real-data --max-spots 600

# run both and build one linked summary report
python benchmarks/run_all.py --toy
```

Add `--no-external` to consider only internal methods.

## Methods

**Bulk** — TissueResolve (flat + hierarchical), internal NNLS baseline; optional
external: MuSiC, Bisque, DWLS, CIBERSORTx (input export only; the web tool is
never run automatically), and BayesPrism-compatible export.

**Spatial** — TissueResolve (flat + hierarchical), internal spot-level NNLS
baseline; optional external: RCTD/spacexr, cell2location, stereoscope,
SPOTlight, Tangram, DestVI, CARD.

Each method declares its capabilities (normalization, protocol, reference type,
hierarchical support) and returns a status of `success`, `skipped`, `failed`,
or `exported_not_run`.

## What is measured

- **Bulk (ground truth):** Pearson, Spearman, RMSE, MAE, signed bias,
  per-cell-type RMSE, dominant-type accuracy, calibration, family- vs fine-level
  metrics, runtime, gene overlap, warnings, flat-vs-hierarchical.
- **Real spatial (no ground truth):** cross-method concordance, pairwise
  correlation, Moran's I, entropy, near-zero fraction, runtime, failure modes.
  **No absolute accuracy is claimed without ground truth.**
- **Synthetic spatial (ground truth):** accuracy metrics as for bulk.
- **Robustness:** normalization decisions, protocol/library detection and
  compatibility, scRNA vs snRNA vs mixed reference, batch confounding and
  batch-aware vs naive marker stability.

## Outputs

Per modality under `benchmarks/outputs/<bulk|spatial>/`: `method_status.tsv`,
`accuracy_metrics.tsv` (when ground truth), `method_concordance.tsv` and
`spatial_structure.tsv` (spatial, no ground truth), `method_compatibility.tsv`,
`benchmark_metadata.json`, figures (each with a `.data.tsv`), and an HTML
report. `benchmarks/run_all.py` writes a single
`benchmark_summary_report.html` linking to both. Benchmark outputs are
git-ignored and never committed.

## Interpretation caveats

- Accuracy claims apply only where ground truth exists (bulk pseudobulk,
  synthetic spatial). For real Visium, the report is about concordance,
  spatial structure, stability, runtime, and failure modes.
- Batch-aware marker selection can improve robustness, but aggressive batch
  correction may remove real biological signal — the default is *diagnose, not
  correct*.
- scRNA and snRNA references differ systematically (gene detection, RNA
  composition); mixed references require library-type confounding checks.
- TissueResolve is not claimed to be universally more accurate; the goal is to
  make uncertainty, spillover, protocol/normalization/library/batch effects, and
  resolution limits explicit.

## Real external-tool benchmark

The internal benchmark compares TissueResolve against internal baselines. To
compare against **published external tools** (bulk: MuSiC, BayesPrism, BisqueRNA;
spatial: cell2location, CARD, SPOTlight), use the real external benchmark.

External tools live in **separate environments** (see
`benchmarks/envs/benchmark_installation.md`) so they cannot break the main
`.venv`. They are included **only when installed or when their results are
imported** — tools that are merely *exported* or *skipped* are never counted as
benchmarked or scored.

```bash
python benchmarks/run_real_external_benchmark.py --dry-run            # plan + availability
python benchmarks/run_real_external_benchmark.py --prepare-inputs --use-existing-real-data
python benchmarks/run_real_external_benchmark.py --install-tools      # attempt installs (slow; logs)
python benchmarks/run_real_external_benchmark.py --run-all --use-existing-real-data --fast
```

Import results produced elsewhere (counted as *executed (imported)*):

```bash
python benchmarks/import_external_results.py --method cell2location --modality spatial \
    --predictions preds.tsv --spot-id-col spot
```

Outputs: `tool_installation_status.tsv`, `composite_scores.tsv`,
`real_external_method_status.tsv`, and `real_external_benchmark_report.html`
(executive summary, method status, composite score, best-tool-by-scenario,
metric explanations, limitations). The **composite score** weights accuracy,
robustness, usability, interpretability, resolution-awareness, and
runtime/resource (config: `benchmarks/configs/composite_score_weights.yaml`);
only executed/imported tools are scored, and real-Visium runs exclude accuracy
(no ground truth).
