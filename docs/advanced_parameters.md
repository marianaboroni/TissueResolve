# Advanced Parameters

## Presets

TissueResolve supports named presets that tune plotting, bootstrap, and
spatial defaults.

Available presets:

- `quick`
- `standard`
- `publication`
- `diagnostic`

Use a preset to avoid manually tuning many options.

## Top-level run options

```bash
tissueresolve run --reference reference.h5ad --query query --out results \
  --mode bulk --preset standard --dry-run
```

- `--reference` — single-cell reference `.h5ad` or saved reference directory.
- `--query` — bulk counts table or Visium `.h5ad` or folder.
- `--out` — output directory for the run bundle.
- `--mode` — `auto`, `bulk`, or `spatial`.
- `--preset` — `quick`, `standard`, `publication`, or `diagnostic`.
- `--dry-run` — write the analysis plan without executing the pipeline.

## Spatial run options

```bash
tissueresolve spatial run --visium visium.h5ad --reference reference.h5ad \
  --output results/spatial --cell-type-col cell_type --lambda-spatial 0.1
```

- `--visium` — Visium `.h5ad` file or Space Ranger output folder.
- `--reference` — saved reference directory or `.h5ad` reference.
- `--output` — output directory for spatial results.
- `--config` — optional YAML config file.
- `--cell-type-col` — cell-type column name in reference `obs`.
- `--marker-genes` — optional marker gene list file.
- `--lambda-spatial` — CAR smoothing strength.
- `--max-iter` — max solver iterations.
- `--random-state` — random seed.
- `--min-counts` — minimum UMI per spot.
- `--min-genes` — minimum detected genes per spot.
- `--neighbourhood / --no-neighbourhood` — compute neighborhood statistics.
- `--genome` — genome assembly (`hg38` or `mm10`).

## Reporting options

```bash
tissueresolve report --modality bulk --results-dir results/bulk --out report.html
```

- `--modality` — bulk or spatial report.
- `--results-dir` — folder containing `tables/` and `figures/`.
- `--out` — optional path for the HTML report.

## When to change these values

- Use `quick` for fast exploratory runs.
- Use `standard` for routine analysis.
- Use `publication` for figures and uncertainty output.
- Use `diagnostic` for deeper bootstrap and tuning.
- Override `--lambda-spatial` when spatial smoothing needs manual control.
- Use `--marker-genes` when you have a validated marker set.
- Increase `--max-iter` if the spatial solver does not converge.

## Beginner recommendations

- Start with `--preset standard`.
- Do not override `--lambda-spatial` unless you understand spatial smoothing.
- Keep `--marker-genes` only if you know the marker set is reliable.
- Use `--dry-run` to inspect `analysis_plan.json` before a full run.

## Hierarchical (broad → fine) parameters

Enable with `--resolution-mode hierarchical`. Configurable via the CLI and the
`HierarchicalConfig` section of `TissueResolveConfig`:

| Parameter | CLI flag | Default | Meaning |
|---|---|---|---|
| `broad_cell_type_col` | `--broad-cell-type-col` | `auto` | obs column with broad/compartment labels (`auto` detects a known candidate). |
| `fine_cell_type_col` | `--fine-cell-type-col` | `auto` | obs column with fine/subpopulation labels. |
| (mapping file) | `--cell-type-hierarchy` | — | fine→broad TSV when the reference has only fine labels. |
| `allow_unresolved` | `--allow-unresolved` / `--no-allow-unresolved` | `True` | Keep non-separable families at the broad level as `unresolved_<family>`. |
| `unresolved_threshold` | — | `0.10` | A family is unresolved when mean within-family separability (`1 − Bhattacharyya`) is below this. |
| `min_discriminating_genes` | — | `10` | A family is unresolved when its worst within-family pair has fewer discriminating genes than this. |
| `within_family_spillover_threshold` | — | `0.30` | A family is unresolved when mean within-family spillover (max correlation to a sibling) exceeds this. |
| `within_family_marker_selection` | — | `auto` | `auto` / `all` / `pairwise` within-family marker strategy. |
| `hierarchy_level` | — | `both` | Which level (`fine` / `family` / `both`) to emphasise in outputs. |

A family is split into subtypes only when **all three** gating signals agree it
is separable (separability, discriminating genes, spillover). This defence in
depth prevents a single noisy metric from forcing or blocking a split. All
thresholds are heuristic and recorded in `run_metadata.json`.

Presets adjust the gating strictness: `publication` is stricter
(`min_discriminating_genes=30`, `within_family_spillover_threshold=0.25`),
while `quick`/`diagnostic` are more permissive.
