# QC thresholds and their provenance

TissueResolve's QC thresholds are **heuristic defaults**, not validated
universal cut-offs. They are tunable hyperparameters and should be revisited
per dataset. Every threshold-based flag is recorded as heuristic in the
output, and no result is suppressed because it crosses a threshold — flags add
warnings, they never remove data.

## Bulk QC (`BulkQCConfig`)

| Threshold | Default | Meaning | Status |
|---|---|---|---|
| `r2_warn` | 0.50 | Reconstruction R² below this warns of poor fit | Heuristic |
| `profile_corr_warn` | 0.60 | Bulk-vs-reference profile correlation warning | Heuristic |
| mismatch flag | low/medium/high | Per-sample protocol-mismatch severity | Heuristic |

These originate from internal calibration on protocol-mismatch benchmarks and
are consistent with the order of magnitude reported in bulk-deconvolution
literature (e.g. Finotello et al., 2019). They are **not** reproduced from a
single authoritative table; treat them as starting points.

## Spatial QC (`SpatialQCConfig`)

| Threshold | Default | Meaning | Status |
|---|---|---|---|
| `spatial_resid_threshold` | `None` (disabled) | Per-spot spatial residual flag | Heuristic, off by default |
| entropy / dominant-fraction | reported, not thresholded | Mixedness of a spot | Descriptive |

High spatial residual can reflect genuine biology (sharp tissue boundaries),
so it is **not** flagged unless the user opts in.

## Separability (`PairSeparability`)

| Bhattacharyya coefficient | Risk level |
|---|---|
| > 0.97 | CRITICAL |
| > 0.90 | HIGH |
| > 0.80 | MEDIUM |
| ≤ 0.80 | OK |

Pairs above the HIGH threshold are surfaced as warnings; estimates for those
cell types are reported as unreliable rather than confident.

## Resolvability (Stage 6)

`reference.resolution` classifies each cell-type pair by its separability score
`s = 1 − Bhattacharyya` (`ResolutionConfig`):

| Class | Separability score `s` | Bhattacharyya | Meaning |
|---|---|---|---|
| `resolved` | `s ≥ 0.20` | BC ≤ 0.80 | reliably distinguishable |
| `partially_resolved` | `0.10 ≤ s < 0.20` | 0.80 < BC ≤ 0.90 | usable with caution |
| `poorly_resolved` | `0.03 ≤ s < 0.10` | 0.90 < BC ≤ 0.97 | unreliable apart |
| `unresolved` | `s < 0.03` | BC > 0.97 | effectively indistinguishable |

Pairs with `BC ≥ family_bc_threshold` (default 0.90) are grouped into a
**family**. These thresholds align with the existing `PairSeparability` risk
levels and are tunable via `ResolutionConfig`.

## Spillover and the abstain (unresolved) mode

`benchmark.spillover` deconvolves pure single-type pseudobulks to build a
**spillover matrix** (row = true type, column = predicted type, rows sum to 1).
`spillover_risk = 1 − self_retention` (off-diagonal mass). Defaults
(`ResolutionConfig`): `spillover_threshold = 0.30`, `unresolved_threshold = 0.10`
(family mean separability), `uncertainty_threshold = 0.10` (bootstrap CI width).

`apply_unresolved_mode` collapses a multi-member family into a single
`unresolved_<family>` column **only when several signals agree** (low
separability, plus — when provided — high spillover risk and/or high bootstrap
uncertainty). It is configurable (`allow_unresolved`, default `True`) and
**preserves total mass**. These are conservative defaults, not validated
universal cut-offs.

## Resolution recommendation system (Stage 8)

When the reference has many poorly separable pairs, TissueResolve builds a
**non-separable graph** (edge when `BC ≥ family_bc_threshold`, default 0.90),
takes its connected components as candidate **merge families**, and names them
with label heuristics (e.g. T/NK lymphocytes, myeloid cells, endothelial cells,
epithelial cells, mural cells). It writes a machine-readable
`recommended_merges.tsv` and a `cell_type_families.tsv`. The separability
warning now points to this recommender and to `--resolution-mode` instead of
telling users to "merge manually". Merging is always explicit and recorded.

## Using your own thresholds

All thresholds live in `TissueResolveConfig` (`bulk_qc`, `spatial_qc`) and are
serialisable to YAML (`config.to_yaml(...)`). Override them per run rather than
treating the defaults as ground truth.
