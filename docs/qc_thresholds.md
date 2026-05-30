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

## Using your own thresholds

All thresholds live in `TissueResolveConfig` (`bulk_qc`, `spatial_qc`) and are
serialisable to YAML (`config.to_yaml(...)`). Override them per run rather than
treating the defaults as ground truth.
