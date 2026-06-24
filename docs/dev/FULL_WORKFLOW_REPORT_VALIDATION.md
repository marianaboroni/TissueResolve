# Full-Workflow & QC-First Report Validation (post Resolution Decision Layer)

End-to-end validation that the reorganized TissueResolve pipeline runs and that the
QC-first report communicates reference quality, trusted resolution, unresolved mass,
and modality-specific behaviour. Driver:
`benchmarks/diagnostics/full_workflow_validation.py` (donor- and seed-disjoint;
query/test donors never used for gene selection; **no model defaults changed**;
spatial λ unchanged at 0.1). Outputs (git-ignored) under
`benchmarks/outputs/full_workflow_{breast,lung}/`.

## Task 1 — pipeline runs end-to-end (both tissues)

Order exercised: reference QC → query compatibility → **trusted-resolution
decision** → modality-aware weighting → broad deconvolution → fine deconvolution →
soft gating → QC-first report.

| | breast | lung (HLCA) |
|---|---|---|
| reference donors / query donors | 20 / 21 | 20 / 20 |
| fine types / families | 32 / 8 | 51 / 10 |
| reference suitability | 0.627 (WARNING) | 0.740 (WARNING) |
| modality / prediction unit | bulk / sample | bulk / sample |
| gene-weighting mode | protocol_aware_weighted_nnls | protocol_aware_weighted_nnls |
| hierarchical gating | soft | soft |
| trusted resolution (full/selected/broad_only) | 2 / 0 / 6 | 1 / 0 / 9 |
| mass-conservation error | 2.2e-16 | 2.2e-16 |
| unresolved-mass fraction | 0.829 | 0.810 |
| broad Pearson / RMSE | 0.816 / 0.087 | 0.683 / 0.104 |
| fine Pearson / RMSE | 0.172 / 0.063 | 0.154 / 0.047 |
| report generated | yes | yes |

The decision layer behaves as designed and consistently with all prior evidence:
the broad level is recovered reasonably (Pearson 0.68–0.82), while the
collinear/shared-lineage families are correctly classified **broad_only** (6/8
breast, 9/10 lung), so their unreliable fine splits are routed to unresolved
(≈0.81–0.83 of mass) and marked diagnostic rather than reported as trusted fine
composition. The `full_fine` families here are the single-subtype families
(Adipocyte, Stromal/Fibroblast; Mesothelium) where fine == family. Mass is
conserved exactly. This is the honest, conservative outcome the reorganization was
built to produce — it does not manufacture fine precision the data cannot support.

Required outputs written per tissue: `run_metadata.json`, `predictions.tsv`,
`unresolved_mass.tsv`, `trusted_resolution.tsv`, `qc.tsv`, `bulk_metrics.tsv`,
`reference_suitability.tsv`, `warnings.json`, `summary.json`, `report.html`
(+ figure source data emitted by the report builder).

## Task 2 — QC-first report content (verified programmatically, both tissues)

All checks PASS on both `report.html` files:

| check | breast | lung |
|---|---|---|
| "Trusted resolution by family" section present | ✅ | ✅ |
| trusted-resolution table appears **before** the fine/broad composition tables | ✅ | ✅ |
| broad_only families flagged **diagnostic only** | ✅ | ✅ |
| soft-gating status shown | ✅ | ✅ |
| unresolved mass shown | ✅ | ✅ |
| cell-level AUROC caveat present | ✅ | ✅ |
| RNA-derived proportions, **not** cell fractions | ✅ | ✅ |
| spillover / false-positive caution present | ✅ | ✅ |

Reference quality (suitability score + classification) and per-family
counts/resolvability are written to `reference_suitability.tsv` / `qc.tsv` and the
report's reference + hierarchical sections. Technical matrices (per-family
resolvability table, raw proportion tables) appear in the report's detailed/lower
sections, after the quality and trusted-resolution sections.

Automated regression coverage lives in
`tests/test_resolution_decision_layer.py` (trusted-resolution table precedes fine
predictions; broad_only → not trusted fine; AUROC caveat; RNA-derived not cell
fractions; soft-gating status; unresolved shown; bulk/spatial estimate-type
statements). Bulk and spatial reports are produced by separate generators
(`generate_bulk_report` / `generate_spatial_report`) with modality recorded in run
metadata (`modality`, `prediction_unit`, `gene_weighting_mode`,
`spatial_smoothing_used`, `coordinates_used`, `h_and_e_available`); the spatial
estimate-type statement makes no accuracy claim and never calls values cell counts.

## Limitations

* **Reference-suitability WARNING** on both subsampled atlases (breast 0.627, lung
  0.740) — surfaced honestly in the report; downstream numbers are relative, not
  absolute-quality guarantees.
* **Abstention trade-off (by design).** The conservative resolution decision routes
  most collinear-family mass to unresolved (≈0.81–0.83 here), so broad-level recovery
  is good (Pearson 0.68–0.82) while raw fine correlation is low — TissueResolve trades
  fine recall for trustworthy, low-false-positive resolution. This is the intended
  behaviour, not a defect; see `EXTERNAL_TOOL_BENCHMARK_REPORT.md` §12 for the powered
  precision-vs-recall analysis.
* **Heuristic thresholds.** Resolution-decision/gating cut-offs are heuristic defaults
  (see `docs/qc_thresholds.md`), tunable per dataset.
* **Sparsity diagnostics.** Raw NNLS coefficients before L1-row normalisation are saved
  (`raw_weights_pre_l1norm.tsv`, `raw_weight_sparsity.tsv`) via the solver's opt-in
  `return_raw_weights` flag (default off; does not change any output).

This validation covers the **bulk pseudobulk** path on both tissues (the path with
known ground truth). The synthetic-spatial path reuses the same reference,
hierarchy, resolution decision, soft gating, and report machinery (validated
separately in the spatial benchmark scripts); a fresh full synthetic-spatial run
was not repeated here to avoid the long multi-panel runtimes. Reference suitability
classifies both atlases as WARNING (driven by gene-overlap/marker-stability on
these subsampled atlases), which the report surfaces honestly.
