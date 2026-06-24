# Report Stability Audit

Audit date: 2026-06-02. Scope: HTML report generation. Two report paths exist:

- **Canonical:** `report/unified.py` (+ `components.py`, `style`, `glossary`,
  `figures`) — the QC-first `report.html` + `technical_appendix.html`. Driven by
  the report layer / `07_generate_reports.py`.
- **Deprecated (still wired):** `report/html.py` (+ `sections`, `templates`,
  `assets`) — powers the per-modality `tissueresolve report --modality ...` and
  the `bulk/spatial report` subcommands. Kept until consolidation; do not extend
  (per `docs/reporting.md`).

## Checklist

| Check | Status | Evidence |
|---|---|---|
| `report.html` exists | ✅ | `tissueresolve report` wrote `report.html`; `generate_report` tested in `tests/report/test_results_dir_report.py`. |
| Report is QC-first | ✅ (canonical) | `tests/report/test_report_redesign_ordering.py::test_qc_sections_appear_before_results` enforces reference/signature/input QC before bulk/spatial. |
| Quality sections before prediction sections | ✅ | Same ordering test; deprecated path also renders "Single-cell reference quality" before "Deconvolution predictions". |
| Bulk and spatial predictions separate | ✅ | Distinct sections (canonical) and distinct per-modality report files. |
| Bulk and spatial benchmarks separate | ✅ | `run_all.py` writes separate `bulk/bulk_benchmark_report.html` and `spatial/spatial_benchmark_report.html`; benchmark docs keep them split. |
| Spatial benchmark does not claim accuracy without ground truth | ✅ | `report/glossary.py` ("concordance … high concordance ≠ correctness"); `docs/benchmarking.md`; `docs/output_interpretation.md` §"Spatial benchmark is not accuracy". |
| Warnings visible | ✅ | `test_not_no_warnings_when_separability_exists`; `converged=False` shown; my synthetic run's report rendered "Warnings & limitations (1)". |
| No empty figure cards | ✅ | Figures are embedded only when the figure file is present; missing artifacts render an explanatory card, not an empty plot (`test_missing_bootstrap_shows_card_not_plot`). |
| No generic captions | ✅ (canonical) | `figure_card` carries title/subtitle/caption/how-to-read/variables; estimate-type subtitle per figure. |
| Source-data links exist | ✅ | `figure_card(source_links=...)`; every figure co-writes a `.data.tsv`; `figures/figure_manifest.tsv` maps figure→data. |
| Technical figures in appendix/collapsible | ✅ | `technical_appendix.html` holds heatmaps/spillover networks/spot pies/signature heatmap; main report uses `<details>` collapsibles. |
| Not overloaded with raw tables | ✅ | Full tables live in collapsible "Detailed outputs"/appendix. |

## Issue found (not a report-internal defect)

**Prediction-empty report when run directly on a raw `tissueresolve run`
directory.** The results-dir report reader expects report-shaped inputs
(top-level `bulk_estimated_proportions.tsv` / `spatial_spot_proportions.tsv`,
`tables/`, `figures/`). A raw `run` writes predictions to
`deconv/proportions.tsv`, so `tissueresolve report --results-dir <run-dir>`
emits a structurally complete report whose **predictions section is empty**.
This is the same layout gap recorded in `END_TO_END_WORKFLOW_AUDIT.md`. It is an
end-to-end wiring gap, not a report-builder bug.

**Update:** `tissueresolve run` now renders `report.html` directly from the
in-memory result (populated predictions/QC/methods/warnings), so the primary
user path no longer yields a prediction-empty report. The layout mismatch only
remains if you point `tissueresolve report --results-dir` at a *raw* run
directory; regenerating from a report-shaped directory works as before.

## Fixes applied

None to the report builders — the report itself is QC-first, surfaces warnings,
separates bulk/spatial, gates figures on artifact presence, and carries the
no-ground-truth caveat. The audit confirms the report is stable; the only gap is
the upstream layout mismatch, addressed by documentation correction (Part 5/8)
and the deferred wiring item.

## Tests (existing + added)

Existing: `tests/report/test_report_redesign.py`,
`test_report_redesign_ordering.py`, `test_results_dir_report.py`,
`test_report.py`, `test_report_components.py`.

Added in `tests/test_report_stability.py` (Part 7): report.html generated; QC
sections before predictions; no empty figure cards (no `fig-body` with empty
contents); bulk vs spatial benchmark separation guard; spatial no-ground-truth
caveat present; methods/source-data section present; warnings surfaced (not
"No warnings" when issues exist).
</content>
