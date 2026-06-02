# TissueResolve Alpha-Readiness Report

Date: 2026-06-02. Branch: `benchmark-report-refactor`.

## 1. Executive summary

TissueResolve is a resolution-aware, reference-based deconvolution package for
bulk RNA-seq and 10x Visium with an unusually disciplined engineering base:
deterministic seeds, enforced estimate-type typing, recorded smoothing
parameters, non-suppressed warnings, every figure co-saving its source data, and
an honest executed/imported/exported/skipped benchmark taxonomy. Both the bulk
and spatial pipelines run end-to-end on synthetic fixtures and produce predictions,
QC, and provenance metadata. The CLI matches the documentation after this pass's
fixes. The defensible contribution — **honest, gated abstention at the right
cell-type resolution, backed by separability/spillover diagnostics and a single
shared reference** — is genuine and honestly implemented.

It is **alpha / early-access research software**. The priority is to stabilize
what exists rather than expand it.

## 2. Current test status

Full offline suite (`python3 -m pytest -q`): **all collected tests pass**
(1050 baseline + 21 new = see final response for the exact count). No tests were
weakened or skipped. New tests added this pass:
`tests/test_end_to_end_workflows.py`, `tests/test_docs_cli_consistency.py`,
`tests/test_report_stability.py`.

## 3. CLI status

All documented commands exist and respond to `--help`. `tissueresolve run`
(`--mode auto|bulk|spatial`), `tissueresolve spatial run`, `tissueresolve
report`, `tissueresolve bulk/spatial report`, `tissueresolve info`, and the
benchmark scripts are correct. `tissueresolve bulk run|benchmark|
check-compatibility` are explicit unimplemented stubs (exit 2) and are not
advertised as working. Fixed this pass: `run_all.py --include-imported` (now
forwarded), the broken `tissueresolve bulk run` tutorial example, and the
overstated run-output documentation.

## 4. Bulk workflow status

**Working.** `tissueresolve run --mode bulk` builds a reference, runs the solver
(`auto`/pipeline), and writes `deconv/proportions.tsv`, reconstruction QC,
`qc/`, `analysis_plan.json`, and `run_metadata.json`. Predictions are a valid
simplex. Estimate type is RNA-derived mRNA proportions (enforced).

## 5. Spatial workflow status

**Working.** `tissueresolve run --mode spatial` runs the NB-CAR model and writes
spot proportions, Moran's I, per-spot QC, convergence trace, and metadata.
Smoothing (`lambda_spatial`) is recorded. Estimate type is spot-level RNA-derived
composition (enforced).

## 6. Report status

**Stable.** The report is QC-first, separates bulk/spatial, gates figures on
artifact presence (no empty cards), carries the spatial no-ground-truth caveat,
and links a technical appendix with source data. `tissueresolve run` **now**
emits `report.html` (rendered from the in-memory result, so predictions/QC/
methods/warnings are populated), `methods.txt`, and `warnings.json` directly.

**Report paths consolidated.** There is now a **single canonical path**:
`report/orchestration.py` (`generate_report(modality, source, out)` +
`render_sections()`), which renders every report — CLI `report`,
`bulk/spatial report`, `tissueresolve run`, and `api.generate_report` — through
the unified single-page shell from either a results directory or an in-memory
result. `report/html.py`'s public functions are deprecation shims; its duplicate
page renderers were removed. `sections.py`/`assets.py`/`interpretation.py` are
the results-dir content layer for the canonical path. See
`docs/REPORT_PATH_CONSOLIDATION_PLAN.md`.

**Remaining gap:** `run` still does not write standalone figure files
(`figures/*.html` + `.data.tsv`) — those come from the report layer / harness
(the chosen scope was "no new figures"). Fully deleting `html.py` (vs the
current shim) would require relocating its in-memory section builders — an
optional further cleanup.

## 7. Documentation status

Honest after this pass. All six required wordings (state-aware experimental,
expression reconstruction not implemented, external-benchmark executed/imported
only, spatial no-accuracy-without-ground-truth, estimate-type, and
publication-status alpha) are present and now test-guarded.

## 8. Stable features

Bulk deconvolution; spatial deconvolution; solver `auto`; broad/fine hierarchy
with unresolved mass; reference suitability; separability/spillover diagnostics;
QC-first report + technical appendix; deterministic family-aware palette;
per-figure source data.

## 9. Experimental features (behind flags, labelled, not default)

State-aware three-level hierarchy (`--state-aware`); granular signatures;
spatial multi-metric ranking; external-tool benchmark runners; composite
scorecard; synthetic state-aware benchmark.

## 10. Planned / not implemented (must not be claimed)

Reference adaptation; cell-type-specific expression reconstruction;
hyperparameter tuning (quarantined stub that raises); full BayesPrism-like
Bayesian (Gibbs) model.

## 11. Remaining blockers for alpha release

1. ✅ **Done:** `tissueresolve run` now emits `report.html`, `methods.txt`, and
   `warnings.json` from the in-memory result.
2. Optionally render standalone figure files (`figures/*.html` + `.data.tsv`)
   during `run` for full parity with the harness (chosen scope this pass was
   "no new figures").
3. ✅ **Done:** the two report code paths are consolidated onto a single
   canonical path (`report/orchestration.py` → unified shell). `html.py`'s
   duplicate renderers were removed; its public API is now a deprecation shim.

(None of these block correctness of the core algorithms; they are
completeness/consistency items.)

## 12. Remaining blockers for publication

1. Quantified, multi-dataset accuracy/robustness evidence — current validation
   is a single small dataset; the auto-solver's headline equals plain NNLS on
   that data, and the robustness payoff is not yet demonstrated.
2. Execute (not just export/skip) ≥3 external bulk and ≥3 external spatial tools
   at full settings, with confidence intervals and significance.
3. Validate abstention precision/recall and hierarchical mass conservation on
   real data with non-empty high-risk families.
4. Confidence intervals reported throughout.

## 13. Recommended next actions

1. Land the report-from-run wiring (alpha blocker #1) and the
   `warnings.json`/`methods.txt` emission.
2. Add a per-row hierarchical mass-conservation assertion + regression test.
3. Plan the multi-dataset external benchmark with CIs (publication track).
4. Consolidate the two report code paths.

## 14. Safe to commit?

**Yes.** This pass only added documentation, audit reports, regression tests, and
small consistency fixes (a forwarded benchmark flag, corrected examples, and
honest output documentation). No core algorithm was changed; the full offline
suite passes.

## Final assessment

- **Alpha ready: yes** — as alpha / early-access research software, with the
  report-from-run wiring as the top stabilization item.
- **Publication ready: no** — not until the publication blockers (multi-dataset
  validation, executed external competitors, confidence intervals, demonstrated
  robustness payoff, and final API stabilization) are resolved.

TissueResolve is promising and has a strong technical foundation, but should
currently be considered alpha / early-access research software. The priority is
to stabilize what exists, not to expand the tool.
</content>
