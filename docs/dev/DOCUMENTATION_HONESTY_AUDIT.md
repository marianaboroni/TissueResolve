# Documentation Honesty Audit

Audit date: 2026-06-02. Reviews documentation for overclaims against the
required wordings, and records what was already correct vs. fixed in this pass.

## Required-wording checklist

| # | Claim | Required wording present | Where | Action |
|---|---|---|---|---|
| 1 | **State-aware** experimental | ✅ | README Status; CLI `--state-aware` help; `FEATURE_STATUS.md`; `V0_1_SCOPE.md` | already correct (guarded by `test_v0_1_stabilization.py`) |
| 2 | **Expression reconstruction** not implemented | ✅ | README Status; `FEATURE_STATUS.md` ("not implemented"); `V0_1_SCOPE.md`; `output_interpretation.md` uses "cell-type-level deconvolution" | already correct |
| 3 | **External benchmarks** — only executed/imported ranked | ✅ | `FEATURE_STATUS.md` (verbatim); README Benchmarking; `benchmarking.md` ("exported or skipped … never counted as benchmarked or scored") | already correct |
| 4 | **Spatial benchmark** — no accuracy without ground truth | ✅ | `benchmarking.md`; `output_interpretation.md` §"Spatial benchmark is not accuracy"; `report/glossary.py`; added verbatim to `FEATURE_STATUS.md` | strengthened |
| 5 | **Estimate type** — bulk RNA-derived proportions / spatial spot-level composition | ✅ | README Interpretation; `output_interpretation.md`; added verbatim to `FEATURE_STATUS.md` | strengthened |
| 6 | **Publication status** — alpha / early-access; not publication-ready until benchmarks + multi-dataset validation + API stabilization | ✅ (added) | added to README Status, `FEATURE_STATUS.md`, `V0_1_SCOPE.md` | **added this pass** |

## Overclaims found and fixed

1. **Run output overclaim** (README "Outputs", `tutorial.md` §6/§9): listed
   `tables/`, `figures/`, `report.html`, `warnings.json`, `methods.txt` as
   products of `tissueresolve run`. A run writes only `analysis_plan.json`,
   `run_metadata.json`, `deconv/`, `qc/`. **Fixed:** corrected to the real
   outputs and clarified the separate report step / harness produces the full
   bundle.

2. **Broken bulk-flat example** (`tutorial.md`): used the unimplemented
   `tissueresolve bulk run` stub with non-existent flags. **Fixed:** rewritten to
   `tissueresolve run --mode bulk … --resolution-mode flat`.

3. **Broken benchmark command** (README, `import_external_results.py` help):
   `run_all.py --include-imported` errored. **Fixed:** added/forwarded the flag
   so the documented command works.

4. **Missing publication-status statement** (README): only said "research
   software / pre-release". **Fixed:** added the explicit alpha / early-access /
   not-yet-publication-ready statement.

## Residual (low severity, noted not fixed)

- `output_interpretation.md` references a report message "run with
  `--n-bootstrap > 0`", but `run` has no `--n-bootstrap` flag — bootstrap is
  controlled by `--preset publication|diagnostic`. The string is emitted by the
  report builder; left unchanged to avoid touching report internals. Documented
  in `CLI_DOCUMENTATION_AUDIT.md`.

## Claims that are honest and well-supported

- Estimate-type typing is enforced in code (`ESTIMATE_TYPE`,
  `# estimate_type:` headers) — not just documentation.
- Smoothing parameters (`lambda_spatial`, `alpha`, `n_smooth`) are recorded, per
  `output_interpretation.md` and the spatial run metadata.
- Hierarchical "unresolved family mass" is mass-preserving and explicitly
  reported; never a silent merge.
- The README comparison table is explicitly framed as "a high-level orientation,
  not a claim that TissueResolve outperforms these tools."
</content>
