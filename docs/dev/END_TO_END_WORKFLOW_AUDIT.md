# End-to-End Workflow Audit

Audit date: 2026-06-02. Verified with small deterministic synthetic fixtures
only (no downloads): a 90-cell / 60-gene single-cell reference (3 cell types), a
4-sample bulk counts table, and a 36-spot / 60-gene Visium-like `.h5ad`. Each
workflow was executed through the real CLI.

## 1. Bulk workflow — `tissueresolve run --mode bulk`

**Runs:** ✅ completes; 3 cell types × 4 samples; resolution `auto → flat`
(no broad/fine labels in the toy reference).

**Outputs actually produced:**

```
analysis_plan.json          ✅ metadata JSON
run_metadata.json           ✅ metadata JSON
deconv/proportions.tsv      ✅ prediction table (rows sum to 1)
deconv/coverage_r2.tsv      ✅ reconstruction QC
deconv/gene_panel.txt       ✅ selected genes
deconv/metadata.json
qc/metadata.json            ✅ QC + warnings (embedded)
qc/recommendations.txt      ✅ QC recommendations
```

**Outputs missing vs. documentation:** `tables/`, `figures/`, `report.html`,
`warnings.json`, `methods.txt` are **not** written by `run`. (`report.html`
comes from the separate `tissueresolve report` step; the full
`tables/`+`figures/`+`methods.txt`+`warnings.json` bundle comes from the
validation-harness script `07_generate_reports.py`.)

## 2. Spatial workflow — `tissueresolve run --mode spatial`

**Runs:** ✅ completes; 3 cell types × 36 spots; NB-CAR model.

**Outputs actually produced:**

```
analysis_plan.json              ✅ metadata JSON
run_metadata.json               ✅ metadata JSON (records lambda_spatial)
deconv/proportions.tsv          ✅ prediction table (spot_rna_composition)
deconv/marker_genes.txt
deconv/mismatch_factors.npy
deconv/convergence_trace.json   ✅ convergence status
deconv/metadata.json
qc/morans_i.tsv                 ✅ Moran's I
qc/spot_qc.tsv                  ✅ per-spot QC
qc/metadata.json                ✅ QC + warnings (embedded)
```

**Outputs missing vs. documentation:** same as bulk — `tables/`, `figures/`,
`report.html`, `warnings.json`, `methods.txt` are not written by `run`.

## 3. Report generation — `tissueresolve report --modality bulk|spatial`

**Runs:** ✅ writes `report.html`. Verified the command succeeds against a run
output directory and (via existing `tests/report/test_results_dir_report.py`)
against a report-shaped results directory.

**Caveat (layout mismatch):** the results-dir report reader expects top-level
files (`bulk_estimated_proportions.tsv`, `tables/`, `figures/`,
`*_spot_proportions.tsv`, `morans_i.tsv`, …). A raw `tissueresolve run` output
puts predictions in `deconv/proportions.tsv` instead, so running `report`
directly on a `run` output produces a structurally complete but
**prediction-empty** report. The harness script `07_generate_reports.py` is the
component that materialises the report-shaped layout from a run.

## Summary matrix

| Required output | Bulk run | Spatial run | Report step | Notes |
|---|---|---|---|---|
| prediction table | ✅ `deconv/proportions.tsv` | ✅ `deconv/proportions.tsv` | reads report-shaped layout | layout differs from report reader |
| QC / warning table | ✅ `qc/*` | ✅ `qc/*` | ✅ surfaces warnings | warnings embedded in `qc/metadata.json`, not a standalone `warnings.json` |
| ≥1 figure | ❌ | ❌ | ✅ when figures present | `run` does not render figures |
| `report.html` | ❌ | ❌ | ✅ | separate step |
| `warnings.json` | ❌ (in `qc/metadata.json`) | ❌ | depends on input dir | not a top-level file from `run` |
| `methods.txt` | ❌ | ❌ | depends on input dir | `methods_text` module exists but is not invoked by `run` |
| metadata JSON | ✅ `run_metadata.json`, `analysis_plan.json` | ✅ | reads `run_metadata.json` | ✅ |

## Findings

- **Which workflows run:** both bulk and spatial `run` complete; the report
  step runs. The core algorithms execute end-to-end on toy data.
- **Which outputs are produced:** prediction table, QC tables, and metadata JSON
  are produced by every `run`. Figures, `report.html`, `methods.txt`, and a
  top-level `warnings.json` are **not** produced by `run`.
- **Are warnings visible:** yes within `qc/metadata.json` and
  `qc/recommendations.txt`, and the generated `report.html` surfaces a Warnings
  section. They are not silently dropped.
- **Is methods.txt generated:** not by `run`. A `methods_text` module exists and
  is used by the report layer/harness.
- **Output folder structure consistency:** **inconsistent** between `run`
  output (`deconv/`, `qc/`) and what the report reader and the documentation
  describe (`tables/`, `figures/`, top-level prediction files). This is the main
  end-to-end gap.

## Decision (per CLAUDE.md / Part 8 — minimal fixes, prefer doc correction)

The output-structure gap is a **documentation** defect (the docs overstate what
`run` emits), not a bug in the core algorithms. Action taken:

1. Corrected README/tutorial/reporting docs to describe `run`'s real outputs and
   the separate report step (Part 5/8).
2. Added honest regression tests (`tests/test_end_to_end_workflows.py`) that
   assert the **real** `run` contract (predictions + QC + metadata JSON) and that
   the report step produces `report.html` from a report-shaped directory.

**Recommended (deferred, not done here to avoid scope creep):** wire report
generation (and `methods.txt` / `warnings.json` emission) into `tissueresolve
run` so a single command produces the full bundle the docs advertise. This is an
alpha-stabilization priority, tracked in `ALPHA_READINESS_REPORT.md`.

## Update — report bundle now wired into `run`

`tissueresolve run` now writes `methods.txt`, `warnings.json`, and `report.html`
in addition to `deconv/`, `qc/`, and the JSON plans. `report.html` is rendered
from the **in-memory result**, so predictions, QC, methods, and warnings are
populated (verified on the synthetic bulk and spatial fixtures). `warnings.json`
always records the estimate-type caveat plus QC recommendations, non-convergence
(spatial), and protocol risk. Report-bundle generation is wrapped so a failure
is reported but never aborts a successful deconvolution.

**Still not produced by `run`:** standalone figure files (`figures/*.html` +
`.data.tsv`) and a `tables/` directory — those remain the responsibility of the
report layer / validation harness (the chosen scope was "methods.txt +
warnings.json + report.html, no new figures"). Regression tests in
`tests/test_end_to_end_workflows.py` assert the new bundle outputs.
</content>
