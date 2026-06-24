# Pipeline Reorganization Plan — Evidence-Based Decision Order

Reorganize TissueResolve so it **decides what resolution is statistically
supportable before interpreting fine-subpopulation predictions**, instead of
"deconvolve everything first, judge later". This is an organizational/decision-layer
change, **not** a new algorithm and **not** a rewrite. Soft gating remains the
validated default and the final uncertainty layer.

## Task 1 — Audit of the current flow

### Current execution order (bulk + spatial hierarchical)

```
tissueresolve run  →  _run_top_level (cli.py:152)
  → _execute_bulk / _execute_spatial (cli.py:491 / :563)
     → deconv_bulk / deconv_spatial (api.py:73 / :249)
        → run_hierarchical_bulk (bulk/hierarchical.py:77)
           1. aggregate_reference_by_family            (reference/hierarchy.py:525)
           2. BulkPipeline.run(bulk, family_ref)       BROAD deconvolution   (hierarchical.py:126)
              ├─ check_gene_overlap (query compat)      (pipeline.py)
              └─ GeneSelector.select + compute_gene_weights  GENE WEIGHTING  (pipeline.py:194/209)
           3. BulkPipeline.run(bulk, fine_ref)         FINE deconvolution    (hierarchical.py:128)
           4. assemble_hierarchical_estimates          (reference/hierarchy.py:950)
              ├─ evaluate_within_family_resolvability  RESOLVABILITY  (hierarchy.py:675)
              ├─ compute_within_family_subtype_confidence            (hierarchy.py:795)
              ├─ gating dispatch (soft default)         SOFT GATING   (hierarchy.py:1000)
              └─ unresolved mass                        UNRESOLVED    (hierarchy.py:1034)
  → _write_run_report_bundle (cli.py:646) → report (orchestration.py:52, result_sections.py:232)
```

### Problems identified

1. **"Deconvolve first, judge later."** Fine deconvolution (`hierarchical.py:128`) runs
   **before** within-family resolvability is evaluated (`assemble_hierarchical_estimates`,
   step 4). The pipeline never decides *whether fine resolution is supportable* before
   attempting it. Fine estimates always exist; gating only reallocates their mass.
2. **No explicit trusted-resolution status.** There is no first-class per-family
   classification (`broad_only` / `selected_fine` / `full_fine`). The `resolvable`
   boolean in the resolvability table is binary and internal; the report reads it for
   display but it is not a formal decision object recorded in metadata.
3. **Reference QC / suitability is not on the critical path.** `compute_reference_suitability_score`
   (suitability.py:222) exists but is **not** automatically computed before deconvolution;
   it does not feed the resolution decision.
4. **Identifiability levels are conflated.** Cell-level classifiability, pairwise mixture
   recovery, and full-panel deconvolution reliability are not separated in the decision;
   the report shows AUROC-style evidence without making clear only full-panel reliability
   should drive trusted fine status.
5. **Decision is implicit/visual, not a recorded algorithmic object.** The numeric output
   (`combined_fine`) already reflects soft gating (good — gating IS algorithmic), but the
   *trusted-resolution status per family* is not computed, recorded in metadata, or placed
   before fine predictions in the report.
6. **Modality differences are implicit.** Bulk uses protocol-aware weighted NNLS; spatial
   uses unweighted markers + NB-CAR. This is correct but not recorded as an explicit
   "modality-aware gene weighting" mode in metadata, and the report does not separate
   shared vs bulk-specific vs spatial-specific compatibility.

### What is already correct (do not change)

* Soft gating is the validated default and is applied exactly once, last (after fine).
* Soft gating IS algorithmic: `combined_fine` already routes low-confidence/unresolvable
  mass to `unresolved_<family>`; the report displays, it does not re-gate.
* Broad-family mass is conserved; total mass is conserved.
* `reference/resolution.py` provides reference-level merge diagnostics (separate concern).

## Proposed flow (after)

```
1. Build/load reference
2. Reference QC (suitability)                         ── feeds decision
3. Query QC + reference-query compatibility           ── feeds decision
4. RESOLUTION DECISION LAYER  →  per-family {broad_only|selected_fine|full_fine}
5. Modality-aware gene weighting (recorded mode)
6. Broad deconvolution
7. Fine deconvolution only where supported (broad_only → diagnostic only)
8. Soft gating + unresolved mass (final layer, once)
9. QC-first report: quality + trusted-resolution table BEFORE fine predictions
```

The decision is **recorded algorithmically** in `HierarchicalEstimates.metadata`
(`trusted_resolution`) and the per-family `qc` table, and consumed by the report.

## Implementation (minimal footprint)

* **New file `src/tissueresolve/resolution.py`** — the decision layer:
  `ResolutionDecision` dataclass, `ResolutionDecisionConfig`,
  `decide_trusted_resolution(...)`, and `decisions_from_estimates(...)` (derives
  per-family decisions from the already-computed resolvability + confidence, so no
  extra deconvolution is needed), plus `compute_modality_aware_gene_weights(...)`
  (thin modality-aware wrapper over the existing gene selection + a recorded mode).
* **`reference/hierarchy.py`** — record the per-family trusted-resolution status in
  `assemble_hierarchical_estimates` metadata + qc (additive; numeric outputs unchanged).
* **`bulk/hierarchical.py`, `spatial/hierarchical.py`** — record modality metadata
  (modality, prediction unit, gene-weighting mode, spatial smoothing used, H&E
  available, whether fine predictions are trusted vs diagnostic). Spatial smoothing
  default unchanged.
* **`report/result_sections.py`** — add a "Trusted resolution by family" table BEFORE
  the fine composition; mark `broad_only` families' fine as **diagnostic, not trusted**;
  keep the existing soft-gating / spillover / false-positive cautions.
* **Docs** — README/method overview, `output_interpretation.md`, `reporting.md`,
  `FEATURE_STATUS.md`, this plan.
* **Tests** — `tests/test_resolution_decision_layer.py` (one compact file).

### Expected behaviour changes

* **No change to numeric estimates by default** — soft gating already routes
  unsupported mass to unresolved; backward-compatible. The decision layer adds an
  explicit *status* and reorders/relabels the report so quality + trusted resolution
  appear before fine predictions and `broad_only` fine is labelled diagnostic.
* New metadata fields; new report table/cards.

### Risks

* Touching `assemble_hierarchical_estimates` metadata could affect tests that inspect
  metadata — mitigated by additive-only changes.
* Report-text tests may assert specific strings — update them alongside.
* Must not change spatial λ default, broad prediction, or soft-gating default.

### Tests to add

Pipeline order, resolution behaviour (broad_only → no trusted fine; selected_fine →
supported subtypes only; full_fine gated; missing evidence → conservative; AUROC alone
cannot trigger full_fine), mass conservation, modality metadata (bulk no smoothing;
spatial records coordinates/smoothing), and report behaviour (QC before predictions;
trusted-resolution table before fine; RNA-derived not cell fractions; spatial no
accuracy claim without truth).
