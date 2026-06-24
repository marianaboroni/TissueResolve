# Bulk prediction sparsity audit

This document describes the TCGA sparsity audit used to investigate why bulk
deconvolution predictions appear to include only ~4–5 detected populations
per sample in some runs.

How to run

1. Ensure predictions are available in:

   - `benchmarks/outputs/bulk/predictions/*.tsv`

2. Run the audit script:

   ```bash
   python -m benchmarks.bulk.tcga_sparsity_audit
   ```

Outputs

- `benchmarks/outputs/tcga_prediction_complexity_audit_by_sample.tsv`
- `benchmarks/outputs/tcga_prediction_complexity_audit_by_method.tsv`

What is computed

- total reference populations
- number of populations with raw coefficient > 0
- number with proportion > 0.0001, 0.001, 0.005, 0.01
- effective number of populations (exp Shannon entropy)
- Gini coefficient
- Shannon entropy
- dominant fraction
- cumulative fraction in top 1/3/5
- unresolved mass (columns starting with `unresolved_`)
- fraction merged into `Other`
- fraction removed by small-threshold (default 0.001)

Notes

- The audit script reads every TSV under `benchmarks/outputs/bulk/predictions`. If
  your TissueResolve run outputs are elsewhere, copy them into that directory (without
  overwriting existing files) or run the script manually pointed at your files.
# Bulk Prediction Sparsity Audit (PART 1)

**Status:** Diagnostic complete. No core algorithm was modified (per `CLAUDE.md`
rule 10 and the study brief: "Do not change TissueResolve algorithms before
completing the diagnostic audit").

**Question.** Bulk TCGA-TNBC deconvolution appears to detect only ~4–5 cell
populations per sample, far fewer than the 32 populations in the reference. Is
this (a) solver behaviour, (b) post-processing/thresholding, (c) hierarchy
gating, (d) report/plot display filtering, or (e) true biological concentration?

**Answer (one line).** It is **solver behaviour** — unregularised NNLS and
weighted-NNLS on a near-collinear 32-population reference produce *sparse vertex
solutions* of ~3–4 active populations **inside the solver, before any
normalisation, threshold, aggregation, or display step**. Ridge regularisation
(which `auto` selects on TCGA) removes the effect entirely (≈23 populations).
The hierarchical path produces a *separate* ~4-population result by a different
mechanism (within-family resolvability gating). Post-processing, "Other"
aggregation, top-N display, and zero-thresholding are **not** responsible.

---

## 1. Method

Reproduced with `benchmarks/audit/tcga_sparsity_audit.py` (deterministic, seed
0), run on the real cohort:

- **Data:** `examples/real_breast_cancer/data/bulk_tcga_tnbc/tcga_tnbc_counts.tsv`
  — 40 TCGA-BRCA TNBC bulk samples × 59 427 genes; **4 982 genes shared** with
  the reference.
- **Reference:** `examples/real_breast_cancer/outputs/reference/breast_cancer_reference`
  — **32 fine populations**, 8 broad families.
- **Public API only** (`tissueresolve.api.deconv_bulk`), per the codebase
  discovery rule. Flat runs use `resolution_mode="none"`, i.e. the *raw
  normalised solver output* with **every optional filter at its default**
  (no min-abundance threshold, no top-N, no "Other" merge).

Per (method, sample) we compute the PART-17 panel: non-zero counts at
thresholds {0, 1e-4, 1e-3, 5e-3, 1e-2}; Shannon entropy; effective number of
populations `exp(H)`; Gini; dominant fraction; top-1/3/5 cumulative mass;
unresolved mass. Full tables:

- `benchmarks/outputs/tcga_prediction_complexity_audit.tsv` (320 rows = 8 methods × 40 samples)
- `benchmarks/outputs/tcga_method_complexity_summary.tsv`
- raw per-method proportions preserved under `benchmarks/outputs/raw/`.

## 2. Result — the source is the solver

Identical input, identical pipeline, only the solver changes:

| Solver (`resolution_mode="none"`) | mean non-zero | median non-zero | effective-N | dominant frac | top-5 cum |
|---|---|---|---|---|---|
| `nnls`            | 2.8  | **3**  | 2.1  | 0.64 | 1.00 |
| `weighted_nnls`   | 3.9  | **4**  | 2.9  | 0.56 | 1.00 |
| `marker_nnls`     | 11.4 | 12     | 7.0  | 0.32 | 0.83 |
| `ridge_nnls`      | 22.9 | 23     | 17.5 | 0.13 | 0.45 |
| `auto` (→ ridge)  | 22.9 | 23     | 17.5 | 0.13 | 0.45 |

The "4–5 populations" observation is reproduced **exactly** by `nnls` (median 3)
and `weighted_nnls` (median 4). The same TCGA data, the same 4 982 genes, and
the same post-processing pipeline yield **23 populations** under ridge. The
result is therefore determined by the **solver's regularisation**, spanning
2.1 → 17.5 effective populations on identical input.

**Mechanism.** With the non-negativity constraint, NNLS on a design matrix whose
columns are highly correlated (the 32 populations include many near-collinear
immune/endothelial/epithelial subtypes; reference condition number ≈ 486)
returns a solution on a low-dimensional *face of the feasible cone* — only a
handful of components are active and the rest are exactly zero. This is the same
implicit sparsity that makes constrained/`L1` problems sparse; it is a property
of the geometry, **not** a thresholding step. Ridge adds an `L2` penalty that
spreads mass and yields a dense solution.

## 3. What is NOT responsible (ruled out with evidence)

Every candidate from the brief was checked against the code and the data:

| Candidate cause | Verdict | Evidence |
|---|---|---|
| Negative clipping | Not a filter | `nnls` enforces θ≥0 *as the model* (`bulk/solver.py:364`); same constraint under ridge, which is dense. |
| Row normalisation | Not the cause | Only operation in flat mode (`bulk/solver.py:365`); preserves zeros; ridge normalised identically stays dense. |
| Min-abundance / zero threshold | Not applied | No data-level abundance threshold in the flat path; `n_gt_0 == n_gt_1e-4` within rounding. |
| "Other" aggregation | Display only, default off | `plotting/bulk_plots.py:284` `_collapse_to_top_n`, default `top_n=None`. Does not mutate stored proportions. |
| Top-N report display filter | Not applied | `report/templates.py` caps only warnings (5) and table rows; no population cap. |
| Post-processing zeroing | Ruled out | `flat_auto` proportions are **identical per-sample** to `flat_ridge_nnls` (23 non-zero) through the *same* pipeline — nothing downstream is removing populations. |
| Hierarchy gating | Separate mechanism, see §4 | Flat modes never touch the hierarchy. |
| Insufficient gene overlap | Contributing factor (conditioning), not a filter | 4 982 shared genes is adequate in count but the *signatures are collinear*; this is what makes NNLS sparse and ridge necessary. |

## 4. The hierarchical path is sparse for a different reason

Hierarchical mode (`resolution_mode="hierarchical"`) yields ~4 populations by an
independent mechanism — **within-family resolvability gating**
(`reference/hierarchy.py:675–897`):

- Broad families: 8 total, ~3 non-zero per sample (family-level NNLS sparsity).
- Fine subtypes: **2.8 non-zero of 32**; effective-N 1.8.
- **All 6 multi-member families are marked *unresolved*** and their subtype mass
  is zeroed into `unresolved_<family>` (`add_unresolved_family_mass`,
  `hierarchy.py:613–661`); mean **unresolved mass = 0.27**.

The gate fails predominantly on **separability** and **within-family spillover**,
not on gene count:

| Family | n_subtypes | mean separability (1−BC) | min disc. genes | mean spillover (Pearson) | resolvable? |
|---|---|---|---|---|---|
| T/NK | 11 | 0.048 | 48 | **0.86** | False |
| Myeloid | 7 | 0.081 | 171 | **0.81** | False |
| Endothelial | 5 | 0.121 | 127 | **0.67** | False |
| Epithelial | 3 | 0.069 | 809 | **0.76** | False |
| B/plasma | 2 | 0.300 | 1503 | **0.80** | False |
| Mural | 2 | 0.029 | 362 | **0.92** | False |

Thresholds: separability ≥ 0.10, disc. genes ≥ 10, spillover ≤ 0.30. Families
pass disc-genes easily (48–1503 ≫ 10) but fail on spillover (0.67–0.92 ≫ 0.30)
and usually separability — i.e. the within-family CPM signatures are
**near-collinear**, the same conditioning problem that drives flat-NNLS sparsity.
Here the behaviour is *intentional and conservative* (TissueResolve abstains
rather than fabricate subtype splits it cannot support), but it is the source of
the fine-level collapse and should be reported as such.

## 5. Important discrepancy to flag

The **currently shipped TCGA output**
(`examples/real_breast_cancer/outputs/bulk_tcga/tcga_bulk_estimated_proportions.tsv`)
was produced by `auto`, which **on TCGA selects ridge** → it actually contains
~23 populations (diffuse), *not* 4–5. The low cross-solver concordance already
recorded there (`nnls` vs `auto` Pearson = 0.44 in
`tcga_bulk_method_benchmark.tsv`) is now fully explained: `nnls` is sparse (≈3),
`auto`/ridge is dense (≈23). A user seeing "4–5 populations" is therefore
looking at one of: a non-auto solver (`nnls`/`weighted_nnls`), the hierarchical
report (broad families / fine collapse), or the dominant-type table — **not** the
default `auto` flat output.

> Note: `auto`'s selection is data-dependent. On the *pseudobulk* benchmark it
> selected `nnls` (`benchmarks/outputs/bulk/selected_solver.json`); on *TCGA* it
> selected ridge. The auto-selector's behaviour across scenarios is itself worth
> benchmarking (Phase 2).

## 6. Is 4–5 populations biologically correct?

**Undetermined — and not assertable without ground truth.** The audit proves the
result is overwhelmingly a function of *solver regularisation and signature
conditioning* (2.1 vs 17.5 effective populations on identical input), not a
display or threshold artefact. Whether the true tumour composition is
concentrated (≈ NNLS) or diffuse (≈ ridge) can only be decided by the
ground-truth pseudobulk benchmark (Study Part 2 / Phase 2). Neither extreme is
self-evidently right: unregularised NNLS is known to *under*-detect on collinear
references; ridge can *over*-disperse small mass across all populations.

## 7. Recommendation (do not act yet)

Per the brief and `CLAUDE.md` rule 10, **no defaults are changed in this stage.**
Candidate actions, to be validated against ground truth in Phase 2 before any
change:

1. Decide the default flat solver on *measured* pseudobulk accuracy +
   composition-complexity preservation (Parts 5/10), not on conditioning alone.
2. Surface the effective-population / dominant-fraction / unresolved-mass
   diagnostics in the bulk report so sparsity is visible to users.
3. Make `auto`'s solver choice and its sparsity consequence explicit in the report.
4. Re-examine the within-family spillover threshold (0.30) against ground-truth
   resolvability once Phase 2 fine-level metrics exist (Part 7).

## 8. Figures (to generate in Phase 1 plotting step)

The per-sample table supports all PART-1 figures:
detected-populations-vs-threshold; raw-vs-post-processed (here: flat-vs-flat,
confirming no removal); entropy-by-sample; dominant & top-5; fraction
removed/aggregated (≈0 in flat, 0.27 unresolved in hierarchical); broad-vs-fine
counts; per-population prevalence; rare-population detection frequency. These are
deferred to the Phase-1 plotting deliverable (every figure exports
PNG/SVG + `.data.tsv` + caption per `CLAUDE.md` plotting rules).

---

*Generated by `benchmarks/audit/tcga_sparsity_audit.py`. Raw outputs preserved
under `benchmarks/outputs/raw/`. No core algorithm modified.*
