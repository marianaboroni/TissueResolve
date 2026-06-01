# TissueResolve — Audit Executive Summary

## What it is
TissueResolve is a reference-based cell-type/cell-state **deconvolution** package for bulk RNA-seq and 10x Visium spatial transcriptomics, built from two legacy codebases (CHIMERA bulk, SpatCAR spatial). Its distinctive framing is **resolution-awareness**: rather than forcing fine subtypes, it decides per cell-family whether subtypes are separable and otherwise reports `unresolved_<family>` mass. It ships a unified reference layer, separability/spillover diagnostics, protocol-awareness, uncertainty estimation, an opt-in real-data validation harness, and HTML reports.

## Defensible core claim
**Honest, gated abstention at the right cell-type resolution, backed by separability/spillover diagnostics and a single reference shared across bulk and spatial.** This is sound and honestly implemented in code (explicit mRNA-proportion typing, recorded smoothing parameters, non-suppressed warnings). It is **not** defensible as an accuracy-competitive method: on the package's own benchmark the best result (`TissueResolve_auto`, Pearson 0.7682150016166125) **equals plain NNLS to 16 digits by construction** (the auto-selector chose the NNLS backbone), the flagship hierarchical mode trails badly on fine accuracy (0.062), and the "combines competitive solvers" language is an overclaim (auto = single backbone, no blending).

## Top 10 BLOCKERs
1. Unified report prints "No warnings recorded" while 30+ real QC/spillover warnings exist (CLAUDE.md rule-2 violation).
2. Composite-score modality leakage: spatial methods (CARD, cell2location) ranked in the bulk leaderboard with meaningless 0.295 scores.
3. License contradiction: on-disk LICENSE is **MIT**, pyproject declares **BSD-3-Clause**, README claims no LICENSE file exists.
4. `TissueResolve_auto` is literally NNLS on this dataset; the headline "best method" is a self-tie presented as a win.
5. Validation breadth: a single tiny dataset (12 bulk samples, 32/41 types, 1915 cells capped 60/type; synthetic-truth spatial), **no confidence intervals anywhere**.
6. The robustness payoff that justifies sacrificing clean-data accuracy is **never demonstrated** (no protocol-mismatch/noisy ground-truth experiment).
7. `tuning/` module fabricates metrics (hardcoded scores, `"Tuning ran (stub)"`) — a scientific-honesty hazard.
8. Figures are invisible in the unified report (22/22 figure cards have empty bodies).
9. External tools not truly benchmarked: only MuSiC + Bisque executed for bulk; spatial cell2location ran degenerate; CIBERSORTx exported-only; BayesPrism failed at timeout.
10. Half-finished report-system migration uncommitted on `main`; HEAD does not build/import the new modules and their only consumer is also uncommitted.

## Top 10 HIGH tasks
1. Route the real protocol-aware `BulkPipeline` into a benchmarked TissueResolve method (currently bypassed by `_auto`).
2. Fix preset→config wiring (`bootstrap:False` ignored → 200 bootstraps; `auto_tune`/`plots`/`spatial` silently dead).
3. Surface and quantify the silent 41→32 cell-type intersection that drops 9 true types' mass.
4. Add a per-row hierarchical mass-conservation assertion + regression test (consistent but unguarded).
5. Validate abstention precision/recall (currently degenerate vs an empty `high_risk_families` set).
6. Replace boilerplate report interpretation with the data-driven module the sub-reports already use.
7. Replace flattened global Pearson with per-cell-type metrics; relabel composite as an opinionated scorecard and show accuracy-only alongside.
8. Add offline end-to-end tests for example scripts 03/04/07/09.
9. Reconcile doc/CLI mismatches (`bulk run` stub, non-existent `--n-bootstrap`, broken `--include-imported`).
10. Execute ≥3 external bulk + ≥3 spatial tools at full settings with CIs and significance.

## Publication-readiness score: **30 / 100**
Justification: the architecture and engineering are unusually disciplined (deterministic seeds, scoped warnings, enforced estimate-type typing, every figure co-saves source data, honest executed/imported/exported/skipped taxonomy) — this earns real credit. But publishability is gated by science and integrity issues: the quantitative case is a tie with NNLS, the flagship mode trails on accuracy, validation is a single tiny dataset with no CIs, competitors are mostly not executed, and there are concrete integrity violations (hidden warnings, fabricated tuning metrics, contradictory license). The defensible contribution (diagnostics + abstention) is genuine but under-validated. Score reflects "strong engineering, unproven scientific claim, fixable integrity defects."

## Is it currently publishable, and where?
**Not yet** — even the most forgiving venue requires fixing the BLOCKERs first. The realistic ceiling on current evidence (after Stages 0–2 of the release plan) is **Bioinformatics Application Note, NAR Genomics & Bioinformatics, GigaScience/Scientific Data, or Genome Biology**, framed as a **resolution-aware deconvolution + honest-abstention diagnostics framework** — not as an accuracy-superior method. **Nature Methods / Nature Biotech are out of reach** without a quantified, multi-dataset, statistically-significant advantage that does not currently exist. Nature Communications is possible only after multi-dataset validation with executed competitors.

## Single recommended next action
**Execute Stage 0 stabilization first — and within it, fix the report "No warnings" BLOCKER and remove the `tuning` fabricated-metrics stub immediately**, because these are scientific-integrity defects that contradict the project's own non-negotiable rules and undermine every other claim. In the same pass, commit-or-revert the uncommitted report system so `main` builds, and reconcile the license. Nothing else (benchmarks, framing, paper) is worth doing until the repo is honest and buildable.