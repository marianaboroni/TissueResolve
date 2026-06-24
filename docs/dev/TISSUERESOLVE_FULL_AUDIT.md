# TissueResolve — Full Audit

_Read-only audit (9-agent review across architecture, scientific utility, algorithms, code quality, tests, benchmark, reports/UX, documentation, and synthesis). Severity: BLOCKER / HIGH / MEDIUM / LOW / QUESTION. Grounded in inspected source._

## PART 1 — Repository & Architecture Audit

Scope inspected: `src/tissueresolve/` (all subpackages), `benchmarks/`, `examples/real_breast_cancer/`, `docs/`, `tests/`, CLI, `scripts/`, legacy folders (`chimera_v1/`, `spatcar/`), scratch (`demo_demo/`, `demo_wiz/`, `data/`), `.gitignore`, `pyproject.toml`. Everything below is grounded in files I read.

### Method / what I actually ran
- `git ls-files` directory histogram + full file listing.
- `git status --porcelain --ignored` and `git check-ignore` to separate tracked / untracked / ignored.
- Read: `src/tissueresolve/__init__.py`, `api.py`, `cli.py` (headers + command map), `report/__init__.py`, `report/html.py` (tail), `solver/__init__.py`, `bulk/solver.py` (head), `tuning/__init__.py`, `validation/__init__.py`, `pyproject.toml`, `.gitignore`, `docs/project_structure_audit.md`, plus module/script headers across all subpackages.

### Repository hygiene findings (git)
- **Legacy folders are correctly excluded.** `chimera_v1/` and `spatcar/` exist on disk but are **not tracked** and are git-ignored (`.gitignore:26-27`). Good — matches CLAUDE.md "do not modify legacy folders". *(QUESTION/LOW: they remain on disk at 100KB+/24KB; intended for reference only.)*
- **Scratch dirs ignored:** `demo_demo/`, `demo_wiz/`, `data/` (46 MB), `.DS_Store`, `.venv/`, `.pytest_cache/`, `src/tissueresolve.egg-info/`, `benchmarks/envs/{Rlib,c2l_venv}/`, `benchmarks/{logs,outputs}/` — all ignored and untracked (verified). `demo_demo`/`demo_wiz` are referenced by **no** tracked source/test/doc. **LOW**.
- **Uncommitted work on `main` (HIGH):** `git status` shows 12 modified tracked files and 5 new untracked files that are clearly part of the codebase, not scratch:
  - Untracked source: `src/tissueresolve/report/{components.py,figures.py,glossary.py,style.py}`, `tests/report/test_report_components.py`.
  - Modified: `src/tissueresolve/report/unified.py`, `examples/real_breast_cancer/scripts/07_generate_reports.py`, `src/tissueresolve/api.py`, `src/tissueresolve/cli.py`, several `benchmarks/*`, `tests/cli/test_cli.py`, `tests/test_api.py`.
  - The committed HEAD therefore does **not** build/import the new report design-system; its only consumer (`07_generate_reports.py`) is also uncommitted. CI as committed does not exercise the new modules.
- **Two genuinely-untracked, NON-ignored env files (LOW):** `benchmarks/envs/Makevars`, `benchmarks/envs/_install_music_patched.R` (`git check-ignore` returns nonzero). Decide intentionally whether to track or ignore.
- **Empty `scripts/` dir (LOW):** top-level `scripts/` is empty and unused (git can't track empty dirs); likely leftover. Distinct from `examples/real_breast_cancer/scripts/`.

---

### Table 1 — Core modules (`src/tissueresolve/`)

| Module | Purpose | Status | Redundancy risk | Recommendation |
|---|---|---|---|---|
| `api.py` | Top-level `deconv_bulk/spatial`, `build_reference`, `generate_report`, `plot_results` | Active, public | Medium — overlaps CLI routing; `_solver_bulk_result` reimplements reconstruction R² | Stop calling them "thin wrappers" (docstring L4-5) — they route resolution modes; centralize resolution logic shared with CLI |
| `cli.py` (921 LOC) | Click CLI: `run`, `bulk {run,check-compatibility,benchmark,report}`, `spatial {run,report,benchmark,info}`, top `report`, `info` | Active | Medium — duplicates resolution/solver gating done in `api.py` | Factor shared routing out of cli+api |
| `config.py`, `presets.py`, `results.py`, `utils.py` | Config dataclasses, presets, result containers (`ReferenceSignature`, `BulkDeconvResult`, `SpatialDeconvResult`, `QCReport`, `SeparabilityReport`, `BenchmarkResult`), helpers | Active, central | Low | Keep |
| `io/` (bulk, spatial, reference, validation, autodetect) | Input loading + autodetection | Active | Low | Keep |
| `reference/` (build, markers, pairwise_markers, gene_filters, separability, resolution, suitability, hierarchy, hierarchical_build) | Reference construction, markers, separability, hierarchy | Active | Medium — `markers.py`+`pairwise_markers.py`, `build.py`+`hierarchical_build.py` are adjacent pairs | Verify pairwise vs single marker paths don't duplicate logic |
| `protocol/` (detect, metadata, mismatch, risk) | Protocol-aware QC (CHIMERA) | Active | Low | Keep |
| `bulk/` (pipeline, solver, qc, hierarchical) | Bulk pipeline + protocol-aware weighted NNLS solver | Active | **See solver note** | Keep; document relationship to `solver/` |
| `solver/` (base, nnls, weighted_nnls, marker_nnls, ridge_nnls, ensemble, auto) | Pluggable solver backbones + auto-selection via gene-masking CV | Active | **High naming overlap** with `bulk/solver.py` (the CHIMERA pipeline solver). Two different "solvers." | Rename or document: `bulk/solver.py` = pipeline default; `solver/` = `--solver` backbones. Confusing as-is |
| `spatial/` (model, graph, pipeline, qc, neighbourhood, auto_params, hierarchical, **benchmark.py** 704 LOC) | SpatCAR NB-CAR model + spatial pipeline + in-package synthetic benchmark | Active | Medium — `spatial/benchmark.py` overlaps the top-level `benchmarks/spatial/` harness conceptually | Clarify scope; consider moving synthetic benchmark out of the importable package |
| `uncertainty/bootstrap.py` | Bootstrap CIs (CHIMERA) | Active | Low | Keep |
| `validation/gene_masking.py` | Gene-masking CV (no-ground-truth model selection) | Active | Low | Keep — well documented |
| `benchmark/` (spillover.py) | Spillover simulation/estimation utils | Active | **Name collision** with `spatial/benchmark.py` and top-level `benchmarks/` | Rename to `spillover/` or fold into `reference/`; "benchmark" is overloaded 3 ways |
| `plotting/` (12 modules) | Bulk/spatial/QC/separability/spillover/benchmark/histology/summary plots + palette/style/captions/export | Active | Low–Medium — `summary_figures.py` vs per-type plot modules may overlap | Keep; audit for duplicate figure builders |
| `report/` (html, templates, sections, assets, interpretation, methods_text, unified, + untracked components/figures/glossary/style) | HTML report generation | **Mid-migration, fragmented** | **HIGH** | See Table 4; consolidate |
| **`tuning/__init__.py`** | "Parameter tuning" | **STUB — fabricates metrics** | Orphaned (imported nowhere in `src`) | **HIGH: delete or implement.** Writes hardcoded `score 0.8`, `lambda=0.05 → 0.9`, `"Tuning ran (stub)"`. If ever surfaced it reports fake numbers — violates CLAUDE.md "do not assume a function is scientifically valid just because it runs" and "transparent QC" |

---

### Table 2 — Scripts (harness + CLI)

| Script | Purpose | Still needed? | Can merge? | User-facing? | Recommendation |
|---|---|---|---|---|---|
| `examples/.../00_download_data.py` | Opt-in download of reference + Visium | Yes | No | Yes (opt-in) | Keep |
| `01_prepare_reference.py` | Build TissueResolve reference | Yes | No | Yes | Keep |
| `02_make_pseudobulk.py` | Pseudobulk mixtures + ground truth | Yes | Could fold into 03 | Yes | Keep (truth generator) |
| `03_run_bulk_validation.py` | Flat bulk deconv vs truth | Yes | Pair w/ 02 | Yes | Keep |
| `04_run_spatial_validation.py` | Spatial deconv + QC | Yes | No | Yes | Keep |
| `05_summarize_results.py` | Markdown/JSON summary | Optional | Superseded by 07 unified report | Yes | Keep as lightweight summary or fold into 07 |
| `06_resolution_spillover_analysis.py` | Separability + spillover matrix + merges | Yes | **Overlaps 08** (both → `outputs/resolution/`) | Yes | **Merge with 08** |
| `07_generate_reports.py` (modified, uncommitted) | HTML reports + unified `report.html` | Yes | — | Yes | Commit; it is the only consumer of the new report design-system |
| `08_resolution_analysis.py` | Merge recommender + family aggregation | Yes | **Overlaps 06** | Yes | **Merge with 06** into one `resolution_analysis` step |
| `09_run_hierarchical_deconvolution.py` | Broad→fine bulk+spatial | Yes | No | Yes | Keep |
| `_download_utils.py`, `_harness.py` | Support (paths/IO/metrics/download) | Yes | No | No (internal) | Keep |
| `benchmarks/run_all.py` | Run both benchmarks + summary report | Yes | No | Yes | Keep |
| `benchmarks/run_real_external_benchmark.py` (mod) | External-tool orchestration, honest skip | Yes | Overlaps `import_external_results.py` partly | Yes | Keep; commit changes |
| `benchmarks/import_external_results.py` | Import externally-run predictions | Yes | No | Yes | Keep |
| `benchmarks/bulk/run_bulk_benchmark.py` (mod) / `spatial/run_spatial_benchmark.py` (mod) | Per-modality benchmark runners | Yes | No | Yes | Keep; commit changes |
| CLI `tissueresolve {run,bulk,spatial,report,info}` | Public entry point | Yes | `bulk benchmark`/`spatial benchmark` subcommands partly duplicate `benchmarks/` harness | Yes | Clarify CLI-benchmark vs `benchmarks/` harness boundary |

---

### Table 3 — Outputs / data tracking

| Path | Should be tracked? | Current status | Recommendation |
|---|---|---|---|
| `examples/real_breast_cancer/data/` | No | Ignored (`.gitignore:30`) | Correct |
| `examples/real_breast_cancer/outputs/` | No | Ignored (`.gitignore:31`) | Correct |
| `benchmarks/outputs/`, `benchmarks/logs/`, `benchmarks/envs/{Rlib,c2l_venv}/` | No | Ignored | Correct (c2l_venv is a full venv with pytorch-lightning — confirmed ignored) |
| `data/` (46 MB incl. `V1_Breast_Cancer_Block_A_Section_1`, `.DS_Store`) | No | Ignored (`.gitignore:49`) | Correct; large local scratch, fine |
| `chimera_v1/`, `spatcar/` (legacy) | No | Ignored | Correct |
| `demo_demo/`, `demo_wiz/` | No | Ignored | Correct; unreferenced scratch — could be deleted locally |
| `src/tissueresolve.egg-info/` | No | Ignored | Correct |
| `docs/assets/*.png` (5 figures) | **Yes** | Tracked | OK — small doc images, intentional |
| `benchmarks/envs/Makevars`, `_install_music_patched.R` | Decide | **Untracked, NOT ignored** | **Action needed:** track (if part of install recipe) or add to `.gitignore` |
| `benchmarks/envs/Rlib/.../NEWS.md` etc. | No | Ignored via `Rlib/` | Correct |
| **No accidentally-tracked large/binary files found** | — | `git ls-files` shows only `.py/.tsv/.txt/.yaml/.R/.md/.png/.sh/.yml/.cff` | Clean — `*.h5ad/*.h5/*.mtx/*.zip/*.tar.gz` all ignored (`.gitignore:38-44`) |

---

### Table 4 — Redundancy candidates

| A | B | Overlap | Recommendation |
|---|---|---|---|
| `report/html.py` (+`templates.py`,`sections.py`) — the **live** results-dir report | `report/unified.py` + untracked `components.py`/`glossary.py`/`figures.py`/`style.py` — newer design-system | Both generate a single `report.html`; the new system duplicates CSS/section/figure-card logic. Only consumer of new system is uncommitted `07_generate_reports.py`. `report/__init__.py` exports only `html.*` and its docstring omits the new modules. | **HIGH: pick one report architecture.** Either migrate the package API to the design-system and retire `html/templates/sections`, or drop `unified/components/...`. Do not ship both. |
| `06_resolution_spillover_analysis.py` | `08_resolution_analysis.py` | Both run `compute_separability`, recommend merges, list unresolved families, write `outputs/resolution/{recommended_merges,unresolved_families,...}` | **MEDIUM: merge into one resolution step** (already noted in `docs/project_structure_audit.md:75-80` but not done) |
| `bulk/solver.py` (CHIMERA weighted-NNLS pipeline solver) | `solver/*` (nnls/weighted/marker/ridge/ensemble/auto backbones) | Two distinct "solver" concepts; `weighted_nnls` exists in both senses | **MEDIUM: rename/document** to remove ambiguity (e.g. `bulk/pipeline_solver.py` vs `solver/backbones/`) |
| `src/tissueresolve/benchmark/` (spillover) | `src/tissueresolve/spatial/benchmark.py` (synthetic Visium) + top-level `benchmarks/` | Three things named "benchmark"; none import each other | **MEDIUM: rename in-package `benchmark/`→`spillover/`**; clarify `spatial/benchmark.py` is synthetic-benchmark tooling |
| `reference/markers.py` | `reference/pairwise_markers.py` | Adjacent marker-selection modules | **QUESTION:** confirm pairwise is a distinct algorithm, not duplicated logic |
| `reference/build.py` | `reference/hierarchical_build.py` | Reference construction | **QUESTION:** confirm hierarchical build reuses (not re-implements) `build.py` |
| `plotting/summary_figures.py` | `plotting/bulk_plots.py` + `spatial_plots.py` | Multi-panel summaries may re-plot the same data | **LOW:** audit for duplicate figure code |

---

### Naming / consistency issues
- **"benchmark" is overloaded 3 ways** (in-package `benchmark/`, `spatial/benchmark.py`, top-level `benchmarks/`) — confusing for contributors.
- **"solver" is overloaded 2 ways** (`bulk/solver.py` vs `solver/`).
- **`report/__init__.py` is stale:** docstring says modules are "implemented in Stage 5" and lists only `html/methods_text/templates/sections/assets`; it never mentions `interpretation.py`, `unified.py`, or the new `components/figures/glossary/style`. Misleads readers about the real surface.
- **`api.py` overstates itself** ("thin wrappers… do not introduce new behaviour") while implementing resolution-mode routing, an auto→flat silent-downgrade-with-warning, and a parallel `_solver_bulk_result` path.

### Overclaiming flags (per the "be critical" mandate)
- `tuning/__init__.py` emits fabricated tuning scores — the single clearest scientific-honesty risk in the tree (currently dormant because unused).
- The new report layer is described in commit message `8c9c87e` ("Prepare TissueResolve documentation and GitHub release") and `47d31ac` ("Redesign reports…") as if complete, yet the redesigned modules are **uncommitted** and not wired into the package's public `generate_report`. Release-readiness is overstated relative to the committed state.

### Net assessment
Architecture is broadly sound and matches `DESIGN_SPEC`/CLAUDE.md intent (clean bulk/spatial separation, gitignored legacy/scratch/outputs, no large binaries tracked). The concrete debt is: (1) a fabricating orphan stub (`tuning/`), (2) a half-finished report-system migration with significant uncommitted work, and (3) three-way "benchmark" / two-way "solver" naming overload plus the documented-but-unmerged 06/08 overlap. None of these are BLOCKERs to the core algorithms, but the `tuning` stub and the uncommitted/duplicated report layer should be resolved before any "release."

## PART 2 — Scientific utility audit

Scope read: `src/tissueresolve/api.py`, `bulk/pipeline.py`, `bulk/hierarchical.py`, `bulk/solver` (via pipeline), `spatial/model.py`, `reference/hierarchy.py`, `reference/separability.py`, `reference/suitability.py`, `solver/{__init__,base,auto,nnls,weighted_nnls}.py`, `README.md`, `docs/accuracy_improvement_audit.md`, `docs/gene_masking_cv.md`, and the committed benchmark outputs under `benchmarks/outputs/`.

### (1) What scientific problem does it solve?

Reference-based **cell-type / cell-state deconvolution** of two modalities from a shared single-cell reference:
- **bulk RNA-seq** → RNA-derived (mRNA) proportions via marker selection + weighted NNLS (`src/tissueresolve/bulk/pipeline.py:217-220`), output explicitly typed `mRNA_proportion`, never silently converted to cell fractions (`bulk/pipeline.py:18-22, 242-245`).
- **10x Visium spatial** → spot-level composition via an NB multiplicative update with optional spatial smoothing (`src/tissueresolve/spatial/model.py:56-245, 368-396`).

The distinctive framing is **resolution-awareness**: instead of always forcing fine subtypes, it decides per broad family whether subtypes are separable and otherwise reports `unresolved_<family>` mass (`reference/hierarchy.py:675-772, 890-1003`). This is a real and honest scientific stance and the cleanest differentiator (see §4).

### (2-3) Better / worse than plain NNLS

The repo's *own* committed benchmark (`benchmarks/outputs/bulk/accuracy_metrics.tsv`, single real breast-cancer pseudobulk: 12 samples, 32 cell types, 5000 genes) is the ground truth here:

| method | fine Pearson | fine RMSE | source |
|---|---|---|---|
| `NNLS_baseline` | 0.7682150016166125 | 0.0328 | accuracy_metrics.tsv |
| `TissueResolve_auto` | **0.7682** (identical) | 0.0328 | best_method_summary.tsv / real_external_method_status.tsv |
| `TissueResolve_flat` | 0.6754 | 0.0450 | accuracy_metrics.tsv |
| `TissueResolve_hierarchical` | 0.0618 | 0.0513 | accuracy_metrics.tsv |
| `MuSiC` (executed) | 0.652 | 0.0356 | external_full_metrics.tsv |
| `BisqueRNA` (executed) | 0.583 | 0.0525 | external_full_metrics.tsv |

**HIGH — the "competitive solver" advantage over NNLS is not demonstrated; it is a tie by construction.** `TissueResolve_auto`'s 0.7682 equals `NNLS_baseline` to >15 significant figures because `AutoSolver.select` chose `nnls` (`benchmarks/outputs/bulk/solver_comparison.tsv`: nnls objective 0.95288 vs weighted_nnls 0.95273 — a 0.0002 margin) and then literally runs the NNLS backbone (`solver/auto.py:72-80` → `solver/nnls.py:15-23`). So "auto matches NNLS" is correct but the matching *is* NNLS.

What it does **worse** than plain NNLS (documented honestly in `docs/accuracy_improvement_audit.md:14-37`):
- its own default panel + weighting **lose** ~0.09 Pearson on clean pseudobulk (0.768 → 0.705 from 500-gene restriction, → 0.675 from wNNLS weighting). The doc correctly frames this as a robustness/clean-data trade-off — but **no benchmark in the repo demonstrates the claimed robustness benefit** (no protocol-mismatch / noisy / cross-platform experiment exists with ground truth). So the trade-off justification is asserted, not shown.
- hierarchical mode's headline fine accuracy (0.062) is far below NNLS; this is an artifact of scoring zeroed abstained columns (`accuracy_improvement_audit.md:69-75`), but it means the flagship mode cannot be advertised on accuracy.

What it plausibly does **better** (but only qualitatively): abstention/uncertainty, separability/spillover diagnostics, unified reference, and reporting. These are genuine, but none is a quantified accuracy win.

### (4) Clearest unique value proposition

**Honest, gated abstention at the right resolution.** The within-family resolvability gate (`reference/hierarchy.py:675-772`) combines three independent signals (mean Bhattacharyya-based separability, min discriminating genes, within-family spillover) and only splits a family into subtypes when all pass; otherwise mass stays at `unresolved_<family>` (`add_unresolved_family_mass` 613-661; partial resolution 829-849). Combined with a **single reference abstraction reused for bulk and spatial** (`api.py:73-236`; both modalities call `build_cell_type_hierarchy` + `assemble_hierarchical_estimates`), this is the defensible, differentiated contribution. The `reference/suitability.py` scorer (PASS/CAUTION/WARNING/FAIL with explicit UNKNOWN for unevaluable components, worst-component override at 251-257) is a good, conservative companion idea.

### (5) Is it doing too much?

**MEDIUM — yes.** 80+ source modules (full `src/tissueresolve` tree) span: 6 solvers + auto + ensemble, two hierarchy code paths (legacy `_BROAD_KEYWORDS` at `hierarchy.py:58-69` *and* `_FAMILY_KEYWORDS` at 79-108 — duplicated, divergent keyword tables), protocol detection/risk/mismatch, spatial graph + NB-CAR + neighbourhood, an 8-component suitability scorer, 12 plotting modules, a multi-file report layer, and a large benchmark harness with import/export adapters for ~13 external tools. The genuinely novel science is a small fraction of this; much (NNLS, ridge, marker NNLS, MuSiC/Bisque-style ideas) reimplements established methods without a demonstrated edge. The breadth dilutes the validation budget (one tiny dataset has to cover everything).

### (6) Essential features

- Reference builder + `ReferenceSignature` abstraction (`api.py:239-287`).
- Bulk wNNLS pipeline with explicit mRNA-proportion typing and gene-overlap guard (`bulk/pipeline.py:160-272`).
- Separability diagnostics (`reference/separability.py:64-179`) — the scientific safety net the package is built around.
- Hierarchy gating + unresolved-mass arithmetic (`reference/hierarchy.py:675-1003`).
- Spatial NB update with always-recorded smoothing λ/α (`spatial/model.py:99-114, 207-214`).
- Reference suitability scoring (`reference/suitability.py`).

### (7) Optional / experimental

- `solver=auto` / ensemble / gene-masking CV (`solver/auto.py`) — useful diagnostically but, on the only benchmark, selects NNLS by a 2e-4 margin and adds runtime (auto 13.9 s vs NNLS 1.7 s, `real_external_method_status.tsv`) for no accuracy gain.
- Partial hierarchical resolution (`hierarchy.py:829-849, 930-951`) — new, plausible, but unvalidated against ground truth at subtype level.
- Heuristic keyword family inference (`infer_broad_cell_type_family`, two competing tables) — explicitly a fallback; risky for non-breast/non-immune tissues and untested outside this domain.
- mRNA-content correction, bootstrap CIs, neighbourhood stats — sound but peripheral to the core claim.

### (8) Deferrable

- External-tool import/export benchmark adapters (most external methods are `skipped`/`failed`/`exported_not_run` in `composite_scores.tsv`; BayesPrism `failed` after 2180 s, CIBERSORTx `exported_not_run`).
- The 12 plotting modules and multi-file HTML report stack — publication polish, not science.
- Ensemble solver and the second (legacy) hierarchy keyword table (dead/duplicated path).

---

### Is the core claim defensible?

Claim under audit: *"TissueResolve is a reference-based, resolution-aware deconvolution framework that combines competitive solvers with broad-to-fine hierarchical inference, reference suitability diagnostics, separability/spillover awareness, and publication-ready reports."*

| Sub-claim | Verdict | Evidence |
|---|---|---|
| reference-based | **Defensible** | `api.py:239-287`; `ReferenceSignature` reused bulk+spatial |
| resolution-aware | **Defensible (design), unvalidated (benefit)** | gating `hierarchy.py:675-772`; but hierarchical fine Pearson 0.062, accuracy_metrics.tsv |
| combines **competitive** solvers | **Overclaim** | auto selects single backbone (no blending in auto path, `solver/auto.py:72-80`); ties NNLS by construction; "competitive" vs external is shown only against MuSiC/Bisque, both of which it beats but which are weaker than plain NNLS here |
| broad-to-fine hierarchical inference | **Defensible (implemented), not shown to help accuracy** | `bulk/hierarchical.py`, `hierarchy.py` |
| reference suitability diagnostics | **Defensible** | `reference/suitability.py` (conservative, UNKNOWN-aware) |
| separability/spillover awareness | **Defensible** | `reference/separability.py:64-179`; within-family spillover gate |
| publication-ready reports | **Plausible (not audited for science here)** | report/plotting modules present |

**Net:** The *architecture* claim is largely true. The *performance/competitiveness* implications a reader will infer are **not** supported: on the package's own data the best result equals plain NNLS, the flagship hierarchical mode trails badly on accuracy, and the robustness rationale for the accuracy sacrifice is undemonstrated.

### What is missing before publication (BLOCKER/HIGH)

1. **BLOCKER — Validation breadth.** One dataset (12 bulk samples / 32 types / 5000 genes / 1915 cells capped 60/type per `input_summary.tsv`; spatial = 36 synthetic-truth spots, no real ground truth per `spatial/accuracy_metrics.tsv` + README). Need multiple references, multiple tissues, and ≥ several dozen mixtures with significance testing. With n=12 the 0.768 vs 0.675 gap has no confidence interval anywhere in the repo.
2. **BLOCKER — Demonstrate the robustness payoff.** The entire justification for losing clean-data accuracy (`accuracy_improvement_audit.md:21-37`) is robustness to protocol mismatch/spillover. There is **no** ground-truth experiment showing TissueResolve_flat/auto beating NNLS under mismatch/noise. Without it, the weighting/panel design is a net negative on every benchmark present.
3. **HIGH — Validate hierarchical mode as *output*, not just abstention.** Report subtype-level accuracy on resolvable families with ground truth (the `fine_metrics_resolvable_only` view exists per audit §6 but is not foregrounded), and quantify the precision/recall of the unresolved decision (does it abstain on truly non-separable families and resolve genuinely separable ones?).
4. **HIGH — Drop or fix overclaiming language.** "combines competitive solvers" (auto = single backbone), and any framing implying superiority over RCTD/cell2location/CIBERSORTx etc. (all `skipped`/`exported`/`failed` in `composite_scores.tsv`). The README table's caveat is good but the prose and the audited claim still imply more than the data shows.
5. **MEDIUM — Composite-scoring transparency.** `composite_scores.tsv` ranks `TissueResolve_auto` #1 partly via subjective `interpretability/usability/resolution_awareness` columns (0.8-1.0 for TissueResolve, 0.1-0.4 for baselines) and `runtime_resource_score=0.0` for all TissueResolve modes — these hand-assigned weights, not measured science, drive the ranking. This must be disclosed or removed before publication.
6. **MEDIUM — Consolidate duplicated hierarchy heuristics** (`_BROAD_KEYWORDS` vs `_FAMILY_KEYWORDS`) and test family inference outside breast/immune tissue, or restrict the package to explicit user mappings for publication claims.

### Bottom line
The science that is *unique* (gated abstention, unified reference, suitability/separability diagnostics) is sound and honestly implemented, with non-negotiable rules visibly respected in code (explicit mRNA-proportion typing, recorded smoothing, non-suppressed warnings). But the package's quantitative case is currently a tie with plain NNLS plus a flagship mode that trails on accuracy, all on a single tiny dataset. It is publishable as a *diagnostics/abstention framework* with honest framing; it is **not** publishable as an *accuracy-competitive deconvolution method* on the present evidence.

## PART 3 — Algorithmic Audit

Scope: files actually read — `src/tissueresolve/reference/{build,hierarchy,separability,suitability}.py`, `src/tissueresolve/solver/{base,nnls,weighted_nnls,marker_nnls,ridge_nnls,auto,ensemble}.py`, `src/tissueresolve/validation/gene_masking.py`, `src/tissueresolve/spatial/{model,auto_params}.py`, `src/tissueresolve/bulk/{solver,pipeline}.py`, `benchmarks/shared/{metrics,composite_score}.py`, `benchmarks/bulk/methods/{tissueresolve,nnls_baseline}.py`, `benchmarks/run_real_external_benchmark.py`, `benchmarks/spatial/run_spatial_benchmark.py`.

### Per-method table

| method | scientific rationale | implementation status | assumptions | failure modes | validation status | recommendation |
|---|---|---|---|---|---|---|
| `ReferenceBuilder._aggregate` (build.py) | Donor-aware mean-of-per-donor-means + cross-donor CV; CPM + L1-norm `phi`; spatial `R_cpm/R_log` | Implemented, validated inputs | Raw counts in; ≥2 donors for CV; per-cell-type CPM is a valid signature | Equal-donor-weight ignores donor cell counts; CPM normalization of an already-mean profile loses absolute mRNA-content (correct by design, but downstream `phi` cannot recover fractions) | Unit-level only (no metric here) | OK. Document that donor averaging is unweighted by donor size (LOW). |
| `_estimate_overdispersion` (build.py) | MoM NB dispersion `φ=μ²/(σ²−μ)`, clip [0.5,100] | Implemented | Genes independent; pooled across ALL cells (not per-type) | **Pools dispersion over the whole matrix, so between-cell-type mean differences inflate variance → underestimates `φ` (overstates dispersion)**; sub-Poisson genes forced to `_PHI_G_MAX` (Poisson) | Not validated against known NB | MEDIUM: dispersion is estimated on the mixed population, not within cell type; this is a biased estimator for the per-type NB used in the spatial model. Flag in docs. |
| `compute_separability` (separability.py) | Bhattacharyya coeff on L1-normalized CPM, Jeffreys, Pearson, discriminating-gene count; warns + optional raise | Implemented, warns properly | BC on CPM profiles reflects identifiability | BC computed on full-gene CPM (housekeeping genes dominate → BC biased high/optimistic about similarity but the warning direction is conservative) | Self-consistent; thresholds (0.90/0.97) are heuristic, uncited | OK, but thresholds are arbitrary constants (LOW). |
| `merge_nonseparable_types` (separability.py) | Union-find grouping; cell-count-weighted CPM merge; re-derive `phi`; copy global `phi_g` | Implemented; documents the SpatCAR `phi_g` bug-fix | Transitive merging via union-find can chain weakly-linked types | **Transitivity: A~B and B~C (each BC>thr) merges A+B+C even if BC(A,C) is low** — can over-merge | Bug-fix asserted but the fix is for a *different* (legacy) bug; no regression test cited here | MEDIUM: union-find merging is transitive and can collapse a chain of distinct types; consider connected-component review or warn on chain length. |
| `infer_broad_cell_type_family` / hierarchy keyword maps | Heuristic fine→broad family from substring keywords | Implemented, warns when inferred | English keyword substrings; first-match-wins ordering | Silent mis-mapping for non-standard labels; `dc`/`caf`/`nk` substrings can false-match (e.g. "dc" in arbitrary IDs) | Warns + `Other` fallback | OK (it warns and prefers explicit mapping). Substring matching is fragile (LOW). |
| Hierarchy mass arithmetic (`combine_*`, `add_unresolved_family_mass`, `estimate_partial_subtype_resolution`) | `fine=family×P(sub\|family)`; unresolved mass parked in `unresolved_<fam>` | Implemented; claims mass preservation | family_props rows sum to a fixed total; conditional sums to 1 within family | See Q6 — **partial-resolution path uses confident-subtype conditional sum, which can leave residual ≠ proper complement if confident subset conditional > family**; generally consistent but unverified at code level | No numeric mass-conservation assertion in these functions | MEDIUM: add an explicit per-row mass-conservation check/test (see Q6). |
| `evaluate_within_family_resolvability` (hierarchy.py) | 3-signal gate (mean 1−BC, min discriminating genes, mean spillover=max Pearson) | Implemented | Three thresholds (0.10, 10 genes, 0.30) chosen heuristically | All three thresholds uncited/arbitrary; `mean_spillover` uses **max positive Pearson per pair averaged** — conflates magnitude with sign | Not externally validated | MEDIUM: thresholds are magic numbers; document provenance, expose as config. |
| `compute_reference_suitability_score` (suitability.py) | Weighted mean of 8 component scores; worst-FAIL override; UNKNOWN never =PASS | Implemented, honest about UNKNOWN | Component-to-score mappings are heuristic (e.g. risk low/med/high→1/0.6/0.3) | Score is a weighted average of mostly heuristic 0..1 maps; gene_overlap "≥60%→full marks" is arbitrary | No calibration; purely advisory | OK as advisory; **must not be read as a probability** (LOW/QUESTION). |
| `NNLSSolver` (nnls.py) | Plain NNLS on `R_cpm.T`, per-sample, normalize | Implemented | Linear mixing in CPM space; library-size invariance via post-normalize | Identical to `NNLS_baseline` (see Q3) | Used as auto candidate | OK. |
| `WeightedNNLSSolver` (weighted_nnls.py) | Row-weight `√max(row-normalized R)`; specificity weighting | Implemented | Specificity weight = per-gene max relative expression | Weight derived **only from reference**, ignores query noise/protocol; `cond` reported on unweighted `R` (mislabeled) | Auto candidate | LOW: reported `condition_number` is for unweighted R, not the solved weighted system. |
| `MarkerNNLSSolver` (marker_nnls.py) | Top-N one-vs-rest log-CPM markers, then NNLS | Implemented | Marker = top (log µ_k − mean log others) | Marker selection uses the **whole** reference; in CV the ref is subset to train genes so it is leak-safe (verified) | Auto candidate | OK. |
| `RidgeNNLSSolver` (ridge_nnls.py) | Ridge closed-form via pinv + clip≥0 + normalize | Implemented | L2 on weights; non-negativity by clipping (not constrained solve) | **Clip-then-normalize is not a true NNLS**; can produce biased solutions vs constrained QP; `λ` scaled by trace/K | Auto candidate | LOW: "ridge_nnls" is clip-projected ridge, not constrained NNLS — name is mildly misleading. |
| `AutoSolver` (auto.py) | Pick backbone by gene-masking CV minus conditioning penalty | Implemented | CV masked-gene reconstruction predicts proportion accuracy | See Q2 — reconstruction-quality proxy; tiny conditioning penalty (0.05·tanh) rarely changes ranking | No ground-truth validation of the selection criterion | HIGH (see Q2). |
| `EnsembleSolver` (ensemble.py) | CV-weighted convex combo of backbones | Implemented | Same as Auto; weights ∝ clipped CV Pearson | Same circularity as Q2; weights can collapse to uniform when all CV≤0 | Not validated | HIGH (inherits Q2). |
| `WNNLSSolver` (bulk/solver.py) + `BulkPipeline` | Protocol/marker-aware weighted NNLS on L1-normalized `phi`; bootstrap; mRNA-correction separate | Implemented; honors mRNA/cell-fraction separation | Gene selection + composite weights add value over plain NNLS | This is the **real** TissueResolve bulk method — but `TissueResolve_auto` in benchmarks does NOT use it (Q3/Q4) | R² coverage reported; bootstrap available | OK; but mislabeled in benchmark (Q4). |
| `SpatCARModel` (spatial/model.py) | NB multiplicative-update MAP + spatial proximal mixing + mismatch factors | Implemented; smoothing always recorded; non-convergence warns | NB per-gene `φ_g`; `α=min(0.5,λ)` mixing with row-normalized adjacency | Multiplicative update **claims monotone NB-NLL decrease, but the spatial mixing + mismatch updates break that guarantee** (the docstring's "non-increasing per step" applies only to the bare NB step); `_nnls_init` uniform `φ_g=1` fallback | Convergence trace recorded; no likelihood reported per iter | MEDIUM: docstring overstates monotonicity for the *combined* iteration (Q misc). |
| `select_lambda_spatial` (auto_params.py) | Sweep λ; objective = reconstruction − w·over-smoothing | Implemented; explicitly disclaims accuracy | Over-smoothing penalty = fractional variance loss vs λ=0 | Reconstruction favors λ=0 (no smoothing always reconstructs observed better); penalty also references λ=0 → objective can be near-degenerate; `oversmooth_weight=0.3` arbitrary | No ground truth (correctly disclaimed) | MEDIUM/QUESTION (Q8). |
| `benchmarks/shared/metrics.py` | Accuracy (Pearson/Spearman/RMSE/MAE/bias), per-type, dominant, family/unresolved-aware | Implemented; flatten-then-correlate | Flattened obs×type Pearson is meaningful | **Flattened global Pearson is inflated by between-cell-type mean differences** (Q10); dominant-accuracy is asymmetric (top-1 est ∈ top-k true) | These are the scoring functions | HIGH (Q10). |
| `compute_composite_scores` (composite_score.py) | Weighted multi-dimension score; non-GT drops accuracy | Implemented | Per-method `_PRIORS` are valid stand-ins for robustness/usability/etc. | **Non-accuracy dims are hardcoded author priors that favor TissueResolve** (Q9); accuracy fed as raw Pearson (can be negative, unclipped) | No external validation of priors | HIGH (Q9). |

---

### The 10 questions

**Q1 — Are the reference signatures scientifically sound?**
Mostly. CPM + L1-`phi` for bulk and `R_cpm/R_log` for spatial are standard. Two caveats: (a) donor averaging is *unweighted* by per-donor cell count (`build.py:322-326`), so a 5-cell donor counts as much as a 5000-cell donor; (b) `_estimate_overdispersion` pools across the entire matrix (`build.py:430-448`), not within cell type, biasing the NB `φ_g` used by the spatial model. Both are defensible defaults but should be documented. **Severity: MEDIUM (overdispersion), LOW (donor weighting).**

**Q2 — Does gene-masking CV risk selecting a solver that reconstructs genes well but estimates proportions poorly? — YES (HIGH).**
The CV criterion (`gene_masking.py` + `auto.py`) scores a solver by how well `proportions · R_ref[masked]` reconstructs the *held-out observed* expression, scaled per-sample. This is a **reconstruction objective, not a proportion-accuracy objective**, and the two are not equivalent:
- The reconstruction uses the *same* reference `R_cpm` the solver fits against (`reconstruct_masked_genes` calls `ref.subset_genes(masked).as_R_cpm()`), so a solver that drives proportions toward whatever linear combination best matches the bulk's marginal gene profile is rewarded — even if that combination splits two collinear cell types arbitrarily (which is exactly the non-separability failure mode the package elsewhere warns about). High masked-gene Pearson is achievable with badly wrong but mutually-compensating proportions.
- The conditioning penalty is negligible: `0.05·tanh(log10(cond+1)/3)` ∈ [0, ~0.05), far smaller than typical Pearson differences, so it almost never overrides reconstruction.
- Per-sample relative scaling (`obs/obs.sum`) removes library size but not the dominant between-gene mean structure, so the Pearson is further inflated and insensitive to proportion errors.
**Conclusion:** the proxy can and will prefer solvers that reconstruct genes while mis-estimating proportions, especially for collinear/non-separable references. This selection criterion is unvalidated against any ground-truth proportion accuracy. Recommend adding a synthetic-mixture validation showing CV-rank correlates with true proportion RMSE before relying on `auto`.

**Q3 — Are NNLS and `TissueResolve_auto` essentially the same for bulk? — Effectively YES in the common case (HIGH).**
`NNLS_baseline` (`benchmarks/bulk/methods/nnls_baseline.py`) = plain `scipy.nnls` on `ref.as_R_cpm().T`, normalized. `AutoSolver`'s candidate set includes `NNLSSolver` (`solver/nnls.py`) which is **byte-for-byte the same computation** (`R = sub.as_R_cpm().T`; `nnls(R, B[:,j])`; `_normalize`). When auto's CV picks `nnls` (a frequent outcome on clean references where plain NNLS reconstructs best), `TissueResolve_auto` ≡ `NNLS_baseline`. In `run_real_external_benchmark.py:126-129` both are even instantiated side-by-side (`AutoSolver` vs `NNLSSolver`). The only difference is auto *might* pick weighted/marker/ridge — but all four candidates operate on raw `R_cpm` with no protocol-awareness, marker selection, or gene weighting.

**Q4 — What does TissueResolve add beyond NNLS? — In the benchmarked `_auto` path, almost nothing; the real value lives in a code path the benchmark does NOT exercise (HIGH).**
The scientifically substantive bulk method is `BulkPipeline`/`WNNLSSolver` (`bulk/pipeline.py`, `bulk/solver.py`): protocol-risk gene exclusion, marker selection, composite gene weights, L1-normalized `phi`, bootstrap CIs, R² coverage, explicit mRNA-proportion-vs-cell-fraction separation. **But `TissueResolve_auto` bypasses all of it** — it calls `AutoSolver` over the lightweight `solver/*.py` backbones on raw CPM with no marker selection and no weighting. So the headline benchmarked "improved" method does not include TissueResolve's actual differentiators. `TissueResolve_flat`/`_hierarchical` *do* call `tr.deconv_bulk` (the real pipeline), so the package's value is real — it is just **not** what `_auto` measures. This is a labeling/benchmark-design problem: `_auto` over-claims to be "improved TissueResolve" while being a thin NNLS-family selector. **Recommend either routing `_auto` through `BulkPipeline` or renaming it to reflect that it is a solver-backbone selector, not the protocol-aware pipeline.**

**Q5 — Is the separability/merge logic sound?**
BC-based screening is reasonable and the warnings are non-suppressible (good). Two issues: (1) BC is computed on *all* genes' CPM, so ubiquitously-expressed genes dominate and the coefficient is an optimistic measure of similarity for rare-marker-distinguished types; (2) `merge_nonseparable_types` uses union-find, which is **transitive** — a chain A~B~C each above threshold merges all three even when A and C are well separated (`separability.py:231-233`). The documented "bug fix" addresses a *legacy* `phi_g` loop, not these issues. **Severity: MEDIUM.**

**Q6 — Is unresolved mass mathematically consistent (mass preservation)?**
Largely yes, with one unverified branch. The simple path (`combine_family_and_conditional_estimates` + `add_unresolved_family_mass`) preserves mass: `fine = family×conditional`, conditional sums to 1 within family (`hierarchy.py:570-583`, with uniform split when family signal is zero), so resolved-fine + unresolved equals the family total. The **partial-resolution path** (`estimate_partial_subtype_resolution`, `hierarchy.py:829-849`) computes `residual = family_mass·(1 − Σ_confident conditional)`; this is mass-consistent *only because* conditional sums to 1 over all members, so the residual is exactly the non-confident members' conditional mass. That holds, but: (a) there is **no runtime assertion** that combined rows sum to the family total (or to 1); (b) the docstrings on `HierarchicalEstimates` assert "fine + unresolved sum to 1 per row" but this is never checked in code and depends on `family_props` rows summing to 1 (true for solver output, but not enforced for arbitrary inputs). **Conclusion: consistent by construction, but unguarded — add a per-row mass-conservation assertion/test. Severity: MEDIUM.**

**Q7 — Is the spatial NB/CAR model sound?**
The NB multiplicative update (`_nb_multiplicative_update`) is a standard MM-style non-negative update and is per-element non-negative + simplex-renormalized. **However the module docstring claims the step is "non-increasing in NB-NLL per step" — that guarantee (if it holds for the bare NB step) does not survive the spatial proximal mixing (`(1−α)π_NB + α·Aπ`) and the periodic mismatch-factor updates**, which are not part of the same descent. So the model is a reasonable heuristic alternating scheme, but the stated monotonic-descent property is **overclaimed** for the actual combined iteration. Convergence is monitored and non-convergence warns (good). The `φ_g=1` uniform fallback when overdispersion was not estimated is honestly warned. **Severity: MEDIUM (docstring overclaim).**

**Q8 — Are spatial metrics appropriate without ground truth? — Mostly yes, with one degenerate objective (MEDIUM/QUESTION).**
The benchmark correctly **refuses to claim accuracy without ground truth** (`run_spatial_benchmark.py` only computes accuracy when `has_gt`; `composite_score` drops the accuracy dimension when `not has_ground_truth`). No-GT metrics used — cross-method Pearson concordance, proportion entropy, near-zero fraction — are appropriate descriptive/stability measures and are not dressed up as accuracy. Moran's I exists in `spatial/qc.py` (referenced) and is a legitimate spatial-autocorrelation descriptor. The weak point is `select_lambda_spatial`: its objective `reconstruction − 0.3·over_smoothing_penalty` is **partly self-referential** — both the reconstruction (best at λ=0, which reproduces observed) and the penalty are measured relative to the λ=0 baseline, so the procedure is biased toward small λ and the `0.3` weight is an unjustified constant. It is honestly disclaimed as "accuracy not claimed," so it is acceptable as a heuristic but should not be presented as principled tuning. **Concordance ≠ correctness** (all methods can agree on a wrong answer) and the report should keep stating that.

**Q9 — Are composite scores defensible? — Partly; the non-accuracy dimensions are not (HIGH).**
The accuracy dimension and the renormalization-when-no-GT logic are defensible. The problem is the other five dimensions: `compute_composite_scores` fills `robustness/usability/interpretability/resolution_awareness` from a **hardcoded `_PRIORS` dict** (`composite_score.py:24-35`) when the per-method table doesn't supply them — and `run_real_external_benchmark.py` does **not** supply `robustness`, so those dimensions are entirely author-assigned constants. Those constants systematically favor TissueResolve (e.g. `TissueResolve_hierarchical` gets interpretability=1.0, resolution_awareness=1.0, robustness=0.85; `NNLS_baseline` gets 0.4/0.1/0.5). With 65% of the default weight on non-accuracy dimensions (0.20+0.15+0.15+0.10+0.05), **a method's rank is largely predetermined by its prior**, independent of measured performance. This is a self-fulfilling benchmark. Defensible only if the priors are explicitly disclosed as subjective and the accuracy-only ranking is shown alongside. **Recommend: report accuracy-only and composite separately, and label priors as subjective.**

**Q10 — Any misleading metrics? — Yes (HIGH).**
1. **Flattened global Pearson** (`accuracy_metrics`, `metrics.py:39-52`): correlating the entire flattened obs×cell-type matrix is inflated by the large between-cell-type mean structure (a method that just predicts each type's grand-mean proportion scores high Pearson). It overstates "accuracy"; per-cell-type Pearson (already present) is the honest view. RMSE/MAE are fine.
2. **Composite accuracy = raw Pearson**: in `run_real_external_benchmark.py:131` the Pearson r is fed directly as the 0..1 `accuracy_score` with no clipping; a negative correlation would enter the weighted sum as a negative contribution, and even valid r∈(0,1) is not a calibrated accuracy.
3. **`unresolved_precision/recall` against `set()`**: in both run scripts `hierarchical_fair_metrics(..., high_risk_families=set())` is called with an **empty** high-risk set, so precision/recall are NaN/degenerate — the "fair abstention" reward is wired but not actually evaluated.
4. **Reported `condition_number`** in `weighted_nnls`/`ridge_nnls` is computed on the *unweighted* `R`, not the system actually solved — mildly misleading diagnostic.
5. **`reconstruction_pearson`/masked-gene Pearson** are presented as quality signals but, as in Q2, do not track proportion accuracy.

---

### Summary of severities
- **HIGH:** Q2 (CV reconstruction proxy ≠ proportion accuracy, unvalidated), Q3/Q4 (`TissueResolve_auto` ≈ NNLS and bypasses the real pipeline; benchmark over-labels it), Q9 (composite non-accuracy dims are author priors favoring TissueResolve), Q10 (flattened Pearson + raw-Pearson-as-accuracy + empty high-risk set).
- **MEDIUM:** overdispersion pooled across types (Q1), union-find transitive over-merge (Q5), unguarded hierarchical mass conservation (Q6), spatial NB monotonicity overclaim (Q7), self-referential λ selection (Q8), arbitrary resolvability thresholds.
- **LOW:** unweighted donor averaging, ridge "NNLS" is clip-projected, mislabeled condition numbers, fragile keyword family inference, advisory-only suitability score read as a number.

No code was modified (read-only audit).

## PART 4 — Code quality audit

Scope: sample-read `api.py`, `cli.py`, `config.py`, `results.py`, `solver/{base,auto}.py`, `bulk/solver.py`, `spatial/model.py`, `plotting/{palette,style,export}.py`, `utils.py`, `presets.py`, plus targeted greps across all of `src/tissueresolve/` and `benchmarks/`. The `benchmarks/envs/` tree is a vendored virtualenv (site-packages) and was excluded from quality assessment.

Overall the code is unusually disciplined for a scientific package: deterministic RNG via `np.random.default_rng(seed)` everywhere, scoped `warnings.catch_warnings()` (no global suppression), `logging` used consistently (~109 hits, no stray `print`/`plt.show`), enforced estimate-type constants, optional deps lazily imported with actionable errors, and "every figure ships its source data" enforced in `plotting/style.py` and `plotting/export.py`. Findings below are mostly correctness/consistency issues at the orchestration layer (CLI presets), not in the core algorithms.

### BLOCKER

None. No bare `except:`, no silent gene removal, no hidden warnings, no missing seeds in the core paths.

### HIGH

| Severity | Area | File:Line | Issue |
|---|---|---|---|
| HIGH | Preset → config wiring | `cli.py:223-233`, `presets.py:5-49` | `_configure_from_preset` only consumes `bootstrap`, `n_bootstrap`, and `hierarchical`. The preset keys `auto_tune`, `plots`, and `spatial` (e.g. `{"lambda": "auto", "n_neighbors": 8/12}`) are **defined but never read** on the `tissueresolve run` path. So `--preset publication`/`diagnostic` advertise `auto_tune: True` and per-preset spatial λ/neighbour settings that have no effect. This is a silent assumption that the preset does what it says. |
| HIGH | Bootstrap gating bug | `cli.py:223-233` + `cli.py:428` + `config.py:204` | Default `BootstrapConfig.n_bootstrap = 200`. `_configure_from_preset` only *lowers* nothing — when a preset has `bootstrap: False` (`quick`, `standard`) it does **not** set `n_bootstrap=0`, so `_execute_bulk` calls `deconv_bulk(..., n_bootstrap=cfg.bootstrap.n_bootstrap)` = 200 and the pipeline runs 200 bootstrap iterations despite the preset declaring `bootstrap: False`. The `standard` preset therefore silently runs full bootstrap (slow) and contradicts its own declaration. Should set `cfg.bootstrap.n_bootstrap = 0` when `not preset_params.get("bootstrap")`. |

### MEDIUM

| Severity | Area | File:Line | Issue |
|---|---|---|---|
| MEDIUM | Global RNG mutation | `utils.py:119-132`, `spatial/model.py:162` | `set_random_state` calls `np.random.seed(seed)` (mutates the legacy global state) in addition to returning a `Generator`. `SpatCARModel.fit` calls it on every fit, so fitting the spatial model has a global side effect on any caller using `np.random.*`. Reproducible but a hidden global mutation; prefer returning the `Generator` only. |
| MEDIUM | Silent `except Exception: pass` in orchestration | `cli.py:114-115`, `cli.py:202-203`, `cli.py:362-363` | Protocol detection, `resolution_mode` metadata stamping, and h5ad annotation probing swallow all exceptions. The CLAUDE.md rule is "never fail silently". These are non-fatal best-effort paths, but a broad `except Exception: pass` can mask real bugs (e.g. a malformed h5ad). At minimum log at DEBUG. `_detect_hierarchy_availability` (`cli.py:263`) returning `(False, None)` on any exception can silently downgrade a hierarchical run to flat. |
| MEDIUM | Duplicated R²/reconstruction logic | `api.py:55-59`, `bulk/solver.py:381-396`, `spatial/auto_params.py` | Per-sample R² is hand-rolled in `_solver_bulk_result` (`api.py`) and again in `bulk/solver.py::_r2_per_sample`, with slightly different NaN/zero handling (`api.py` returns NaN when `ss_tot==0`, solver returns `0.0`). One shared helper would avoid divergence. |
| MEDIUM | Duplicated TSV-with-comments writer | `results.py:49-54`, `plotting/style.py:258-263`, `plotting/export.py:83-87` | Three near-identical `_write_tsv` implementations. Low risk but a maintenance smell; consolidate into `utils` or `io`. |
| MEDIUM | Untyped orchestration signatures | `cli.py:236,268,313,386,399,439,838`; `api.py` (`cfg`, `ref`, `bulk` params) | Many CLI/api helpers take `cfg`, `ref`, `result` with no annotations and return un-annotated objects, weakening the "clear API" goal. `deconv_bulk`/`deconv_spatial` accept `bulk`, `ref`, `Y` positionally with no type hints and `**kwargs` forwarding (hard to discover valid kwargs without reading the pipeline). |
| MEDIUM | `__post_init__` dead branch | `results.py:657-661` | `SpatialDeconvResult.__post_init__` has an `if self.n_smooth is not None and self.n_smooth > 0:` block whose body is only `pass` with a comment — dead code that reads as if validation is enforced when it is not. |

### LOW / quick wins

| Severity | Item | File:Line |
|---|---|---|
| LOW | `deconv_bulk` resolution-mode validation message lists `flat` but the canonical set check excludes it after remapping — message says "auto, hierarchical, flat, none, suggest" while internal set is `{none,suggest,auto,hierarchical}` (works, but the error text/`Choice` lists differ between `api.py:115`, `cli.py:41`, `cli.py:506`). | `api.py:111-116,189-194` |
| LOW | `_solver_bulk_result` is dense one-liners with inline `;` statements (`api.py:53`, `46-52`) — hard to read; split for clarity. | `api.py:33-70` |
| LOW | `MRNAContentCorrector.from_reference_counts` per-cell-type loop builds index sets in Python (`valid = [c for c in cell_ids if c in total_umi.index]`) — O(cells²)-ish on large refs; vectorize via groupby. | `bulk/solver.py:310-320` |
| LOW | `config.heuristic_thresholds()` omits the `bootstrap_cv_*`, `marker_recall`, and all spatial-QC thresholds it documents as heuristic, so `run_metadata.json` under-reports heuristic provenance vs. the docstring claim. | `config.py:425-448` |
| LOW | `infer_cell_type_family` keyword list (`palette.py:71-87`) is order-sensitive and English-only; `"basal"`→epithelial could misclassify `basal` immune contexts. Deterministic and documented, but a fragile heuristic worth a docstring caveat. | `palette.py:71-101` |
| LOW | Report embeds figures as base64 PNG inline (`report/html.py:526-527`) — fine for a few panels, but for `--preset diagnostic` ("plots: all") this can produce very large self-contained HTML. No size cap / external-asset option. | `report/html.py:526` |

### Files needing docstrings / type hints
- `cli.py` — internal helpers (`_run_top_level`, `_execute_bulk/spatial`, `_load_reference_signature`) lack return-type annotations and parameter types for `cfg`/`ref`.
- `api.py` — public `deconv_bulk`/`deconv_spatial`/`build_reference` have good docstrings but untyped `bulk`/`ref`/`Y` and opaque `**kwargs`; the accepted kwargs should be enumerated or typed.
- `solver/base.py` — `BaseSolver.solve` and `align_query_to_reference` use `pd.DataFrame`/`ReferenceSignature` but `_normalize`/`_frame` are fine; minor.

### Files needing simplification
- `api.py::_solver_bulk_result` (compressed multi-statement lines, inline reconstruction loop).
- `cli.py` (921 lines, 25+ functions; the `run` flow spreads resolution-mode logic across `_select_resolution_mode`, `_detect_hierarchy_availability`, `_resolve_hierarchy_mapping`, `_configure_from_preset`, and `_execute_*` with overlapping h5ad-probing — consolidate the hierarchy-detection paths).
- `reference/hierarchy.py` (1003 lines, 23 funcs) and `spatial/benchmark.py` (704 lines, 23 funcs) are the largest modules; not read in full but flagged as candidates for module-splitting given function density.

### Checklist answers

- **Does every public function validate inputs?** Mostly. Core entry points validate well: `WNNLSSolver.solve` (empty-overlap raise + <20-gene warning), `BulkPipeline.run`/`SpatialPipeline.run` call `validate_counts_matrix`/`check_gene_overlap` and raise on empty overlap, `MRNAContentCorrector` rejects non-positive/uniform content, `ReferenceSignature.validate()` checks all shapes, `build_reference`/`generate_report`/`plot_results` raise clear `TypeError`/`ValueError` on bad source/result types, `resolution_mode` is validated against an allowed set. Gap: CLI helpers and `deconv_*` wrappers do little of their own validation (they delegate), and `_solver_bulk_result` does no shape validation before matrix multiply.
- **Are warnings actionable?** Yes — consistently. Warnings name the cause and the fix (e.g. gene-naming hint in `WNNLSSolver`, `estimate_overdispersion=True` hint in `SpatCARModel.fit`, "provide a fine→broad mapping" in `deconv_bulk`, kaleido install hint in `export_figure`).
- **Are errors user-friendly?** Yes. `ImportError` messages point at the exact extra (`pip install 'tissueresolve[report]'`), hierarchy errors give exact column/flag remedies, file-not-found and empty-table errors are explicit.
- **Outputs reproducible?** Yes. Every result object writes `metadata.json` with versions/estimate-type/heuristic flags; figures always co-write source-data TSVs; estimate-type warnings embedded in every TSV header. Caveat: heuristic-threshold provenance in `metadata` is incomplete (see LOW on `heuristic_thresholds()`).
- **Seeds controlled?** Yes. All RNG via `np.random.default_rng(seed)` with seeds from `BulkSolverConfig`/`SpatialSolverConfig`/`BootstrapConfig` (default 42/0). Only caveat: `set_random_state` additionally mutates the global `np.random.seed` (MEDIUM).
- **Color maps deterministic?** Yes. `palette.QUALITATIVE_BASE`/`FAMILY_RAMPS` are fixed; `assign_family_palette` and `build_hierarchical_color_map` assign by sorted order and can lock to an existing map for cross-run stability; `style.PALETTE`/`color_sequence` cycle a fixed list. No random color assignment.
- **Reports reproducible?** Yes — report reads saved tables/figures and embeds them; figures are reproducible from saved data + metadata. Caveat: no size cap on base64-embedded figures for `plots: all` (LOW).
- **Large outputs avoided by default?** Partially. Spatial bootstrap CIs are off unless requested; figure formats are bounded. **But** the bulk `run` path does not honor the preset's `bootstrap: False`, so `standard`/`quick` runs perform 200 bootstrap iterations by default (HIGH) — a large/slow output produced when the preset says it should be skipped.

## PART 5 — Test Audit

Scope read: `pyproject.toml` `[tool.pytest.ini_options]`; `tests/conftest.py`; `examples/real_breast_cancer/tests/conftest.py`; all test files enumerated under `tests/`, `examples/real_breast_cancer/tests/`, `benchmarks/tests/`; plus source under `src/tissueresolve/solver/`, `src/tissueresolve/spatial/model.py`, and the numbered scripts in `examples/real_breast_cancer/scripts/`.

### Test inventory

- **Total project tests:** 882 (`grep -rc "def test_"` over `tests`, `examples/real_breast_cancer/tests`, `benchmarks/tests`, excluding `benchmarks/envs/` venv site-packages, which polluted the naive count).
- **pytest config** (`pyproject.toml:87-89`): `testpaths = ["tests", "examples/real_breast_cancer/tests", "benchmarks/tests"]`, `addopts = "-v --tb=short"`. No custom markers registered, no `--run-real-data` flag wired into pytest itself (the flag is a CLI arg on `00_download_data.py`, not a pytest option).

#### Largest test files (sanity of distribution)
| File | tests |
|---|---|
| `tests/shared/test_protocol.py` | 67 |
| `tests/test_results.py` | 52 |
| `tests/shared/test_reference.py` | 43 |
| `examples/real_breast_cancer/tests/test_harness.py` | 38 |
| `benchmarks/tests/test_benchmark.py` | 33 |
| `tests/shared/test_io.py` | 33 |

### Offline / network discipline — GOOD
The non-negotiable "default tests must be offline" rule (CLAUDE.md) is genuinely enforced:
- `examples/real_breast_cancer/tests/conftest.py:42-69` builds tiny synthetic AnnData; no downloads.
- `test_scripts.py:35-77` exercises the download path in **dry-run only**, deletes the env var (`monkeypatch.delenv(harness.REAL_DATA_ENV)`), monkeypatches `download_reference`/`download_spatial` to raise if called, and even blocks `socket.socket.connect` to prove offline (`test_download_dry_run_is_offline`).
- Skips are import-guards (`pytest.importorskip("matplotlib"/"plotly")`) or documented data-availability skips (`tests/shared/test_protocol.py:346,363`) — no failing tests masked as skipped.

### Coverage by required area
| Area | Current coverage | Missing tests | Priority |
|---|---|---|---|
| **Bulk real-data validation** | Indirect only. `test_harness.py` covers underlying `_harness` functions: `prepare_reference` (`:116-137`), pseudobulk generation/determinism/simplex (`:138-187`), metrics perfect-recovery & known-error (`:337-364`), gene orientation/overlap (`:381-429`). Script `03_run_bulk_validation.py` gets ONLY `test_script_imports` (`test_scripts.py:25`). | No offline end-to-end run of `03` (load synthetic ref + pseudobulk → run → assert the documented output files: `bulk_estimated_proportions.tsv`, `bulk_validation_metrics.tsv`, `bulk_per_celltype_metrics.tsv`, `bulk_qc.tsv`, `bulk_warnings.json`). | **HIGH** |
| **Spatial real-data validation** | `04_run_spatial_validation.py` tested ONLY for missing-input error path (`test_scripts.py:79-84`). Spatial model invariants well covered in `tests/spatial/test_model.py` (NB multiplicative update, convergence trace decreasing, λ/α recording, non-convergence warning) and `tests/spatial/test_pipeline.py` (28). | No offline e2e of script `04` producing documented outputs (`spatial_spot_proportions.tsv`, `morans_i.tsv`, `spatial_warnings.json`, `spatial_run_metadata.json`). | **HIGH** |
| **Hierarchical mode** | Strong: `tests/shared/test_hierarchical.py` (27) covers hierarchy build/inference/conflict, mass-preservation (`test_aggregate_predictions_by_family_preserves_mass`, `test_combine_and_unresolved_preserve_total_mass`), conditional sums-to-one, resolvability flagging of non-separable families, partial resolution, bulk+spatial hierarchical workflows. Benchmark side: `test_benchmark.py:222` `test_hierarchical_not_judged_only_by_fine_level`. | Script `09_run_hierarchical_deconvolution.py` has **zero** test references (not even import). | MEDIUM |
| **Solver auto-selection** | Good: `tests/test_solver.py:43` `test_auto_selects_and_records_reason`, `:54` `test_auto_not_chosen_by_in_sample_reconstruction_alone` (verifies CV-recon-minus-conditioning objective, not in-sample fit), all solvers run + simplex. `AutoSolver` in `src/tissueresolve/solver/auto.py`. | No test of tie-breaking / penalty weighting edge cases or `candidate_solvers()` composition; only `n_splits=2`. | LOW |
| **Gene masking CV** | `tests/test_solver.py:64` `test_gene_masking_reproducible`, `:73` `test_gene_masking_scores_reconstruction`. Backs `mask_fraction`/`n_splits`/`seed` in `auto.py`/`ensemble.py`. | No test that masking handles tiny gene panels / panel smaller than fold count, or degenerate splits; no boundary test on `mask_fraction`. | MEDIUM |
| **Benchmark framework** | Strong (33 in `test_benchmark.py`): metrics, counts/CPM/log detection, library-type detection, batch-confounding, internal bulk/spatial baselines run on toy, capability matrix, unified report links, family/unresolved-aware metrics, dry-run offline. | Largely smoke/contract tests; no assertion that benchmark *rankings* are scientifically correct on a controlled ground-truth case beyond `perfect`/`known-error`. | LOW |
| **Report generation** | Structural only. `tests/report/test_report.py` (16), `test_report_components.py` (17), `test_report_redesign.py` (7), `test_results_dir_report.py` (7), `tests/test_cli_report.py` (5). Asserts HTML substrings (`"<html"`, `"Methods"`, `"mRNA"`/`"cell fractions"` warning surfaced, λ recorded, non-convergence shown, "no ground truth" stated, figures have title/caption/source-link). | No visual regression (see below). No test that report renders with REAL pipeline outputs rather than hand-built fixtures. Scripts `07_generate_reports.py` untested. | MEDIUM |
| **External benchmark runners** | Good contract coverage (`test_external_benchmark.py`, 18): MuSiC/RCTD/CARD/CIBERSORTx marked skipped/imported/exported when absent (not failed), composite excludes non-executed, imported-result CLI, Bisque R runner records "skipped" without package, `test_at_least_four_non_tissueresolve_executed`. | All external tools are absent in CI, so the actual import/parse of REAL external outputs is exercised only via tiny synthetic TSVs, never a real tool's output format. | LOW |
| **Documentation examples** | Existence/keyword only: `test_benchmark.py:189` `test_docs_pages_exist`, README keyword greps (`:183`, `:318`, `:393`), `test_compliance.py` reads CLAUDE.md/DESIGN_SPEC. `docs/` has `quickstart.md`, `tutorial.md`, `gene_masking_cv.md`, etc. | **No doctest, no executed examples.** Code snippets in `docs/quickstart.md`, `docs/tutorial.md`, `chimera_v1/TUTORIAL.md`, `examples/real_breast_cancer/README.md` can drift from the real API undetected. | MEDIUM |

### Targeted assessment questions

**Are tests too synthetic?** Yes, by design and largely defensibly. `tests/conftest.py` uses 40 genes / 4 cell types / Dirichlet proportions / exponential φ; the harness fixture uses 240 cells, 40 block-structured genes. These are excellent for *invariants* (simplex, mass conservation, non-negativity, reproducibility, error paths) but say little about behavior on realistic gene-panel sizes, sparsity, donor variability, or batch structure. The CLAUDE.md rule expects real-data fidelity to live in the opt-in harness — but that harness's user-facing scripts (03/04/07/09) are themselves under-tested (above), so the "real workflow" layer is the weakest link.

**Do tests reflect real workflows?** Partially. Algorithm/IO/QC/results layers are well covered. The end-user workflows — the numbered validation scripts and the report-generation/hierarchical scripts — are mostly verified via their *helper functions* (`_harness`) rather than the script `main()` entry points. `test_script_imports` (`test_scripts.py:25`) only checks `hasattr(mod, "main")`. Scripts `06`/`08` do have e2e tests (`test_resolution_script.py`, `test_resolution_analysis_script.py`), proving the pattern is feasible and just not applied to 03/07/09.

**Are reports tested structurally but not visually?** Confirmed. Every report test asserts on HTML string content or `path.exists()`. There is **no** `mpl_image_compare`, baseline image, PNG hash, or pixel comparison anywhere (grep for `imread|pixel|baseline_image|mpl_image_compare|compare_images|perceptual` returned nothing). Figure tests (`tests/plotting/test_plots.py:36-39`) assert the PDF/SVG/PNG and the accompanying data TSV *exist* — satisfying the "save underlying data" rule — but never validate figure correctness. Plotly publication figures (`test_publication.py`, `test_summary_and_histology.py`) are `importorskip`-guarded, so they silently don't run without plotly/matplotlib installed.

**Is CLI / examples adequately tested?** CLI: reasonably — `tests/cli/` (autodetect 3, cli 10, presets 1, protocol 1, wizard 1), `tests/spatial/test_cli.py` (8), `tests/test_cli_report.py` (5). `presets`/`protocol`/`wizard` are thin at 1 test each. Examples: the *harness library* is well tested, but the numbered example scripts are not consistently e2e-tested (03/07/09 gaps).

### What is NOT covered
- Offline e2e of `03_run_bulk_validation.py`, `04_run_spatial_validation.py` producing their documented output files. **(HIGH)**
- Any test touching `07_generate_reports.py`, `09_run_hierarchical_deconvolution.py`. **(HIGH)**
- Visual/figure regression for any plot or report. **(MEDIUM)**
- Executable documentation / doctests. **(MEDIUM)**
- Real external-tool output parsing (only synthetic TSV stand-ins). **(LOW)**
- Realistic-scale / sparse / donor-structured inputs in unit tests (relegated to opt-in harness). **(LOW/QUESTION)**

### Caveat on confidence
I did not execute the suite (read-only audit), so the 882 figure is a static `def test_` count, not a passed-test count; parametrized tests (e.g. `test_script_imports` over 6 scripts, `test_detect_cell_type_col_auto_tokens`) expand at runtime, so the *executed* test count is higher than 882. Coverage percentages are not asserted because no coverage run was performed; `[tool.coverage.run]` targets `src/tissueresolve` but the example scripts and benchmarks are outside that source path and thus excluded from any coverage metric even if run.

## PART 6 — Benchmark Audit

Scope: `benchmarks/run_real_external_benchmark.py`, `run_all.py`, `shared/{analysis,composite_score,metrics,imported}.py`, `shared/prepare_external_inputs.py`, `bulk/methods/*`, `spatial/methods/*`, and the committed outputs under `benchmarks/outputs/`. All findings are grounded in files I read.

### Executive summary of findings

| Severity | Finding |
|---|---|
| BLOCKER | Composite-score table labels spatial methods (CARD, cell2location) as `modality="bulk"` and ranks them in the bulk leaderboard with `final_score=0.295` — spatial-vs-bulk methods are scored in one table (`benchmarks/outputs/composite_scores.tsv` lines 12-13). |
| BLOCKER | `TissueResolve_auto` and `NNLS_baseline` are literally the same estimator on this dataset (AutoSolver selected `nnls`), so the headline "best method" win is a self-tie, not an external win (`selected_solver.json`, `real_external_method_status.tsv`). |
| HIGH | Truth has 41 cell types, reference 32; `metrics.align()` silently intersects to 32, dropping 9 true cell types' mass from every accuracy number. Violates CLAUDE.md "never silently remove" spirit and inflates all Pearson values. |
| HIGH | Benchmark scale is tiny: 12 pseudobulk samples, reference subsampled to 60 cells/type (1915 cells total), 600 Visium spots. Not statistically convincing; no replicates / CIs / significance tests. |
| HIGH | cell2location result is degenerate (near-uniform, all 32 types ~0.03 mean, `near_zero_fraction=0.0`); it still appears as `executed` with a composite score equal to CARD's, with no degeneracy flag in the score table. |
| MEDIUM | Composite non-accuracy dimensions are hand-set priors favoring TissueResolve (`_PRIORS` in `composite_score.py`), not measured. The "best overall composite = TissueResolve_auto" conclusion is partly built on these priors. |
| MEDIUM | Two parallel, inconsistent output sets: `outputs/bulk/executive_summary.tsv` shows MuSiC/Bisque `skipped`, while `outputs/real_external_method_status.tsv` shows them `executed`. Reports can disagree about what ran. |
| MEDIUM | `runtime_resource_score` ≈ 0.0 for nearly every executed method (CorrelationMatcher at 0.007s dominates the min), making the runtime dimension contribute almost nothing and distorting comparisons. |
| LOW | SPOTlight is `skipped (not installed)`, not "blocked"; BayesPrism `failed` after ~2180s (timeout). Framing in the task as "blocked/timed out" is roughly right but the report should state the install-vs-runtime cause precisely. |
| QUESTION | `requires_external_install` tools (MuSiC, Bisque, cell2location) have Python/R wrappers that always `raise NotImplementedError` and only export inputs; the *actual* external runs came from separate R/py runner scripts. Two execution paths for the same tool is confusing and a maintenance hazard. |

---

### Answers to the 10 questions

**1. Is the benchmark scientifically convincing?**
Not yet. It is a careful, honest *engineering* harness but a weak *scientific* benchmark:
- n=12 pseudobulk samples, one Visium section, reference downsampled to ≤60 cells/type (`prepared_inputs/input_summary.json`: 1915 cells, 5000 genes, 32 types). No replicate datasets, no bootstrap CIs on the leaderboard, no significance testing of method differences.
- The single accuracy number per method is a flattened Pearson over (samples × 32 shared types). With the auto-solver collapsing to NNLS, the "winner" is an internal baseline tie.
- Real spatial has no ground truth (correctly acknowledged), so the spatial half contributes structure/concordance only.

**2. Is the harness honest about what ran?**
Mostly yes — the status taxonomy (`executed`/`imported`/`exported_only`/`skipped`/`failed`) in `analysis.build_status_table` and `imported.py` is a genuine strength and the report repeatedly states exported-only tools are NOT benchmarked. But there are two divergent committed status tables (see MEDIUM above) and the NotImplementedError wrappers vs R-runner duality undercuts the clarity.

**3. Are the external runs trustworthy?**
Partially. MuSiC/Bisque R runners (`run_music.R`, `run_bisque.R`) use the same harmonised inputs as TissueResolve — good. But cell2location ran CPU/fast and produced a degenerate near-uniform output (verified directly); CARD ran (luminal-dominant, plausible). BayesPrism failed at the 1800s subprocess timeout (`_invoke` timeout=1800; status shows runtime_seconds=2180 from an earlier run). These caveats are partly in metadata but not surfaced in the composite ranking.

**4. Is input preparation fair across methods?**
Largely fair (shared signature/counts, same gene set, recorded subsampling in `prepare_external_inputs.py`). The 41→32 cell-type intersection (HIGH) is the main unfairness: methods are only ever scored on the 32 types the reference contains, and the 9 truth-only types are dropped silently for everyone.

**5. Are metrics fair to hierarchical methods?**
This is the best part of the design. `metrics.hierarchical_fair_metrics`, `fine_metrics_resolvable_only`, `family_level_metrics`, and `unresolved_aware_metrics` correctly avoid penalizing honest abstention, reward abstaining on `high_risk_families`, and report family-level + resolvable-fine views. `analysis.hierarchical_rankings` gives a multi-criterion view. Caveats:
- `unresolved_precision`/`recall` reward abstaining on families the *same reference* flagged as non-separable (`_high_risk_families` from `evaluate_within_family_resolvability`) — somewhat circular (the method is graded against a diagnostic computed from its own reference).
- The composite score still injects fine-level `accuracy` as the single accuracy term (weight 0.35), so hierarchical mode's intentionally lower fine Pearson (0.399 vs 0.768) drags its composite down despite the "fair" tables existing elsewhere. The fairness lives in the per-modality report, not in the composite ranking.

**6. Is it fair/meaningful that TissueResolve_auto selects NNLS then ties NNLS — how to frame it?**
This is the central credibility issue. `selected_solver.json`: AutoSolver picked `nnls` by gene-masking CV (masked-gene Pearson=0.986, condition≈234.6). Consequently `TissueResolve_auto` accuracy = `NNLS_baseline` accuracy = 0.7682150016166125 *to 16 digits* (`real_external_method_status.tsv`). So "TissueResolve_auto is the best method" reduces to "NNLS is the best method, and our auto-selector correctly chose NNLS."
- Fair framing: present it as *model selection working correctly* ("the auto-solver recovered the best backbone for this reference, which happened to be NNLS"), NOT as TissueResolve beating NNLS on accuracy. The current `benchmark_metadata.json` `best_bulk_fine_pearson: TissueResolve_auto` and the composite rank-1 are misleading without that caveat.
- The genuine TissueResolve value-add here is diagnostics (spillover, separability, abstention, protocol/normalization awareness, uncertainty), not a raw-accuracy win. The report's conclusion text says this, but the headline tables do not. Recommend: report the auto tie explicitly, and make the accuracy column show "= NNLS (selected)" rather than a duplicated number presented as an independent win.

**7. Is the report understandable / best-tool-by-scenario / no overclaiming / spatial-no-GT / enough plots / composite justified?**
- Understandable: yes, well-structured (`report.build_modality_report`, numbered sections, intro `note`/`warn` boxes).
- Best-tool-by-scenario: yes — `analysis.recommendation_by_use_case` and `best_method_summary` give per-criterion picks. But several picks are hard-coded to TissueResolve (`publication_report → TissueResolve_flat`, `broad_family → TissueResolve_hierarchical`) regardless of measured performance.
- Overclaiming: mostly avoided for spatial (explicit no-GT language in `report.py` and `analysis.conclusion_text`), but the composite-score rank-1 and `best_bulk_fine_pearson` overclaim for bulk given the NNLS tie.
- Spatial no-ground-truth: handled well and stated repeatedly.
- Enough plots: leaderboards, runtime, capability heatmap, concordance heatmap, near-zero bars — adequate for internal use, thin for publication (no per-cell-type scatter, no calibration plots, no CI error bars).
- Composite scores justified? Partially. Weights are config-driven and renormalized when GT absent (good, `composite_score.py`), and dimensions are written out (transparent). But 5 of 6 dimensions for non-TissueResolve methods come from a flat `_DEFAULT_PRIOR`, and TissueResolve gets favorable hard-coded priors. The composite is therefore not an objective ranking; it should be labeled as an opinionated scorecard, not a benchmark result.

**8-9. What benchmark + external methods + datasets are needed for publication?**
External methods to actually execute (not export):
- Bulk: CIBERSORTx, MuSiC, BisqueRNA, BayesPrism (fix timeout), DWLS, SCDC, EPIC, quanTIseq, Scaden, Bisque + DWLS. Currently only MuSiC + Bisque truly ran.
- Spatial: cell2location (GPU, full epochs — not the degenerate fast run), RCTD/spacexr, SPOTlight (install), stereoscope, Tangram, DestVI, CARD. Currently only CARD + a degenerate cell2location ran.
Datasets:
- Multiple independent references and tissues (not one breast-cancer reference), several Visium sections, plus a held-out reference (cross-dataset, to test reference-mismatch robustness).
- Bulk: many pseudobulk replicates per difficulty tier with full-depth references, AND a true bulk RNA-seq cohort with orthogonal proportions (flow/IHC) for mRNA-content discussion.
- Spatial: at least one platform with quasi-ground-truth (e.g. Visium + matched snRNA, or simulated spatial with known proportions; ideally a higher-resolution platform like Xenium/MERFISH aggregated to spots for partial ground truth).
Statistics: bootstrap/permutation CIs on every metric, per-cell-type metrics, calibration, paired significance tests between methods, runtime/memory on standardized hardware.

**Minimum benchmark by release tier:**
- **Internal release:** current harness is acceptable *if* (a) the composite bulk/spatial modality bug is fixed, (b) the auto=NNLS tie is labeled, (c) the 41→32 type drop is reported. It demonstrates the pipeline runs end-to-end against ≥2 real external tools.
- **Preprint:** ≥3-4 *executed* external bulk tools + ≥3 *executed* spatial tools, ≥2 datasets per modality, bootstrap CIs and per-cell-type metrics, full-depth (not 60-cell) reference, full-epoch cell2location, BayesPrism completing, and the headline reframed around diagnostics + correct model selection rather than a self-tie.
- **Nature-level:** multi-dataset, multi-tissue, multi-platform, cross-reference robustness, orthogonal ground truth (flow/IHC for bulk; matched high-plex imaging for spatial), all major competitors executed at full settings on common hardware, statistical significance with multiple-testing control, ablations of every TissueResolve component (protocol-awareness, weighting, abstention, spatial CAR, uncertainty), and demonstrated added value on a real biological question — not just correlation tables.

**10. Evaluate the composite score specifically.**
- Transparent and config-driven (strength).
- BLOCKER bug: spatial methods carried into the bulk modality scoring (`compute_composite_scores` is called once in `run_real_external_benchmark._write_report` with `modality="bulk"` default and `has_ground_truth=True` even though the status frame mixes spatial rows). Result: CARD/cell2location get bulk composite scores with `accuracy=NaN` → only priors/runtime → both exactly 0.295, ranked 11. This is meaningless and should never appear.
- Priors bias the ranking toward TissueResolve; should be removed or measured.
- Runtime dimension is effectively dead (min runtime ~0.007s dominates).
- `final_score` for `has_ground_truth=True` mixes a measured accuracy with five mostly-prior dimensions — not defensible as an objective leaderboard.

---

### Strengths worth preserving
- Honest executed/imported/exported/skipped/failed taxonomy (`analysis.build_status_table`, `imported.py`); never counts exported-only as benchmarked.
- Fair hierarchical metrics (`metrics.hierarchical_fair_metrics` and friends) — genuinely thoughtful handling of abstention.
- Harmonised inputs across tools (`prepare_external_inputs.py`), recorded subsampling, reproducibility metadata (`benchmark_metadata.json` with package versions).
- Spatial no-ground-truth framing is correct and repeated.
- Real R/py runner scripts (`run_music.R`, `run_bisque.R`, `run_cell2location.py`) use the same inputs and fail gracefully.

### Specific fixes before any release
1. Fix composite-score modality leakage (score bulk and spatial in separate calls; exclude accuracy for spatial). `composite_score.py` / `run_real_external_benchmark._write_report`.
2. Label the `TissueResolve_auto == NNLS` selection explicitly; do not present a duplicated accuracy as an independent win.
3. Surface the 41→32 cell-type intersection in every report and quantify dropped true mass (`metrics.align`).
4. Reconcile the two status tables (`outputs/bulk/executive_summary.tsv` vs `outputs/real_external_method_status.tsv`).
5. Flag cell2location's degenerate output; do not score a near-uniform prediction as a successful comparison.
6. Replace or measure the hand-set `_PRIORS`; relabel the composite as an opinionated scorecard, not a benchmark ranking.

## PART 7 — Report & UX Audit

Scope read in full: `src/tissueresolve/report/{unified.py, components.py, style.py, glossary.py, figures.py, html.py, sections.py, templates.py}`, `examples/real_breast_cancer/scripts/07_generate_reports.py`, plus the actually-rendered `examples/real_breast_cancer/outputs/report.html` (969 lines / 77 KB), `figure_manifest.tsv`, `bulk/bulk_warnings.json`, `validation_summary/warnings.json`, and the bulk/spatial sub-reports.

### Direct answers to the audit questions

| Question | Verdict |
|---|---|
| Would a first-time user understand what to do? | **Partially.** `index.html` + `README.md` + the executive `interpretation_guide` give a clear "start here". But the headline numbers section is generic and the figures are not visible inline. |
| Would a reviewer understand how results were generated? | **Partially.** Methodology boxes exist per section, but methods text is `<pre>` dumped, the figure manifest is incomplete (empty format paths), and the unified report hides the real QC warnings. |
| Are figures publication-ready? | **No, as presented.** In the unified report figures are NOT shown — 22/22 figure cards have empty bodies; only download links. Static PDF/SVG/PNG do exist on disk so the *assets* can be publication-grade, but the report itself shows no figure. |
| Are raw tables hidden appropriately? | **Yes.** `collapsible_table` / `_df_collapsible` keep full matrices closed; only top rows shown. This is well done. |
| One clear entry point? | **Yes.** `outputs/report.html`, reinforced by `index.html` and `README.md` ("start here"). Good. |
| Are bulk/spatial/hierarchical/benchmark reports coherently linked? | **Weakly.** They are linked (source-data links per section), but the sub-reports use a *different, older* visual style, and the benchmark link only surfaces a 2 KB summary file. |

### The OLD-style sub-report gap (explicitly requested)

Confirmed: `outputs/bulk/report.html` and `outputs/spatial/report.html` are generated via `sections.py` → `templates.py` (`CSS` with `#4C72B0`, `max-width:1080px`, simple `.toc`, no sidebar). `grep -c "report-container\|sidebar-nav"` returns **0** for both; `#4C72B0` appears in both. The unified report uses `style.py` (`--accent:#2c6e8f`, sticky dark sidebar). So a user who clicks "detailed bulk report" (`unified.py` link wiring, script lines 682–684 / 704–706) lands on a page that looks like a different product. Ironically the sub-reports are *more useful* (data-driven interpretation, summary cards, per-sample QC tables, severity-ranked warnings) but look older. This is a real coherence regression.

### Severity-tagged findings

- **BLOCKER — Warnings are silently dropped in the main report.** `07_generate_reports.py:768–786` reads `validation_summary/warnings.json`; on disk that file is `{"bulk": [], "spatial": []}`. Meanwhile `bulk/bulk_warnings.json` contains ~30 real `qc_recommendations` (high spillover for T-cell subsets, macrophages, endothelial; `mismatch_flag=medium`). The unified report therefore prints "No warnings recorded for this run." This directly violates CLAUDE.md rule #2 ("Never hide warnings") and rule #6 (separability warnings). Root cause: `main()` passes `bulk_warns`/`spatial_warns` (only figure-generation exceptions, almost always empty) into `_write_report_bundle_artifacts`, never the QC/spillover/separability warnings. The structured, data-driven warnings that `interpretation.collect_structured_warnings` produces (used by sub-reports via `sections._build_cards_and_warnings`) are not surfaced in the unified report at all.

- **HIGH — Figures are invisible in the unified report.** `_figure_cards` (script line 524–528) calls `C.figure_card(..., body_html="")`. Rendered HTML confirms 22 `class='fig-body'></div>` empty bodies, `<img`=0, `iframe`=1. The redesigned "figure cards" are effectively link lists. For a publication/clinical audience this is the single biggest UX miss — there is no visual at the entry point. Fix: embed the static PNG (`<img src=...>`) or an `<iframe>` to the interactive HTML in `body_html`.

- **HIGH — Generic captions for every figure.** All 22 unified figure cards use `_CAPTIONS["_default"]` ("Visual summary of a TissueResolve output…"). `figure_card`'s rich API (caption / legend / how_to_read / methodology / variables answering "what is plotted / axes / colors") is wired but fed boilerplate. A `bulk_main_summary_figure.caption.txt` exists on disk (written by the plotting layer) but is never read by `_figure_cards` (only `sections._fig` reads `.caption.txt`). So a curated caption is discarded in the main report.

- **HIGH — Unified report's interpretation is non-data-driven.** Every `C.interpretation_guide(...)` and `C.methodology_summary(...)` in `07_generate_reports.py` is a hardcoded constant string; none reflect the actual numbers (no "T cells dominate", no overlap level, no convergence status). The sub-reports DO this correctly via `report/interpretation.py` (`I.interpret_bulk_predictions`, `I.executive_summary_paragraph`, `I.generate_key_findings`). The redesign moved users to a report with *weaker* interpretation than the legacy one.

- **MEDIUM — Figure manifest is incomplete / not fully reproducible.** `figures/figure_manifest.tsv`: all 22 rows have empty `png_path`, `svg_path`, `pdf_path`, `caption_path`, `variables_defined`. `_figure_cards` (script 529–532) only populates `html_path`, `source_data_path`, `methodology`, `status`, even though PNG/SVG/PDF exist next to each figure and are linked in the HTML. The manifest's stated purpose ("every figure traceable to its source data and generation parameters", `figures.py` docstring) is only half met.

- **MEDIUM — Redundant double numbering in navigation.** `navigation_sidebar` (components.py:158) prepends `01`, `02`… while section titles already begin "1.", "2."… (script passes `"1. Executive summary"` etc.). Sidebar reads "01 1. Executive summary". Pick one numbering scheme.

- **MEDIUM — Not responsive / not print-friendly.** `style.py` has **0** `@media` queries. `.report-container{display:flex;max-width:1180px}` + `.sidebar-nav{width:248px;position:sticky;height:100vh}` will overflow/clip on narrow screens and waste a 248 px dark column on print/PDF export. No print stylesheet despite a publication mandate.

- **LOW — No active-section highlight / scrollspy.** No JS and no `:target` styling; the sidebar never indicates the current section in a long single-page report. `scroll-margin-top:16px` is small given content density.

- **LOW — Over-escaping limits caption formatting.** `figure_card` `esc()`s caption, legend, how_to_read, methodology, subtitle (components.py:119–131). Safe, but a caption can never use `<sub>`/italics/superscripts (e.g. "Moran's *I*", "log₂"), which a publication figure caption usually needs. Note the inconsistency: `estimate_note`, `interpretation_guide`, `method_box` pass raw HTML, so the escaping policy is uneven across the component set.

- **QUESTION — Benchmark linkage thin.** Section 8 inlines `best_method_summary.tsv` and links `benchmark_summary_report.html` (2 KB). Is that 2 KB file the intended full benchmark report, or a stub? If the latter, "coherently linked benchmark" overstates what a reviewer actually reaches.

- **QUESTION — `metric_grid` shows "—"/"generated" placeholders.** Executive cards include `"Report":"generated"` and several `—` values (samples/spots when context probing fails). Cosmetic, but a first-time user sees a dash where they expect a count.

### Top 10 report UX issues
1. **BLOCKER** Unified Warnings section says "No warnings" while 30+ real QC/spillover warnings exist (warnings hidden).
2. **HIGH** Figures not shown inline anywhere in the unified report (empty figure-card bodies).
3. **HIGH** Every figure uses the same generic `_default` caption; curated `.caption.txt` ignored.
4. **HIGH** Interpretation/methodology text is static boilerplate, not data-driven (sub-reports are better).
5. **MEDIUM** Sub-reports use the old visual style → jarring inconsistency from the main report.
6. **MEDIUM** Double section numbering in the sidebar ("01 1. …").
7. **MEDIUM** No responsive/print CSS; sidebar breaks layout on narrow screens and wastes space in PDF.
8. **LOW** No scrollspy/active-section indicator in a long single-page document.
9. **LOW** Methods rendered as a raw `<pre>` dump (script 790) — hard to read for non-bioinformaticians.
10. **LOW** Placeholder values ("generated", "—") in headline metric cards reduce trust.

### Top 10 figure improvements
1. Embed the static PNG (or interactive iframe) inside each figure card body; do not ship a "figure section" with no figures.
2. Read and use the per-figure `.caption.txt` (it exists on disk) instead of `_default`.
3. Add figure-specific axes/color/legend text (`figure_card` already supports `legend=` and `variables=`, both unused in the script).
4. Curate human titles instead of `stem.replace("_"," ").capitalize()` ("Bulk qc summary", "He abundance capillary endothelial cell" look auto-generated).
5. Populate the figure manifest's png/svg/pdf/caption columns (assets exist) for true traceability.
6. De-duplicate `spillover_network` (it appears in both bulk and spatial figure dirs → two identical cards).
7. Allow limited HTML/markdown in captions for units and gene/italic formatting (relax `esc()` for trusted caption text, or pass through a small sanitizer).
8. Add per-figure "what to check" tied to the data (e.g. flag the dominated/high-spillover types) rather than the generic how-to-read.
9. Provide a print/export-safe figure layout (current `max-width:100%` inside a flex+sidebar can shrink figures awkwardly on PDF export).
10. Surface separability/spillover *visually* in-report (heatmaps exist as files) rather than only as collapsed tables/links.

### Must-fix before release
1. **BLOCKER:** Wire real QC/separability/spillover warnings into the unified report (reuse `interpretation.collect_structured_warnings` and `T.severity_warning_box`, or populate `validation_summary/warnings.json` from the QC layer). "No warnings recorded" while spillover warnings exist is a scientific-integrity violation per CLAUDE.md.
2. **HIGH:** Render at least one image/iframe per figure card; an entry-point report with zero visible figures is not release-ready.
3. **HIGH:** Replace boilerplate interpretation with the data-driven `interpretation` module already used by the sub-reports, or the main report regresses on usefulness.
4. **MEDIUM:** Unify the sub-report styling with the new design system (or stop calling them "detailed" links that look like a different product).
5. **MEDIUM:** Complete the figure manifest (format/caption columns) so the "reproducible figure manifest" claim holds.
6. **MEDIUM:** Add responsive + print CSS before this is used to export publication PDFs.

### What is genuinely good (do not regress)
- Single clear entry point with `index.html` + `README.md` ("start here").
- Raw tables correctly hidden in `collapsible`/`_df_collapsible`; only top-N rows shown inline.
- Estimate-type framing is consistently correct and prominent (mRNA-derived proportions / spot-level composition, "no ground truth for real Visium") in `estimate_note`, `limitation_box`, glossary, and `methods_text.estimate_type_statement` — this honors the bulk/spatial scientific rules and avoids over-claiming.
- Glossary (`glossary.py`) is plain-language and explicitly cautious ("higher concordance ≠ correctness", "diagnostic, not accuracy"), and is actually wired into sections via `G.subset(...)`.
- `warning_box`/`limitation_box` return `""` on empty input by design (components.py:59–71) so they never print a false "No warnings" — good; the BLOCKER above is the *unified script* substituting an explicit "No warnings recorded" string, not the component.

## PART 8 — Documentation Audit

Scope: `README.md`, `docs/{tutorial,benchmarking,input_formats,output_interpretation,advanced_parameters,gene_masking_cv,project_structure_audit,accuracy_improvement_audit}.md`, `examples/real_breast_cancer/README.md`, `CLAUDE.md`. Cross-checked against `src/tissueresolve/cli.py`, `presets.py`, `api.py`, `bulk/solver.py`, `pyproject.toml`, `benchmarks/`, `examples/real_breast_cancer/scripts/`, and the on-disk `LICENSE`.

Overall the documentation is unusually detailed, honest about estimate types (mRNA proportions vs cell fractions), and careful not to overclaim accuracy. Most file paths, scripts, and APIs it references actually exist. The problems are concentrated in (1) a contradictory/incorrect license story, (2) several documented commands/flags that do not match the real CLI, and (3) a few stale references.

### Documentation issues table

| # | Severity | File / lines | Issue | Evidence |
|---|---|---|---|---|
| 1 | **BLOCKER** | `LICENSE`; `pyproject.toml:10`; `README.md:306-309` | Three-way license contradiction. On-disk `LICENSE` is **MIT** ("MIT License / Copyright (c) 2026"). `pyproject.toml` declares `license = { text = "BSD-3-Clause" }` and classifier "BSD Software License". README says "License metadata is declared as BSD-3-Clause" **and** "A top-level `LICENSE` file is not currently present" — but the file *does* exist and is MIT. A reader cannot determine the actual license. | `head LICENSE` → MIT; `pyproject.toml:10,22`; `README.md:306-309` |
| 2 | **HIGH** | `docs/tutorial.md:304-311` | Documents a working `tissueresolve bulk run --reference … --bulk bulk_counts.tsv --cell-type-col sub_cell_type --resolution-mode none --out …` command. The real `bulk run` (cli.py:568-572) is a stub: no options, prints "not yet implemented", `sys.exit(2)`. The flags `--bulk`, `--cell-type-col`, `--resolution-mode`, `--out` do not exist on it. README:170-171 *does* correctly warn `bulk run` is unimplemented, so the docs contradict each other. | `cli.py:568-572`; `README.md:170-171` |
| 3 | **HIGH** | `docs/output_interpretation.md:135,169`; `examples/real_breast_cancer/README.md:169` | Instruct users to "run with `--n-bootstrap > 0`". No such CLI flag exists. Bootstrap is only enabled by `--preset publication` or `--preset diagnostic` (presets.py:25-45 set `bootstrap: True`); `quick`/`standard` disable it. `n_bootstrap` is set from the preset, never from a CLI option (cli.py:227-228). | `cli.py` has no `--n-bootstrap`; `presets.py:5-46` |
| 4 | **HIGH** | `README.md:262`; `benchmarks/import_external_results.py:17,82` | The command `python benchmarks/run_all.py --use-existing-real-data --include-imported` will error. `run_all.py` only defines `--toy`, `--use-existing-real-data`, `--no-external` (run_all.py:28-30). `--include-imported` exists only on `benchmarks/bulk/run_bulk_benchmark.py:103` and `.../spatial/run_spatial_benchmark.py:97`. So the documented "fair comparison" workflow as written fails. | `run_all.py:28-30` vs `run_bulk_benchmark.py:103` |
| 5 | **MEDIUM** | `README.md:84-98` | Install story is internally inconsistent with the extras. README first shows `pip install -e ".[all]"`, then "for finer control" `".[spatial,report,realdata,benchmark]"`. But `[all]` = `spatial,report,plots,dev` only (pyproject.toml:74-76) — it excludes `realdata` and `benchmark`. So `[all]` does **not** give you the realdata/benchmark features the README and example then exercise. The example README:43-50 partially clarifies this for `realdata`, but the top-level README does not. | `pyproject.toml:74-76` |
| 6 | **MEDIUM** | `src/tissueresolve/cli.py:1-13` (module docstring) | The CLI module docstring lists only `bulk`, `spatial`, `info` commands and omits the primary `run` and `report` commands that all docs steer users to. README:162-169 lists the right commands, but the source-level help/docstring is stale and misleading. | `cli.py:1-13` vs registered commands `run`, `report`, `spatial run/report/benchmark/info`, `info` |
| 7 | **MEDIUM** | `README.md:162-169`; `docs/advanced_parameters.md` | Working commands omitted from "Supported CLI commands": `tissueresolve info` (cli.py:894) and `tissueresolve spatial info` (cli.py:814), and `tissueresolve spatial benchmark` (cli.py:770). These run; they are simply undocumented at the command-list level. | `cli.py:770,814,894` |
| 8 | **MEDIUM** | `docs/output_interpretation.md:84-89` | Says `suggest` is the **default** resolution mode ("`suggest` *(default)*"). The actual CLI default is `auto` (cli.py:40-46, 504-507), which resolves to `hierarchical` when broad/fine labels exist else `none`. README:13-17 correctly states `auto` is the default. So output_interpretation.md misstates the default. | `cli.py:42` `default="auto"` |
| 9 | **LOW** | `docs/tutorial.md:18` | Typo in install instructions: `dcd TissueResolve` (should be `cd TissueResolve`). Copy-paste of the block fails. | `tutorial.md:18` |
| 10 | **LOW** | `docs/tutorial.md:234-235` | Lists harness scripts in an order that runs `09_run_hierarchical_deconvolution.py` *before* `07_generate_reports.py`. Functionally fine, but inconsistent with the example README run order (07 then 09 optional) and the numbered convention; minor confusion. | `tutorial.md:225-236` vs `examples/.../README.md:53-101` |
| 11 | **LOW** | `docs/output_interpretation.md:52` ("Stage 6"), `:75` etc.; `examples/.../README.md` "(Stage 7/8/9)" | User-facing docs leak internal "Stage N" development labels. Harmless but unpolished for a public README/docs set. | `output_interpretation.md:51`; example README headers |
| 12 | **QUESTION** | `docs/input_formats.md:24`; `tutorial.md:46-47` | Bulk table orientation is stated as "genes in rows, samples in columns", but the CLI loader `_read_counts_table` (cli.py:374-383) only `pd.read_csv(index_col=0)` — it does not transpose or validate orientation. Whether downstream `deconv_bulk` expects genes-as-index is not confirmed in this audit; docs assert an orientation the loader does not enforce. Worth a one-line "the loader assumes genes are the row index" note and a validation check. | `cli.py:374-383` |
| 13 | **QUESTION** | `docs/gene_masking_cv.md:27-28`; `docs/accuracy_improvement_audit.md` | Specific accuracy numbers ("auto selects NNLS, masked-gene Pearson ≈ 0.99, fine-level Pearson 0.768", "hierarchical 0.749 vs flat 0.661 vs NNLS 0.844") are quoted with no pointer to a reproducible artifact (script, output TSV, seed). Not necessarily wrong, but unverifiable from the docs; benchmark outputs are git-ignored so a reader cannot reproduce them. Recommend linking the exact script + committed metadata, or labeling as illustrative. | `gene_masking_cv.md:27-28`; `accuracy_improvement_audit.md:8-37` |

### Verified-correct claims (no action needed)

- All `examples/real_breast_cancer/scripts/00–09`, `_harness.py`, `_download_utils.py` exist as documented.
- All benchmark entry points exist: `benchmarks/{run_all,import_external_results,run_real_external_benchmark}.py`, `bulk/run_bulk_benchmark.py`, `spatial/run_spatial_benchmark.py`, `envs/benchmark_installation.md`, `configs/composite_score_weights.yaml`.
- `MRNAContentCorrector` and the `deconv_bulk(..., mrna_corrector=...)` kwarg and separate `cell_fractions` field exist (solver.py:153, api.py:81, results.py:453-468) — the output_interpretation.md mRNA-correction snippet is accurate.
- Hierarchical preset gating values (`min_discriminating_genes` 5/10/30, spillover thresholds) in advanced_parameters.md match `presets.py:5-46`.
- `--solver` choices and `auto`/`pipeline` behavior in advanced_parameters.md / gene_masking_cv.md match `cli.py:59-64` and `api.py:75` (`solver` param exists).
- Hierarchical-mode requirement of broad/fine labels OR a mapping file, with a clear error otherwise, is real (`cli.py:_resolve_hierarchy_mapping`, 313-371) and well documented.

### Overclaiming statements to revise (quoted)

The docs are notably restrained about accuracy. The overclaims are mostly **capability/correctness** claims that the CLI does not back, plus license overstatement:

1. README:308-309 — *"License metadata is declared as BSD-3-Clause in `pyproject.toml`. A top-level `LICENSE` file is not currently present in this repository."* Both halves are problematic: a `LICENSE` file **is** present and it is **MIT**, contradicting the BSD-3-Clause metadata. (BLOCKER #1.)
2. tutorial.md:307-311 — the entire `tissueresolve bulk run … --bulk … --resolution-mode none --out …` example implies a working command; it is an unimplemented stub. (HIGH #2.)
3. output_interpretation.md:135 — *"if not computed, the report shows a message card ('run with `--n-bootstrap > 0`')"* and example README:169 *"Bootstrap uncertainty, if not computed, appears as a message card (run with `--n-bootstrap > 0`)."* The suggested remediation flag does not exist. (HIGH #3.)
4. README:259-262 — the *"run an external tool yourself, then import its predictions for a fair comparison"* block ending in `run_all.py … --include-imported` presents a complete workflow whose final command fails. (HIGH #4.)
5. gene_masking_cv.md:27-28 — *"auto selects full-gene **NNLS** (masked-gene Pearson ≈ 0.99) and reaches fine-level Pearson **0.768**, matching the strongest baseline"* — a concrete benchmarking claim with no reproducible, committed evidence (outputs are git-ignored). Either link the producing script + committed summary or mark as illustrative. (QUESTION #13.)

Conversely, the README "How TissueResolve differs" table (README:55-80) is appropriately hedged ("partial"/"varies", explicit "not a claim that TissueResolve outperforms these tools") and the benchmarking caveats (benchmarking.md:70-82) are honest about no-ground-truth spatial. These are good and should be kept.

### Missing sections

- **A correct, single source of truth for licensing.** Reconcile LICENSE/pyproject/README to one license (MIT or BSD-3-Clause), then state it once.
- **An accurate "extras matrix."** A short table mapping each extra (`spatial`, `report`, `plots`, `dev`, `realdata`, `benchmark`, `all`) to what it enables, making clear `[all]` excludes `realdata`/`benchmark`.
- **A complete CLI command reference** (currently scattered): all real commands — `run`, `report`, `info`, `spatial run`, `spatial report`, `spatial benchmark`, `spatial info` — and an explicit "not yet implemented" list (`bulk run`, `bulk check-compatibility`, `bulk benchmark`).
- **External-tool dependency clarity.** benchmarking.md mentions `benchmarks/envs/benchmark_installation.md` and "separate environments" but the top-level docs never state which external tools are R vs Python, that they are *not* pip-installed by any TissueResolve extra, and that CIBERSORTx/BayesPrism are export-only (the benchmarking.md does say this for CIBERSORTx, but it is buried).
- **Bulk input orientation contract** (genes-as-rows) stated as an enforced requirement, ideally with a validation note (#12).
- **Reproducibility pointer for the quoted accuracy numbers** (script + seed + committed metric file), since benchmark outputs are git-ignored.

### Recommended README structure

1. Title + one-line description
2. What it does / who it's for
3. **License** (single authoritative statement, matching the LICENSE file)
4. Installation
   - Base install
   - **Extras matrix** (table; note `[all]` ⊄ realdata/benchmark)
   - Python version + external-tool note (R tools, export-only tools)
5. Quickstart (bulk, spatial) — only commands that run today
6. **CLI command reference** (implemented vs not-yet-implemented, explicit)
7. Hierarchical broad→fine (required reference annotations; mapping-file alternative)
8. Input formats (incl. bulk orientation contract)
9. Outputs & "Where are my results?"
10. Interpretation (mRNA proportions vs cell fractions; spot composition; unresolved mass)
11. Benchmarking (internal baselines vs external tools; ground-truth caveats; correct commands)
12. Real-data validation harness (opt-in, offline-by-default)
13. Reproducibility / how to regenerate the quoted accuracy numbers
14. Status / contributing

### CLAUDE.md note

`CLAUDE.md` is a development-instruction file, not user docs; nothing to fix there. It is worth noting that several findings above are CLAUDE.md rule violations in the *documentation*, not the code: Rule "Never overstate" / honesty principles are breached by the license claim (#1), the non-existent `--n-bootstrap` remediation (#3), and the documented-but-unimplemented `bulk run` command (#2). Fixing the docs to match the real CLI brings them back in line.

Confirmed: license contradiction (MIT file vs BSD-3-Clause metadata) and substantial uncommitted report-system work. Both facts are load-bearing for my synthesis. I have enough verified ground to produce the deliverables.

# PART 9 — Publication-Readiness Audit

## 9.1 Readiness by journal tier

For each: novelty / benchmark / maturity / docs / use-case rated, plus likely reviewer criticisms.

| Journal tier | Readiness | Missing pieces | Realistic timeline |
|---|---|---|---|
| **Nature Methods** | **Not ready (≈15%)** | A demonstrated, statistically-significant accuracy or robustness advantage over the field's leading tools (RCTD, cell2location, CIBERSORTx, BayesPrism) on **multiple independent datasets**; orthogonal ground truth (flow/IHC for bulk, high-plex imaging for spatial); full ablations of every component; method novelty beyond "honest abstention." Currently the headline result is a **tie with plain NNLS by construction** (Part 2/3/6), competitors are mostly `skipped`/`exported`/`failed` (Part 6), and the flagship hierarchical mode trails on fine accuracy (0.062). | **18–30 months**, and only if a real, quantified advantage emerges. High risk it never reaches NM. |
| **Nature Biotechnology** | **Not ready (≈10%)** | NBT wants a tool with biological/translational impact and adoption potential. Needs everything NM needs **plus** a compelling real-world biological application (not a correlation table) and external user uptake. Out of scope for current evidence. | **24–36 months**; lowest fit of the Nature-family for a methods/diagnostics package. |
| **Nature Communications** | **Possible (≈30%)** | Multi-dataset validation (≥3 references, ≥2 tissues, ≥2 Visium sections), executed (not exported) external comparators (≥3 bulk, ≥3 spatial) with CIs and significance, validated abstention precision/recall, and reframing the contribution as **resolution-aware diagnostics + honest abstention** rather than accuracy superiority. Fix BLOCKERs first. | **12–18 months** with focused effort. The most realistic Nature-family target. |
| **Genome Biology** | **Possible (≈35%)** | Similar to Nat Commun but GB tolerates a strong methods+benchmarking story with rigorous evaluation even without a raw-accuracy win, **if** the diagnostics/abstention framework is validated and benchmarked honestly against several executed tools on several datasets. Composite-score and metric BLOCKERs must be fixed. | **10–16 months**. Good fit for the "robust, honest deconvolution + diagnostics" framing. |
| **Bioinformatics (Application Note / full paper)** | **Near-ready for App Note (≈55%); not ready for full paper** | App Note: stabilize the codebase (commit report system, remove `tuning` stub), fix the license, fix the benchmark BLOCKERs, ship working CLI/docs, and one honest benchmark table. Full Bioinformatics paper would still need multi-dataset benchmarking with significance. | App Note: **3–5 months**. Full paper: **9–14 months**. |
| **NAR Genomics & Bioinformatics** | **Approaching (≈45%)** | NAR GAB accepts solid, well-engineered tools with honest, reproducible benchmarking and clear utility. Needs the BLOCKERs fixed, ≥3 executed external tools per modality, ≥2 datasets, CIs, and the abstention contribution validated. More forgiving than Nat-family on requiring a raw-accuracy win. | **8–14 months**. Strong realistic target alongside GB. |
| **GigaScience / Scientific Data** | **Software/Technical-Note possible (≈50%)** | These value reproducibility, data/software resources, and engineering rigor over methodological novelty. A reproducible deconvolution+diagnostics framework with an opt-in real-data harness fits well. Needs codebase stabilization, license fix, executable docs, and a citable, reproducible benchmark artifact (outputs are currently git-ignored, Part 8 #13). GigaScience would want the benchmark to be genuinely reproducible end-to-end. | **5–9 months**. Strong fit for a software/resource paper if novelty is not overclaimed. |

### Per-tier dimension grid

| Tier | Novelty | Benchmark | Maturity | Docs | Use-case |
|---|---|---|---|---|---|
| Nature Methods | Insufficient | Far insufficient | Mid | Mid | Strong (concept) |
| Nature Biotech | Insufficient | Far insufficient | Mid | Mid | Weak fit |
| Nature Commun | Borderline | Insufficient | Mid | Mid | Strong |
| Genome Biology | Borderline+ | Insufficient | Mid | Mid | Strong |
| Bioinformatics | Adequate (App Note) | Borderline | Mid | Mid–good | Strong |
| NAR GAB | Adequate | Borderline | Mid | Mid–good | Strong |
| GigaScience/Sci Data | Adequate (resource) | Borderline | Mid | Mid | Strong |

## 9.2 Responses / fixes to the listed potential criticisms

| # | Criticism | Response / fix |
|---|---|---|
| 1 | **NNLS ties auto (the "win" is automatic)** | Honest reframe: this is **model selection working correctly** — the gene-masking CV selector recovered NNLS as the best backbone for a clean reference (margin 0.0002 over weighted_nnls, `solver/auto.py`). Do **not** present `TissueResolve_auto` accuracy (0.7682150016166125) as an independent win over `NNLS_baseline` (same 16 digits). Fix: in benchmark tables show `= NNLS (auto-selected)`; route the *real* contribution (protocol-aware `BulkPipeline`/`WNNLSSolver`) into a benchmarked method; and **demonstrate the robustness payoff** (the only justification for sacrificing clean-data accuracy) on a protocol-mismatch/noisy ground-truth experiment. Until that experiment exists, the weighting/panel choices are a net negative on every benchmark present (Part 2 BLOCKER, Part 3 Q3/Q4). |
| 2 | **External tools not truly benchmarked** | Only MuSiC + Bisque genuinely executed for bulk; CARD + a *degenerate* cell2location for spatial; CIBERSORTx `exported_not_run`, BayesPrism `failed` at timeout, SPOTlight not installed (Part 6). Fix: actually execute ≥3 bulk (CIBERSORTx, MuSiC, Bisque, DWLS/BayesPrism) and ≥3 spatial (RCTD, cell2location full-epoch, CARD, SPOTlight/stereoscope) at full settings on standard hardware, with recorded versions. Reconcile the two contradictory status tables (`outputs/bulk/executive_summary.tsv` vs `outputs/real_external_method_status.tsv`). Never count exported-only as benchmarked (the harness taxonomy already enforces this — keep it). |
| 3 | **Spatial has no ground truth** | Keep the (already correct) explicit no-ground-truth framing — do not claim spatial accuracy. For publication, add at least one quasi-ground-truth source: simulated spatial with known proportions, matched snRNA-seq, or aggregated Xenium/MERFISH-to-spot. Report concordance/entropy/Moran's I as descriptive only and state "concordance ≠ correctness" (glossary already does, Part 7). |
| 4 | **Reports not polished** | Fix the BLOCKER: the unified report prints "No warnings recorded" while 30+ real QC/spillover warnings exist in `bulk/bulk_warnings.json` (Part 7) — a CLAUDE.md rule-2 violation. Then embed figures inline (22/22 cards currently empty), use the curated `.caption.txt`, wire data-driven interpretation (the sub-reports already do this), unify sub-report styling, complete the figure manifest, add responsive/print CSS. Pick **one** report architecture (retire either `html/templates/sections` or `unified/components/...`, Part 1 Table 4). |
| 5 | **Too many partial features** | 80+ modules; the novel science is a small fraction (Part 2 §5). Fix: ship a focused v0.1 (core bulk + core spatial + diagnostics + one report) and demote everything experimental behind clear "experimental" flags or remove it (`tuning` stub, ensemble solver, second legacy hierarchy keyword table). See Part 10/12. |
| 6 | **Unclear core contribution** | State it in one sentence: *honest, gated abstention at the right cell-type resolution, backed by separability/spillover diagnostics and a single reference shared across bulk and spatial*. Drop "combines competitive solvers" (auto = single backbone, no blending; Part 2/3). The defensible claim is **diagnostics + abstention**, not accuracy superiority. |
| 7 | **Hierarchy / unresolved-mass validation** | Mass arithmetic is consistent by construction but **unguarded** (Part 3 Q6). Fix: add a per-row mass-conservation assertion + regression test. Validate the abstention *decision* (precision/recall of resolving genuinely separable families and abstaining on non-separable ones) — currently `high_risk_families=set()` makes unresolved precision/recall degenerate (Part 3 Q10, Part 6). Report subtype accuracy on resolvable families (the `fine_metrics_resolvable_only` view exists but is not foregrounded). |
| 8 | **Limited / no independent datasets** | One reference (12 bulk samples, 32/41 types, 1915 cells capped 60/type), one synthetic-truth spatial set, no CIs anywhere (Part 2/6 BLOCKER). Fix: ≥3 references, ≥2 tissues, ≥2 Visium sections, dozens of mixtures across difficulty tiers, bootstrap/permutation CIs and paired significance tests on every metric. Also surface and quantify the silent 41→32 cell-type intersection that drops 9 true types' mass (Part 6 HIGH). |
| 9 | **External install complexity** | R tools (MuSiC, Bisque, CARD, SPOTlight) and Python GPU tools (cell2location) are not pip-installed by any extra; CIBERSORTx/BayesPrism are export-only. Fix: a clear extras matrix (note `[all]` excludes `realdata`/`benchmark`, Part 8 #5), an external-tool dependency table (R vs Python, export-only), pinned env recipes (decide tracking of `benchmarks/envs/Makevars` + `_install_music_patched.R`), and graceful skip-with-instructions (already largely present — keep). |

---

# PART 10 — Product Definition

Classification: **A. Core for v0.1** / **B. Optional but useful** / **C. Experimental** / **D. Defer-or-remove**.

| Feature | Class | Rationale |
|---|---|---|
| **Core bulk deconvolution** (`BulkPipeline`/`WNNLSSolver`: marker selection, protocol-aware weighting, L1 `phi`, explicit mRNA-proportion typing, gene-overlap guard) | **A** | The real scientific bulk method; honors non-negotiable rules. Must be the benchmarked TissueResolve, not `auto`. |
| **Core spatial deconvolution** (SpatCAR NB multiplicative update + spatial smoothing, λ/α always recorded, non-convergence warns) | **A** | The spatial method; sound heuristic. Fix the overclaimed monotonicity docstring (Part 3 Q7). |
| **Separability / spillover diagnostics** (`reference/separability.py`) | **A** | The scientific safety net the whole package is built around; differentiator. |
| **Reference builder + `ReferenceSignature`** (unified bulk+spatial abstraction) | **A** | Foundational; the "single reference for both modalities" claim rests here. |
| **Hierarchical broad→fine + unresolved-mass abstention** (`reference/hierarchy.py`, `bulk/hierarchical.py`) | **A** (as **diagnostic/abstention**), **C** (as an *accuracy* mode) | The defensible core contribution as abstention; but fine accuracy trails badly, so do not ship it as an accuracy mode. Add mass-conservation guard + abstention validation. |
| **Bootstrap uncertainty (bulk)** | **B** | Sound CHIMERA idea; peripheral to core claim. Fix preset wiring so `bootstrap:False` actually disables it (Part 4 HIGH). |
| **mRNA-content correction** | **B** | Scientifically important and correctly gated/explicit; optional step. |
| **Reference suitability scoring** (`reference/suitability.py`) | **B** | Good conservative advisory; must be labeled advisory, not a probability (Part 3 Q1). |
| **Reports (unified HTML)** | **A** (one architecture, fixed) | Needed for usability/publication, but must fix the warnings BLOCKER, show figures, and de-duplicate the two report stacks. |
| **Plotting modules (12)** | **B** | Publication polish; co-saving source data is good. Audit for duplicate figure builders. |
| **Internal benchmarking baselines (NNLS etc.)** | **A** | Needed for honest self-comparison; keep but fix metric BLOCKERs (flattened Pearson, composite priors, modality leakage). |
| **External-tool benchmark harness** (import/export adapters, status taxonomy) | **B** | Honest engineering; required to *execute* (not export) competitors for publication. Fix composite modality-leakage BLOCKER. |
| **Solver `auto` selection + gene-masking CV** | **C** | Selects NNLS by 2e-4 margin, adds runtime, no accuracy gain; CV reconstruction proxy is unvalidated against proportion accuracy (Part 3 Q2). Keep as a diagnostic, do not headline. |
| **Ensemble solver** | **D** | Inherits the unvalidated CV proxy; no demonstrated benefit. Defer or remove. |
| **Ridge-NNLS / marker-NNLS backbones** | **C** | Reimplement known ideas without demonstrated edge; ridge is clip-projected (not true NNLS, misnamed). Keep experimental or remove. |
| **Expression reconstruction (R²/reconstruction metrics)** | **B** | Useful QC signal, but must not be presented as proportion accuracy (Part 3 Q2/Q10). Consolidate the duplicated R² helpers (Part 4). |
| **Spatial in-package synthetic benchmark (`spatial/benchmark.py`)** | **C** | Useful tooling; move out of the importable package and disambiguate the 3-way "benchmark" naming. |
| **H&E QC** (histology QC plots) | **C** | Present but not validated/foregrounded; experimental. |
| **H&E modeling** (using histology in the model) | **D** | Not a demonstrated, validated capability; defer until there is a real use-case + validation. |
| **Reference-free deconvolution** | **D** | Not implemented as a validated capability; out of scope for v0.1. Remove from any roadmap claims. |
| **`tuning/` module** | **D — remove now** | Orphaned stub that fabricates metrics (hardcoded `score 0.8`, `"Tuning ran (stub)"`); a scientific-honesty hazard if ever surfaced (Part 1). |
| **Heuristic keyword family inference (two competing tables)** | **C** | Fragile outside breast/immune; consolidate `_BROAD_KEYWORDS`/`_FAMILY_KEYWORDS`, prefer explicit user mappings (Part 2/3). |

---

# PART 11 — Release Plan

### Stage 0 — Stabilization
- **Tasks:** Commit/wire the uncommitted report system (or revert it) so HEAD builds with its only consumer (`07_generate_reports.py`); remove the `tuning` stub; fix the license (MIT-vs-BSD); fix preset→config wiring (`auto_tune`/`plots`/`spatial` ignored; `bootstrap:False` not honored, Part 4 HIGH); fix the report warnings BLOCKER (Part 7); decide tracking of `benchmarks/envs/{Makevars,_install_music_patched.R}`; reconcile docs with the real CLI (Part 8 #2–#8).
- **Success criteria:** `git status` clean; full test suite green; `tissueresolve run` + report produce a report that surfaces real warnings and shows figures; one authoritative license; no fabricated-metric code paths.
- **Effort:** 2–3 weeks.
- **Risks:** Choosing the wrong report architecture to keep; preset fix changing default behavior in tests.

### Stage 1 — v0.1 (focused, honest tool)
- **Tasks:** Lock Part-10 Class-A features as the product; demote Class-C behind explicit "experimental" flags; remove Class-D; disambiguate "solver"/"benchmark" naming; consolidate scripts 06/08; add mass-conservation assertion + abstention validation tests; add offline e2e tests for scripts 03/04/07/09 (Part 5 HIGH); write the extras matrix + complete CLI reference docs.
- **Success criteria:** A clearly-scoped package whose README claims match the CLI; no overclaiming language ("competitive solvers" removed); core paths fully tested incl. e2e; reproducible single-command demo.
- **Effort:** 4–6 weeks.
- **Risks:** Scope creep; underestimating the report consolidation.

### Stage 2 — Benchmark (credible evaluation)
- **Tasks:** Fix composite-score modality leakage + flattened-Pearson + raw-Pearson-as-accuracy + empty-high-risk-set (Part 3/6 BLOCKER/HIGH); label `auto==NNLS`; route real `BulkPipeline` into a benchmarked method; surface/quantify the 41→32 type drop; flag degenerate outputs (cell2location); **execute** ≥3 external bulk + ≥3 spatial tools at full settings; add bootstrap CIs, per-cell-type metrics, paired significance; run on ≥2 datasets per modality; full-depth (not 60-cell) reference; commit a reproducible benchmark artifact (or make outputs reproducible from a script + seed).
- **Success criteria:** A defensible leaderboard with CIs and significance where competitors are actually run; composite labeled as opinionated scorecard with accuracy-only shown alongside; the robustness-vs-clean-accuracy trade-off demonstrated on a mismatch experiment.
- **Effort:** 8–14 weeks (external-tool installs and GPU runs dominate).
- **Risks:** External tools fail/timeout (BayesPrism already does); GPU access for cell2location; the robustness advantage may fail to materialize, which would force a stronger reframe.

### Stage 3 — Preprint
- **Tasks:** Multi-dataset/multi-tissue validation (≥3 references, ≥2 tissues, ≥2 Visium sections); validated abstention precision/recall; subtype accuracy on resolvable families; honest figures with CIs; methods text; reframed contribution (diagnostics + abstention); publication-quality + reproducible reports/figures.
- **Success criteria:** A bioRxiv preprint whose central claim is fully supported by committed, reproducible evidence; no claim a reviewer can refute from the repo.
- **Effort:** 8–12 weeks after Stage 2.
- **Risks:** Reviewer-anticipated "no accuracy win" — mitigated only if the robustness/abstention case is quantified.

### Stage 4 — Nature-level
- **Tasks:** Multi-platform (Visium + Xenium/MERFISH-to-spot quasi-truth), cross-reference robustness, orthogonal ground truth (flow/IHC for bulk), all major competitors at full settings on common hardware with multiple-testing control, full ablations of every component, and a real biological application demonstrating added value.
- **Success criteria:** A quantified, statistically-significant, multi-dataset advantage (accuracy and/or robustness/abstention) plus biological impact.
- **Effort:** 9–18 months.
- **Risks:** High — may never reach a defensible raw-accuracy advantage; Nat Commun / GB / NAR GAB are the realistic ceilings on current trajectory.

---

# PART 12 — Actionable Prioritized Task List

Priority 1 = do first. Severity per audit. Effort: S(<1d), M(1–3d), L(1–2wk), XL(weeks). Pub impact = how much it gates publishability.

| # | Pri | Sev | Task | Why | Files likely affected | Effort | Pub impact |
|---|---|---|---|---|---|---|---|
| 1 | 1 | BLOCKER | Wire real QC/spillover/separability warnings into unified report (stop printing "No warnings") | CLAUDE.md rule-2 violation; integrity (Part 7) | `examples/.../07_generate_reports.py`, `report/unified.py`, `report/interpretation.py` | M | High |
| 2 | 1 | HIGH | Remove `tuning/` stub (fabricated metrics) | Scientific-honesty hazard (Part 1) | `src/tissueresolve/tuning/` | S | High |
| 3 | 1 | BLOCKER | Resolve license: make pyproject + README match the MIT LICENSE file | Reader cannot determine license (Part 8 #1) | `pyproject.toml`, `README.md` | S | High |
| 4 | 1 | HIGH | Commit/wire or revert the uncommitted report system so HEAD builds | HEAD doesn't import new modules; only consumer uncommitted (Part 1) | report/*, `07_generate_reports.py`, tests | M | High |
| 5 | 1 | BLOCKER | Fix composite-score modality leakage (spatial scored in bulk leaderboard) | Meaningless 0.295 scores (Part 6) | `benchmarks/shared/composite_score.py`, `run_real_external_benchmark.py` | M | High |
| 6 | 1 | HIGH | Label `TissueResolve_auto == NNLS`; stop presenting duplicated accuracy as a win | Central credibility issue (Part 2/3/6) | benchmark report/analysis modules, README, docs | M | High |
| 7 | 2 | HIGH | Embed figures inline in unified report (22/22 cards empty) | Entry-point report shows no figures (Part 7) | `07_generate_reports.py`, `report/components.py` | M | High |
| 8 | 2 | HIGH | Fix preset wiring: honor `bootstrap:False`, `auto_tune`, `plots`, `spatial` | Presets silently do nothing / run 200 bootstraps (Part 4) | `cli.py`, `presets.py`, `config.py` | M | Med |
| 9 | 2 | HIGH | Surface + quantify the silent 41→32 cell-type intersection | Drops 9 true types' mass, inflates Pearson (Part 6) | `benchmarks/shared/metrics.py`, reports | M | High |
| 10 | 2 | HIGH | Route real `BulkPipeline`/`WNNLSSolver` into a benchmarked TissueResolve method | `_auto` bypasses the actual differentiators (Part 3 Q4) | benchmark bulk methods, `run_real_external_benchmark.py` | M | High |
| 11 | 2 | HIGH | Fix doc/CLI mismatches: `bulk run` stub, `--n-bootstrap`, `--include-imported` | Documented commands fail (Part 8 #2–#4) | `docs/*.md`, `examples/.../README.md` | M | Med |
| 12 | 2 | HIGH | Replace boilerplate report interpretation with data-driven module | Regresses below legacy sub-reports (Part 7) | `07_generate_reports.py`, `report/interpretation.py` | M | Med |
| 13 | 2 | HIGH | Use curated `.caption.txt` + figure-specific captions in unified report | Generic `_default` caption for all 22 (Part 7) | `07_generate_reports.py`, `report/figures.py` | S | Med |
| 14 | 3 | HIGH | Add per-row hierarchical mass-conservation assertion + regression test | Consistent but unguarded (Part 3 Q6) | `reference/hierarchy.py`, `tests/shared/test_hierarchical.py` | M | High |
| 15 | 3 | HIGH | Validate abstention (set real `high_risk_families`, compute precision/recall) | Currently degenerate vs `set()` (Part 3 Q10, Part 6) | benchmark run scripts, `metrics.py` | M | High |
| 16 | 3 | HIGH | Demonstrate robustness payoff under protocol-mismatch/noise with ground truth | Sole justification for clean-accuracy loss is undemonstrated (Part 2 BLOCKER) | new benchmark experiment | XL | High |
| 17 | 3 | HIGH | Add offline e2e tests for scripts 03/04/07/09 | Only smoke-tested via helpers (Part 5) | `examples/.../tests/test_scripts.py` | M | Med |
| 18 | 3 | MEDIUM | Replace/measure composite `_PRIORS`; show accuracy-only alongside; relabel as scorecard | Priors predetermine TissueResolve rank (Part 3 Q9, Part 6) | `composite_score.py`, reports | M | High |
| 19 | 3 | MEDIUM | Replace flattened global Pearson headline with per-cell-type metrics | Inflated by between-type mean structure (Part 3 Q10) | `benchmarks/shared/metrics.py`, reports | M | High |
| 20 | 3 | MEDIUM | Execute ≥3 external bulk + ≥3 spatial tools at full settings | Most competitors skipped/exported/failed (Part 6) | `benchmarks/envs/*`, runners | XL | High |
| 21 | 4 | BLOCKER (validation) | Add ≥2 datasets/modality, full-depth ref, bootstrap CIs, significance | n=12, one dataset, no CIs (Part 2/6) | benchmark harness, new data scripts | XL | High |
| 22 | 4 | MEDIUM | Add spatial quasi-ground-truth (simulation or imaging-to-spot) | Cannot claim spatial accuracy otherwise (Part 6) | spatial benchmark | XL | High |
| 23 | 4 | MEDIUM | Reconcile the two contradictory benchmark status tables | Reports disagree on what ran (Part 6) | `benchmarks/outputs/*`, analysis | S | Med |
| 24 | 4 | MEDIUM | Flag degenerate competitor outputs (cell2location near-uniform) | Scored as success despite degeneracy (Part 6) | `composite_score.py`, analysis | S | Med |
| 25 | 4 | MEDIUM | Fix overdispersion: estimate within cell type, not pooled | Biased NB φ for spatial (Part 3 Q1) | `reference/build.py` | M | Med |
| 26 | 4 | MEDIUM | Fix union-find transitive over-merge in separability | Chains A~B~C collapse distinct types (Part 3 Q5) | `reference/separability.py` | M | Med |
| 27 | 5 | MEDIUM | Disambiguate naming: `benchmark`→`spillover`; `bulk/solver.py` vs `solver/` | 3-way "benchmark", 2-way "solver" (Part 1) | rename modules + imports | M | Low |
| 28 | 5 | MEDIUM | Merge scripts 06 + 08 into one resolution-analysis step | Documented overlap unfixed (Part 1/2) | `examples/.../scripts/06,08` | M | Low |
| 29 | 5 | MEDIUM | Consolidate `_BROAD_KEYWORDS`/`_FAMILY_KEYWORDS`; prefer explicit mapping | Duplicated, fragile outside breast/immune (Part 2/3) | `reference/hierarchy.py` | M | Med |
| 30 | 5 | MEDIUM | Fix spatial NB monotonicity docstring overclaim | Guarantee doesn't survive mixing/mismatch (Part 3 Q7) | `spatial/model.py` | S | Low |
| 31 | 5 | MEDIUM | Add responsive + print CSS to report; fix double sidebar numbering | No `@media`; breaks PDF export (Part 7) | `report/style.py`, `components.py` | M | Low |
| 32 | 5 | MEDIUM | Unify sub-report styling with the new design system (or stop calling them "detailed") | Looks like a different product (Part 7) | `report/sections.py`, `templates.py` | M | Low |
| 33 | 6 | MEDIUM | Defer/remove Ensemble solver; mark ridge/marker/auto as experimental | Unvalidated CV proxy, no demonstrated edge (Part 3 Q2, Part 10) | `solver/ensemble.py`, `solver/auto.py`, docs | M | Med |
| 34 | 6 | MEDIUM | Remove reference-free / H&E-modeling from roadmap claims (defer) | Not validated capabilities (Part 10) | docs, README | S | Low |
| 35 | 6 | LOW | Fix mislabeled `condition_number` (reported on unweighted R) | Misleading diagnostic (Part 3 Q10) | `solver/weighted_nnls.py`, `ridge_nnls.py` | S | Low |
| 36 | 6 | LOW | Stop calling `api.py` functions "thin wrappers"; consolidate shared CLI/api routing | They route resolution modes; duplication (Part 1/4) | `api.py`, `cli.py` | M | Low |
| 37 | 6 | LOW | Make `set_random_state` stop mutating global `np.random.seed` | Hidden global side effect on every spatial fit (Part 4) | `utils.py`, `spatial/model.py` | S | Low |
| 38 | 6 | LOW | Complete figure manifest (png/svg/pdf/caption columns) | "Reproducible manifest" claim half-met (Part 7) | `report/figures.py`, `07_generate_reports.py` | S | Low |
| 39 | 6 | LOW | Add doctests / executable doc examples; fix `dcd` typo, Stage-N leakage | Docs can drift from API (Part 5/8) | `docs/*.md`, test config | M | Low |
| 40 | 6 | LOW | Document/enforce bulk input orientation (genes-as-rows) | Loader doesn't enforce documented contract (Part 8 #12) | `cli.py`, `io/bulk`, docs | S | Low |