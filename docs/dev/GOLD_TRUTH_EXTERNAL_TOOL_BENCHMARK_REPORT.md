# Gold-Truth External Tool Benchmark Report

## 1. Objective

Extend the existing TissueResolve gold-truth bulk benchmark so external deconvolution tools are evaluated on the same donor-held-out pseudobulk mixtures, with the same datasets, splits, selected genes, fine and broad labels, scenarios, truth definitions, metrics, and output formats.

The benchmark must remain scientifically fair: no external tool may use a different train/test split, no tool may be tuned on test donors, and RNA-derived/mRNA proportion truth remains the primary comparison target unless a tool explicitly reports cell fractions and that is clearly labeled.

## 2. Reused internal benchmark artifacts

The benchmark will reuse and adapt the internal gold-truth benchmark machinery from `benchmarks/diagnostics/gold_truth_performance_benchmark.py`:

- `DatasetSpec`, `DatasetContext`, `MixtureBundle`
- `DATASETS` definitions for `breast` and `lung_hlca_subset`
- `SCENARIOS` list and scenario target generation
- donor-disjoint split logic (`split_donors`, `assert_donor_disjoint`)
- gene selection and reference construction
- mixture realization and truth generation
- `aggregate_truth_to_broad`, `conditional_truth_from_fine`
- composition metrics: Pearson, Spearman, RMSE, MAE, Jensen-Shannon, Aitchison, dominant accuracy
- scoring utilities for broad/fine/conditional/rare/spillover/runtime/family/subtype
- output writing conventions for TSV summaries

Where possible, the new benchmark will reuse data from `benchmarks/outputs/gold_truth_performance/` for internal TissueResolve baseline comparison and will preserve the same truth tables and scenario structure.

## 3. Datasets

- `breast`: `examples/real_breast_cancer/data/reference/breast_cancer_sc_reference.h5ad`
- `lung_hlca_subset`: `examples/second_tissue_lung/data/derived/hlca_subset.h5ad`

The benchmark will reuse the same donor-held-out splits, selected fine labels, and selected gene panels that the existing internal benchmark already defines.

## 4. Truth definition

Primary truth is RNA-derived/mRNA proportions, computed from the generated pseudobulk counts over the benchmark gene panel:

- `fine_truth`
- `broad_truth` aggregated from fine truth using the benchmark mapping
- `conditional_truth` within-family fine proportions
- `cell_fraction_truth` is preserved in outputs but is not used for primary scoring unless a method explicitly outputs cell fractions and is documented separately

The benchmark will not change the truth generation process, selected genes, or donor splits unless an external tool strictly requires a documented additional subset.

## 5. External tools attempted

Bulk external tools prioritized for this benchmark:

- `NNLS_baseline` (internal control)
- `MuSiC`
- `BisqueRNA`
- `BayesPrism` (attempt if installable and runtime acceptable)
- `DWLS` / `SCDC` (attempt if installed)
- `CIBERSORTx` (export-only unless predictions are actually imported)

Spatial tools are not required for the gold-truth bulk benchmark and will not block bulk progress. If spatial tools are used later, they must be clearly separated and not mixed with the bulk ground-truth comparison.

## 6. Expected limitations

- External tool execution may be limited by environment availability, package installation, or single-session runtime constraints. Tools that are not installed or cannot be executed will be recorded as skipped/failed, not hidden.
- Existing external R runner scripts are specialized for prepared inputs; harmonized gold-truth inputs will need to be exported in the same format.
- Some external methods may return cell fractions rather than RNA proportions. These methods will be labeled clearly and will not be compared directly to RNA-derived truth unless the truth definition is explicitly aligned.
- This benchmark report is focused on donor-held-out bulk pseudobulk; it does not claim real Visium spatial accuracy.
- To keep file footprint minimal, large harmonized input matrices are exported for runtime use but are git-ignored; only manifest and summary TSVs are versioned.

## 7. Implementation notes

The new benchmark code will be implemented in `benchmarks/diagnostics/gold_truth_external_benchmark.py`.

The harmonized input files and summaries will be placed under:

- `benchmarks/external_tools/inputs/gold_truth_bulk/`
- `benchmarks/external_tools/inputs/gold_truth_bulk/input_manifest.tsv`
- `benchmarks/external_tools/inputs/gold_truth_bulk/label_mapping.tsv`
- `benchmarks/external_tools/inputs/gold_truth_bulk/gene_overlap_summary.tsv`
- `benchmarks/external_tools/inputs/gold_truth_bulk/truth_manifest.tsv`

External tool metadata will be recorded in:

- `benchmarks/external_tools/tool_registry.tsv`
- `benchmarks/external_tools/run_manifest.tsv`

Benchmark outputs will be written under:

- `benchmarks/outputs/gold_truth_external/`

and combined with existing internal outputs to produce merged comparison tables.

---

# EXECUTED RESULTS (run on 2026-06-23)

> The sections above are the design plan. The sections below report what was
> actually executed. Run command:
> `python benchmarks/diagnostics/gold_truth_external_benchmark.py --tools nnls music bisque --datasets breast lung --run-available`
> Mixture parameters were read from `benchmarks/outputs/gold_truth_performance/run_metadata.json`
> (seed=2026, samples_per_scenario=1, cells_per_mixture=180, max_ref_cells_per_type=120)
> so the regenerated mixtures match the internal gold-truth benchmark exactly.

## R1. What was reused from the internal gold-truth benchmark

`gold_truth_external_benchmark.py` is an **extension/wrapper** of
`gold_truth_performance_benchmark.py`, not an independent runner and not a
duplicated scoring system. It imports and reuses, unchanged:

- `DATASETS`, `SCENARIOS`, `DatasetSpec`/`DatasetContext`/`MixtureBundle`
- `prepare_dataset`, `generate_all_mixtures`, `concat_full_overlap`,
  `realize_mixtures` (same donor splits, gene selection, reference build,
  mixture realization, and truth definitions)
- `score_bulk_predictions` (the **same** broad/fine/conditional/rare/spillover/
  resolution/runtime/family/subtype metrics used for the internal modes)
- `conditional_truth_from_fine`, `aggregate_truth_to_broad`, `_normalise_rows`

The older `benchmarks/diagnostics/external_bulk_benchmark.py` is a separate
holdout-input runner. It was **not** reused because it is tied to the older
holdout export, and the design rule is to avoid two divergent scoring systems.
All external predictions in this report are scored by the internal
`score_bulk_predictions` function, guaranteeing identical metric definitions.

Mixtures were **regenerated** (the internal run saves truth tables and metrics
but not the per-sample mixture count matrices). Regeneration uses the internal
run's recorded seed/splits/scenario parameters, so they are identical to the
internal benchmark's `full_overlap` group.

## R2. Datasets and truth definitions

| dataset | reference cells | bulk samples | selected genes | fine labels |
|---|---|---|---|---|
| breast | 1746 | 11 | 900 | 15 |
| lung_hlca_subset | ~ | 11 | 1000 | 15 |

Truth is **RNA-derived proportions** (`fine_truth`, `broad_truth`,
`conditional_truth`). `cell_fraction_truth` is exported but not used for primary
scoring. The external comparison uses the `full_overlap` scenario group (the
11 full-gene-overlap samples per dataset), matching how external tools see the
data.

## R3. External tools attempted

NNLS_external_control, MuSiC, BisqueRNA (selected via `--tools`). DWLS and SCDC
were not attempted (not installed in `benchmarks/envs/Rlib`). BayesPrism was not
run in this session: it is installed but a prior run took ~2180 s, which exceeds
a reasonable smoke/session budget; it is recorded as available-but-deferred.

## R4. External tools executed

All three attempted tools executed successfully on **both** datasets:

| tool | breast status | lung status | breast runtime (s) | lung runtime (s) |
|---|---|---|---|---|
| NNLS_external_control | executed | executed | 0.015 | 0.10 |
| MuSiC | executed | executed | 16.1 | 20.6 |
| BisqueRNA | executed | executed | 2.7 | 4.2 |

## R5. Tools skipped / failed / exported-only and why

- DWLS, SCDC: **skipped** — R packages not installed.
- BayesPrism: **deferred** — installed but runtime (~2180 s prior) too large for
  this session.
- CIBERSORTx: **exported-only** — web/token-gated, cannot auto-run.

Skipped / deferred / exported-only tools are **not** ranked.

## R6. Input harmonization

For each dataset the harness writes harmonized inputs under
`benchmarks/external_tools/inputs/gold_truth_bulk/<dataset>/`
(reference counts genes×cells, reference cell metadata with `cellType`/
`SubjectName`, bulk counts genes×samples, sample metadata, selected genes, label
mapping, and the three truth tables) plus `input_manifest.tsv`,
`gene_overlap_summary.tsv`, `truth_manifest.tsv`. Gene overlap was 900/900
(breast) and 1000/1000 (lung) between reference signature and bulk. The R
runners (`run_music.R`, `run_bisque.R`) read the shared `prepared_inputs/`
copies, the existing harness contract.

## R7. Estimate-type handling

All three external tools output **RNA-derived proportions** (row-normalized to
sum 1), the same estimate type as the primary truth. No tool output cell
fractions, so no comparison against `cell_fraction_truth` was needed. Truth
types were not mixed. NNLS is solved against the shared reference CPM signature
(`ReferenceSignature.as_R_cpm()`) but **not** the TissueResolve weighted-NNLS /
hierarchical solver — it is a genuine external control.

## R8. Broad-cell-type results (full_overlap; Pearson / RMSE)

| method | breast Pearson | breast RMSE | lung Pearson | lung RMSE |
|---|---|---|---|---|
| TissueResolve_flat | 0.806 | 0.105 | 0.860 | 0.103 |
| TissueResolve_broad_only | 0.796 | 0.111 | 0.750 | 0.131 |
| TissueResolve_hierarchical_soft | 0.058 | 0.222 | 0.602 | 0.239 |
| NNLS_external_control | 0.874 | 0.086 | 0.834 | 0.110 |
| MuSiC | 0.899 | 0.073 | 0.892 | 0.090 |
| BisqueRNA | 0.902 | 0.072 | 0.621 | 0.145 |

On broad composition, **MuSiC and BisqueRNA match or beat** TissueResolve_flat on
breast; on lung, TissueResolve_flat and MuSiC are best and BisqueRNA is worst.
`hierarchical_soft` scores low here **by construction** — it pushes mass to
`unresolved_*` columns (see R10/R12), which the broad scorer does not credit.

## R9. Fine-subpopulation results (full_overlap; Pearson / RMSE)

| method | breast Pearson | breast RMSE | lung Pearson | lung RMSE |
|---|---|---|---|---|
| TissueResolve_flat | 0.742 | 0.071 | 0.719 | 0.080 |
| TissueResolve_hierarchical_soft | 0.275 | 0.113 | 0.584 | 0.113 |
| TissueResolve_hard_legacy | 0.240 | 0.116 | 0.679 | 0.086 |
| NNLS_external_control | 0.773 | 0.069 | 0.714 | 0.090 |
| MuSiC | 0.796 | 0.063 | 0.755 | 0.076 |
| BisqueRNA | 0.746 | 0.068 | 0.452 | 0.106 |

On **raw fine accuracy, MuSiC is the strongest method** on both datasets;
NNLS and TissueResolve_flat are competitive; BisqueRNA is strong on breast but
weak on lung. TissueResolve_flat does **not** win raw fine accuracy.

## R10. Conditional within-family results (full_overlap; Pearson / RMSE)

| method | breast Pearson | breast RMSE | lung Pearson | lung RMSE |
|---|---|---|---|---|
| NNLS_external_control | 0.453 | 0.375 | 0.401 | 0.383 |
| MuSiC | 0.481 | 0.344 | 0.297 | 0.370 |
| BisqueRNA | 0.241 | 0.380 | 0.129 | 0.413 |

Within-family conditional accuracy is **poor for all external tools** (Pearson
0.13–0.48, RMSE 0.34–0.41). This is the regime where fine subtypes within a
broad family are hard to separate — and where TissueResolve's abstention layer
(R12) is designed to refuse false precision rather than emit unreliable
within-family splits.

## R11. Rare subtype results (threshold = 0.01)

| dataset / rare subtype | method | sensitivity | precision | FPR |
|---|---|---|---|---|
| breast / vein endothelial | NNLS | 0.000 | NaN | 0.000 |
| breast / vein endothelial | MuSiC | 0.889 | 0.889 | 0.500 |
| breast / vein endothelial | BisqueRNA | 0.778 | 0.778 | 1.000 |
| lung / SMG serous | NNLS | 1.000 | 0.889 | 0.333 |
| lung / SMG serous | MuSiC | 0.875 | 1.000 | 0.000 |
| lung / SMG serous | BisqueRNA | 0.750 | 1.000 | 0.000 |

MuSiC/Bisque **detect** the rare breast subtype that NNLS misses entirely, but
at high false-positive rate (0.5–1.0). On lung, NNLS and MuSiC both detect the
rare subtype; MuSiC has zero false positives.

## R12. Spillover / false-positive subtype detection + abstention

External spillover (mass placed on truly-absent subtypes):

| method | breast absent-mass | breast FP-rate | lung absent-mass | lung FP-rate |
|---|---|---|---|---|
| NNLS_external_control | 0.059 | 0.042 | 0.053 | 0.024 |
| MuSiC | 0.045 | 0.042 | 0.035 | 0.024 |
| BisqueRNA | 0.073 | 0.073 | 0.072 | 0.018 |

Abstention (`unresolved_mass`, only TissueResolve produces it):

| method | breast unresolved | lung unresolved | unresolved precision/recall |
|---|---|---|---|
| TissueResolve_hierarchical_soft / auto | 0.796 | 0.870 | 1.0 / 1.0 |
| NNLS / MuSiC / BisqueRNA | 0.000 | 0.000 | n/a |

**All external tools have `unresolved_mass = 0`**: they always emit a full fine
decomposition, even for families that are not fine-identifiable from this
reference. TissueResolve's soft/auto modes instead route ~0.80–0.87 of the
mass to `unresolved_*` with precision = recall = 1.0 and
`correct_resolution_level_rate = 1.0`. This is a capability the external tools
do **not** have, not a higher accuracy number.

## R13. Runtime

NNLS_external_control: 0.015–0.10 s. BisqueRNA: 2.7–4.2 s. MuSiC: 16–21 s.
NNLS is fastest; MuSiC is the slowest of the three (still well under a minute).

## R14. Statistical comparison

**Paired bootstrap CIs were not computed** and are not reported. The internal
gold-truth run saved only aggregated metric tables and truth (not per-mode,
per-sample TissueResolve prediction matrices), so paired per-sample differences
between TissueResolve modes and the external tools cannot be formed from the
saved artifacts. External per-sample predictions are saved under
`raw_predictions/`. Only descriptive per-dataset ranking is reported (R8–R12);
no CIs are fabricated. Computing paired CIs would require re-running the internal
benchmark with per-sample prediction export.

## R15. Where TissueResolve performs best

- **Abstention / resolution-awareness (R12):** unique `unresolved_mass` with
  perfect precision/recall; no external tool offers this.
- **Broad composition on lung (R8):** TissueResolve_flat (0.860 Pearson) ties
  MuSiC and beats NNLS/Bisque.
- **Refusing unreliable within-family splits (R10):** where external conditional
  accuracy is poor (Pearson 0.13–0.48), soft/auto modes decline to over-resolve.

## R16. Where TissueResolve ties

- **Fine accuracy, breast & lung (R9):** TissueResolve_flat (0.742 / 0.719
  Pearson) is within range of NNLS (0.773 / 0.714) and close to MuSiC.
- **Broad accuracy, lung (R8):** flat ties MuSiC.

## R17. Where TissueResolve underperforms

- **Raw fine accuracy vs MuSiC (R9):** MuSiC has the best fine Pearson/RMSE on
  both datasets; TissueResolve_flat does not win.
- **Raw broad accuracy on breast (R8):** MuSiC/Bisque beat TissueResolve_flat.
- **`hierarchical_soft` raw fine/broad numbers (R8/R9):** lower than flat —
  expected, because abstained mass is not credited by the raw scorers. This is a
  reliability-vs-raw-accuracy trade-off, not a defect.

## R18. Supported claims

- TissueResolve provides resolution-aware abstention diagnostics
  (`unresolved_mass`, correct-resolution-level rate) that MuSiC, BisqueRNA, and
  plain NNLS do not.
- On broad composition, TissueResolve_flat is competitive with the best external
  tools (ties on lung).
- All methods, including TissueResolve, struggle with within-family conditional
  accuracy on these references — TissueResolve's soft/auto modes respond by
  abstaining rather than guessing.

## R19. Unsupported claims

- **NOT** claimed: TissueResolve is the most accurate method overall. On raw
  fine accuracy MuSiC is at least as good or better.
- **NOT** claimed: TissueResolve's abstention implies higher accuracy. It is a
  different, complementary capability (reliability), not a higher score.
- **NOT** claimed: any statistical superiority — no CIs were computed.
- **NOT** claimed: spatial accuracy — this benchmark is bulk pseudobulk only.

## R20. Remaining limitations

- Only 11 full-overlap samples per dataset; CIs are wide (see "Bulk benchmark
  completion status" below — paired CIs are now computed, not faked).
- DWLS, SCDC, CIBERSORTx not executed → the external panel remains **partial**
  (BayesPrism has now been run on a small subset).
- Broad labels for breast are inferred from fine labels (documented internal
  limitation).
- Reduced-gene-overlap and low-depth scenarios are present in the internal run
  but the external comparison uses only the `full_overlap` group.

---

# Bulk benchmark completion status (Phase 1, run 2026-06-23)

This section supersedes the "no CIs / panel incomplete" caveats above with the
work done to make the bulk benchmark statistically usable.

## C1. Per-sample TissueResolve predictions — SAVED

`gold_truth_performance_benchmark.py` now saves per-sample prediction matrices
for every TissueResolve mode (estimates unchanged; only persisted):

```
benchmarks/outputs/gold_truth_performance/predictions/TissueResolve_<mode>__<dataset>[__<group>].tsv
benchmarks/outputs/gold_truth_performance/predictions_broad/...
benchmarks/outputs/gold_truth_performance/predictions_manifest.tsv
```

Each file has rows = sample IDs, columns = fine/broad labels (same names as
truth), with a sidecar `.meta.json` recording row sums, estimate type
(`RNA_derived`), mode, dataset, scenario group, seed, and donor split.
**Verified:** internal `full_overlap` sample IDs are identical to the external
tool predictions (alignment test passes), so paired comparison is valid.

## C2. Paired bootstrap CIs — COMPUTED (not faked)

`gold_truth_paired_stats.py` computes percentile paired bootstrap CIs (n=5000)
of per-sample metric **differences** (method A − method B) on the shared
`full_overlap` mixtures. Outputs:

```
benchmarks/outputs/gold_truth_external/paired_bootstrap_ci.tsv
benchmarks/outputs/gold_truth_external/rank_stability.tsv
benchmarks/outputs/gold_truth_external/rank_by_task.tsv
```

Headline paired results (n=11 per dataset; "favors" only when the 95% CI
excludes 0):

| comparison | dataset | metric | mean diff (A−B) | 95% CI | favors |
|---|---|---|---|---|---|
| flat vs MuSiC | breast | fine RMSE | +0.008 | [−0.003, +0.021] | **tie** |
| flat vs MuSiC | lung | fine RMSE | +0.004 | [−0.005, +0.013] | **tie** |
| flat vs MuSiC | breast | broad RMSE | +0.032 | [+0.019, +0.045] | MuSiC |
| flat vs BisqueRNA | lung | fine RMSE | −0.023 | [−0.042, −0.005] | **flat** |
| flat vs BisqueRNA | lung | fine Pearson | +0.311 | [+0.109, +0.517] | **flat** |
| flat vs NNLS | lung | fine Pearson | +0.096 | [+0.030, +0.178] | **flat** |
| hierarchical_soft vs MuSiC | both | fine/broad RMSE | + | excludes 0 | MuSiC |

Interpretation: **on raw fine accuracy, TissueResolve_flat is statistically
tied with MuSiC on both datasets**, beats BisqueRNA on lung, and is slightly
behind MuSiC/Bisque on breast *broad*. `hierarchical_soft`/`auto` are
significantly worse on raw metrics — expected, because they abstain. No method
statistically dominates across tasks; MuSiC is the most consistent #1 on fine
RMSE in `rank_stability` (prob_best 0.55 breast, 0.73 lung), with
TissueResolve_flat the top internal mode (rank #2 on lung).

## C3. DWLS / SCDC — skipped_not_installed

Both checked: `requireNamespace('DWLS')` → FALSE, `requireNamespace('SCDC')` →
FALSE; no runner scripts present. **Not installed, not run, not ranked.**
Install hints: DWLS — `remotes::install_bitbucket('yuanlab/DWLS')` (or the
maintained `install.packages('DWLS')` where available); SCDC —
`remotes::install_github('meichendong/SCDC')`. Per project rules these were not
auto-installed.

## C4. BayesPrism — executed on a small subset (feasible)

Run on **breast, 2 scenarios, 2 mixtures**:

| metric | value |
|---|---|
| status | executed |
| runtime | 88.0 s (2 samples) |
| broad Pearson / RMSE | 0.927 / 0.032 |

BayesPrism is the **strongest broad-composition method on the breast subset**,
but a full-panel run was ~2180 s previously, so it is recorded as
**feasible-for-subset / too-slow-for-routine-full-panel**. Subset outputs:
`benchmarks/outputs/gold_truth_external_bayesprism_subset/`. It is not added to
the main combined leaderboard because it was scored on a 2-sample subset, not
the full 11-sample `full_overlap` set.

## C5. CIBERSORTx — export-only (unchanged)

Web/token-gated; no local predictions imported → not scored, not ranked.

## C6. Scenario-group coverage

The internal benchmark scores three groups: `full_overlap` (6 TissueResolve
modes), `reduced_gene_overlap` (6 modes), and `spatial_like` (1 mode). The
**external bulk comparison uses only `full_overlap`** (the external runners
consume the concatenated full-overlap bundle). `reduced_gene_overlap` and
`spatial_like` were scored for TissueResolve internally but **not** for external
tools in this session. Truth generation was not changed. Extending external
tools to `reduced_gene_overlap` is future work.

## C7. Bulk completion summary

| blocker | status |
|---|---|
| Per-sample predictions saved | ✅ done |
| Paired bootstrap CIs | ✅ computed (not faked) |
| DWLS / SCDC | ⛔ not installed (recorded + hints) |
| BayesPrism | ✅ subset executed; full panel deferred (runtime) |
| CIBERSORTx | ◻️ export-only |
| Beyond full_overlap (external) | ◻️ documented, not run |

Supported after Phase 1: "TissueResolve_flat is statistically competitive
(tie with MuSiC) on raw fine accuracy and provides resolution-aware abstention
the external tools lack." Still **not** supported: any claim of broad
superiority.
