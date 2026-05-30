# TissueResolve — Audit and Migration Plan

---

## 1. High-Level Summary of Each Legacy Package

### CHIMERA (chimera_v1)

CHIMERA is a protocol-aware bulk RNA-seq deconvolution package. Its core scientific contribution is recognising that bulk sequencing protocols (polyA, ribodepleted) and single-cell reference protocols (scRNA, snRNA, 10x 3′ capture) create systematic gene-level biases that distort proportion estimates if uncorrected.

The pipeline is:
1. Build a reference matrix `φ` (G × K, L1-normalised probability vectors) from sc/snRNA-seq, with optional donor-aware averaging and cross-donor CV computation.
2. Assess protocol compatibility and compute per-gene risk scores using curated gene lists (intronic-dominant, dissociation-stress, length-biased, hypervariable).
3. Select a protocol-safe gene panel using a composite log-additive score (specificity × stability × protocol-safety × bulk-concordance).
4. Optionally filter discordant genes using a paired bulk↔pseudobulk BH-corrected t-test (Tier 3).
5. Solve per-sample weighted NNLS.
6. Bootstrap gene-panel CIs.
7. Compute per-sample QC (R², profile correlation, mismatch flag) and per-cell-type QC (marker recall, spillover risk, condition number).
8. Emit HTML + TSV + JSON report.

A standalone `MRNAContentCorrector` exists (but is not in the public API) to convert RNA proportions into cell fractions by dividing by per-type mRNA content. The output is correctly labelled as RNA proportions throughout.

Modules: `config`, `protocol`, `reference`, `gene_selection`, `wnnls`, `bootstrap`, `qc`, `discordance`, `mrna_correction`, `report`, `benchmark`, `plots`, `paper_plots`, `cli`.

---

### SpatCAR (spatcar v1.1.0)

SpatCAR deconvolves cell-type composition from 10x Visium spots using a negative-binomial count model with a conditional autoregressive (CAR) spatial prior. The key idea is that each Visium spot is a spatial mixture, and neighbouring spots share biology, so proportions should vary smoothly across the hexagonal array while still fitting the NB count model per spot.

The pipeline is:
1. Load Visium data from SpaceRanger output or h5ad. Parse array coordinates.
2. Build a pseudobulk reference by streaming sc/snRNA-seq counts, estimating per-gene NB overdispersion using Welford's online algorithm.
3. Select marker genes per cell type by log2FC + Visium overlap, optionally augmented with pairwise discriminability genes.
4. Check pairwise cell-type separability (Bhattacharyya coefficient + Jeffreys divergence + Pearson r). Warn loudly on poorly separable pairs. Optionally merge.
5. Build a hexagonal spot graph (row-normalised CSR adjacency matrix from Visium array coordinates). Optionally weight edges by expression similarity.
6. Initialise proportions per spot via NNLS warm-start.
7. Fit via block coordinate descent: multiplicative NB update → spatial proximal mixing (α-weighted blend with neighbourhood mean) → periodic mismatch scale-factor update (closed-form ratio estimator).
8. Compute per-spot QC (entropy, NB log-likelihood, dominant type, spatial residual), Moran's I per cell type.
9. Compute neighbourhood co-occurrence statistics via spatial permutation test with BH-FDR correction.
10. Detect spatial niches via Laplacian-smoothed K-means.
11. Optionally compute parametric bootstrap CIs by resampling from the fitted NB model.
12. Save results into h5ad, emit HTML report.

Modules: `io`, `reference`, `gene_selection`, `spatial_graph`, `deconvolution`, `protocol`, `qc`, `separability`, `neighborhood_stats`, `utils`, `benchmark`, `plots`, `paper_plots`, `report`, `cli`.

---

## 2. Function-by-Function Migration Map

### CHIMERA functions → TissueResolve destination

| CHIMERA source | Function / class | Destination in TissueResolve | Action |
|---|---|---|---|
| `config.py` | `ChimeraConfig` | `config.py → TissueResolveConfig` (bulk sub-config) | Rewrite |
| `config.py` | `ReferenceConfig` | `config.py → ReferenceConfig` (shared) | Preserve & generalise |
| `config.py` | `GeneConfig` | `config.py → GeneConfig` (shared) | Preserve & generalise |
| `config.py` | `DiscordanceConfig` | `config.py → DiscordanceConfig` (bulk only) | Preserve |
| `config.py` | `SolverConfig` | `config.py → BulkSolverConfig` | Preserve |
| `config.py` | `BootstrapConfig` | `config.py → BootstrapConfig` (shared) | Preserve |
| `config.py` | `QCConfig` | `config.py → QCConfig` (shared base, modality-specific sub-configs) | Rewrite |
| `protocol.py` | `BulkProtocol`, `RefModality`, `RefCapture`, `RefCounting` | `protocol/metadata.py` | Preserve |
| `protocol.py` | `ProtocolMetadata` | `protocol/metadata.py` | Preserve |
| `protocol.py` | `ProtocolRiskReport` | `protocol/risk.py` | Preserve |
| `protocol.py` | `ProtocolRiskAssessor` | `protocol/risk.py` | Preserve |
| `reference.py` | `ReferenceBundle` | `results.py → ReferenceBundle` (shared base class) | Rewrite; merge with SpatCAR `ReferenceData` |
| `reference.py` | `ReferenceBuilder` | `reference/build.py` | Merge with SpatCAR builder |
| `gene_selection.py` | `GeneSelector.select` | `reference/markers.py` | Preserve; merge pairwise discriminability from SpatCAR |
| `gene_selection.py` | `GeneSelector.compute_gene_weights` | `reference/markers.py` | Preserve |
| `gene_selection.py` | `_pairwise_max_fc`, `_max_log2fc` | `reference/markers.py` (private) | Preserve |
| `gene_selection.py` | `_load_gene_set_from_tsv`, `_load_prefix_blacklist` | `reference/gene_filters.py` | Preserve; remove global singletons |
| `gene_selection.py` | `_get_protein_coding`, `_get_hypervariable`, `_get_blacklist_prefixes` | `reference/gene_filters.py` | Rewrite (eliminate module-level singletons) |
| `wnnls.py` | `DeconvResult` | `results.py → BulkDeconvResult` | Preserve |
| `wnnls.py` | `WNNLSSolver` | `bulk/solver.py` | Preserve |
| `bootstrap.py` | `BootstrapCI` | `uncertainty/bootstrap.py` | Preserve; create shared interface for spatial bootstrap |
| `qc.py` | `QCMetrics` | `results.py → BulkQCMetrics` | Preserve; share base class with spatial |
| `qc.py` | `QCReport` | `bulk/qc.py` | Preserve |
| `discordance.py` | `DiscordanceFilter` | `bulk/pipeline.py` | Preserve; promote to documented feature |
| `discordance.py` | `_bh_correction` | `utils.py` | Merge with SpatCAR `_bh_fdr` |
| `mrna_correction.py` | `MRNAContentCorrector` | `bulk/solver.py` (or `bulk/mrna.py`) | **Promote to public API** |
| `mrna_correction.py` | `from_reference`, `from_file`, `from_bundle` | same | Preserve |
| `report.py` | `ChimeraRunReport` | `report/html.py` | Rewrite; share template with spatial |
| `benchmark.py` | All scenarios | `bulk/benchmark.py` | Preserve |
| `cli.py` | `run`, `check-compatibility`, `benchmark` | `cli.py → tissueresolve bulk ...` | Rewrite under unified CLI |
| `plots.py` | All | `plotting/bulk_plots.py` | Preserve style; adapt to shared `style.py` |
| `paper_plots.py` | All | `plotting/benchmark_plots.py` | Preserve |

---

### SpatCAR functions → TissueResolve destination

| SpatCAR source | Function / class | Destination in TissueResolve | Action |
|---|---|---|---|
| `io.py` | `load_visium` | `io/spatial.py` | Preserve |
| `io.py` | `load_reference` | `io/reference.py` | Merge with CHIMERA's `build_from_h5ad` |
| `io.py` | `save_results` | `io/spatial.py` | Preserve |
| `io.py` | `check_gene_overlap` | `io/validation.py` | Preserve |
| `reference.py` | `ReferenceData` | `results.py → ReferenceData` (shared base) | Merge with CHIMERA `ReferenceBundle` |
| `reference.py` | `build_pseudobulk_reference` | `reference/build.py` | Merge with CHIMERA `ReferenceBuilder` |
| `reference.py` | `estimate_overdispersion` | `reference/build.py` | Preserve (spatial-specific) |
| `gene_selection.py` | `select_marker_genes` | `reference/markers.py` | Merge with CHIMERA `GeneSelector` |
| `gene_selection.py` | `validate_gene_overlap` | `reference/markers.py` | Preserve |
| `gene_selection.py` | `get_marker_col_indices` | `reference/markers.py` | Preserve |
| `gene_selection.py` | `compute_pairwise_discriminability` | `reference/markers.py` | Merge into CHIMERA's pairwise FC logic |
| `gene_selection.py` | `select_marker_genes_with_separability` | `reference/markers.py` | Preserve; make shared |
| `spatial_graph.py` | `SpatialGraph` | `spatial/graph.py` | Preserve |
| `spatial_graph.py` | `build_hex_graph` | `spatial/graph.py` | Preserve |
| `spatial_graph.py` | `build_hex_graph_from_arrays` | `spatial/graph.py` | Preserve |
| `spatial_graph.py` | `build_expression_weighted_graph` | `spatial/graph.py` | Preserve |
| `spatial_graph.py` | `compute_spatial_batches` | `spatial/graph.py` | Preserve |
| `deconvolution.py` | `SpatCARModel` | `spatial/model.py` | Preserve |
| `deconvolution.py` | `_nb_multiplicative_update` | `spatial/model.py` | Preserve |
| `deconvolution.py` | `_update_mismatch_inline` | `spatial/model.py` | Preserve |
| `deconvolution.py` | `_store_fit_state` | `spatial/model.py` | Preserve |
| `protocol.py` | `ProtocolMismatch` | `protocol/mismatch.py` | Preserve; clarify naming (spatial mismatch vs bulk risk) |
| `protocol.py` | `compute_discordance` | `protocol/mismatch.py` | Preserve |
| `protocol.py` | `update_mismatch_factors` | `protocol/mismatch.py` | Preserve |
| `qc.py` | `compute_spot_qc` | `spatial/qc.py` | Preserve |
| `qc.py` | `compute_morans_i` | `spatial/qc.py` | Preserve |
| `qc.py` | `flag_low_quality_spots` | `spatial/qc.py` | Preserve |
| `qc.py` | `compute_model_qc` | `spatial/qc.py` | Preserve |
| `qc.py` | `bootstrap_proportions` | `uncertainty/bootstrap.py` | Move; align interface with CHIMERA bootstrap |
| `qc.py` | `boundary_sharpness` | `spatial/benchmark.py` | Move to benchmark |
| `separability.py` | `PairSeparability`, `SeparabilityReport`, `SeparabilityWarning` | `reference/separability.py` | Preserve; expose for bulk too |
| `separability.py` | `compute_separability` | `reference/separability.py` | Preserve; make shared |
| `separability.py` | `merge_nonseparable_types` | `reference/separability.py` | **Bug fix required**, then preserve |
| `separability.py` | `separability_heatmap_data` | `plotting/qc_plots.py` | Move |
| `neighborhood_stats.py` | `NeighbourhoodStats` | `spatial/neighbourhood.py` | Preserve |
| `neighborhood_stats.py` | `compute_neighbourhood_stats` | `spatial/neighbourhood.py` | Preserve |
| `neighborhood_stats.py` | `detect_spatial_niches` | `spatial/neighbourhood.py` | Preserve; **add smoothing param to output** |
| `neighborhood_stats.py` | `distance_decay_cooccurrence` | `spatial/neighbourhood.py` | Preserve |
| `neighborhood_stats.py` | `_bh_fdr` | `utils.py` | Merge with CHIMERA `_bh_correction` |
| `utils.py` | `project_simplex`, `project_simplex_batch` | `utils.py` | Preserve |
| `utils.py` | `sparse_column_subset`, `dense_subset` | `utils.py` | Preserve |
| `utils.py` | `set_random_state` | `utils.py` | Preserve |
| `utils.py` | `MemoryTracker`, `timer` | `utils.py` | Preserve |
| `utils.py` | `safe_log`, `nb_loglik_batch` | `utils.py` | Preserve |
| `utils.py` | `setup_logging` | `utils.py` | Merge with CHIMERA's logging setup |
| `benchmark.py` | `BenchmarkDataset`, `simulate_visium` | `spatial/benchmark.py` | Preserve |
| `benchmark.py` | Scenario runners | `spatial/benchmark.py` | Preserve |
| `cli.py` | `run`, `report`, `benchmark`, `info` | `cli.py → tissueresolve spatial ...` | Rewrite under unified CLI |
| `plots.py` | All | `plotting/spatial_plots.py` | Preserve style; adapt to shared `style.py` |
| `paper_plots.py` | All | `plotting/benchmark_plots.py` | Merge with CHIMERA paper plots |
| `report.py` | All | `report/html.py` | Rewrite; share template with CHIMERA |

---

## 3. Module Preservation Classification

### Preserve as-is (into TissueResolve)

These modules have solid scientific implementations and clear separation of concerns. Copy logic directly, adjusting only imports and signatures as needed for the new package structure.

- `chimera/protocol.py` → `protocol/metadata.py` + `protocol/risk.py`
- `chimera/wnnls.py` → `bulk/solver.py`
- `chimera/bootstrap.py` → `uncertainty/bootstrap.py`
- `chimera/discordance.py` (logic) → `bulk/pipeline.py`
- `chimera/mrna_correction.py` → `bulk/solver.py` or `bulk/mrna.py`
- `spatcar/spatial_graph.py` → `spatial/graph.py`
- `spatcar/deconvolution.py` → `spatial/model.py`
- `spatcar/protocol.py` → `protocol/mismatch.py`
- `spatcar/qc.py` (spot metrics, Moran's I, flagging, model QC) → `spatial/qc.py`
- `spatcar/separability.py` (after bug fix) → `reference/separability.py`
- `spatcar/neighborhood_stats.py` (after fix) → `spatial/neighbourhood.py`
- `spatcar/utils.py` → `utils.py`
- `spatcar/io.py` → `io/spatial.py` + `io/reference.py`

---

### Rewrite

These modules are correct scientifically but need architectural changes for TissueResolve.

- **`config.py` (both packages)** → Rewrite as a single `TissueResolveConfig` with `bulk:` and `spatial:` sub-sections. Keep the YAML round-trip.
- **`chimera/reference.py` + `spatcar/reference.py`** → Rewrite into `reference/build.py` with a unified `ReferenceData` object. See Section 6 (P0-1) for design constraints.
- **`report.py` (both packages)** → Rewrite as `report/html.py` with a shared Jinja2 (or f-string) template base.
- **`cli.py` (both packages)** → Rewrite as a single `tissueresolve` CLI with `bulk` and `spatial` sub-command groups.
- **`chimera/gene_selection.py` + `spatcar/gene_selection.py`** → Rewrite into `reference/markers.py` that incorporates both CHIMERA's composite scoring and SpatCAR's pairwise discriminability. Remove hardcoded test gene prefixes from CHIMERA's protein-coding mask.

---

### Merge into shared TissueResolve modules

These functions exist in parallel in both packages and must be unified.

| Functionality | CHIMERA | SpatCAR | Merged destination |
|---|---|---|---|
| Reference construction | `ReferenceBuilder` | `build_pseudobulk_reference` | `reference/build.py` |
| Gene selection | `GeneSelector` | `select_marker_genes` | `reference/markers.py` |
| Separability diagnostics | `_spillover_risk` in `qc.py` | `separability.py` | `reference/separability.py` |
| BH-FDR correction | `_bh_correction` in `discordance.py` | `_bh_fdr` in `neighborhood_stats.py` | `utils.py → bh_fdr()` |
| HTML report | `ChimeraRunReport` | SpatCAR `report.py` | `report/html.py` |
| Benchmark metrics | CHIMERA metrics in `benchmark.py` | SpatCAR metrics in `benchmark.py` | `benchmark/metrics.py` |
| Logging setup | CHIMERA `cli.py` | `spatcar/utils.py setup_logging` | `utils.py` |
| Result containers | `DeconvResult` | `ProtocolMismatch`, model attributes | `results.py` |
| Paper/publication plots | `paper_plots.py` (both) | same | `plotting/benchmark_plots.py` |

---

## 4. Bulk-Specific Modules

Modules that contain logic specific to bulk RNA-seq and must NOT be shared with spatial:

| Module | Content | Destination |
|---|---|---|
| `bulk/solver.py` | `WNNLSSolver`, weighted NNLS, L1-normalisation, per-sample R² | bulk-only |
| `bulk/pipeline.py` | Full bulk pipeline including DiscordanceFilter (Tier 3), mRNA correction | bulk-only |
| `bulk/qc.py` | Per-sample R², profile correlation, mismatch flag; per-cell-type marker recall, condition number | bulk-only |
| `bulk/benchmark.py` | Scenario definitions for protocol-mismatch benchmarking | bulk-only |
| `protocol/risk.py` | `ProtocolRiskAssessor`, gene list loading, risk scoring for bulk-reference combinations | bulk-specific; gene lists can be shared with spatial |
| Data gene lists | `intronic_dominant`, `dissociation_stress`, `length_biased_3prime`, `safe_universal`, `hypervariable_inflammatory` | bulk gene selection primarily; see P1-6 |

The `protocol/metadata.py` enumerations (`BulkProtocol`, `RefModality`, `RefCapture`, `RefCounting`) need a `VisiumProtocol` or `SpatialProtocol` enum added for the spatial side, but the existing bulk enums remain valid.

---

## 5. Spatial-Specific Modules

Modules that contain logic specific to spatial transcriptomics and must NOT be shared with bulk:

| Module | Content | Destination |
|---|---|---|
| `spatial/graph.py` | `SpatialGraph`, hexagonal graph construction, expression-weighted edges, spatial batching | spatial-only |
| `spatial/model.py` | `SpatCARModel`, NB multiplicative updates, CAR proximal mixing, mismatch factor coordinate descent | spatial-only |
| `spatial/qc.py` | Spot entropy, NB log-likelihood, spatial residual, Moran's I, bootstrap CIs, boundary sharpness | spatial-only |
| `spatial/neighbourhood.py` | `NeighbourhoodStats`, permutation co-occurrence tests, niche detection, distance decay | spatial-only |
| `spatial/benchmark.py` | Synthetic Visium simulation, DWLS/Spatial-NNLS/NNLS comparison | spatial-only |
| `io/spatial.py` | `load_visium`, `save_results`, SpaceRanger parsing, tissue positions CSV | spatial-only |
| `protocol/mismatch.py` | `ProtocolMismatch`, per-gene scale factors `d_g`, discordance scoring, update step | spatial-only |
| `plotting/spatial_plots.py` | Spatial proportion maps, convergence traces, mismatch factor plots | spatial-only |

---

## 6. P0 — Critical Issues

Must be resolved before any migration commit is made.

**P0-1: Incompatible reference object formats.**

CHIMERA produces `ReferenceBundle` (L1-normalised φ, G × K, donor_cv). SpatCAR produces `ReferenceData` (CPM matrix K × G + log1p, phi_g dispersion per gene). These orientations differ (G × K vs K × G), the normalisation differs (L1 probability vs CPM), and SpatCAR carries NB dispersion which CHIMERA does not need.

The unified `TissueResolve` reference object must carry all of:
- L1-normalised `phi` (G × K) for bulk deconvolution
- `R_cpm` (K × G) and `R_log` for SpatCAR's NB model
- `phi_g` (G,) NB dispersion for spatial
- `donor_cv` (G × K) cross-donor CV for bulk gene selection
- `n_cells_per_type` for provenance

This requires a single `ReferenceData` class with a `.as_phi()` convenience method and a factory that produces both representations during construction. Do not store duplicate data — build on demand.

**P0-2: `MRNAContentCorrector` is absent from CHIMERA's public API.**

`chimera/mrna_correction.py` is a complete, well-documented implementation of the correction that converts RNA proportions to cell fractions, but it is not imported in `chimera/__init__.py` and is never called from the CLI or the pipeline module. The DESIGN_SPEC mandates that TissueResolve "always communicate what its estimates represent." This module must be promoted to the public API, exposed in the CLI (`--mrna-content` flag), and its results written to a dedicated output file.

**P0-3: `merge_nonseparable_types` in SpatCAR has a computation bug.**

In `spatcar/separability.py`, the dispersion merge step:
```python
phi_new = np.maximum(phi_new, ref_data.phi_g)  # global max
```
This computes the global maximum across ALL original cell types' dispersions, not just the types being merged into the current group. The result is that `phi_new` becomes `max(ref_data.phi_g)` for every new type regardless of group membership — every merged cell type will have the same dispersion equal to the most overdispersed gene in the entire reference. This must be corrected before migration and a regression test added.

**P0-4: `gene_selection_orig.py` exists alongside `gene_selection.py` in CHIMERA.**

There is a legacy file `chimera/gene_selection_orig.py` that has not been removed. Before migration, audit whether it is dead code (delete it) or still imported by tests (migrate those tests to the current `gene_selection.py`).

**P0-5: Global singletons for gene lists in CHIMERA's `gene_selection.py`.**

Module-level mutable singletons `_PROTEIN_CODING`, `_HYPERVARIABLE`, `_BL_PREFIXES` are populated lazily on first call and cached as module globals. This causes test contamination (a test that loads hg38 lists will contaminate a subsequent test expecting mm10) and makes the genome parameter effectively non-functional on second call. All gene list loading must be moved to instance state or factory functions.

---

## 7. P1 — Important Issues

Must be resolved before the first stable release.

**P1-1: SpatCAR's `compute_spot_qc` never receives model arrays in the `save_results` path.**

`spatcar/io.py:save_results` calls `compute_spot_qc(model.proportions_, graph)` without passing `Y_marker`, `R_d`, `phi_g`, or `lib_sizes`. This means `nb_loglik` is always `NaN` in the output `prop_qc`. The call in `save_results` must be updated to extract and pass the relevant arrays from the model, since the NB log-likelihood is one of the most informative per-spot QC metrics.

**P1-2: SpatCAR's bootstrap CI under-coverage is not surfaced to the user.**

`spatcar/qc.py:bootstrap_proportions` documents (in a code comment only) that empirical coverage is ~88–93% for 95% nominal with `n_iter_per_boot=30`. In TissueResolve, bootstrap CI output must include a `coverage_note` field or metadata entry documenting this limitation, consistent with the DESIGN_SPEC rule on transparent uncertainty.

**P1-3: CHIMERA's `from_bundle` in `mrna_correction.py` silently disables correction.**

`from_bundle()` warns that L1-normalised phi cannot be used for mRNA content estimation, then returns a `MRNAContentCorrector` with all weights = 1.0 (no correction). In TissueResolve's pipeline, if mRNA correction is requested but raw counts are unavailable, the pipeline must raise an explicit error or at minimum emit a prominent non-suppressible warning in the report.

**P1-4: CHIMERA benchmark scenarios may be partially unimplemented.**

The `benchmark.py` docstring describes five scenarios (`SANITY_PSEUDOBULK`, `TISSUE_MISMATCH_CONTEXT`, `PROTOCOL_MISMATCH_CENTRAL`, `EXTERNAL_REFERENCE`, `PAIRED_CALIBRATION`) but the `EXTERNAL_REFERENCE` scenario in particular requires real-world data adapters. Before migration, the actual implementation extent must be audited and stubs must raise `NotImplementedError` rather than silently returning empty results.

**P1-5: Two independent BH-FDR implementations with subtly different edge-case handling.**

CHIMERA's `_bh_correction` and SpatCAR's `_bh_fdr` are both correct in principle but differ in implementation. A single implementation in `utils.py` must be verified against a known-good reference implementation, documented with a test against scipy if available, and used throughout.

**P1-6: Protocol risk gene lists are CHIMERA-only; spatial deconvolution benefits from the same filters.**

The intronic-dominant and dissociation-stress gene lists in `chimera/data/gene_lists/` are equally relevant to Visium deconvolution when the reference is snRNA-seq. SpatCAR currently does not apply any protocol filter. In TissueResolve, the gene filter logic in `reference/gene_filters.py` must be available to both modalities' marker selection pipelines.

**P1-7: SpatCAR's `detect_spatial_niches` smooths proportions but does not record the parameter.**

The function applies `n_smooth` rounds of Laplacian smoothing to Pi before clustering, but the return value is only a label array — the smoothing parameter is not stored in the output. This violates the DESIGN_SPEC rule "Never smooth spatial estimates without reporting the smoothing parameter." The function must return or expose `n_smooth` in the output, and `save_results` must persist it.

---

## 8. P2 — Improvements

Quality enhancements for post-v1.0.

**P2-1: Unify pairwise log2FC computation.**

CHIMERA uses `_pairwise_max_fc` (max across all pairs) and `_max_log2fc` (vs. mean of others) on L1-normalised phi. SpatCAR uses one-vs-mean-others on log1p-CPM. These are not equivalent. For marker selection in TissueResolve, one approach should be chosen and justified; the other removed.

**P2-2: CHIMERA has hardcoded test gene prefixes in the protein-coding filter.**

In `chimera/gene_selection.py`, test data gene prefixes (`"M0"`, `"M1"`, `"SH_"`, `"g0"`, `"MAC"`, etc.) are hardcoded into the production protein-coding filter. The fix is to make the protein-coding filter configurable (allowlist or skip-flag) and handle test data via `GeneConfig` or test fixtures.

**P2-3: CHIMERA's `QCConfig` thresholds need richer provenance.**

Thresholds like `r2_warn=0.50` and `profile_corr_warn=0.60` are documented as "empirical" from "Finotello 2019 + BAL 2026" but the exact figures are not traceable without cross-referencing. TissueResolve should either reproduce the derivation in `docs/qc_thresholds.md` or clearly mark them as unvalidated defaults that must be treated as tunable hyperparameters rather than universal thresholds.

**P2-4: CHIMERA CLI lacks a `--genome` flag in some commands.**

The genome parameter (hg38 vs mm10) affects which gene lists are loaded in `ProtocolRiskAssessor` and `GeneSelector`. The unified CLI must expose `--genome [hg38|mm10|auto]` as a top-level option.

**P2-5: SpatCAR's `build_expression_weighted_graph` duplicates `_HEX_OFFSETS`.**

`spatcar/spatial_graph.py` defines `_HEX_OFFSETS` as a module-level constant, but `build_expression_weighted_graph` re-defines it locally. The local definition should be removed and the module constant used.

**P2-6: No `__all__` in several modules.**

`chimera/mrna_correction.py`, `chimera/discordance.py`, `spatcar/utils.py`, `spatcar/separability.py` do not declare `__all__`. All TissueResolve public modules must declare `__all__` to define a stable API surface.

**P2-7: CHIMERA `DeconvResult` column names do not encode the output type.**

The comment `# CHIMERA v1 output: mRNA proportions (not cell fractions)` in the TSV header is ignored by any TSV reader. The DESIGN_SPEC says to never report bulk estimates as absolute cell fractions. Column names should encode the estimate type, e.g. `{cell_type}_mRNA_prop`, or a machine-readable metadata field should be added.

---

## 9. Proposed Final Architecture for TissueResolve

```
TissueResolve/
├── src/
│   └── tissueresolve/
│       ├── __init__.py              # version + top-level re-exports
│       ├── config.py                # TissueResolveConfig, BulkConfig, SpatialConfig
│       │                              # ReferenceConfig, GeneConfig, BootstrapConfig,
│       │                              # QCConfig, DiscordanceConfig — all with YAML I/O
│       ├── cli.py                   # tissueresolve bulk / spatial / info
│       ├── results.py               # Shared result containers: ReferenceData, BulkDeconvResult,
│       │                              # SpatialDeconvResult, BulkQCMetrics, SpatialQCMetrics
│       ├── utils.py                 # bh_fdr, project_simplex, dense_subset, nb_loglik_batch,
│       │                              # safe_log, set_random_state, setup_logging, MemoryTracker
│       ├── io/
│       │   ├── bulk.py              # read_bulk_counts (TSV/CSV/h5ad)
│       │   ├── spatial.py           # load_visium, save_results (SpaceRanger + h5ad)
│       │   ├── reference.py         # load_reference_h5ad, load_reference_csv
│       │   └── validation.py        # check_gene_overlap, validate_array_coords
│       ├── reference/
│       │   ├── build.py             # unified ReferenceBuilder: donor-aware + streaming
│       │   │                          # builds ReferenceData with phi + R_cpm + phi_g + donor_cv
│       │   ├── markers.py           # GeneSelector (composite weights + pairwise discriminability)
│       │   ├── gene_filters.py      # gene list loading (intronic, stress, length, hypervariable,
│       │   │                          # blacklist, protein-coding) — instance-based, no singletons
│       │   └── separability.py      # compute_separability, merge_nonseparable_types,
│       │                              # separability_heatmap_data, SeparabilityReport
│       ├── protocol/
│       │   ├── metadata.py          # BulkProtocol, RefModality, RefCapture, ProtocolMetadata
│       │   ├── risk.py              # ProtocolRiskAssessor, ProtocolRiskReport — bulk-focused
│       │   └── mismatch.py          # ProtocolMismatch (spatial d_g factors), compute_discordance,
│       │                              # update_mismatch_factors
│       ├── bulk/
│       │   ├── solver.py            # WNNLSSolver, MRNAContentCorrector (now public)
│       │   ├── pipeline.py          # ChimeraPipeline (orchestrates all bulk steps + CLI entry)
│       │   │                          # integrates DiscordanceFilter
│       │   ├── qc.py                # BulkQCReport: R², profile_corr, marker_recall,
│       │   │                          # spillover_risk, condition_number
│       │   └── benchmark.py         # All 5 benchmark scenarios
│       ├── spatial/
│       │   ├── model.py             # SpatCARModel (NB-MAP + CAR proximal mixing)
│       │   ├── graph.py             # SpatialGraph, build_hex_graph, expression-weighted graph
│       │   ├── pipeline.py          # SpatialPipeline (orchestrates all spatial steps + CLI entry)
│       │   ├── qc.py                # compute_spot_qc, compute_morans_i, flag_low_quality_spots,
│       │   │                          # compute_model_qc, boundary_sharpness
│       │   ├── neighbourhood.py     # NeighbourhoodStats, compute_neighbourhood_stats,
│       │   │                          # detect_spatial_niches, distance_decay_cooccurrence
│       │   └── benchmark.py         # BenchmarkDataset, simulate_visium, method runners
│       ├── uncertainty/
│       │   ├── bootstrap.py         # BootstrapCI (bulk gene-panel bootstrap) +
│       │   │                          # bootstrap_proportions (spatial parametric bootstrap)
│       │   │                          # Shared interface: returns (ci_lo, ci_hi, metadata)
│       │   └── stability.py         # Cross-run stability checks (future)
│       ├── plotting/
│       │   ├── style.py             # Shared matplotlib rcParams, colour palettes, figure factory
│       │   ├── bulk_plots.py        # Proportion heatmaps, R² scatter, CI bar charts
│       │   ├── spatial_plots.py     # Spatial proportion maps, convergence traces,
│       │   │                          # mismatch factor plots, niche maps
│       │   ├── qc_plots.py          # Separability heatmap, marker recall plots, Moran's I
│       │   ├── benchmark_plots.py   # Method comparison figures (merges paper_plots.py from both)
│       │   └── captions.py          # Auto-generated caption strings with parameter provenance
│       ├── report/
│       │   ├── html.py              # Unified HTMLReport (bulk and spatial variants)
│       │   ├── methods_text.py      # Auto-generated methods section text
│       │   └── templates/           # Jinja2 or f-string templates
│       └── data/
│           └── gene_lists/          # All gene list TSVs from CHIMERA (shared)
├── tests/
│   ├── shared/                      # Tests for utils, reference, gene_filters, separability
│   ├── bulk/                        # Tests ported from chimera_v1/tests/
│   └── spatial/                     # Tests ported from spatcar/tests/
├── examples/
├── benchmark/
├── docs/
│   └── qc_thresholds.md             # Threshold provenance documentation
├── scripts/
├── pyproject.toml
├── README.md
├── DESIGN_SPEC.md
└── CLAUDE.md
```

---

## 10. Staged Implementation Plan

The plan has five stages. Each stage ends with a passing test suite and a reviewable diff.

---

### Stage 0: Scaffold and audits (no scientific code yet)

**Goal:** Establish the package skeleton, resolve pre-migration blockers, and confirm legacy test suites pass unchanged.

Tasks:
1. Create the directory skeleton above. Add `pyproject.toml` with dependencies (numpy, scipy, pandas, anndata, scanpy, click, pyyaml, joblib, tqdm, scikit-learn, jinja2).
2. Run both legacy test suites (`pytest chimera_v1/tests/ spatcar/tests/`) and record baseline pass/fail.
3. Audit `chimera/gene_selection_orig.py`: if unused, delete; if referenced by tests, migrate tests. **Closes P0-4.**
4. Fix `merge_nonseparable_types` dispersion bug in SpatCAR (P0-3). Add regression test.
5. Add `__all__` to all legacy modules that are missing it.
6. Write a single `utils.py` in TissueResolve with `bh_fdr` (unified, tested against scipy), `project_simplex_batch`, `dense_subset`, `nb_loglik_batch`, `safe_log`, `set_random_state`, `setup_logging`.

Deliverable: Empty package with working imports + fixed bugs + unified `utils.py`.

---

### Stage 1: Shared reference layer

**Goal:** Unified reference construction that serves both modalities. This is the architectural foundation that everything else depends on.

Tasks:
1. Design and implement `results.py:ReferenceData` with both phi (G × K, L1-normalised) and R_cpm (K × G, CPM) representations. Add `.as_phi()`, `.subset_genes()`, `.save()`, `.load()`.
2. Implement `reference/build.py:ReferenceBuilder` that:
   - Accepts `build_from_df`, `build_from_adata`, `build_from_h5ad`, `build_from_csv`.
   - Performs donor-aware aggregation (from CHIMERA) and streaming aggregation (from SpatCAR).
   - Estimates NB overdispersion (from SpatCAR) when requested.
   - Computes cross-donor CV (from CHIMERA) when ≥2 donors present.
   - Emits `ReferenceData` with all fields.
3. Implement `reference/gene_filters.py` with instance-based gene list loading. Remove singletons. **Closes P0-5.**
4. Implement `reference/separability.py` (direct port from SpatCAR, with bug fix already applied).
5. Implement `reference/markers.py:GeneSelector` (CHIMERA composite scoring + SpatCAR pairwise discriminability). Remove hardcoded test prefixes from protein-coding mask. **Closes P2-2.**
6. Implement `io/reference.py` and `io/validation.py`.
7. Port all reference and gene-selection tests from both legacy packages. Add new unified reference tests.

Deliverable: `reference/` and `results.py` fully tested.

> **Gate:** Do not proceed to Stage 2 until the `ReferenceData` design is reviewed and approved. The data layout, normalisation conventions, and save/load format determined here govern all downstream modules.

---

### Stage 2: Protocol layer

**Goal:** Unified protocol metadata and risk/mismatch handling.

Tasks:
1. Implement `protocol/metadata.py` (enums from CHIMERA + Visium-specific enum).
2. Implement `protocol/risk.py` (direct port of `ProtocolRiskAssessor`). Copy gene list data files into `src/tissueresolve/data/gene_lists/`. Verify the data path resolver works in an installed package.
3. Implement `protocol/mismatch.py` (direct port of SpatCAR `protocol.py`). Clarify that this is "Visium-vs-reference spatial mismatch" distinct from the bulk protocol risk.
4. Port all protocol tests from both legacy packages.

Deliverable: `protocol/` fully tested.

---

### Stage 3: Bulk workflow

**Goal:** Complete bulk deconvolution pipeline with all CHIMERA features.

Tasks:
1. Implement `bulk/solver.py`: `WNNLSSolver` (direct port), `MRNAContentCorrector` (direct port). **Promote `MRNAContentCorrector` to public API. Closes P0-2.**
2. Implement `uncertainty/bootstrap.py:BootstrapCI` (direct port from CHIMERA bootstrap).
3. Implement `bulk/qc.py:BulkQCReport` (direct port). Define `BulkQCMetrics` in `results.py`.
4. Implement `bulk/pipeline.py:ChimeraPipeline` that orchestrates all bulk steps. Wire in DiscordanceFilter as optional Tier 3 step. Wire in MRNAContentCorrector output as a separate result file (not silently omitted). **Closes P1-3.**
5. Implement `bulk/benchmark.py` (port from CHIMERA benchmark; mark unimplemented scenarios explicitly with `NotImplementedError`). **Closes P1-4.**
6. Implement `io/bulk.py`.
7. Implement `plotting/bulk_plots.py`.
8. Implement `report/html.py` (bulk variant).
9. Wire `tissueresolve bulk` CLI: `run`, `check-compatibility`, `benchmark`. Add `--genome` flag. **Closes P2-4.**
10. Port all CHIMERA tests. Add integration test running the full pipeline on synthetic data.

Deliverable: `tissueresolve bulk run` end-to-end working.

---

### Stage 4: Spatial workflow

**Goal:** Complete spatial deconvolution pipeline with all SpatCAR features.

Tasks:
1. Implement `spatial/graph.py` (direct port; remove `_HEX_OFFSETS` duplication). **Closes P2-5.**
2. Implement `spatial/model.py:SpatCARModel` (direct port from SpatCAR deconvolution).
3. Implement `spatial/qc.py` (direct port from SpatCAR qc). Fix the `save_results` path to pass model arrays to `compute_spot_qc`. **Closes P1-1.**
4. Implement `spatial/neighbourhood.py` (direct port). Add `n_smooth` to `detect_spatial_niches` output. **Closes P1-7.**
5. Implement `uncertainty/bootstrap.py:bootstrap_proportions` (port from SpatCAR qc.py). Add coverage metadata to return value. **Closes P1-2.**
6. Implement `spatial/benchmark.py` (direct port from SpatCAR benchmark).
7. Implement `io/spatial.py` (direct port from SpatCAR io).
8. Implement `plotting/spatial_plots.py`.
9. Extend `report/html.py` for spatial variant.
10. Wire `tissueresolve spatial` CLI: `run`, `report`, `benchmark`, `info`.
11. Port all SpatCAR tests. Add integration test on synthetic Visium.

Deliverable: `tissueresolve spatial run` end-to-end working.

---

### Stage 5: Unification, QC review, publication outputs

**Goal:** Verify design-spec compliance, complete shared documentation, clean up all remaining issues.

Tasks:
1. Implement `plotting/style.py` (unified matplotlib config). Apply to all plot modules.
2. Implement `plotting/qc_plots.py` (separability heatmap, marker recall, Moran's I, protocol risk summary). Ensure every plot saves underlying data alongside the figure.
3. Implement `plotting/benchmark_plots.py` (merge both `paper_plots.py` files).
4. Implement `plotting/captions.py` (auto-generated captions with parameter provenance).
5. Implement `report/methods_text.py` (auto-generated methods text for papers).
6. Validate `DeconvResult` column naming to encode output type. **Closes P2-7.**
7. Document QC threshold provenance in `docs/qc_thresholds.md`. **Closes P2-3.**
8. Write full end-to-end compliance tests:
   - Bulk: snRNA reference → high intronic risk flag, correct gene exclusion.
   - Bulk: non-separable types → spillover warning in QC report.
   - Spatial: poorly separable types → `SeparabilityWarning` emitted.
   - Spatial: spatial smoothing recorded in output metadata.
   - Spatial: `nb_loglik` is non-NaN in `prop_qc` after `save_results`.
9. Final review against every non-negotiable rule in `CLAUDE.md`. Confirm:
   - No silent gene removal anywhere.
   - No hidden warnings.
   - No bulk estimates reported as absolute cell fractions.
   - No spatial smoothing without recording λ and n_smooth.
   - No tests removed.
   - No plots without saved underlying data.
   - No confident estimates for non-separable types without warnings.
   - No mixed bulk/spatial model assumptions.

Deliverable: TissueResolve v1.0.0-rc1 ready for internal review.

---

### Stage summary

| Stage | Scope | Key risk |
|---|---|---|
| 0 — Scaffold + audits | Package skeleton, bug fixes, unified utils | Low |
| 1 — Shared reference layer | Unified `ReferenceData`, gene selection, separability | **High** — architectural decision governs all downstream |
| 2 — Protocol layer | Metadata enums, risk assessor, mismatch factors | Low |
| 3 — Bulk workflow | NNLS solver, bootstrap, QC, pipeline, CLI | Medium |
| 4 — Spatial workflow | NB-CAR model, graph, QC, neighbourhood, CLI | Medium |
| 5 — Unification + QC | Plotting, reports, compliance tests | Low–Medium |
