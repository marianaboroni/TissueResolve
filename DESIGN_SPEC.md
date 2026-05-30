# TissueResolve — Design Specification

## 1. Purpose

TissueResolve is a unified, reproducible, and scientifically robust package for cell-type and cell-state deconvolution across two complementary transcriptomic settings:

1. **Bulk RNA-seq deconvolution**, based on the strengths of CHIMERA.
2. **Spatial transcriptomics deconvolution**, especially 10x Genomics Visium, based on the strengths of SpatCAR.

The goal is not simply to merge two packages. The goal is to create a single stable tool with a shared scientific core, consistent input validation, transparent quality control, uncertainty estimation, reproducible outputs, and publication-quality figures.

TissueResolve should always communicate what its estimates represent. For bulk RNA-seq, outputs should be described as **RNA-derived proportions**, not necessarily absolute cell fractions. For spatial transcriptomics, outputs should be described as **spot-level RNA-derived cellular composition estimates**, not direct single-cell counts unless explicitly calibrated.

---

## 2. Legacy package summary

### 2.1 CHIMERA functionality to preserve

CHIMERA is a protocol-aware bulk RNA-seq deconvolution package. Its strongest elements are:

- Weighted non-negative least squares for bulk deconvolution.
- Protocol-aware gene selection.
- Explicit distinction between RNA proportions and absolute cell fractions.
- Gene-panel construction from a single-cell reference.
- Compatibility checks between bulk and reference protocols.
- Protocol-risk assessment for combinations such as polyA bulk, ribodepleted bulk, scRNA reference, snRNA reference, and 10x 3-prime capture.
- Exclusion or down-weighting of genes affected by known biases, including intronic-retention genes, length-biased 3-prime genes, dissociation-stress genes, hypervariable inflammatory genes, and blacklist gene prefixes.
- Bootstrap confidence intervals.
- Per-sample quality control using reconstruction R², profile correlation, and mismatch flags.
- Per-cell-type quality control using marker recall and spillover risk.
- Benchmarking against baseline deconvolution.
- Self-contained HTML report.
- CLI commands for running deconvolution, checking compatibility, building references, and benchmarking.

These elements should be migrated into TissueResolve’s bulk module and generalized where useful.

### 2.2 SpatCAR functionality to preserve

SpatCAR is a spatial deconvolution package for 10x Visium. Its strongest elements are:

- Loading of SpaceRanger/Visium outputs and h5ad files.
- Pseudobulk reference construction from scRNA-seq or snRNA-seq.
- Marker-gene selection with gene-overlap validation.
- Spot-level NNLS initialization.
- Negative-binomial count model for spatial deconvolution.
- CAR-style spatial smoothing using a hexagonal spot graph.
- Optional protocol mismatch correction using per-gene scale factors.
- Sparse and memory-conscious computation.
- Spatial graph construction from array coordinates.
- Spot-level quality control.
- Moran’s I spatial autocorrelation.
- Boundary sharpness and spatial residual metrics.
- Pairwise cell-type separability analysis using similarity/discriminability metrics.
- Optional merging of non-separable cell types.
- Neighbourhood enrichment and spatial niche detection.
- Synthetic Visium benchmark with spatial patterns.
- Publication-style plotting functions for spatial maps, convergence, QC, mismatch factors, and benchmark comparisons.
- HTML report generation.
- CLI commands for run, report, benchmark, and info.

These elements should be migrated into TissueResolve’s spatial module and integrated with the shared QC, plotting, and reporting system.

---

## 3. Design principle

TissueResolve should have one unified package but two clearly separated scientific workflows:

1. `tissueresolve bulk` for bulk RNA-seq.
2. `tissueresolve spatial` for 10x Visium and related spatial transcriptomics data.

The two workflows should share:

- Reference construction.
- Gene harmonization.
- Marker selection.
- Protocol-aware gene filtering.
- Separability diagnostics.
- QC framework.
- Uncertainty estimation.
- Benchmark metrics.
- Plot style.
- Report generation.
- Result objects.
- Logging and reproducibility metadata.

The two workflows should not share algorithmic assumptions when those assumptions are modality-specific. Bulk deconvolution and spatial deconvolution must remain separate modules.

---

## 4. Required package structure

```text
TissueResolve/
├── src/
│   └── tissueresolve/
│       ├── __init__.py
│       ├── config.py
│       ├── cli.py
│       ├── io/
│       │   ├── bulk.py
│       │   ├── spatial.py
│       │   ├── reference.py
│       │   └── validation.py
│       ├── reference/
│       │   ├── build.py
│       │   ├── markers.py
│       │   ├── gene_matching.py
│       │   ├── gene_filters.py
│       │   └── separability.py
│       ├── protocol/
│       │   ├── metadata.py
│       │   ├── risk.py
│       │   └── mismatch.py
│       ├── bulk/
│       │   ├── solver.py
│       │   ├── pipeline.py
│       │   ├── qc.py
│       │   └── benchmark.py
│       ├── spatial/
│       │   ├── model.py
│       │   ├── graph.py
│       │   ├── pipeline.py
│       │   ├── qc.py
│       │   ├── neighbourhood.py
│       │   └── benchmark.py
│       ├── uncertainty/
│       │   ├── bootstrap.py
│       │   └── stability.py
│       ├── plotting/
│       │   ├── style.py
│       │   ├── bulk_plots.py
│       │   ├── spatial_plots.py
│       │   ├── qc_plots.py
│       │   ├── benchmark_plots.py
│       │   └── captions.py
│       ├── report/
│       │   ├── html.py
│       │   ├── methods_text.py
│       │   └── templates/
│       ├── results.py
│       └── utils.py
├── tests/
├── examples/
├── benchmark/
├── docs/
├── scripts/
├── pyproject.toml
├── README.md
├── DESIGN_SPEC.md
└── CLAUDE.md