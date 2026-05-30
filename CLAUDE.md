You are helping develop TissueResolve, a scientific package for robust cell-type and cell-state deconvolution in bulk RNA-seq and spatial transcriptomics.

TissueResolve is being built from two legacy packages:

1. CHIMERA: protocol-aware bulk RNA-seq deconvolution.
2. SpatCAR: spatial deconvolution for 10x Visium.

Always prioritize:
1. Scientific validity.
2. Reproducibility.
3. Clear software architecture.
4. Explicit assumptions.
5. Tests.
6. Robust error handling.
7. Publication-quality outputs.
8. Transparent uncertainty and QC.

Do not simply merge code mechanically.

Preserve the best scientific ideas from both packages:
- CHIMERA’s protocol-aware gene selection, weighted NNLS, bootstrap confidence intervals, compatibility QC, and explicit RNA-proportion warning.
- SpatCAR’s Visium input handling, spatial graph, negative-binomial CAR model, mismatch correction, spatial QC, separability diagnostics, neighbourhood statistics, and benchmark framework.

Non-negotiable rules:
- Never silently remove genes.
- Never hide warnings.
- Never report bulk estimates as absolute cell fractions.
- Never smooth spatial estimates without reporting the smoothing parameter.
- Never remove tests to make the package pass.
- Never generate plots without saving the underlying data.
- Never report confident estimates for non-separable cell types without warnings.
- Never mix bulk and spatial assumptions in a way that makes the model unclear.

Before editing code:
1. Read DESIGN_SPEC.md.
2. Audit the legacy packages.
3. Create a migration map.
4. Identify P0, P1, and P2 issues.
5. Propose a staged implementation plan.