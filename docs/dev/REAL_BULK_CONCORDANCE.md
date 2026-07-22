# Real-bulk (TCGA-TNBC) validation — concordance + marker coherence

**Status: truth-free real-bulk gate; NO accuracy claim; NO default changed.** TCGA-TNBC
bulk (40 primary tumors) has no cell-composition ground truth, so this evaluates broad-level
**marker coherence** (does each estimated broad family track its canonical markers across
tumors) and **concordance with Rectangle** (a strong external tool). Script:
`benchmarks/dev/real_bulk_concordance.py` (real-data-gated; gitignored outputs).

## Results

Mean marker coherence (Spearman, higher = fraction tracks its markers):

| method | mean coherence | concordance w/ Rectangle |
|---|---|---|
| wNNLS (current default) | 0.154 | 0.250 |
| Poisson (all genes) | 0.285 | **0.891** |
| **Poisson + donor_de** | **0.457** | (not re-run; ~0.89 expected) |
| Rectangle | 0.421 | — |

Per-family, `donor_de` fixes the negative-coherence families of all-genes Poisson:
**T/NK −0.04 → +0.51**, Endothelial −0.38 → +0.41, Mural −0.22 → +0.42, Epithelial
0.41 → 0.51, Myeloid 0.74 → 0.80 — at a cost to the already-excellent B/Plasma
(0.94 → 0.64) and Stromal (0.90 → 0.65, still strong). Net mean 0.285 → **0.457**.

*(First run, before the align fix below, `Poisson_donorDE` was identical to `Poisson`
because the stored panel was ignored by the solver backbone; the numbers above are after
the fix. Rectangle mean coherence in the first run was 0.421, concordance Poisson↔Rectangle
0.891, wNNLS↔Rectangle 0.250.)*

Per-family highlights: B/Plasma (Poisson 0.94 > Rectangle 0.81), Stromal/Fibroblast (0.90 >
0.86), Myeloid (0.74 ≈ 0.75, both ≫ wNNLS 0.53), Epithelial (Poisson 0.41 vs wNNLS **−0.13**).
The wNNLS default returned **degenerate (near-zero/constant) fractions** for Adipocyte,
B/Plasma, Mural, Stromal (coherence undefined) — Poisson estimated them non-degenerately.

## Findings

1. **Poisson beats the wNNLS default on real bulk** (broad level): ~2× marker coherence,
   non-degenerate family estimates, and **high concordance with Rectangle (0.89)** while the
   default is an outlier (0.25). Supportive real-bulk evidence for promoting Poisson.
2. **Rectangle still leads mean coherence (0.42)** — Poisson closes most of the gap; the
   residual is the signature, consistent with the pseudobulk benchmark.
3. **Caveat — T/NK**: Poisson coherence −0.04 vs wNNLS 0.46 / Rectangle 0.44. On the
   collinear T/NK family the all-genes flat Poisson misallocates. A real weakness to watch.
4. **Endothelial / Adipocyte / Mural**: poor for all methods (rare/hard in bulk TNBC).

## Integration nuance (follow-up)

`Poisson` and `Poisson_donorDE` were identical because **`ref.selected_genes` (the donor_de
panel) is honored by the `BulkPipeline` path (default-wNNLS flat, hierarchical,
pipeline-Poisson) but NOT by the `solver="poisson"` backbone flat path** (`_solver_bulk_result`
→ `PoissonGLMSolver`, which uses all genes). So donor_de's real-bulk effect is untested here.
Fix: make the solver backbones default `genes=ref.selected_genes` when `genes is None` and a
panel is stored (or route flat GLM through the pipeline). Deferred.

## Verdict

- **Poisson → stronger `candidate_for_default`**: it now has broad-level real-bulk support
  (coherence + Rectangle concordance) on top of the pseudobulk §28 gates. A cautious flip
  would still address the **T/NK coherence regression** first and confirm on a second real
  bulk cohort. No fine/rare claim on TCGA (no truth).
- **donor_de**: real-bulk effect **not yet tested** (integration nuance above) — remains
  pseudobulk-validated only.
