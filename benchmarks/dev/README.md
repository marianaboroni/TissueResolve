# benchmarks/dev — exploratory & negative-result benchmark scripts

These are **exploratory / negative-result** benchmark scripts moved out of the core
`benchmarks/diagnostics/` surface for tidiness. They are kept (not deleted) for
reproducibility and provenance — including the experiments that **did not** work,
which are deliberately preserved.

They import shared infrastructure by absolute path (e.g.
`benchmarks.spatial.run_weak_smoothing_grid`, `benchmarks.shared.*`), so they still
run from the repository root with `--run-real-data`; only their location changed.

Contents (negative / exploratory):

- Spatial smoothing / state regularization (negative results):
  `combined_spatial_smoothing_benchmark.py`, `in_solver_edge_aware_benchmark.py`,
  `in_solver_state_regularization_benchmark.py`, `state_similarity_regularization_benchmark.py`,
  `adaptive_method_benchmark.py`
- Granularity / fine refinement (negative): `fine_refiner_benchmark.py`,
  `high_granularity_benchmark.py`
- Ceilings / diagnostics: `identifiability_ceiling.py`, `mixture_recovery_ceiling.py`,
  `hierarchical_oracle.py`
- Development phases / one-offs: `phase2a_soft_gating.py`, `phase2b_stage1_validation.py`,
  `phase2b_stage2_joint.py`, `phase_taskD_level_spatial.py`, `full_workflow_validation.py`,
  `external_bulk_benchmark.py` (superseded by `diagnostics/external_bulk_comparison.py`)

**Core / publishable benchmarks remain in [`benchmarks/diagnostics/`](../diagnostics/):**
`nb_bulk_solver_benchmark.py`, `external_bulk_comparison.py`,
`score_external_bulk_generic.py`, `rare_detection_calibration.py`, and the
`gold_truth_*` scripts (the last are also imported by the test suite). Shared
infra stays in `benchmarks/shared/`; pipeline runners in `benchmarks/bulk/` and
`benchmarks/spatial/`.
