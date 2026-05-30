"""
Spatial deconvolution workflow for 10x Visium (SpatCAR algorithm).

This subpackage must not import from ``bulk`` or mix spatial and bulk
model assumptions.

Public API (Stage 4)
--------------------
``graph.SpatialGraph``
    Hexagonal adjacency graph for a Visium section.
``model.SpatCARModel``
    NB-CAR deconvolution model (array-based interface).
``qc.*``
    Spot-level QC, Moran's I, model-level QC, boundary sharpness.
``neighbourhood.*``
    Co-occurrence statistics, niche detection, distance decay.
``pipeline.SpatialPipeline``
    Full orchestration → SpatialPipelineResult.
``benchmark.*``
    Synthetic Visium simulation and method comparison (SpatCAR/DWLS/
    Spatial-NNLS/NNLS).
"""
from tissueresolve.spatial.graph import (
    SpatialGraph,
    build_hex_graph,
    build_hex_graph_from_arrays,
    build_expression_weighted_graph,
    compute_spatial_batches,
)
from tissueresolve.spatial.model import SpatCARModel
from tissueresolve.spatial.qc import (
    compute_spot_qc,
    compute_morans_i,
    flag_low_quality_spots,
    compute_model_qc,
    boundary_sharpness,
)
from tissueresolve.spatial.neighbourhood import (
    NeighbourhoodStats,
    compute_neighbourhood_stats,
    detect_spatial_niches,
    distance_decay_cooccurrence,
)
from tissueresolve.spatial.pipeline import SpatialPipeline, SpatialPipelineResult
from tissueresolve.spatial.benchmark import (
    BenchmarkDataset,
    simulate_visium,
    apply_mismatch,
    shuffle_coordinates,
    run_baseline_nnls,
    run_baseline_dwls,
    run_baseline_spatial_nnls,
    run_spatcar,
    compute_benchmark_metrics,
    run_benchmark_suite,
    benchmark_results_to_frame,
)

__all__ = [
    "SpatialGraph",
    "build_hex_graph",
    "build_hex_graph_from_arrays",
    "build_expression_weighted_graph",
    "compute_spatial_batches",
    "SpatCARModel",
    "compute_spot_qc",
    "compute_morans_i",
    "flag_low_quality_spots",
    "compute_model_qc",
    "boundary_sharpness",
    "NeighbourhoodStats",
    "compute_neighbourhood_stats",
    "detect_spatial_niches",
    "distance_decay_cooccurrence",
    "SpatialPipeline",
    "SpatialPipelineResult",
    "BenchmarkDataset",
    "simulate_visium",
    "apply_mismatch",
    "shuffle_coordinates",
    "run_baseline_nnls",
    "run_baseline_dwls",
    "run_baseline_spatial_nnls",
    "run_spatcar",
    "compute_benchmark_metrics",
    "run_benchmark_suite",
    "benchmark_results_to_frame",
]
