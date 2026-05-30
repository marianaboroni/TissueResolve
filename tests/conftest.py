"""
Shared pytest fixtures for TissueResolve tests.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import (
    BenchmarkResult,
    BulkDeconvResult,
    PairSeparability,
    QCReport,
    ReferenceSignature,
    SeparabilityReport,
    SpatialDeconvResult,
)


# ---------------------------------------------------------------------------
# Dimensions
# ---------------------------------------------------------------------------

N_GENES = 40
N_CELL_TYPES = 4
N_SAMPLES = 6
N_SPOTS = 20
N_MARKERS = 12

GENE_NAMES = [f"GENE_{i:03d}" for i in range(N_GENES)]
CELL_TYPES = ["CellTypeA", "CellTypeB", "CellTypeC", "CellTypeD"]
SAMPLE_IDS = [f"sample_{i}" for i in range(N_SAMPLES)]
SPOT_IDS = [f"AAACCTGAGCAG-{i}" for i in range(N_SPOTS)]
MARKER_GENES = [f"GENE_{i:03d}" for i in range(N_MARKERS)]


# ---------------------------------------------------------------------------
# ReferenceSignature
# ---------------------------------------------------------------------------

@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(42)


@pytest.fixture
def phi_matrix(rng: np.random.Generator) -> np.ndarray:
    """L1-normalised (G × K) probability matrix."""
    raw = rng.exponential(1.0, (N_GENES, N_CELL_TYPES))
    col_sums = raw.sum(axis=0, keepdims=True)
    return (raw / col_sums).astype(np.float64)


@pytest.fixture
def R_cpm_matrix(phi_matrix: np.ndarray) -> np.ndarray:
    """CPM matrix (K × G) derived from phi."""
    return (phi_matrix.T * 1e6).astype(np.float32)


@pytest.fixture
def phi_g_array(rng: np.random.Generator) -> np.ndarray:
    """Per-gene NB dispersion (G,)."""
    return np.clip(rng.gamma(2.0, 5.0, N_GENES), 0.5, 100.0).astype(np.float32)


@pytest.fixture
def donor_cv_matrix(rng: np.random.Generator) -> np.ndarray:
    """Cross-donor CV (G × K)."""
    return rng.uniform(0.0, 0.5, (N_GENES, N_CELL_TYPES)).astype(np.float64)


@pytest.fixture
def ref_sig_full(
    phi_matrix: np.ndarray,
    R_cpm_matrix: np.ndarray,
    phi_g_array: np.ndarray,
    donor_cv_matrix: np.ndarray,
) -> ReferenceSignature:
    """ReferenceSignature with all arrays populated."""
    return ReferenceSignature(
        gene_names=GENE_NAMES,
        cell_types=CELL_TYPES,
        phi=phi_matrix,
        R_cpm=R_cpm_matrix,
        R_log=np.log1p(R_cpm_matrix).astype(np.float32),
        phi_g=phi_g_array,
        donor_cv=donor_cv_matrix,
        n_cells_per_type={"CellTypeA": 100, "CellTypeB": 80, "CellTypeC": 120, "CellTypeD": 90},
        genome="hg38",
    )


@pytest.fixture
def ref_sig_phi_only(phi_matrix: np.ndarray) -> ReferenceSignature:
    """ReferenceSignature with only phi (bulk-only case)."""
    return ReferenceSignature(
        gene_names=GENE_NAMES,
        cell_types=CELL_TYPES,
        phi=phi_matrix,
    )


@pytest.fixture
def ref_sig_cpm_only(
    R_cpm_matrix: np.ndarray,
    phi_g_array: np.ndarray,
) -> ReferenceSignature:
    """ReferenceSignature with only R_cpm + phi_g (spatial-only case)."""
    return ReferenceSignature(
        gene_names=GENE_NAMES,
        cell_types=CELL_TYPES,
        R_cpm=R_cpm_matrix,
        R_log=np.log1p(R_cpm_matrix).astype(np.float32),
        phi_g=phi_g_array,
    )


# ---------------------------------------------------------------------------
# BulkDeconvResult
# ---------------------------------------------------------------------------

@pytest.fixture
def bulk_proportions(rng: np.random.Generator) -> pd.DataFrame:
    """(N_SAMPLES × N_CELL_TYPES) mRNA proportion DataFrame."""
    raw = rng.dirichlet(np.ones(N_CELL_TYPES), size=N_SAMPLES)
    return pd.DataFrame(raw, index=SAMPLE_IDS, columns=CELL_TYPES)


@pytest.fixture
def bulk_result(bulk_proportions: pd.DataFrame, rng: np.random.Generator) -> BulkDeconvResult:
    r2 = pd.Series(rng.uniform(0.6, 0.98, N_SAMPLES), index=SAMPLE_IDS, name="coverage_r2")
    return BulkDeconvResult(
        proportions=bulk_proportions,
        coverage_r2=r2,
        gene_panel=GENE_NAMES[:20],
        run_metadata={"chimera_version": "0.1.0"},
    )


@pytest.fixture
def bulk_result_with_ci(
    bulk_proportions: pd.DataFrame,
    rng: np.random.Generator,
) -> BulkDeconvResult:
    r2 = pd.Series(rng.uniform(0.6, 0.98, N_SAMPLES), index=SAMPLE_IDS, name="coverage_r2")
    half_width = rng.uniform(0.01, 0.05, (N_SAMPLES, N_CELL_TYPES))
    lo = pd.DataFrame(
        np.clip(bulk_proportions.values - half_width, 0, 1),
        index=SAMPLE_IDS, columns=CELL_TYPES,
    )
    hi = pd.DataFrame(
        np.clip(bulk_proportions.values + half_width, 0, 1),
        index=SAMPLE_IDS, columns=CELL_TYPES,
    )
    return BulkDeconvResult(
        proportions=bulk_proportions,
        coverage_r2=r2,
        gene_panel=GENE_NAMES[:20],
        lower_ci=lo,
        upper_ci=hi,
    )


# ---------------------------------------------------------------------------
# SpatialDeconvResult
# ---------------------------------------------------------------------------

@pytest.fixture
def spatial_proportions(rng: np.random.Generator) -> pd.DataFrame:
    raw = rng.dirichlet(np.ones(N_CELL_TYPES), size=N_SPOTS)
    return pd.DataFrame(raw, index=SPOT_IDS, columns=CELL_TYPES)


@pytest.fixture
def spatial_result(
    spatial_proportions: pd.DataFrame,
    rng: np.random.Generator,
) -> SpatialDeconvResult:
    mismatch = rng.uniform(0.8, 1.2, N_MARKERS).astype(np.float32)
    return SpatialDeconvResult(
        proportions=spatial_proportions,
        cell_types=CELL_TYPES,
        marker_genes=MARKER_GENES,
        n_iter=87,
        converged=True,
        convergence_trace=[0.01 * (0.9 ** i) for i in range(87)],
        lambda_spatial=0.1,
        mismatch_factors=mismatch,
        run_metadata={"spatcar_version": "0.1.0"},
    )


# ---------------------------------------------------------------------------
# QCReport
# ---------------------------------------------------------------------------

@pytest.fixture
def bulk_qc_report(rng: np.random.Generator) -> QCReport:
    return QCReport(
        modality="bulk",
        recommendations=["Check sample s3 — high mismatch flag."],
        recon_r2=pd.Series(rng.uniform(0.5, 0.99, N_SAMPLES), index=SAMPLE_IDS),
        profile_corr=pd.Series(rng.uniform(0.6, 0.99, N_SAMPLES), index=SAMPLE_IDS),
        mismatch_flag=pd.Series(["low"] * (N_SAMPLES - 1) + ["high"], index=SAMPLE_IDS),
        marker_recall=pd.Series(rng.uniform(0.4, 1.0, N_CELL_TYPES), index=CELL_TYPES),
        spillover_risk=pd.Series(rng.uniform(0.1, 0.6, N_CELL_TYPES), index=CELL_TYPES),
        condition_number=320.0,
        protocol_risk_level="low",
        n_genes_excluded_protocol=0,
    )


@pytest.fixture
def spatial_qc_report(rng: np.random.Generator) -> QCReport:
    spot_df = pd.DataFrame({
        "prop_entropy": rng.uniform(0.5, 2.0, N_SPOTS),
        "nb_loglik": rng.uniform(-500, -100, N_SPOTS),
        "dominant_frac": rng.uniform(0.3, 0.9, N_SPOTS),
        "spatial_resid": rng.uniform(0.01, 0.4, N_SPOTS),
    })
    morans = pd.Series(rng.uniform(0.3, 0.9, N_CELL_TYPES), index=CELL_TYPES)
    return QCReport(
        modality="spatial",
        recommendations=[],
        spot_qc=spot_df,
        morans_i=morans,
        model_qc={"n_iter": 87, "converged": True, "d_g_median": 1.03},
    )


# ---------------------------------------------------------------------------
# SeparabilityReport
# ---------------------------------------------------------------------------

@pytest.fixture
def separability_report() -> SeparabilityReport:
    pairs = [
        PairSeparability("CellTypeA", "CellTypeB", 0.72, 1.4, 0.65, 80),
        PairSeparability("CellTypeA", "CellTypeC", 0.55, 2.8, 0.40, 130),
        PairSeparability("CellTypeA", "CellTypeD", 0.91, 0.6, 0.88, 12),
        PairSeparability("CellTypeB", "CellTypeC", 0.60, 2.2, 0.45, 110),
        PairSeparability("CellTypeB", "CellTypeD", 0.50, 3.1, 0.30, 150),
        PairSeparability("CellTypeC", "CellTypeD", 0.98, 0.1, 0.97, 3),
    ]
    pairs.sort(key=lambda p: p.bhattacharyya_coeff, reverse=True)
    return SeparabilityReport(pairs=pairs)


# ---------------------------------------------------------------------------
# BenchmarkResult
# ---------------------------------------------------------------------------

@pytest.fixture
def benchmark_result() -> BenchmarkResult:
    per_ct = pd.DataFrame(
        {"rmse": [0.04, 0.05, 0.03, 0.06], "pcc": [0.92, 0.88, 0.94, 0.85]},
        index=CELL_TYPES,
    )
    return BenchmarkResult(
        method="SpatCAR",
        scenario="basic",
        metrics={"rmse": 0.045, "pcc": 0.90, "jsd": 0.03, "morans_i_mean": 0.71},
        per_celltype_metrics=per_ct,
        n_samples_or_spots=N_SPOTS,
        run_time_s=12.4,
        run_metadata={"lambda_spatial": 0.1, "n_marker_genes": 100},
    )
