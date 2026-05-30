"""
Synthetic Visium benchmark framework for the TissueResolve spatial workflow.

This module ports SpatCAR's benchmark framework onto TissueResolve's
array-based :class:`~tissueresolve.spatial.model.SpatCARModel` interface and
the unified :class:`~tissueresolve.results.BenchmarkResult` container.

Four deconvolution methods are compared:

==============  ===========================================================
SpatCAR         NB-MAP + CAR spatial prior (this package)
DWLS            Dampened Weighted Least Squares (RCTD-like per-spot WLS)
Spatial-NNLS    NNLS + post-hoc neighbour smoothing (naive-spatial / CARD-like)
NNLS            Plain per-spot non-negative least squares (no model, no spatial)
==============  ===========================================================

Scenarios:

==========  ==================================================================
basic       spatially-structured proportions, matched reference
mismatch    in-silico protocol mismatch (a fraction of reference genes rescaled)
ablation    SpatCAR with ``lambda_spatial=0`` (no spatial prior)
shuffle     spot coordinates permuted — negative control for spatial structure
==========  ==================================================================

Metrics: RMSE, Jensen–Shannon divergence, Pearson r (per cell type), Moran's I.

Scientific notes
----------------
* All synthetic counts are drawn from the same negative-binomial generative
  model assumed by :class:`SpatCARModel`, so this is a *self-consistency*
  benchmark, not an independent validation.  Real-data validation lives under
  ``examples/`` and must never weaken these algorithms.
* The ``shuffle`` negative control exists precisely so that a method which
  "improves" RMSE only by exploiting (now meaningless) spatial structure is
  exposed: a well-behaved spatial method must *not* beat plain NNLS once the
  coordinates are randomised.
* Every run is seeded.  ``run_benchmark_suite`` writes the underlying tidy
  results table to disk before returning, so any downstream plot has saved
  data behind it.
"""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy.optimize as opt
import scipy.sparse as sp

from tissueresolve.results import BenchmarkResult, ReferenceSignature
from tissueresolve.spatial.graph import SpatialGraph, build_hex_graph_from_arrays
from tissueresolve.spatial.qc import compute_morans_i

__all__ = [
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

logger = logging.getLogger("tissueresolve.spatial.benchmark")

SCENARIOS = ("basic", "mismatch", "ablation", "shuffle")
METHODS = ("SpatCAR", "DWLS", "Spatial-NNLS", "NNLS")


# ---------------------------------------------------------------------------
# Synthetic dataset container
# ---------------------------------------------------------------------------


@dataclass
class BenchmarkDataset:
    """A synthetic Visium section with known ground-truth proportions.

    Attributes
    ----------
    Y_sparse:
        Integer count matrix, shape ``(N, G)``, CSR sparse.
    array_row / array_col:
        Integer Visium array coordinates, shape ``(N,)``.
    Pi_true:
        Ground-truth proportion matrix, shape ``(N, K)``, rows sum to 1.
    R_cpm:
        Reference CPM matrix, shape ``(K, G)`` float32.
    phi_g:
        Per-gene NB dispersion, shape ``(G,)`` float32.
    cell_types / gene_names:
        Ordered labels.
    scenario:
        Scenario tag (``"basic"``, ``"mismatch"``, ``"shuffle"``).
    lib_sizes:
        Per-spot library sizes, shape ``(N,)``.  Derived from ``Y_sparse``
        row sums when not supplied.
    """

    Y_sparse: sp.csr_matrix
    array_row: np.ndarray
    array_col: np.ndarray
    Pi_true: np.ndarray
    R_cpm: np.ndarray
    phi_g: np.ndarray
    cell_types: list[str]
    gene_names: list[str]
    scenario: str = "basic"
    lib_sizes: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))

    def __post_init__(self) -> None:
        if len(self.lib_sizes) == 0:
            self.lib_sizes = np.asarray(
                self.Y_sparse.sum(axis=1), dtype=np.float32
            ).ravel()

    @property
    def n_spots(self) -> int:
        return self.Y_sparse.shape[0]

    @property
    def n_types(self) -> int:
        return len(self.cell_types)

    def dense_counts(self) -> np.ndarray:
        """Return the dense ``(N, G)`` count matrix as float32."""
        return np.asarray(self.Y_sparse.todense(), dtype=np.float32)

    def to_reference(self) -> ReferenceSignature:
        """Build the matched :class:`ReferenceSignature` for the model."""
        return ReferenceSignature(
            gene_names=list(self.gene_names),
            cell_types=list(self.cell_types),
            R_cpm=self.R_cpm.astype(np.float32),
            R_log=np.log1p(self.R_cpm).astype(np.float32),
            phi_g=self.phi_g.astype(np.float32),
            n_cells_per_type={ct: 100 for ct in self.cell_types},
        )

    def to_graph(self, *, batch_size: int = 200) -> SpatialGraph:
        """Build the hexagonal spatial graph from this dataset's coordinates."""
        spot_ids = np.array([f"s{i}" for i in range(self.n_spots)])
        return build_hex_graph_from_arrays(
            self.array_row, self.array_col, spot_ids=spot_ids, batch_size=batch_size
        )


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def simulate_visium(
    n_spots: int = 500,
    n_types: int = 5,
    n_genes: int = 200,
    *,
    mean_lib_size: float = 3000.0,
    phi_mean: float = 10.0,
    spatial_smoothness: float = 0.7,
    seed: int = 0,
) -> BenchmarkDataset:
    """Generate a synthetic Visium section with spatially-structured NB counts.

    Each cell type owns a disjoint block of high-expression marker genes.
    Proportions vary smoothly across the hexagonal array (Gaussian-kernel
    interpolation of random anchors), and counts are drawn from the same
    NB model the :class:`SpatCARModel` assumes.

    Parameters
    ----------
    n_spots, n_types, n_genes:
        Section dimensions.
    mean_lib_size:
        Mean per-spot total counts (gamma-distributed around this).
    phi_mean:
        Mean per-gene NB dispersion (``phi_g``).
    spatial_smoothness:
        Larger → smoother spatial proportion fields.
    seed:
        RNG seed (reproducible).

    Returns
    -------
    BenchmarkDataset
    """
    rng = np.random.default_rng(seed)

    # Reference: each type has a block of high-expression markers.
    R_raw = rng.exponential(1.0, (n_types, n_genes)).astype(np.float64)
    for k in range(n_types):
        top = rng.choice(n_genes, max(1, n_genes // n_types), replace=False)
        R_raw[k, top] *= 20.0
    R_cpm = (R_raw / R_raw.sum(1, keepdims=True) * 1e6 + 1.0).astype(np.float32)

    phi_g = np.clip(
        rng.gamma(2.0, phi_mean / 2.0, n_genes), 0.5, 100.0
    ).astype(np.float32)

    array_row, array_col = _make_hex_grid(n_spots, rng)
    N = len(array_row)
    Pi_true = _make_smooth_proportions(
        array_row, array_col, n_types, spatial_smoothness, rng
    )
    lib_sizes = rng.gamma(5.0, mean_lib_size / 5.0, N).astype(np.float32)
    Y = _sample_nb_counts(Pi_true, R_cpm, phi_g, lib_sizes, rng)

    return BenchmarkDataset(
        Y_sparse=sp.csr_matrix(Y.astype(np.int32)),
        array_row=array_row,
        array_col=array_col,
        Pi_true=Pi_true,
        R_cpm=R_cpm,
        phi_g=phi_g,
        cell_types=[f"CellType_{k}" for k in range(n_types)],
        gene_names=[f"Gene_{g}" for g in range(n_genes)],
        scenario="basic",
        lib_sizes=lib_sizes,
    )


def apply_mismatch(
    dataset: BenchmarkDataset,
    *,
    mismatch_frac: float = 0.2,
    mismatch_log_sd: float = 1.5,
    seed: int = 1,
) -> tuple[BenchmarkDataset, np.ndarray]:
    """Introduce an in-silico protocol mismatch into the *reference*.

    A fraction of reference genes are rescaled by a log-normal factor while
    the observed counts ``Y`` stay fixed — this is exactly the situation the
    model's per-gene mismatch factors ``d_g`` are designed to correct.  The
    returned ``true_d`` is the ground-truth scale vector for evaluating the
    recovered factors.

    Parameters
    ----------
    dataset:
        Source dataset (counts are shared, reference is rescaled).
    mismatch_frac:
        Fraction of genes affected.
    mismatch_log_sd:
        Log-scale standard deviation of the rescaling factors.
    seed:
        RNG seed.

    Returns
    -------
    (BenchmarkDataset, np.ndarray)
        New dataset (``scenario="mismatch"``) and the true per-gene factors.
    """
    rng = np.random.default_rng(seed)
    G = dataset.R_cpm.shape[1]
    n = max(1, int(G * mismatch_frac))
    idx = rng.choice(G, n, replace=False)
    true_d = np.ones(G, np.float32)
    mults = np.exp(rng.normal(0, mismatch_log_sd, n)).astype(np.float32)
    true_d[idx] = mults

    R_mis = dataset.R_cpm.copy()
    R_mis[:, idx] /= mults
    R_mis = np.maximum(R_mis, 1.0)

    return (
        BenchmarkDataset(
            Y_sparse=dataset.Y_sparse,
            array_row=dataset.array_row,
            array_col=dataset.array_col,
            Pi_true=dataset.Pi_true,
            R_cpm=R_mis,
            phi_g=dataset.phi_g,
            cell_types=dataset.cell_types,
            gene_names=dataset.gene_names,
            scenario="mismatch",
            lib_sizes=dataset.lib_sizes,
        ),
        true_d,
    )


def shuffle_coordinates(dataset: BenchmarkDataset, seed: int = 2) -> BenchmarkDataset:
    """Randomly permute spot coordinates — spatial negative control.

    Counts, reference and ground-truth proportions are unchanged; only the
    ``(array_row, array_col)`` assignment is permuted, destroying the spatial
    correlation between neighbouring spots.  Spatial methods must not gain an
    advantage over plain NNLS here.
    """
    rng = np.random.default_rng(seed)
    p = rng.permutation(dataset.n_spots)
    return BenchmarkDataset(
        Y_sparse=dataset.Y_sparse,
        array_row=dataset.array_row[p],
        array_col=dataset.array_col[p],
        Pi_true=dataset.Pi_true,
        R_cpm=dataset.R_cpm,
        phi_g=dataset.phi_g,
        cell_types=dataset.cell_types,
        gene_names=dataset.gene_names,
        scenario="shuffle",
        lib_sizes=dataset.lib_sizes,
    )


# ---------------------------------------------------------------------------
# Baseline methods
# ---------------------------------------------------------------------------


def run_baseline_nnls(dataset: BenchmarkDataset) -> np.ndarray:
    """Plain per-spot NNLS — no likelihood model, no spatial structure.

    Solves ``argmin_{π≥0} ‖y_s/lib_s − R·π‖₂²`` independently per spot,
    then renormalises to the simplex.  The simplest deconvolution baseline.
    """
    Y = dataset.dense_counts().astype(np.float64)
    R = dataset.R_cpm.T.astype(np.float64) / 1e6   # (G, K)
    lib = dataset.lib_sizes.astype(np.float64)

    Pi = np.zeros((Y.shape[0], R.shape[1]), dtype=np.float32)
    for s in range(Y.shape[0]):
        pi, _ = opt.nnls(R, Y[s] / max(lib[s], 1.0))
        sm = pi.sum()
        Pi[s] = (pi / max(sm, 1e-10)).astype(np.float32)
    return Pi


def run_baseline_dwls(
    dataset: BenchmarkDataset,
    *,
    n_iter: int = 10,
    damping: float = 0.01,
) -> np.ndarray:
    """Dampened Weighted Least Squares — RCTD-like per-spot inference.

    Iteratively reweights genes by ``w_g ← 1/(fitted_g + η)`` so that highly
    expressed genes receive lower weight.  No spatial component.
    """
    Y = dataset.dense_counts().astype(np.float64)
    R = dataset.R_cpm.T.astype(np.float64) / 1e6   # (G, K)
    lib = dataset.lib_sizes.astype(np.float64)
    G, K = R.shape

    Pi = np.zeros((Y.shape[0], K), dtype=np.float32)
    for s in range(Y.shape[0]):
        y_n = Y[s] / max(lib[s], 1.0)
        w = np.ones(G)
        pi = np.zeros(K)
        for _ in range(n_iter):
            sw = np.sqrt(w)
            pi, _ = opt.nnls(R * sw[:, None], y_n * sw)
            w = 1.0 / (R @ pi + damping)
        sm = pi.sum()
        Pi[s] = (pi / max(sm, 1e-10)).astype(np.float32)
    return Pi


def run_baseline_spatial_nnls(
    dataset: BenchmarkDataset,
    graph: SpatialGraph,
    *,
    alpha: float = 0.3,
    n_smooth: int = 3,
) -> np.ndarray:
    """NNLS + post-hoc neighbour smoothing — naive / CARD-like spatial.

    Runs plain NNLS, then applies ``n_smooth`` passes of spatial averaging
    ``Π ← (1−α)·Π + α·A·Π`` (renormalised each pass).  Spatial information is
    injected *after* inference rather than inside it.

    The smoothing parameters ``alpha`` and ``n_smooth`` are explicit
    arguments and are recorded by :func:`run_benchmark_suite` in the run
    metadata — spatial estimates are never smoothed without reporting how.
    """
    Pi = run_baseline_nnls(dataset).astype(np.float64)
    A = graph.A   # row-normalised CSR
    for _ in range(n_smooth):
        Pi = (1.0 - alpha) * Pi + alpha * (A @ Pi)
        Pi = np.maximum(Pi, 1e-10)
        Pi /= Pi.sum(1, keepdims=True)
    return Pi.astype(np.float32)


def run_spatcar(
    dataset: BenchmarkDataset,
    graph: SpatialGraph,
    *,
    lambda_spatial: float = 0.1,
    max_iter: int = 100,
    tol: float = 1e-4,
    random_state: int = 0,
) -> np.ndarray:
    """Run the TissueResolve :class:`SpatCARModel` on a benchmark dataset.

    All synthetic genes act as marker genes, so the full count matrix and
    full reference are passed to the array-based model interface.
    """
    from tissueresolve.spatial.model import SpatCARModel

    Y_marker = dataset.dense_counts()
    ref_marker = dataset.to_reference()
    model = SpatCARModel(
        lambda_spatial=lambda_spatial,
        max_iter=max_iter,
        tol=tol,
        random_state=random_state,
        verbose=False,
    )
    model.fit(Y_marker, ref_marker, graph, dataset.lib_sizes)
    return model.proportions_


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def compute_benchmark_metrics(
    Pi_pred: np.ndarray, Pi_true: np.ndarray
) -> dict[str, float]:
    """Aggregate accuracy metrics comparing predicted vs true proportions.

    Returns a flat dict with ``rmse``, ``jsd`` (mean Jensen–Shannon
    divergence per spot), ``pcc`` (mean per-cell-type Pearson r) and
    ``pcc_std``.  Cell types with no variance in the truth are excluded from
    the Pearson average (their correlation is undefined).
    """
    rmse = float(np.sqrt(np.mean((Pi_pred - Pi_true) ** 2)))
    jsd = float(np.mean(_jsd_batch(Pi_pred, Pi_true)))
    pccs = _per_celltype_pcc(Pi_pred, Pi_true)
    valid = [r for r in pccs if np.isfinite(r)]
    return {
        "rmse": rmse,
        "jsd": jsd,
        "pcc": float(np.mean(valid)) if valid else float("nan"),
        "pcc_std": float(np.std(valid)) if valid else float("nan"),
    }


def _per_celltype_metrics_frame(
    Pi_pred: np.ndarray, Pi_true: np.ndarray, cell_types: list[str]
) -> pd.DataFrame:
    """Per-cell-type RMSE and Pearson r as a tidy DataFrame."""
    rmse = np.sqrt(np.mean((Pi_pred - Pi_true) ** 2, axis=0))
    pcc = _per_celltype_pcc(Pi_pred, Pi_true)
    return pd.DataFrame(
        {"rmse": rmse.astype(np.float64), "pcc": pcc},
        index=cell_types,
    )


def _per_celltype_pcc(Pi_pred: np.ndarray, Pi_true: np.ndarray) -> np.ndarray:
    K = Pi_pred.shape[1]
    out = np.full(K, np.nan, dtype=np.float64)
    for k in range(K):
        if Pi_true[:, k].std() < 1e-8 or Pi_pred[:, k].std() < 1e-8:
            continue
        r = np.corrcoef(Pi_pred[:, k], Pi_true[:, k])[0, 1]
        if np.isfinite(r):
            out[k] = float(r)
    return out


def _jsd_batch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    eps = 1e-10
    P = np.maximum(P, eps)
    P = P / P.sum(1, keepdims=True)
    Q = np.maximum(Q, eps)
    Q = Q / Q.sum(1, keepdims=True)
    M = 0.5 * (P + Q)
    return 0.5 * (
        np.sum(P * np.log(P / M), 1) + np.sum(Q * np.log(Q / M), 1)
    )


# ---------------------------------------------------------------------------
# Full suite
# ---------------------------------------------------------------------------


def run_benchmark_suite(
    n_spots: int = 500,
    n_types: int = 5,
    n_genes: int = 200,
    *,
    scenarios: tuple[str, ...] | str = "all",
    methods: tuple[str, ...] = METHODS,
    output_dir: Optional[Path | str] = None,
    seed: int = 0,
    spatcar_max_iter: int = 100,
    spatcar_tol: float = 1e-4,
    snnls_alpha: float = 0.3,
    snnls_n_smooth: int = 3,
) -> list[BenchmarkResult]:
    """Run the full benchmark: each scenario × each method.

    Parameters
    ----------
    n_spots, n_types, n_genes:
        Synthetic section dimensions.
    scenarios:
        ``"all"`` or a subset of :data:`SCENARIOS`.
    methods:
        Subset of :data:`METHODS`.
    output_dir:
        When given, the tidy results table is written to
        ``<output_dir>/benchmark_results.csv`` *before* returning, so that
        any later figure has its underlying data saved.
    seed:
        Master RNG seed; per-scenario seeds are offset deterministically.
    spatcar_max_iter, spatcar_tol:
        SpatCAR solver settings.
    snnls_alpha, snnls_n_smooth:
        Spatial-NNLS post-hoc smoothing parameters (recorded in metadata).

    Returns
    -------
    list[BenchmarkResult]
        One :class:`~tissueresolve.results.BenchmarkResult` per
        (scenario, method) pair.

    Raises
    ------
    ValueError
        If an unknown scenario or method is requested.
    """
    run_sc = SCENARIOS if scenarios == "all" else _validate_choices(
        scenarios, SCENARIOS, "scenario"
    )
    run_methods = _validate_choices(methods, METHODS, "method")

    base = simulate_visium(n_spots, n_types, n_genes, seed=seed)
    results: list[BenchmarkResult] = []

    for scenario in run_sc:
        logger.info("── Scenario: %s", scenario)

        if scenario == "basic":
            ds, lam = base, 0.1
        elif scenario == "mismatch":
            ds, _ = apply_mismatch(base, seed=seed + 10)
            lam = 0.1
        elif scenario == "ablation":
            ds, lam = base, 0.0
        elif scenario == "shuffle":
            ds, lam = shuffle_coordinates(base, seed=seed + 20), 0.1
        else:  # pragma: no cover — guarded by _validate_choices
            raise ValueError(f"Unknown scenario {scenario!r}.")

        graph = ds.to_graph()

        for method in run_methods:
            t0 = time.perf_counter()
            meta: dict = {"scenario_lambda_spatial": lam}

            if method == "SpatCAR":
                Pi = run_spatcar(
                    ds, graph, lambda_spatial=lam,
                    max_iter=spatcar_max_iter, tol=spatcar_tol,
                    random_state=seed,
                )
                meta.update(lambda_spatial=lam, max_iter=spatcar_max_iter, tol=spatcar_tol)
            elif method == "DWLS":
                Pi = run_baseline_dwls(ds)
            elif method == "Spatial-NNLS":
                Pi = run_baseline_spatial_nnls(
                    ds, graph, alpha=snnls_alpha, n_smooth=snnls_n_smooth
                )
                meta.update(smoothing_alpha=snnls_alpha, n_smooth=snnls_n_smooth)
            elif method == "NNLS":
                Pi = run_baseline_nnls(ds)
            else:  # pragma: no cover — guarded by _validate_choices
                raise ValueError(f"Unknown method {method!r}.")

            elapsed = time.perf_counter() - t0

            metrics = compute_benchmark_metrics(Pi, ds.Pi_true)
            mi = compute_morans_i(Pi, graph, ds.cell_types)
            mi_true = compute_morans_i(ds.Pi_true, graph, ds.cell_types)
            metrics["morans_i_mean"] = float(mi.mean())
            metrics["morans_i_true_mean"] = float(mi_true.mean())

            results.append(
                BenchmarkResult(
                    method=method,
                    scenario=scenario,
                    metrics=metrics,
                    per_celltype_metrics=_per_celltype_metrics_frame(
                        Pi, ds.Pi_true, ds.cell_types
                    ),
                    n_samples_or_spots=ds.n_spots,
                    run_time_s=round(elapsed, 3),
                    run_metadata={"seed": seed, "n_genes": n_genes, **meta},
                )
            )
            logger.info(
                "  %-13s RMSE=%.4f  PCC=%.3f  Moran=%.3f  (%.2fs)",
                method, metrics["rmse"], metrics["pcc"],
                metrics["morans_i_mean"], elapsed,
            )

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        frame = benchmark_results_to_frame(results)
        frame.to_csv(output_dir / "benchmark_results.csv", index=False)
        logger.info("Saved benchmark table → %s", output_dir / "benchmark_results.csv")

    return results


def benchmark_results_to_frame(results: list[BenchmarkResult]) -> pd.DataFrame:
    """Flatten a list of :class:`BenchmarkResult` into a tidy DataFrame.

    One row per (scenario, method); metric keys become columns.
    """
    rows = []
    for r in results:
        row = {"scenario": r.scenario, "method": r.method,
               "n_spots": r.n_samples_or_spots, "wall_time_s": r.run_time_s}
        row.update({k: round(float(v), 4) for k, v in r.metrics.items()})
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _validate_choices(
    requested: tuple[str, ...] | str, allowed: tuple[str, ...], kind: str
) -> tuple[str, ...]:
    if isinstance(requested, str):
        requested = tuple(s.strip() for s in requested.split(",") if s.strip())
    unknown = [r for r in requested if r not in allowed]
    if unknown:
        raise ValueError(
            f"Unknown {kind}(s) {unknown}.  Allowed: {list(allowed)}."
        )
    return tuple(requested)


def _make_hex_grid(n_spots: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n_rows = max(2, int(np.ceil(np.sqrt(n_spots / 0.8))))
    n_cols = max(2, int(np.ceil(n_spots / n_rows)))
    rs, cs = [], []
    for r in range(n_rows):
        c0 = 0 if r % 2 == 0 else 1
        for c in range(n_cols):
            rs.append(r)
            cs.append(c0 + c * 2)
    ar = np.array(rs, np.int32)
    ac = np.array(cs, np.int32)
    if len(ar) > n_spots:
        idx = np.sort(rng.choice(len(ar), n_spots, replace=False))
        ar, ac = ar[idx], ac[idx]
    return ar, ac


def _make_smooth_proportions(
    ar: np.ndarray, ac: np.ndarray, K: int, smoothness: float, rng: np.random.Generator
) -> np.ndarray:
    N = len(ar)
    coords = np.stack([ar, ac], 1).astype(np.float32)
    n_a = max(K, int(N ** 0.5))
    a_i = rng.choice(N, min(n_a, N), replace=False)
    av = rng.standard_normal((len(a_i), K)).astype(np.float32)
    sc = max(smoothness * float(np.std(coords)), 1e-3)
    diff = coords[:, None, :] - coords[a_i][None, :, :]
    w = np.exp(-(diff ** 2).sum(2) / (2 * sc ** 2))
    w /= w.sum(1, keepdims=True) + 1e-10
    logits = w @ av
    logits -= logits.max(1, keepdims=True)
    Pi = np.exp(logits)
    Pi /= Pi.sum(1, keepdims=True)
    return Pi.astype(np.float32)


def _sample_nb_counts(
    Pi: np.ndarray,
    R_cpm: np.ndarray,
    phi_g: np.ndarray,
    lib_sizes: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    Mu = lib_sizes[:, None] * (Pi @ R_cpm) / 1e6
    N, G = Mu.shape
    Y = np.zeros((N, G), np.int32)
    for g in range(G):
        p = phi_g[g] / (phi_g[g] + Mu[:, g] + 1e-8)
        p = np.clip(p, 1e-6, 1 - 1e-6)
        Y[:, g] = rng.negative_binomial(phi_g[g], p).astype(np.int32)
    return Y
