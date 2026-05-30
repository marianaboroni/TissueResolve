"""
Neighbourhood co-occurrence statistics for spatial deconvolution.

Statistical framework
---------------------
For each ordered pair of cell types (k1, k2)::

    C[k1, k2] = (1/N) Σ_{(s,t)∈E} w_{st} · π_{s,k1} · π_{t,k2}
               = (1/N) Π^T (A Π) [k1, k2]

Under the spatial-independence null (spot compositions are spatially
randomly distributed), the expected value is π̄_{k1} · π̄_{k2}.

The null is built by permuting spot indices of Π and recomputing C.
BH-FDR correction uses the shared ``tissueresolve.utils.bh_fdr``.

Limitations (explicit)
----------------------
- The permutation null assumes exchangeable spots, which is false (nearby
  spots are correlated).  P-values are conservative.  Use z-scores for
  ranking and FDR values as soft guides.
- For N < 200 the permutation distribution is noisy.
- The test operates at spot resolution (Visium ~55 μm diameter).

``detect_spatial_niches`` always records ``n_smooth`` in the returned metadata
so that smoothing is never hidden.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np

from tissueresolve.utils import bh_fdr

__all__ = [
    "NeighbourhoodStats",
    "compute_neighbourhood_stats",
    "detect_spatial_niches",
    "distance_decay_cooccurrence",
]

logger = logging.getLogger("tissueresolve.spatial.neighbourhood")


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class NeighbourhoodStats:
    """Results of the neighbourhood co-occurrence analysis.

    Attributes
    ----------
    observed:
        Observed weighted co-occurrence matrix ``(K, K)`` float32.
    expected:
        Expected co-occurrence (mean over permutations) ``(K, K)`` float32.
    zscore:
        (observed − mean_null) / std_null ``(K, K)`` float32.
        Positive = enriched; negative = depleted.
    pvalue_two_sided:
        Empirical two-sided p-values ``(K, K)`` float32.
    pvalue_fdr:
        BH-FDR-corrected q-values (over all K² tests) ``(K, K)`` float32.
    cell_types:
        Ordered cell-type names.
    n_edges:
        Number of edges used (= nnz of adjacency matrix).
    n_permutations:
        Permutations used for the null.
    """

    observed: np.ndarray
    expected: np.ndarray
    zscore: np.ndarray
    pvalue_two_sided: np.ndarray
    pvalue_fdr: np.ndarray
    cell_types: list[str]
    n_edges: int
    n_permutations: int

    @property
    def K(self) -> int:
        return len(self.cell_types)

    def top_enriched_pairs(self, n: int = 10) -> list[tuple[str, str, float, float]]:
        """Return top-n enriched (type_a, type_b, z, fdr_q) pairs."""
        results = [
            (self.cell_types[i], self.cell_types[j],
             float(self.zscore[i, j]), float(self.pvalue_fdr[i, j]))
            for i in range(self.K) for j in range(self.K) if i != j
        ]
        results.sort(key=lambda x: x[2], reverse=True)
        return results[:n]

    def top_depleted_pairs(self, n: int = 10) -> list[tuple[str, str, float, float]]:
        """Return top-n depleted (type_a, type_b, z, fdr_q) pairs."""
        results = [
            (self.cell_types[i], self.cell_types[j],
             float(self.zscore[i, j]), float(self.pvalue_fdr[i, j]))
            for i in range(self.K) for j in range(self.K) if i != j
        ]
        results.sort(key=lambda x: x[2])
        return results[:n]


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------


def compute_neighbourhood_stats(
    Pi: np.ndarray,
    graph: "SpatialGraph",  # type: ignore[name-defined]
    cell_types: list[str],
    *,
    n_permutations: int = 999,
    seed: int = 0,
    min_spots: int = 50,
) -> NeighbourhoodStats:
    """Compute pairwise neighbourhood co-occurrence statistics.

    Parameters
    ----------
    Pi:
        Proportion matrix ``(N, K)``.  Each row sums to 1.
    graph:
        :class:`~tissueresolve.spatial.graph.SpatialGraph`.
    cell_types:
        Names of the K cell types.
    n_permutations:
        Spatial permutations for the null (default 999).
    seed:
        Random seed.
    min_spots:
        Warn if fewer spots than this.

    Returns
    -------
    NeighbourhoodStats
    """
    N, K = Pi.shape
    rng = np.random.default_rng(seed)

    if N < min_spots:
        logger.warning(
            "Only %d spots — neighbourhood statistics will be noisy (min=%d).",
            N, min_spots,
        )

    A = graph.A
    n_edges = int(A.nnz)

    # Observed co-occurrence
    Pi_neigh = A @ Pi                       # (N, K)
    C_obs = (Pi.T @ Pi_neigh) / N           # (K, K)

    # Null by permutation
    C_null = np.zeros((n_permutations, K, K), dtype=np.float64)
    for b in range(n_permutations):
        perm = rng.permutation(N)
        Pi_perm = Pi[perm]
        C_null[b] = (Pi_perm.T @ (A @ Pi_perm)) / N

    null_mean = C_null.mean(axis=0)
    null_std = C_null.std(axis=0) + 1e-10
    zscore = (C_obs - null_mean) / null_std

    obs_dev = np.abs(C_obs - null_mean)
    null_dev = np.abs(C_null - null_mean[np.newaxis])
    pvalue = (null_dev >= obs_dev[np.newaxis]).mean(axis=0)
    pvalue = np.clip(pvalue, 1.0 / (n_permutations + 1), 1.0)

    # BH-FDR using shared tissueresolve.utils.bh_fdr (not a local copy)
    pvalue_fdr = bh_fdr(pvalue.ravel()).reshape(K, K)

    stats = NeighbourhoodStats(
        observed=C_obs.astype(np.float32),
        expected=null_mean.astype(np.float32),
        zscore=zscore.astype(np.float32),
        pvalue_two_sided=pvalue.astype(np.float32),
        pvalue_fdr=pvalue_fdr.astype(np.float32),
        cell_types=list(cell_types),
        n_edges=n_edges,
        n_permutations=n_permutations,
    )

    n_sig = int((pvalue_fdr < 0.05).sum())
    logger.info(
        "Neighbourhood stats: N=%d, %d edges, %d perm, %d sig pairs (FDR<0.05).",
        N, n_edges, n_permutations, n_sig,
    )
    return stats


# ---------------------------------------------------------------------------
# Spatial niche detection
# ---------------------------------------------------------------------------


def detect_spatial_niches(
    Pi: np.ndarray,
    graph: "SpatialGraph",  # type: ignore[name-defined]
    cell_types: list[str],
    *,
    n_niches: int = 5,
    n_smooth: int = 2,
    seed: int = 0,
) -> tuple[np.ndarray, dict]:
    """Cluster spots into spatial niches based on smoothed composition.

    Applies ``n_smooth`` rounds of Laplacian smoothing (always recorded in
    returned metadata), then clusters using K-means.

    Parameters
    ----------
    Pi:
        Proportion matrix ``(N, K)``.
    graph:
        Spatial graph.
    cell_types:
        Cell-type names.
    n_niches:
        Number of niche clusters.
    n_smooth:
        Rounds of neighbour smoothing before clustering.  **Always
        recorded in metadata** so smoothing is never hidden.
    seed:
        Random seed.

    Returns
    -------
    tuple (labels, metadata)
        ``labels``: integer niche labels ``(N,)`` int32.
        ``metadata``: dict including ``n_smooth``, ``n_niches``, ``seed``.
    """
    from sklearn.cluster import KMeans  # lazy import

    A = graph.A
    Pi_smooth = Pi.copy().astype(np.float64)
    for _ in range(n_smooth):
        Pi_smooth = 0.5 * Pi_smooth + 0.5 * (A @ Pi_smooth)
        Pi_smooth /= Pi_smooth.sum(axis=1, keepdims=True) + 1e-10

    km = KMeans(n_clusters=n_niches, random_state=seed, n_init=10)
    labels = km.fit_predict(Pi_smooth).astype(np.int32)

    for niche_id in range(n_niches):
        mask = labels == niche_id
        if not mask.any():
            continue
        mean_comp = Pi[mask].mean(axis=0)
        top_k = mean_comp.argsort()[::-1][:3]
        logger.info(
            "Niche %d (n=%d): %s",
            niche_id, int(mask.sum()),
            ", ".join(f"{cell_types[k]}={mean_comp[k]:.2f}" for k in top_k),
        )

    metadata: dict = {
        "n_smooth": n_smooth,      # always recorded — smoothing never hidden
        "n_niches": n_niches,
        "seed": seed,
        "n_spots": int(len(labels)),
    }
    logger.info(
        "detect_spatial_niches: %d niches, n_smooth=%d.", n_niches, n_smooth
    )
    return labels, metadata


# ---------------------------------------------------------------------------
# Distance-decay co-occurrence
# ---------------------------------------------------------------------------


def distance_decay_cooccurrence(
    Pi: np.ndarray,
    array_row: np.ndarray,
    array_col: np.ndarray,
    k1: int,
    k2: int,
    *,
    max_distance: float = 10.0,
    n_bins: int = 10,
    n_sample_anchors: int = 200,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Co-occurrence π_{s,k1} · π_{t,k2} as a function of spatial distance.

    Useful for assessing whether co-occurrence decays with distance (genuine
    biological interaction) or is flat (global correlation, not local).

    Parameters
    ----------
    Pi:
        Proportion matrix ``(N, K)``.
    array_row / array_col:
        Integer Visium coordinates ``(N,)``.
    k1, k2:
        Cell-type indices.
    max_distance:
        Maximum grid distance.
    n_bins:
        Distance bins.
    n_sample_anchors:
        Max anchor spots sampled (for efficiency).
    seed:
        Random seed.

    Returns
    -------
    tuple (bin_centres, mean_cooccurrence)
        Both shape ``(n_bins,)`` float32.
    """
    N = len(array_row)
    rng = np.random.default_rng(seed)
    n_sample = min(n_sample_anchors, N)
    anchor_idx = rng.choice(N, n_sample, replace=False)

    bin_edges = np.linspace(0, max_distance, n_bins + 1)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    co_sums = np.zeros(n_bins, dtype=np.float64)
    co_counts = np.zeros(n_bins, dtype=np.int64)

    coords = np.stack([array_row, array_col], axis=1).astype(np.float32)
    for a in anchor_idx:
        dists = np.linalg.norm(coords - coords[a], axis=1)
        co = Pi[a, k1] * Pi[:, k2]
        for b in range(n_bins):
            mask = (dists >= bin_edges[b]) & (dists < bin_edges[b + 1])
            if mask.any():
                co_sums[b] += co[mask].sum()
                co_counts[b] += int(mask.sum())

    mean_co = np.where(co_counts > 0, co_sums / co_counts, np.nan)
    return bin_centres, mean_co.astype(np.float32)
