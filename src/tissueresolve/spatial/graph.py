"""
Spatial graph construction for Visium deconvolution.

Builds the hexagonal neighbourhood graph from Visium array coordinates
(``array_row``, ``array_col``) using a hash-map lookup.

10x Visium array layout
-----------------------
Each interior spot has exactly **6** neighbours at offsets:
``(0, ±2)`` (same row) and ``(±1, ±1)`` (adjacent rows).
This offset pattern is universal across all standard Visium products.

The duplicated ``_HEX_OFFSETS`` redefinition present in the SpatCAR v1
``build_expression_weighted_graph`` has been removed; a single module-level
constant is shared across all builders.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.sparse as sp

__all__ = [
    "SpatialGraph",
    "build_hex_graph",
    "build_hex_graph_from_arrays",
    "build_expression_weighted_graph",
    "compute_spatial_batches",
]

logger = logging.getLogger("tissueresolve.spatial.graph")

# Hexagonal offsets (Δrow, Δcol) — universal across Visium products.
_HEX_OFFSETS: list[tuple[int, int]] = [
    (0, -2), (0,  2),   # same row, ±2 cols
    (-1, -1), (-1, 1),  # row above, ±1 col
    (1,  -1), (1,  1),  # row below, ±1 col
]


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------


@dataclass
class SpatialGraph:
    """Hexagonal adjacency graph for a single Visium section.

    Attributes
    ----------
    A:
        Row-normalised weighted adjacency matrix ``(N, N)``, CSR float32.
        Each non-zero ``A[i, j]`` is the normalised weight from spot *i* to *j*.
        Row sums equal 1 for all spots with at least one neighbour.
    L:
        Graph Laplacian ``I − A``, same shape, CSR float32.
        Used for the spatial regularisation term.
    degree:
        Integer degree (number of neighbours) per spot, shape ``(N,)``.
    n_spots:
        Total number of spots.
    spot_ids:
        Optional barcode/label array, length N.
    batches:
        Spatial tile partition produced by :func:`compute_spatial_batches`.
        Populated during construction.
    """

    A: sp.csr_matrix
    L: sp.csr_matrix
    degree: np.ndarray
    n_spots: int
    spot_ids: Optional[np.ndarray] = None
    batches: list[np.ndarray] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.A.shape != (self.n_spots, self.n_spots):
            raise ValueError(
                f"A shape {self.A.shape} does not match n_spots={self.n_spots}."
            )
        if self.L.shape != (self.n_spots, self.n_spots):
            raise ValueError(
                f"L shape {self.L.shape} does not match n_spots={self.n_spots}."
            )
        if len(self.degree) != self.n_spots:
            raise ValueError(
                f"len(degree)={len(self.degree)} does not match n_spots={self.n_spots}."
            )


# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------


def build_hex_graph(
    adata: "anndata.AnnData",  # type: ignore[name-defined]
    *,
    batch_size: int = 500,
) -> SpatialGraph:
    """Build the hexagonal neighbourhood graph from a Visium AnnData.

    Reads ``adata.obs["array_row"]`` and ``adata.obs["array_col"]``.

    Parameters
    ----------
    adata:
        Visium AnnData.  Must have ``array_row`` and ``array_col`` in obs.
    batch_size:
        Approximate spots per spatial tile for :func:`compute_spatial_batches`.

    Returns
    -------
    SpatialGraph
    """
    for col in ("array_row", "array_col"):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column {col!r} not found in adata.obs.  "
                "Load the AnnData with load_visium() or attach array coordinates."
            )

    array_row = adata.obs["array_row"].values.astype(np.int32)
    array_col = adata.obs["array_col"].values.astype(np.int32)
    spot_ids = adata.obs_names.values

    _validate_coordinates(array_row, array_col)

    graph = _build_graph_from_arrays(array_row, array_col, spot_ids)
    graph.batches = compute_spatial_batches(array_row, array_col, batch_size=batch_size)

    logger.info(
        "Spatial graph (from AnnData): %d spots, mean degree %.2f, "
        "boundary spots: %d, %d spatial batches.",
        graph.n_spots,
        float(graph.degree.mean()),
        int((graph.degree < 6).sum()),
        len(graph.batches),
    )
    return graph


def build_hex_graph_from_arrays(
    array_row: np.ndarray,
    array_col: np.ndarray,
    spot_ids: Optional[np.ndarray] = None,
    *,
    batch_size: int = 500,
) -> SpatialGraph:
    """Build a :class:`SpatialGraph` directly from coordinate arrays.

    Useful for synthetic data and benchmarks where no AnnData is available.

    Parameters
    ----------
    array_row:
        Integer array-row indices, shape ``(N,)``.
    array_col:
        Integer array-col indices, shape ``(N,)``.
    spot_ids:
        Optional string labels, shape ``(N,)``.
    batch_size:
        Spatial tile size for batching.

    Returns
    -------
    SpatialGraph

    Raises
    ------
    ValueError
        If *array_row* and *array_col* have different lengths or contain
        non-integer-like values, or if no spots are provided.
    """
    array_row = np.asarray(array_row, dtype=np.int32)
    array_col = np.asarray(array_col, dtype=np.int32)
    _validate_coordinates(array_row, array_col)

    graph = _build_graph_from_arrays(array_row, array_col, spot_ids)
    graph.batches = compute_spatial_batches(array_row, array_col, batch_size=batch_size)

    logger.info(
        "Spatial graph: %d spots, mean degree %.2f, boundary spots: %d.",
        graph.n_spots,
        float(graph.degree.mean()),
        int((graph.degree < 6).sum()),
    )
    return graph


def build_expression_weighted_graph(
    Y_marker: np.ndarray,
    array_row: np.ndarray,
    array_col: np.ndarray,
    lib_sizes: np.ndarray,
    spot_ids: Optional[np.ndarray] = None,
    *,
    gamma: float = 1.0,
    batch_size: int = 500,
) -> SpatialGraph:
    """Build a spatial graph with expression-similarity edge weights.

    Edges across tissue boundaries (spots with very different marker-gene
    expression) receive lower weights, reducing over-smoothing at sharp
    tissue transitions.

    Edge weight from spot *s* to neighbour *t*::

        w_{st} = exp(−γ · ‖(Y_s / lib_s) − (Y_t / lib_t)‖₁)

    Parameters
    ----------
    Y_marker:
        Dense marker-gene count matrix ``(N, G_m)`` float32.
    array_row / array_col:
        Integer Visium coordinates, shape ``(N,)``.
    lib_sizes:
        Per-spot library sizes, shape ``(N,)``.
    spot_ids:
        Optional barcode labels.
    gamma:
        Decay rate.  Higher → sharper boundary preservation.
    batch_size:
        Spatial tile size.

    Returns
    -------
    SpatialGraph
    """
    array_row = np.asarray(array_row, dtype=np.int32)
    array_col = np.asarray(array_col, dtype=np.int32)
    _validate_coordinates(array_row, array_col)

    N = len(array_row)
    lib_safe = np.maximum(lib_sizes, 1.0).astype(np.float32)
    Y_norm = Y_marker.astype(np.float32) / lib_safe[:, None]

    coord_to_idx: dict[tuple[int, int], int] = {
        (int(r), int(c)): i
        for i, (r, c) in enumerate(zip(array_row, array_col))
    }

    rows_coo, cols_coo, data_coo = [], [], []
    degree_arr = np.zeros(N, dtype=np.int32)

    for i, (r, c) in enumerate(zip(array_row.tolist(), array_col.tolist())):
        for dr, dc in _HEX_OFFSETS:
            j = coord_to_idx.get((r + dr, c + dc))
            if j is not None:
                expr_dist = float(np.abs(Y_norm[i] - Y_norm[j]).sum())
                w = float(np.exp(-gamma * expr_dist))
                rows_coo.append(i)
                cols_coo.append(j)
                data_coo.append(w)
                degree_arr[i] += 1

    if not rows_coo:
        warnings.warn(
            "build_expression_weighted_graph: no edges found.  "
            "Check array coordinates.",
            stacklevel=2,
        )

    A_raw = sp.coo_matrix(
        (np.array(data_coo, dtype=np.float32),
         (np.array(rows_coo, np.int32), np.array(cols_coo, np.int32))),
        shape=(N, N), dtype=np.float32,
    ).tocsr()

    row_sums = np.asarray(A_raw.sum(axis=1), dtype=np.float32).ravel()
    inv = np.where(row_sums > 0, 1.0 / row_sums, 0.0)
    A_norm = sp.diags(inv.astype(np.float32), format="csr") @ A_raw
    L = sp.eye(N, format="csr", dtype=np.float32) - A_norm

    graph = SpatialGraph(
        A=A_norm.tocsr(),
        L=L.tocsr(),
        degree=degree_arr,
        n_spots=N,
        spot_ids=spot_ids,
    )
    graph.batches = compute_spatial_batches(array_row, array_col, batch_size=batch_size)

    logger.info(
        "Expression-weighted graph: %d spots, mean degree %.2f, "
        "mean edge weight %.4f.",
        N, float(degree_arr.mean()),
        float(np.asarray(A_raw.data).mean()) if len(A_raw.data) else 0.0,
    )
    return graph


# ---------------------------------------------------------------------------
# Spatial tiling
# ---------------------------------------------------------------------------


def compute_spatial_batches(
    array_row: np.ndarray,
    array_col: np.ndarray,
    *,
    batch_size: int = 500,
) -> list[np.ndarray]:
    """Partition spots into spatial tiles of approximately *batch_size* spots.

    Tiles are rectangular blocks in ``(array_row, array_col)`` space.
    Spots within the same tile are spatially adjacent, improving cache
    locality when accessing the proportion matrix.

    Parameters
    ----------
    array_row / array_col:
        Integer coordinates, shape ``(N,)``.
    batch_size:
        Target spots per tile.

    Returns
    -------
    list[np.ndarray]
        Integer index arrays, one per tile.  Union covers all N spots
        exactly once.
    """
    N = len(array_row)
    if N == 0:
        return []

    row_min, row_max = int(array_row.min()), int(array_row.max())
    col_min, col_max = int(array_col.min()), int(array_col.max())

    n_tiles_target = max(1, round(N / batch_size))
    n_row_tiles = max(1, round(n_tiles_target ** 0.5))
    n_col_tiles = max(1, round(n_tiles_target / n_row_tiles))

    row_edges = np.linspace(row_min - 0.5, row_max + 0.5, n_row_tiles + 1)
    col_edges = np.linspace(col_min - 0.5, col_max + 0.5, n_col_tiles + 1)

    batches: list[np.ndarray] = []
    for r_lo, r_hi in zip(row_edges[:-1], row_edges[1:]):
        for c_lo, c_hi in zip(col_edges[:-1], col_edges[1:]):
            mask = (
                (array_row >= r_lo) & (array_row < r_hi) &
                (array_col >= c_lo) & (array_col < c_hi)
            )
            idx = np.where(mask)[0].astype(np.int32)
            if len(idx) > 0:
                batches.append(idx)

    total = sum(len(b) for b in batches)
    if total != N:
        raise RuntimeError(
            f"Batch partition incomplete: {total} != {N} spots.  "
            "This is a bug — please report it."
        )
    return batches


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_coordinates(
    array_row: np.ndarray,
    array_col: np.ndarray,
) -> None:
    """Raise ValueError on invalid or inconsistent coordinates."""
    if len(array_row) == 0:
        raise ValueError("array_row is empty — no spots provided.")
    if len(array_row) != len(array_col):
        raise ValueError(
            f"array_row (len={len(array_row)}) and array_col (len={len(array_col)}) "
            "must have the same length."
        )
    if np.any(array_row < 0) or np.any(array_col < 0):
        raise ValueError(
            "Negative coordinates found in array_row / array_col.  "
            "Visium array coordinates are non-negative integers."
        )


def _build_graph_from_arrays(
    array_row: np.ndarray,
    array_col: np.ndarray,
    spot_ids: Optional[np.ndarray],
) -> SpatialGraph:
    """Core graph construction shared by all builders."""
    N = len(array_row)

    coord_to_idx: dict[tuple[int, int], int] = {
        (int(r), int(c)): i
        for i, (r, c) in enumerate(zip(array_row, array_col))
    }

    rows_coo: list[int] = []
    cols_coo: list[int] = []
    degree_arr = np.zeros(N, dtype=np.int32)

    for i, (r, c) in enumerate(zip(array_row.tolist(), array_col.tolist())):
        for dr, dc in _HEX_OFFSETS:
            j = coord_to_idx.get((r + dr, c + dc))
            if j is not None:
                rows_coo.append(i)
                cols_coo.append(j)
                degree_arr[i] += 1

    if not rows_coo:
        warnings.warn(
            "No edges found in spatial graph.  "
            "Check that array_row / array_col are valid Visium coordinates.",
            stacklevel=3,
        )

    data = np.ones(len(rows_coo), dtype=np.float32)
    A_raw = sp.coo_matrix(
        (data,
         (np.array(rows_coo, np.int32), np.array(cols_coo, np.int32))),
        shape=(N, N), dtype=np.float32,
    ).tocsr()

    inv_degree = np.where(degree_arr > 0, 1.0 / degree_arr.astype(np.float32), 0.0)
    D_inv = sp.diags(inv_degree.astype(np.float32), format="csr")
    A_norm = D_inv @ A_raw

    I = sp.eye(N, format="csr", dtype=np.float32)
    L = I - A_norm

    return SpatialGraph(
        A=A_norm.tocsr(),
        L=L.tocsr(),
        degree=degree_arr,
        n_spots=N,
        spot_ids=spot_ids,
    )
