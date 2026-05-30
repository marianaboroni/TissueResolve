"""
I/O routines for spatial (Visium) deconvolution in TissueResolve.

Functions
---------
load_visium
    Unified loader: accepts a SpaceRanger output directory or an ``.h5ad``
    file.  Always attaches ``array_row``, ``array_col`` and validates them.
load_visium_h5ad
    Load directly from an ``.h5ad`` file (without scanpy).
make_synthetic_visium
    Create a minimal AnnData-compatible dict for testing without real data.

Coordinate convention
---------------------
``array_row`` and ``array_col`` in ``adata.obs`` are integer Visium array
coordinates, **not** pixel coordinates.  The spatial graph builder
(:func:`~tissueresolve.spatial.graph.build_hex_graph`) requires them.
Pixel coordinates are stored in ``adata.obsm["spatial"]`` when available.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

__all__ = [
    "load_visium",
    "load_visium_h5ad",
    "make_synthetic_visium",
]

logger = logging.getLogger("tissueresolve.io.spatial")


# ---------------------------------------------------------------------------
# Public loaders
# ---------------------------------------------------------------------------


def load_visium(
    path: Union[str, Path],
    *,
    min_counts: int = 100,
    min_genes: int = 200,
    backed: bool = False,
    load_images: bool = False,
) -> "anndata.AnnData":  # type: ignore[name-defined]
    """Load 10x Visium data into an :class:`anndata.AnnData`.

    Accepts:

    - A **SpaceRanger output directory** (contains
      ``filtered_feature_bc_matrix.h5`` and ``spatial/``).
    - A pre-processed ``.h5ad`` file.

    Parameters
    ----------
    path:
        SpaceRanger directory or ``.h5ad`` file.
    min_counts:
        Minimum total UMI per spot.  Spots below this threshold are dropped.
    min_genes:
        Minimum detected genes per spot.
    backed:
        Open in memory-mapped mode (only for ``.h5ad``).
    load_images:
        Whether to load the tissue image (SpaceRanger directory only).

    Returns
    -------
    anndata.AnnData
        With ``adata.obs["array_row"]``, ``adata.obs["array_col"]``
        (integer), ``adata.obs["total_counts"]`` (float32), and
        ``adata.obsm["spatial"]`` (pixel coords when available).

    Raises
    ------
    ValueError
        If path is neither a directory nor an ``.h5ad`` file, or if
        array coordinates are missing.
    FileNotFoundError
        If *path* does not exist.
    """
    import anndata as ad
    import scipy.sparse as sp

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Visium data not found: {path}")

    if path.suffix in {".h5ad", ".h5"}:
        adata = load_visium_h5ad(path, backed=backed)
    elif path.is_dir():
        adata = _load_from_spaceranger(path, load_images=load_images)
    else:
        raise ValueError(
            f"Cannot load Visium data from {path!r}.  "
            "Expected a SpaceRanger output directory or a .h5ad file."
        )

    # Quality filters (optional — set min_counts=0 to disable)
    if min_counts > 0 or min_genes > 0:
        try:
            import scanpy as sc
            n_before = adata.n_obs
            if min_counts > 0:
                sc.pp.filter_cells(adata, min_counts=min_counts)
            if min_genes > 0:
                sc.pp.filter_cells(adata, min_genes=min_genes)
            n_after = adata.n_obs
            if n_before != n_after:
                logger.info(
                    "Filtered %d / %d spots → %d remain.",
                    n_before - n_after, n_before, n_after,
                )
        except ImportError:
            logger.warning(
                "scanpy not installed — skipping cell-level filtering "
                "(min_counts=%d, min_genes=%d).", min_counts, min_genes,
            )

    # Library sizes
    if "total_counts" not in adata.obs.columns:
        if sp.issparse(adata.X):
            adata.obs["total_counts"] = np.asarray(
                adata.X.sum(axis=1), dtype=np.float32
            ).ravel()
        else:
            adata.obs["total_counts"] = adata.X.sum(axis=1).astype(np.float32)

    # Validate array coordinates
    _validate_array_coords(adata)

    logger.info(
        "Loaded Visium: %d spots × %d genes "
        "(median lib size %.0f, array coords attached).",
        adata.n_obs, adata.n_vars,
        float(np.median(adata.obs["total_counts"])),
    )
    return adata


def load_visium_h5ad(
    path: Union[str, Path],
    *,
    backed: bool = False,
) -> "anndata.AnnData":  # type: ignore[name-defined]
    """Load a Visium AnnData directly from an ``.h5ad`` file.

    Does not require scanpy.  Uses anndata only.

    Parameters
    ----------
    path:
        Path to the ``.h5ad`` file.
    backed:
        Open in backed (memory-mapped) mode.

    Returns
    -------
    anndata.AnnData

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    import anndata as ad
    import scipy.sparse as sp

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f".h5ad file not found: {path}")

    logger.info("Loading Visium from .h5ad: %s (backed=%s)", path, backed)
    mode = "r" if backed else None
    adata = ad.read_h5ad(path, backed=mode)

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    else:
        adata.X = adata.X.tocsr()

    return adata


def make_synthetic_visium(
    n_spots: int = 40,
    n_genes: int = 80,
    n_cell_types: int = 3,
    seed: int = 42,
) -> "anndata.AnnData":  # type: ignore[name-defined]
    """Create a synthetic Visium AnnData for testing.

    Generates a rectangular grid of spots with:
    - Random NB-distributed counts from known proportions.
    - ``array_row``, ``array_col``, ``total_counts`` in obs.

    Returns
    -------
    anndata.AnnData
        With ``adata.uns["true_proportions"]`` (N × K) for benchmark use.
    """
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)

    # Build a small rectangular hex grid
    rows_per_col = max(2, int(np.ceil(np.sqrt(n_spots))))
    rows_: list[int] = []
    cols_: list[int] = []
    for r in range(rows_per_col):
        for c in range(rows_per_col):
            rows_.append(r)
            cols_.append(c * 2 + (r % 2))  # offset columns for hex layout
    array_row = np.array(rows_[:n_spots], dtype=np.int32)
    array_col = np.array(cols_[:n_spots], dtype=np.int32)
    N = len(array_row)

    # Reference profiles (block-diagonal for separability)
    block = n_genes // n_cell_types
    R_lin = np.zeros((n_cell_types, n_genes), dtype=np.float32)
    for k in range(n_cell_types):
        R_lin[k, k * block:(k + 1) * block] = rng.uniform(1e-4, 3e-4, block)
    R_lin += rng.uniform(0, 1e-5, R_lin.shape)
    R_lin /= R_lin.sum(axis=1, keepdims=True)

    # True proportions
    true_props = rng.dirichlet(np.ones(n_cell_types), size=N).astype(np.float32)

    # Synthetic counts
    lib_sizes = rng.integers(1000, 5000, size=N).astype(np.float32)
    Mu = lib_sizes[:, None] * (true_props @ R_lin)  # (N, G)
    phi = 5.0
    p_nb = phi / (phi + Mu + 1e-8)
    p_nb = np.clip(p_nb, 1e-6, 1 - 1e-6)
    counts = rng.negative_binomial(
        phi * np.ones_like(Mu), p_nb
    ).astype(np.float32)

    obs = pd.DataFrame({
        "array_row": array_row,
        "array_col": array_col,
        "total_counts": counts.sum(axis=1),
    }, index=[f"spot_{i}" for i in range(N)])

    var = pd.DataFrame(index=[f"GENE_{i:04d}" for i in range(n_genes)])

    adata = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=obs,
        var=var,
    )
    adata.uns["true_proportions"] = true_props
    adata.uns["n_cell_types"] = n_cell_types
    adata.uns["R_lin_true"] = R_lin
    return adata


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _load_from_spaceranger(path: Path, load_images: bool) -> "anndata.AnnData":
    """Load from a SpaceRanger directory using scanpy."""
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError(
            "scanpy is required to load from SpaceRanger directories.  "
            "Install it with: pip install scanpy"
        ) from exc

    logger.info("Loading Visium from SpaceRanger directory: %s", path)
    adata = sc.read_visium(
        path,
        count_file="filtered_feature_bc_matrix.h5",
        load_images=load_images,
    )
    import scipy.sparse as sp
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    else:
        adata.X = adata.X.tocsr()

    if "array_row" not in adata.obs.columns:
        _attach_array_coords_from_csv(adata, path)
    return adata


def _attach_array_coords_from_csv(adata: "anndata.AnnData", spaceranger_path: Path) -> None:
    """Parse tissue_positions CSV and attach array_row / array_col."""
    spatial_dir = spaceranger_path / "spatial"
    candidates = [
        spatial_dir / "tissue_positions.csv",
        spatial_dir / "tissue_positions_list.csv",
    ]
    pos_file: Optional[Path] = None
    for c in candidates:
        if c.exists():
            pos_file = c
            break
    if pos_file is None:
        raise FileNotFoundError(
            f"tissue_positions CSV not found in {spatial_dir}"
        )

    with open(pos_file) as fh:
        first = fh.readline().strip()
    has_header = first.startswith("barcode")

    cols = ["barcode", "in_tissue", "array_row", "array_col",
            "pxl_col_in_fullres", "pxl_row_in_fullres"]
    pos = pd.read_csv(
        pos_file,
        header=0 if has_header else None,
        names=None if has_header else cols,
        index_col=0,
    )
    pos = pos[pos["in_tissue"] == 1]
    shared = adata.obs_names.intersection(pos.index)
    adata.obs["array_row"] = pos.loc[shared, "array_row"].astype(np.int32)
    adata.obs["array_col"] = pos.loc[shared, "array_col"].astype(np.int32)


def _validate_array_coords(adata: "anndata.AnnData") -> None:
    """Raise ValueError if array_row / array_col are missing or invalid."""
    for col in ("array_row", "array_col"):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column {col!r} missing from adata.obs.  "
                "Ensure the AnnData was loaded with load_visium() or "
                "contains valid Visium coordinates."
            )
    adata.obs["array_row"] = adata.obs["array_row"].astype(np.int32)
    adata.obs["array_col"] = adata.obs["array_col"].astype(np.int32)

    if (adata.obs["array_row"] < 0).any() or (adata.obs["array_col"] < 0).any():
        raise ValueError(
            "Negative values found in array_row / array_col.  "
            "Visium array coordinates must be non-negative integers."
        )
