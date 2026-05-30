"""
Reference I/O helpers for TissueResolve.

Loads single-cell / single-nucleus RNA-seq AnnData objects from disk so
that :class:`~tissueresolve.reference.build.ReferenceBuilder` can aggregate
them.  Does NOT perform aggregation — that lives in ``reference/build.py``.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from tissueresolve.io.validation import validate_cell_type_column

__all__ = [
    "load_reference_h5ad",
    "load_reference_csv",
]

logger = logging.getLogger("tissueresolve.io.reference")


def load_reference_h5ad(
    path: Union[str, Path],
    *,
    cell_type_col: str = "cell_type",
    backed: bool = False,
) -> "anndata.AnnData":  # type: ignore[name-defined]
    """Load a single-cell / single-nucleus reference from an ``.h5ad`` file.

    Parameters
    ----------
    path:
        Path to the ``.h5ad`` file.
    cell_type_col:
        Column in ``adata.obs`` holding cell-type labels.
    backed:
        Open in memory-mapped (backed) mode.  Use for very large references
        that do not fit in RAM.  Note: backed mode may slow down streaming
        aggregation.

    Returns
    -------
    anndata.AnnData
        Reference object.  ``adata.obs[cell_type_col]`` is guaranteed present.

    Raises
    ------
    FileNotFoundError
        When the file does not exist.
    ImportError
        When ``anndata`` is not installed.
    KeyError
        When ``cell_type_col`` is not in ``adata.obs``.
    """
    try:
        import anndata as ad
        import scipy.sparse as sp
    except ImportError as exc:
        raise ImportError(
            "anndata is required to load .h5ad files.  "
            "Install it with:  pip install anndata"
        ) from exc

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Reference file not found: {path}")

    logger.info("Loading reference from %s (backed=%s)", path, backed)
    mode: Optional[str] = "r" if backed else None
    adata = ad.read_h5ad(path, backed=mode)

    validate_cell_type_column(adata.obs, cell_type_col)

    # Ensure CSR sparse so that batch slicing is fast.
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    else:
        adata.X = adata.X.tocsr()

    n_types = adata.obs[cell_type_col].nunique()
    logger.info(
        "Reference loaded: %d cells × %d genes, %d cell types.",
        adata.n_obs, adata.n_vars, n_types,
    )
    return adata


def load_reference_csv(
    counts_path: Union[str, Path],
    meta_path: Union[str, Path],
    *,
    cell_type_col: str = "cell_type",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load a reference from TSV/CSV files.

    Parameters
    ----------
    counts_path:
        Genes-as-rows, cells-as-columns count matrix (TSV or CSV).
        First column is treated as the gene index.
    meta_path:
        Cell metadata file.  First column is treated as the cell ID index.
        Must contain *cell_type_col*.
    cell_type_col:
        Name of the cell-type annotation column.

    Returns
    -------
    tuple (counts, metadata)
        ``counts``: genes × cells DataFrame.
        ``metadata``: cells × features DataFrame.

    Raises
    ------
    FileNotFoundError
        When either file does not exist.
    KeyError
        When ``cell_type_col`` is not in the metadata.
    """
    counts_path = Path(counts_path)
    meta_path = Path(meta_path)

    if not counts_path.exists():
        raise FileNotFoundError(f"Counts file not found: {counts_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")

    sep_counts = "\t" if ".tsv" in counts_path.name else ","
    sep_meta = "\t" if ".tsv" in meta_path.name else ","

    counts = pd.read_csv(counts_path, sep=sep_counts, index_col=0)
    metadata = pd.read_csv(meta_path, sep=sep_meta, index_col=0)

    validate_cell_type_column(metadata, cell_type_col)

    logger.info(
        "Loaded reference CSV: %d genes × %d cells; metadata: %d cells × %d cols.",
        counts.shape[0], counts.shape[1],
        metadata.shape[0], metadata.shape[1],
    )
    return counts, metadata
