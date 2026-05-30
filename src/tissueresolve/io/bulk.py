"""
I/O routines for bulk RNA-seq count matrices and sample metadata.

Functions
---------
read_bulk_counts
    Load a genes × samples count matrix from TSV or CSV.
read_bulk_metadata
    Load sample-level metadata (e.g. protocol, tissue type) from TSV or CSV.

Validation policy
-----------------
All functions raise ``ValueError`` with a clear message rather than silently
proceeding with invalid data.  Duplicated gene names are always rejected;
duplicated sample names always emit a warning.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

__all__ = ["read_bulk_counts", "read_bulk_metadata"]

logger = logging.getLogger("tissueresolve.io.bulk")


def read_bulk_counts(
    path: str | Path,
    *,
    sep: Optional[str] = None,
    index_col: int = 0,
    comment: str = "#",
    dtype: type = np.float64,
) -> pd.DataFrame:
    """Load a bulk RNA-seq count matrix from a delimited file.

    The file is expected to be **genes × samples**: rows are genes (first
    column is the gene identifier), subsequent columns are samples.

    Parameters
    ----------
    path:
        Path to a TSV (``.tsv``, ``.txt``) or CSV (``.csv``) file.
    sep:
        Column separator.  Auto-detected from the file extension when ``None``:
        ``,`` for ``.csv``, ``\\t`` otherwise.
    index_col:
        Column index to use as the gene-name index (default ``0``).
    comment:
        Lines beginning with this character are ignored.
    dtype:
        Numeric dtype for the data matrix.

    Returns
    -------
    pd.DataFrame
        Genes × samples DataFrame.  Index name is ``"gene"``.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file has no numeric columns, contains duplicate gene names, or
        contains negative values.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Bulk count file not found: {path}")

    if sep is None:
        sep = "," if path.suffix.lower() == ".csv" else "\t"

    # Detect duplicate column names from the raw header before pandas mangles them
    _check_duplicate_columns(path, sep=sep, comment=comment)

    df = pd.read_csv(path, sep=sep, index_col=index_col, comment=comment)

    if df.empty or df.shape[1] == 0:
        raise ValueError(
            f"Bulk count file '{path}' has no sample columns after reading.  "
            "Check that the gene column is the first column and that the "
            "separator is correct."
        )

    # Cast to numeric
    try:
        df = df.astype(dtype)
    except (ValueError, TypeError) as exc:
        raise ValueError(
            f"Could not convert bulk count matrix to numeric dtype {dtype.__name__}.  "
            f"Non-numeric data found in '{path}'.  Original error: {exc}"
        ) from exc

    # Duplicate gene names
    dup_genes = df.index[df.index.duplicated()].unique().tolist()
    if dup_genes:
        raise ValueError(
            f"Duplicate gene names found in '{path}': {dup_genes[:10]}"
            + ("…" if len(dup_genes) > 10 else "") + ".  "
            "Deduplicate genes (e.g. by summing counts) before loading."
        )

    # Duplicate sample names are detected pre-parse (see _check_duplicate_columns).
    # After pandas reads the file, it has already mangled duplicate names.

    # Negative values
    if (df.values < 0).any():
        raise ValueError(
            f"Negative values found in bulk count matrix '{path}'.  "
            "Count matrices must be non-negative."
        )

    df.index.name = "gene"
    logger.info(
        "Loaded bulk counts: %d genes × %d samples from '%s'.",
        df.shape[0], df.shape[1], path,
    )
    return df


def read_bulk_metadata(
    path: str | Path,
    *,
    sep: Optional[str] = None,
    index_col: int = 0,
    comment: str = "#",
) -> pd.DataFrame:
    """Load sample-level metadata from a delimited file.

    Parameters
    ----------
    path:
        Path to a TSV or CSV file.  First column is the sample identifier.
    sep:
        Column separator.  Auto-detected from extension when ``None``.
    index_col:
        Column to use as the sample index.
    comment:
        Lines beginning with this character are ignored.

    Returns
    -------
    pd.DataFrame
        Samples × metadata columns DataFrame.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If duplicate sample names are found.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Bulk metadata file not found: {path}")

    if sep is None:
        sep = "," if path.suffix.lower() == ".csv" else "\t"

    df = pd.read_csv(path, sep=sep, index_col=index_col, comment=comment)

    dup_samples = df.index[df.index.duplicated()].unique().tolist()
    if dup_samples:
        raise ValueError(
            f"Duplicate sample names in metadata file '{path}': {dup_samples}.  "
            "Each sample must have a unique identifier."
        )

    logger.info(
        "Loaded bulk metadata: %d samples × %d columns from '%s'.",
        df.shape[0], df.shape[1], path,
    )
    return df


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _check_duplicate_columns(
    path: Path,
    *,
    sep: str,
    comment: str,
) -> None:
    """Warn if the file header contains duplicate column names.

    Must be called *before* ``pd.read_csv`` because pandas silently mangles
    duplicate column names (appending ``.1``, ``.2`` etc.) without warning.
    """
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header_line: Optional[str] = None
        for line in fh:
            stripped = line.strip()
            if stripped and not stripped.startswith(comment):
                header_line = stripped
                break
    if header_line is None:
        return
    cols = header_line.split(sep)
    # Drop the index column (first element)
    sample_cols = cols[1:]
    seen: set[str] = set()
    dups: list[str] = []
    for c in sample_cols:
        if c in seen:
            dups.append(c)
        seen.add(c)
    if dups:
        warnings.warn(
            f"Duplicate sample names in '{path}': {dups}.  "
            "Samples will be disambiguated by pandas with a numeric suffix.",
            stacklevel=3,
        )
