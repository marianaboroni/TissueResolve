"""
Input validation utilities for TissueResolve.

All functions raise explicit errors with clear messages.  No silent failure.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

__all__ = [
    "check_gene_overlap",
    "validate_no_duplicate_genes",
    "validate_matrix_shape",
    "validate_cell_type_column",
    "validate_gene_names_match",
]

logger = logging.getLogger("tissueresolve.io.validation")


def check_gene_overlap(
    query_genes: list[str],
    reference_genes: list[str],
    *,
    min_overlap: int = 50,
    context: str = "",
) -> list[str]:
    """Return genes present in both lists; warn when overlap is small.

    Parameters
    ----------
    query_genes:
        Genes from the query dataset (e.g. bulk counts or Visium).
    reference_genes:
        Genes in the reference.
    min_overlap:
        Emit a warning when fewer shared genes are found.
    context:
        Optional label for the warning message (e.g. ``"Visium vs reference"``).

    Returns
    -------
    list[str]
        Sorted intersection.
    """
    shared = sorted(set(query_genes) & set(reference_genes))
    q, r = len(query_genes), len(reference_genes)
    pct = 100.0 * len(shared) / max(r, 1)
    label = f" ({context})" if context else ""
    logger.info(
        "Gene overlap%s: %d / %d reference genes in query (%.1f%%).",
        label, len(shared), r, pct,
    )
    if len(shared) < min_overlap:
        warnings.warn(
            f"Gene overlap{label}: only {len(shared)} shared genes "
            f"(query {q}, reference {r}).  "
            "Check that gene naming conventions match "
            "(e.g. Ensembl IDs vs gene symbols, capitalisation).",
            stacklevel=2,
        )
    return shared


def validate_no_duplicate_genes(
    gene_names: list[str] | pd.Index,
    source: str = "gene list",
) -> None:
    """Raise ``ValueError`` if *gene_names* contains duplicates.

    Parameters
    ----------
    gene_names:
        Sequence of gene identifiers.
    source:
        Label for the error message (e.g. ``"counts matrix index"``).
    """
    if isinstance(gene_names, pd.Index):
        dupes = gene_names[gene_names.duplicated()].tolist()
    else:
        seen: set = set()
        dupes = []
        for g in gene_names:
            if g in seen and g not in dupes:
                dupes.append(g)
            seen.add(g)

    if dupes:
        raise ValueError(
            f"Duplicate gene names in {source}: {dupes[:5]}"
            + ("…" if len(dupes) > 5 else "")
            + f"  ({len(dupes)} total duplicates).  "
            "Resolve duplicates before building a reference."
        )


def validate_matrix_shape(
    matrix: np.ndarray,
    expected_shape: tuple,
    name: str,
) -> None:
    """Raise ``ValueError`` when *matrix* does not have *expected_shape*.

    Parameters
    ----------
    matrix:
        Array to check.
    expected_shape:
        Expected shape tuple.  Use ``-1`` as a wildcard for any dimension.
    name:
        Human-readable matrix name for the error message.
    """
    if matrix.ndim != len(expected_shape):
        raise ValueError(
            f"{name}: expected {len(expected_shape)}-D array, "
            f"got {matrix.ndim}-D array with shape {matrix.shape}."
        )
    for i, (actual, expected) in enumerate(zip(matrix.shape, expected_shape)):
        if expected != -1 and actual != expected:
            raise ValueError(
                f"{name}: dimension {i} expected {expected}, got {actual} "
                f"(full shape {matrix.shape}, expected {expected_shape})."
            )


def validate_cell_type_column(
    metadata: pd.DataFrame,
    celltype_col: str,
) -> None:
    """Raise ``KeyError`` if *celltype_col* is not in *metadata*.

    Parameters
    ----------
    metadata:
        Cell metadata DataFrame.
    celltype_col:
        Name of the cell-type annotation column.
    """
    if celltype_col not in metadata.columns:
        raise KeyError(
            f"Cell-type column '{celltype_col}' not found in metadata.  "
            f"Available columns: {list(metadata.columns)}.  "
            "Set ReferenceConfig.celltype_col to the correct column name."
        )


def validate_gene_names_match(
    gene_names_a: list[str],
    gene_names_b: list[str],
    label_a: str = "A",
    label_b: str = "B",
) -> None:
    """Raise ``ValueError`` if two gene lists differ.

    Parameters
    ----------
    gene_names_a, gene_names_b:
        Gene identifier lists to compare.
    label_a, label_b:
        Labels used in the error message.
    """
    if gene_names_a != gene_names_b:
        set_a, set_b = set(gene_names_a), set(gene_names_b)
        only_a = sorted(set_a - set_b)[:5]
        only_b = sorted(set_b - set_a)[:5]
        raise ValueError(
            f"Gene names in {label_a} and {label_b} do not match.  "
            f"Only in {label_a}: {only_a}…  "
            f"Only in {label_b}: {only_b}…  "
            "Call ref.subset_genes() to align before proceeding."
        )


def validate_counts_matrix(
    counts: pd.DataFrame,
    *,
    allow_float: bool = True,
) -> None:
    """Validate a genes × cells count matrix.

    Checks:
    - Gene index has no duplicates.
    - Cell column index has no duplicates.
    - All values are finite.
    - All values are ≥ 0.

    Parameters
    ----------
    counts:
        Genes-as-rows, cells-as-columns count DataFrame.
    allow_float:
        If False, raise when values are not integer-valued.
    """
    validate_no_duplicate_genes(counts.index, source="counts matrix index (genes)")

    if counts.columns.duplicated().any():
        dupes = counts.columns[counts.columns.duplicated()].tolist()
        raise ValueError(
            f"Duplicate cell IDs in counts columns: {dupes[:5]}…"
        )

    vals = counts.to_numpy(dtype=float)
    if not np.isfinite(vals).all():
        raise ValueError("counts matrix contains non-finite values (NaN or Inf).")
    if (vals < 0).any():
        raise ValueError("counts matrix contains negative values.")

    if not allow_float:
        if not np.allclose(vals, np.round(vals)):
            raise ValueError(
                "counts matrix contains non-integer values.  "
                "Pass raw integer counts for reference construction."
            )
