"""
I/O routines for TissueResolve.

Public API (Stage 1)
--------------------
``load_reference_h5ad``, ``load_reference_csv``
    Load single-cell / single-nucleus reference data.

``check_gene_overlap``, ``validate_no_duplicate_genes``,
``validate_matrix_shape``, ``validate_cell_type_column``,
``validate_gene_names_match``, ``validate_counts_matrix``
    Input validation utilities used throughout the pipeline.

Public API (Stage 3 — bulk)
----------------------------
``read_bulk_counts``
    Load a genes × samples bulk count matrix from TSV or CSV.
``read_bulk_metadata``
    Load sample-level metadata from TSV or CSV.

Future (Stage 4)
----------------
``spatial`` — load 10x Visium data; save spatial results.
"""
from tissueresolve.io.reference import load_reference_csv, load_reference_h5ad
from tissueresolve.io.validation import (
    check_gene_overlap,
    validate_cell_type_column,
    validate_counts_matrix,
    validate_gene_names_match,
    validate_matrix_shape,
    validate_no_duplicate_genes,
)
from tissueresolve.io.bulk import read_bulk_counts, read_bulk_metadata

__all__ = [
    # reference I/O
    "load_reference_h5ad",
    "load_reference_csv",
    # validation
    "check_gene_overlap",
    "validate_no_duplicate_genes",
    "validate_matrix_shape",
    "validate_cell_type_column",
    "validate_gene_names_match",
    "validate_counts_matrix",
    # bulk I/O (Stage 3)
    "read_bulk_counts",
    "read_bulk_metadata",
]
