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
    "detect_broad_cell_type_col",
    "detect_fine_cell_type_col",
    "validate_hierarchical_annotations",
    "summarize_hierarchical_annotations",
]

# Candidate obs columns for broad (compartment / lineage) labels, in priority
# order.  Detection is a *convenience*; hierarchical mode never silently infers
# a hierarchy without telling the user (the CLI prints what it found).
BROAD_COL_CANDIDATES = (
    "broad_cell_type", "major_cell_type", "cell_type_major", "cell_type_broad",
    "compartment", "lineage", "cell_class", "broad_annotation",
)
# Candidate obs columns for fine (subpopulation / cell-state) labels.
FINE_COL_CANDIDATES = (
    "sub_cell_type", "cell_type_fine", "cell_type", "cell_state", "subtype",
    "annotation", "author_cell_type", "fine_annotation",
)

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


# ---------------------------------------------------------------------------
# Hierarchical (broad/fine) annotation detection & validation
# ---------------------------------------------------------------------------


def _detect_col(obs: pd.DataFrame, candidates, *, exclude: Optional[str] = None):
    """Return the first *candidates* column present in *obs* (case-insensitive)."""
    lower = {str(c).lower(): c for c in obs.columns}
    for cand in candidates:
        col = lower.get(cand.lower())
        if col is not None and col != exclude:
            return col
    return None


def detect_broad_cell_type_col(obs: pd.DataFrame) -> Optional[str]:
    """Detect a broad/compartment annotation column in ``adata.obs``.

    Returns the matched column name, or ``None`` if no known candidate is
    present.  Never guesses silently — callers must report what was detected.
    """
    return _detect_col(obs, BROAD_COL_CANDIDATES)


def detect_fine_cell_type_col(
    obs: pd.DataFrame, *, exclude: Optional[str] = None,
) -> Optional[str]:
    """Detect a fine/subpopulation annotation column in ``adata.obs``.

    *exclude* (e.g. the detected broad column) is skipped so the same column is
    never returned for both roles.
    """
    return _detect_col(obs, FINE_COL_CANDIDATES, exclude=exclude)


def validate_hierarchical_annotations(
    obs: pd.DataFrame,
    broad_col: str,
    fine_col: str,
    *,
    max_fine_per_broad: int = 25,
    min_cells_per_fine: int = 10,
) -> dict:
    """Validate broad/fine annotation columns for hierarchical deconvolution.

    Validation rules (errors are raised, soft issues are warnings):

    * both columns must exist in *obs*;
    * neither column may contain missing values;
    * every fine label must map to exactly **one** broad label;
    * (warn) a broad family with > ``max_fine_per_broad`` fine labels;
    * (warn) a fine label with < ``min_cells_per_fine`` cells.

    Returns
    -------
    dict
        ``{"mapping": {fine: broad}, "warnings": [...],
        "n_broad": int, "n_fine": int}``.  The returned mapping is suitable for
        :func:`tissueresolve.reference.hierarchy.build_cell_type_hierarchy`.
    """
    avail = list(obs.columns)
    for role, col in (("broad", broad_col), ("fine", fine_col)):
        if col not in obs.columns:
            raise KeyError(
                f"Hierarchical {role} cell-type column {col!r} not found in "
                f"adata.obs.  Available columns: {avail}.  "
                "Set --broad-cell-type-col / --fine-cell-type-col to existing "
                "columns, or provide --cell-type-hierarchy mapping.tsv."
            )
    if broad_col == fine_col:
        raise ValueError(
            "Hierarchical broad and fine cell-type columns must differ "
            f"(both are {broad_col!r})."
        )

    broad = obs[broad_col].astype("object")
    fine = obs[fine_col].astype("object")
    for role, col, ser in (("broad", broad_col, broad), ("fine", fine_col, fine)):
        n_missing = int(ser.isna().sum() + (ser.astype(str).str.strip() == "").sum())
        if n_missing:
            raise ValueError(
                f"Hierarchical {role} column {col!r} has {n_missing} "
                "missing/blank value(s).  Annotate every cell or subset the "
                "reference before hierarchical deconvolution."
            )

    # every fine label must map to exactly one broad label
    pairs = (
        pd.DataFrame({"fine": fine.astype(str), "broad": broad.astype(str)})
        .drop_duplicates()
    )
    multi = pairs.groupby("fine")["broad"].nunique()
    ambiguous = multi[multi > 1]
    if len(ambiguous):
        examples = {
            f: sorted(pairs.loc[pairs["fine"] == f, "broad"].unique())
            for f in list(ambiguous.index)[:5]
        }
        raise ValueError(
            "Each fine cell type must map to exactly one broad family, but "
            f"{len(ambiguous)} fine label(s) map to multiple broad labels: "
            f"{examples}.  Fix the reference annotations (one broad family per "
            "fine cell type) before hierarchical deconvolution."
        )

    mapping = dict(zip(pairs["fine"], pairs["broad"]))

    issues: list[str] = []
    # broad families with too many fine subtypes
    per_broad = pairs.groupby("broad")["fine"].nunique()
    crowded = per_broad[per_broad > max_fine_per_broad]
    for fam, n in crowded.items():
        msg = (f"broad family {fam!r} contains {n} fine subtypes "
               f"(> {max_fine_per_broad}); within-family separability may be low.")
        issues.append(msg)
        warnings.warn(msg, stacklevel=2)
    # rare fine labels
    fine_counts = fine.astype(str).value_counts()
    rare = fine_counts[fine_counts < min_cells_per_fine]
    for label, n in rare.items():
        msg = (f"fine cell type {label!r} has only {n} cell(s) "
               f"(< {min_cells_per_fine}); its signature may be unreliable.")
        issues.append(msg)
        warnings.warn(msg, stacklevel=2)

    return {
        "mapping": mapping,
        "warnings": issues,
        "n_broad": int(pairs["broad"].nunique()),
        "n_fine": int(pairs["fine"].nunique()),
    }


def summarize_hierarchical_annotations(
    obs: pd.DataFrame,
    broad_col: str,
    fine_col: str,
) -> pd.DataFrame:
    """Per-(broad, fine) summary table: cell counts and subtypes-per-family.

    Returns a tidy DataFrame with columns
    ``broad_cell_type, fine_cell_type, n_cells, n_fine_in_family`` sorted by
    family then descending cell count.  Saved by reference construction as
    ``hierarchy_summary``-style tables.
    """
    df = pd.DataFrame({
        "broad_cell_type": obs[broad_col].astype(str),
        "fine_cell_type": obs[fine_col].astype(str),
    })
    counts = (
        df.groupby(["broad_cell_type", "fine_cell_type"])
        .size().rename("n_cells").reset_index()
    )
    n_fine = (
        counts.groupby("broad_cell_type")["fine_cell_type"]
        .transform("nunique")
    )
    counts["n_fine_in_family"] = n_fine
    return counts.sort_values(
        ["broad_cell_type", "n_cells"], ascending=[True, False]
    ).reset_index(drop=True)
