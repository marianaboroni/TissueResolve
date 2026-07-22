"""Leakage- and consistency guards for the multipanel signature benchmark (Stage 1.5).

Every guard RAISES a clear error on violation — the benchmark must abort rather than
produce a silently-invalid result. These are offline, deterministic, and unit-tested.
The explicit regression target is: a hierarchy from one tissue (e.g. breast) must never
be applied silently to another (e.g. lung).
"""
from __future__ import annotations

from typing import Iterable, Mapping


class BenchmarkValidationError(ValueError):
    """Raised when a benchmark precondition (leakage / consistency) is violated."""


def validate_hierarchy_covers_dataset(mapping: Mapping[str, str],
                                      cell_types: Iterable[str],
                                      *, dataset: str = "") -> None:
    """Every fine cell type in the dataset must map to a broad family.

    Regression guard against applying the wrong tissue's hierarchy: if a large fraction
    of the dataset's types are absent from the mapping, that almost certainly means the
    hierarchy does not belong to this dataset.
    """
    cts = [str(c) for c in cell_types]
    if not cts:
        raise BenchmarkValidationError("no cell types provided")
    missing = [c for c in cts if c not in mapping]
    if missing:
        frac = len(missing) / len(cts)
        raise BenchmarkValidationError(
            f"hierarchy does not cover {len(missing)}/{len(cts)} cell type(s) "
            f"({frac:.0%}) for dataset {dataset!r}: {missing[:8]}"
            + ("…" if len(missing) > 8 else "")
            + (" — this looks like the WRONG tissue's hierarchy applied to this dataset."
               if frac > 0.5 else " — add them to the mapping."))


def validate_donor_disjoint(train_donors: Iterable[str], test_donors: Iterable[str]) -> None:
    tr, te = set(map(str, train_donors)), set(map(str, test_donors))
    overlap = tr & te
    if overlap:
        raise BenchmarkValidationError(
            f"train/test donors overlap ({len(overlap)}): {sorted(overlap)[:8]} — "
            "donor-disjoint evaluation is mandatory (no reference-donor leakage).")
    if not tr or not te:
        raise BenchmarkValidationError("train or test donor set is empty")


def validate_selection_donors(selection_donors: Iterable[str],
                              eval_donors: Iterable[str]) -> None:
    """Marker selection must never see evaluation donors."""
    leak = set(map(str, selection_donors)) & set(map(str, eval_donors))
    if leak:
        raise BenchmarkValidationError(
            f"marker selection used {len(leak)} evaluation donor(s): {sorted(leak)[:8]} "
            "— evaluation-donor leakage into gene selection is forbidden.")


def validate_budget_recorded(strategy_meta: Mapping) -> None:
    if strategy_meta.get("n_genes") in (None, 0):
        raise BenchmarkValidationError(
            f"strategy {strategy_meta.get('strategy','?')!r} did not record a gene count; "
            "every strategy must record its gene budget for a fair comparison.")


def validate_truth_pred_alignment(truth, pred) -> None:
    """Truth and prediction must be reconcilable (shared samples + at least the truth
    cell types representable)."""
    if list(truth.index) and not set(truth.index) & set(pred.index):
        raise BenchmarkValidationError("truth and prediction share no samples (mis-orientation?)")
    missing_types = [c for c in truth.columns if c not in pred.columns]
    # a population present in truth must be explicitly handled (0 is fine; absence of
    # the column entirely is a silent drop)
    if len(missing_types) == len(list(truth.columns)):
        raise BenchmarkValidationError(
            "no truth cell type appears in predictions — likely mis-orientation or label mismatch")


def validate_matrix_orientation(bulk, reference_genes, *, min_overlap: int = 20) -> None:
    """bulk must be genes×samples (index = genes) with enough overlap to the reference."""
    ref = set(map(str, reference_genes))
    idx_overlap = len(set(map(str, bulk.index)) & ref)
    col_overlap = len(set(map(str, bulk.columns)) & ref)
    if idx_overlap < min_overlap and col_overlap >= min_overlap:
        raise BenchmarkValidationError(
            "bulk appears transposed (genes on columns, not index) — expected genes×samples")
    if idx_overlap < min_overlap:
        raise BenchmarkValidationError(
            f"only {idx_overlap} bulk genes overlap the reference (< {min_overlap}); "
            "check gene identifiers / orientation")
