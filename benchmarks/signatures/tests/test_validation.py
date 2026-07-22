"""Offline regression tests for the multipanel-benchmark guards (Stage 1.5)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from benchmarks.signatures import validation as V  # noqa: E402


def test_hierarchy_covers_dataset_ok():
    V.validate_hierarchy_covers_dataset({"A": "fam", "B": "fam"}, ["A", "B"], dataset="x")


def test_wrong_tissue_hierarchy_raises():
    """THE regression: a breast hierarchy applied to lung types must abort."""
    breast_hier = {"T cell": "T/NK", "myeloid cell": "Myeloid", "fibroblast": "Stroma"}
    lung_types = ["Alveolar Type 1", "Alveolar Type 2", "Club", "Ciliated", "Basal"]
    with pytest.raises(V.BenchmarkValidationError) as e:
        V.validate_hierarchy_covers_dataset(breast_hier, lung_types, dataset="lung")
    assert "WRONG tissue" in str(e.value)


def test_partial_missing_raises_without_wrong_tissue_hint():
    with pytest.raises(V.BenchmarkValidationError) as e:
        V.validate_hierarchy_covers_dataset({"A": "f", "B": "f", "C": "f"}, ["A", "B", "C", "D"])
    assert "add them to the mapping" in str(e.value)


def test_donor_disjoint():
    V.validate_donor_disjoint(["d1", "d2"], ["d3", "d4"])
    with pytest.raises(V.BenchmarkValidationError):
        V.validate_donor_disjoint(["d1", "d2"], ["d2", "d3"])
    with pytest.raises(V.BenchmarkValidationError):
        V.validate_donor_disjoint([], ["d3"])


def test_selection_donor_leakage():
    V.validate_selection_donors(["d1", "d2"], ["d3"])
    with pytest.raises(V.BenchmarkValidationError):
        V.validate_selection_donors(["d1", "d2"], ["d2"])


def test_budget_recorded():
    V.validate_budget_recorded({"strategy": "x", "n_genes": 300})
    with pytest.raises(V.BenchmarkValidationError):
        V.validate_budget_recorded({"strategy": "x", "n_genes": 0})


def test_matrix_orientation():
    genes = [f"g{i}" for i in range(50)]
    good = pd.DataFrame(0.0, index=genes, columns=["s1", "s2"])          # genes×samples
    V.validate_matrix_orientation(good, genes)
    with pytest.raises(V.BenchmarkValidationError):
        V.validate_matrix_orientation(good.T, genes)                     # transposed


def test_truth_pred_alignment():
    t = pd.DataFrame(0.0, index=["s1", "s2"], columns=["A", "B"])
    p = pd.DataFrame(0.0, index=["s1", "s2"], columns=["A", "B"])
    V.validate_truth_pred_alignment(t, p)
    with pytest.raises(V.BenchmarkValidationError):
        V.validate_truth_pred_alignment(t, p.rename(columns={"A": "X", "B": "Y"}))
