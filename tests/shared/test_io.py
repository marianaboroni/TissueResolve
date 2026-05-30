"""
Tests: io/validation.py and io/reference.py.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.io.validation import (
    check_gene_overlap,
    validate_cell_type_column,
    validate_counts_matrix,
    validate_gene_names_match,
    validate_matrix_shape,
    validate_no_duplicate_genes,
)


# ---------------------------------------------------------------------------
# check_gene_overlap
# ---------------------------------------------------------------------------

class TestCheckGeneOverlap:

    def test_returns_sorted_intersection(self):
        q = ["C", "A", "B", "D"]
        r = ["A", "B", "E", "F"]
        shared = check_gene_overlap(q, r)
        assert shared == ["A", "B"]

    def test_empty_intersection(self):
        shared = check_gene_overlap(["X", "Y"], ["A", "B"])
        assert shared == []

    def test_full_overlap(self):
        genes = ["A", "B", "C"]
        shared = check_gene_overlap(genes, genes)
        assert shared == sorted(genes)

    def test_warning_below_min_overlap(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            check_gene_overlap(["A"], ["A", "B", "C"], min_overlap=10)
        assert any("overlap" in str(x.message).lower() for x in w)

    def test_no_warning_above_min_overlap(self):
        genes = [f"G{i}" for i in range(100)]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            check_gene_overlap(genes, genes, min_overlap=50)
        overlap_warnings = [x for x in w if "overlap" in str(x.message).lower()]
        assert len(overlap_warnings) == 0


# ---------------------------------------------------------------------------
# validate_no_duplicate_genes
# ---------------------------------------------------------------------------

class TestValidateNoDuplicateGenes:

    def test_no_duplicates_passes(self):
        validate_no_duplicate_genes(["A", "B", "C"])  # should not raise

    def test_duplicate_raises(self):
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            validate_no_duplicate_genes(["A", "B", "A"])

    def test_pandas_index_no_duplicates(self):
        idx = pd.Index(["A", "B", "C"])
        validate_no_duplicate_genes(idx)  # should not raise

    def test_pandas_index_with_duplicates(self):
        idx = pd.Index(["A", "B", "A"])
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            validate_no_duplicate_genes(idx)

    def test_empty_list_passes(self):
        validate_no_duplicate_genes([])


# ---------------------------------------------------------------------------
# validate_matrix_shape
# ---------------------------------------------------------------------------

class TestValidateMatrixShape:

    def test_correct_shape_passes(self):
        M = np.zeros((5, 3))
        validate_matrix_shape(M, (5, 3), "test_matrix")  # no raise

    def test_wrong_shape_raises(self):
        M = np.zeros((5, 3))
        with pytest.raises(ValueError, match="test_matrix"):
            validate_matrix_shape(M, (5, 4), "test_matrix")

    def test_wildcard_dimension(self):
        M = np.zeros((5, 3))
        validate_matrix_shape(M, (-1, 3), "test_matrix")  # -1 matches any

    def test_wrong_ndim_raises(self):
        M = np.zeros((5, 3))
        with pytest.raises(ValueError, match="3-D"):
            validate_matrix_shape(M, (5, 3, 1), "test_matrix")


# ---------------------------------------------------------------------------
# validate_cell_type_column
# ---------------------------------------------------------------------------

class TestValidateCellTypeColumn:

    def test_present_column_passes(self):
        meta = pd.DataFrame({"cell_type": ["A", "B"]})
        validate_cell_type_column(meta, "cell_type")  # no raise

    def test_missing_column_raises(self):
        meta = pd.DataFrame({"other": ["A", "B"]})
        with pytest.raises(KeyError, match="cell_type"):
            validate_cell_type_column(meta, "cell_type")

    def test_custom_column_name(self):
        meta = pd.DataFrame({"annotation": ["A"]})
        validate_cell_type_column(meta, "annotation")  # no raise


# ---------------------------------------------------------------------------
# validate_gene_names_match
# ---------------------------------------------------------------------------

class TestValidateGeneNamesMatch:

    def test_matching_lists_pass(self):
        g = ["A", "B", "C"]
        validate_gene_names_match(g, g, "X", "Y")  # no raise

    def test_different_lists_raise(self):
        with pytest.raises(ValueError, match="do not match"):
            validate_gene_names_match(["A", "B"], ["A", "C"], "X", "Y")

    def test_different_order_raises(self):
        with pytest.raises(ValueError):
            validate_gene_names_match(["A", "B", "C"], ["C", "B", "A"])

    def test_labels_in_error_message(self):
        with pytest.raises(ValueError, match="REF") as exc_info:
            validate_gene_names_match(["A"], ["B"], "REF", "DATA")
        assert "DATA" in str(exc_info.value)


# ---------------------------------------------------------------------------
# validate_counts_matrix
# ---------------------------------------------------------------------------

class TestValidateCountsMatrix:

    def test_valid_matrix_passes(self):
        rng = np.random.default_rng(0)
        X = rng.integers(0, 100, (20, 50)).astype(float)
        df = pd.DataFrame(X, index=[f"G{i}" for i in range(20)],
                          columns=[f"c{i}" for i in range(50)])
        validate_counts_matrix(df)  # no raise

    def test_duplicate_genes_raises(self):
        X = np.ones((3, 5))
        df = pd.DataFrame(X, index=["G1", "G2", "G1"],
                          columns=[f"c{i}" for i in range(5)])
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            validate_counts_matrix(df)

    def test_nan_raises(self):
        X = np.ones((3, 5))
        X[0, 0] = float("nan")
        df = pd.DataFrame(X, index=[f"G{i}" for i in range(3)],
                          columns=[f"c{i}" for i in range(5)])
        with pytest.raises(ValueError, match="non-finite"):
            validate_counts_matrix(df)

    def test_negative_values_raise(self):
        X = np.ones((3, 5))
        X[0, 0] = -1.0
        df = pd.DataFrame(X, index=[f"G{i}" for i in range(3)],
                          columns=[f"c{i}" for i in range(5)])
        with pytest.raises(ValueError, match="negative"):
            validate_counts_matrix(df)

    def test_duplicate_cell_ids_raise(self):
        X = np.ones((3, 5))
        df = pd.DataFrame(X, index=[f"G{i}" for i in range(3)],
                          columns=["c1", "c2", "c3", "c4", "c1"])
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            validate_counts_matrix(df)


# ---------------------------------------------------------------------------
# io/reference.py
# ---------------------------------------------------------------------------

class TestLoadReferenceCsv:

    def test_load_from_csv(self, tmp_path):
        from tissueresolve.io.reference import load_reference_csv

        # Write a small counts TSV
        genes = [f"G{i}" for i in range(5)]
        cells = [f"c{i}" for i in range(4)]
        X = np.arange(20).reshape(5, 4).astype(float)
        counts_path = tmp_path / "counts.tsv"
        meta_path = tmp_path / "meta.tsv"

        pd.DataFrame(X, index=genes, columns=cells).to_csv(counts_path, sep="\t")
        pd.DataFrame(
            {"cell_type": ["A", "B", "A", "B"]}, index=cells
        ).to_csv(meta_path, sep="\t")

        counts, meta = load_reference_csv(counts_path, meta_path)
        assert list(counts.index) == genes
        assert "cell_type" in meta.columns

    def test_missing_counts_file_raises(self, tmp_path):
        from tissueresolve.io.reference import load_reference_csv

        with pytest.raises(FileNotFoundError):
            load_reference_csv(tmp_path / "nonexistent.tsv", tmp_path / "meta.tsv")

    def test_missing_meta_file_raises(self, tmp_path):
        from tissueresolve.io.reference import load_reference_csv

        counts_path = tmp_path / "counts.tsv"
        pd.DataFrame(np.ones((3, 3))).to_csv(counts_path, sep="\t")
        with pytest.raises(FileNotFoundError):
            load_reference_csv(counts_path, tmp_path / "nonexistent.tsv")

    def test_missing_celltype_col_raises(self, tmp_path):
        from tissueresolve.io.reference import load_reference_csv

        counts_path = tmp_path / "counts.tsv"
        meta_path = tmp_path / "meta.tsv"
        pd.DataFrame(np.ones((3, 4)), index=[f"G{i}" for i in range(3)]).to_csv(counts_path, sep="\t")
        pd.DataFrame({"other": ["x", "y", "z", "w"]}).to_csv(meta_path, sep="\t")

        with pytest.raises(KeyError, match="cell_type"):
            load_reference_csv(counts_path, meta_path)


class TestLoadReferenceH5ad:

    def test_load_h5ad_basic(self, tmp_path):
        import anndata as ad
        import scipy.sparse as sp
        from tissueresolve.io.reference import load_reference_h5ad

        rng = np.random.default_rng(0)
        X = rng.integers(0, 50, (20, 10)).astype(np.float32)
        obs = pd.DataFrame({"cell_type": ["A"] * 10 + ["B"] * 10})
        var = pd.DataFrame(index=[f"G{i}" for i in range(10)])
        adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
        h5_path = tmp_path / "ref.h5ad"
        adata.write_h5ad(h5_path)

        loaded = load_reference_h5ad(h5_path, cell_type_col="cell_type")
        assert loaded.n_obs == 20
        assert loaded.n_vars == 10

    def test_missing_file_raises(self, tmp_path):
        from tissueresolve.io.reference import load_reference_h5ad

        with pytest.raises(FileNotFoundError):
            load_reference_h5ad(tmp_path / "nonexistent.h5ad")

    def test_missing_celltype_col_raises(self, tmp_path):
        import anndata as ad
        from tissueresolve.io.reference import load_reference_h5ad

        X = np.ones((4, 5)).astype(np.float32)
        obs = pd.DataFrame({"other": ["x"] * 4})
        var = pd.DataFrame(index=[f"G{i}" for i in range(5)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        h5_path = tmp_path / "bad.h5ad"
        adata.write_h5ad(h5_path)

        with pytest.raises(KeyError, match="cell_type"):
            load_reference_h5ad(h5_path, cell_type_col="cell_type")
