"""
Tests for tissueresolve.io.bulk (read_bulk_counts, read_bulk_metadata).
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tissueresolve.io.bulk import read_bulk_counts, read_bulk_metadata


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t")


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep=",")


def _make_counts(n_genes: int = 20, n_samples: int = 5, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    data = rng.integers(0, 1000, (n_genes, n_samples)).astype(float)
    genes = [f"GENE_{i:03d}" for i in range(n_genes)]
    samples = [f"sample_{i}" for i in range(n_samples)]
    df = pd.DataFrame(data, index=genes, columns=samples)
    df.index.name = "gene"
    return df


# ---------------------------------------------------------------------------
# read_bulk_counts — basic loading
# ---------------------------------------------------------------------------


class TestReadBulkCountsTSV:
    def test_loads_tsv(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert result.shape == df.shape

    def test_loads_csv(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.csv"
        _write_csv(df, path)
        result = read_bulk_counts(path)
        assert result.shape == df.shape

    def test_index_is_gene_names(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert list(result.index) == list(df.index)

    def test_columns_are_sample_names(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert list(result.columns) == list(df.columns)

    def test_values_match(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert np.allclose(result.values, df.values)

    def test_index_name_is_gene(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert result.index.name == "gene"

    def test_comments_skipped(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        with path.open("w") as fh:
            fh.write("# estimate_type: mRNA_proportion\n")
            fh.write("# WARNING: ...\n")
            df.to_csv(fh, sep="\t")
        result = read_bulk_counts(path)
        assert result.shape == df.shape

    def test_auto_sep_tsv(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "data.txt"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert result.shape == df.shape


# ---------------------------------------------------------------------------
# read_bulk_counts — validation errors
# ---------------------------------------------------------------------------


class TestReadBulkCountsValidation:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_bulk_counts(tmp_path / "does_not_exist.tsv")

    def test_duplicate_genes_raise(self, tmp_path):
        df = _make_counts(n_genes=5)
        # Force duplicate
        df = pd.concat([df, df.iloc[:1]])
        path = tmp_path / "dup_genes.tsv"
        _write_tsv(df, path)
        with pytest.raises(ValueError, match="Duplicate gene names"):
            read_bulk_counts(path)

    def test_negative_values_raise(self, tmp_path):
        df = _make_counts()
        df.iloc[0, 0] = -1.0
        path = tmp_path / "neg.tsv"
        _write_tsv(df, path)
        with pytest.raises(ValueError, match="Negative values"):
            read_bulk_counts(path)

    def test_duplicate_sample_names_warn(self, tmp_path):
        df = _make_counts()
        df.columns = pd.Index(["s1", "s1", "s2", "s3", "s4"])
        path = tmp_path / "dup_samples.tsv"
        _write_tsv(df, path)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            read_bulk_counts(path)
        texts = [str(w.message) for w in caught]
        assert any("Duplicate sample" in t for t in texts)

    def test_returns_non_negative_matrix(self, tmp_path):
        df = _make_counts()
        path = tmp_path / "counts.tsv"
        _write_tsv(df, path)
        result = read_bulk_counts(path)
        assert (result.values >= 0).all()


# ---------------------------------------------------------------------------
# read_bulk_metadata
# ---------------------------------------------------------------------------


class TestReadBulkMetadata:
    def _make_meta(self, n_samples: int = 5) -> pd.DataFrame:
        samples = [f"sample_{i}" for i in range(n_samples)]
        return pd.DataFrame({
            "protocol": ["polyA"] * n_samples,
            "tissue": ["PBMC"] * n_samples,
        }, index=pd.Index(samples, name="sample_id"))

    def test_loads_tsv(self, tmp_path):
        meta = self._make_meta()
        path = tmp_path / "meta.tsv"
        meta.to_csv(path, sep="\t")
        result = read_bulk_metadata(path)
        assert result.shape == meta.shape

    def test_loads_csv(self, tmp_path):
        meta = self._make_meta()
        path = tmp_path / "meta.csv"
        meta.to_csv(path, sep=",")
        result = read_bulk_metadata(path)
        assert result.shape == meta.shape

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_bulk_metadata(tmp_path / "missing.tsv")

    def test_duplicate_sample_ids_raise(self, tmp_path):
        meta = self._make_meta()
        meta = pd.concat([meta, meta.iloc[:1]])
        path = tmp_path / "dup.tsv"
        meta.to_csv(path, sep="\t")
        with pytest.raises(ValueError, match="Duplicate sample"):
            read_bulk_metadata(path)

    def test_columns_preserved(self, tmp_path):
        meta = self._make_meta()
        path = tmp_path / "meta.tsv"
        meta.to_csv(path, sep="\t")
        result = read_bulk_metadata(path)
        assert set(result.columns) == {"protocol", "tissue"}
