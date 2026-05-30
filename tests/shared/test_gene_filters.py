"""
Tests: GeneFilterSet — instance-based loading, independence, and fallback.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from tissueresolve.reference.gene_filters import GeneFilterSet


class TestGeneFilterSetNoGlobalState:

    def test_two_instances_are_independent(self):
        """Caches on instance A must not affect instance B."""
        gfs_hg = GeneFilterSet(genome="hg38")
        gfs_mm = GeneFilterSet(genome="mm10")
        # Loading on one instance must not fill the other instance's cache
        _ = gfs_hg.blacklist_prefixes()
        assert gfs_mm._prefix_cache is None or gfs_mm._prefix_cache != gfs_hg._prefix_cache

    def test_genome_attribute_set_correctly(self):
        gfs = GeneFilterSet(genome="mm10")
        assert gfs.genome == "mm10"

    def test_instance_cache_starts_empty(self):
        gfs = GeneFilterSet(genome="hg38")
        assert len(gfs._cache) == 0

    def test_cache_populated_after_access(self):
        gfs = GeneFilterSet(genome="hg38")
        _ = gfs.blacklist_prefixes()
        # prefix cache is separate from _cache; just verify it's populated
        assert gfs._prefix_cache is not None


class TestGeneFilterSetMissingFiles:

    def test_missing_gene_list_returns_empty_frozenset(self, tmp_path):
        """A missing gene list file must return an empty frozenset, not raise."""
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        result = gfs.intronic_dominant()
        assert isinstance(result, frozenset)
        assert len(result) == 0

    def test_missing_file_logs_warning(self, tmp_path, caplog):
        import logging
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        with caplog.at_level(logging.WARNING, logger="tissueresolve.reference.gene_filters"):
            _ = gfs.dissociation_stress()
        assert "not found" in caplog.text.lower() or "disabled" in caplog.text.lower()

    def test_missing_blacklist_file_returns_defaults(self, tmp_path):
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        prefixes = gfs.blacklist_prefixes()
        assert isinstance(prefixes, list)
        assert len(prefixes) > 0
        assert "MT-" in prefixes


class TestGeneFilterSetWithRealFiles:

    def test_load_gene_set_from_real_file(self, tmp_path):
        """Verify loading a gene list file works correctly."""
        gene_list_dir = tmp_path
        test_file = gene_list_dir / "intronic_dominant_hg38.tsv"
        test_file.write_text("GENE_A\nGENE_B\n# comment\nGENE_C\n")
        gfs = GeneFilterSet(genome="hg38", data_dir=gene_list_dir)
        result = gfs.intronic_dominant()
        assert "GENE_A" in result
        assert "GENE_B" in result
        assert "GENE_C" in result
        assert "# comment" not in result
        assert len(result) == 3

    def test_load_blacklist_prefixes_from_file(self, tmp_path):
        prefix_file = tmp_path / "blacklist_prefixes.txt"
        prefix_file.write_text("MT-\nRPL\n# this is a comment\nHSP\n")
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        prefixes = gfs.blacklist_prefixes()
        assert "MT-" in prefixes
        assert "RPL" in prefixes
        assert "HSP" in prefixes
        assert "# this is a comment" not in prefixes

    def test_cached_after_first_load(self, tmp_path):
        test_file = tmp_path / "hypervariable_inflammatory.tsv"
        test_file.write_text("FOS\nJUN\n")
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        result1 = gfs.hypervariable_inflammatory()
        result2 = gfs.hypervariable_inflammatory()
        assert result1 is result2  # same object from cache


class TestBlacklistPrefixFilter:

    def test_filter_removes_mt_genes(self, tmp_path):
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        genes = ["MT-ND1", "GAPDH", "MT-CO1", "ACTB", "RPL3"]
        kept, removed = gfs.filter_by_blacklist_prefixes(genes)
        assert "MT-ND1" in removed
        assert "MT-CO1" in removed
        assert "RPL3" in removed
        assert "GAPDH" in kept
        assert "ACTB" in kept

    def test_no_genes_removed_when_no_blacklist_match(self, tmp_path):
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        genes = ["GAPDH", "ACTB", "TP53"]
        kept, removed = gfs.filter_by_blacklist_prefixes(genes)
        assert set(kept) == set(genes)
        assert removed == []

    def test_empty_input(self, tmp_path):
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        kept, removed = gfs.filter_by_blacklist_prefixes([])
        assert kept == []
        assert removed == []


class TestGenesRemovedReport:

    def test_report_returns_dict_with_categories(self, tmp_path):
        # Set up a protein_coding file
        pc_file = tmp_path / "protein_coding_hg38.tsv"
        pc_file.write_text("GAPDH\nACTB\n")
        hv_file = tmp_path / "hypervariable_inflammatory.tsv"
        hv_file.write_text("FOS\nJUN\n")

        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        genes = ["GAPDH", "ACTB", "FOS", "MT-ND1", "LINC01234"]
        report = gfs.genes_removed_report(genes)

        assert isinstance(report, dict)
        assert "blacklist_prefix" in report
        assert "not_protein_coding" in report
        assert "hypervariable_inflammatory" in report

        assert "MT-ND1" in report["blacklist_prefix"]
        assert "LINC01234" in report["blacklist_prefix"]
        assert "FOS" in report["hypervariable_inflammatory"]

    def test_report_no_silent_removal(self, tmp_path):
        """Every key in the report must have an explicit list (may be empty)."""
        gfs = GeneFilterSet(genome="hg38", data_dir=tmp_path)
        report = gfs.genes_removed_report(["GAPDH"])
        for key, val in report.items():
            assert isinstance(val, list), f"Value for '{key}' must be a list"
