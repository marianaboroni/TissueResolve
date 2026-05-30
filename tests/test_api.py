"""
Tests for the top-level convenience API (tissueresolve.api) — Stage 5.

These are thin wrappers; the tests confirm they dispatch correctly and forward
to the discovered pipeline APIs without changing behaviour.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import tissueresolve as tr
from tissueresolve.results import BulkDeconvResult, SpatialDeconvResult


def test_top_level_exports():
    for name in ("deconv_bulk", "deconv_spatial", "build_reference",
                 "generate_report", "plot_results"):
        assert hasattr(tr, name), name


def _toy_bulk_ref():
    genes = [f"GENE_{i:03d}" for i in range(48)]
    cts = ["A", "B", "C"]
    rng = np.random.default_rng(0)
    phi = rng.exponential(1.0, (len(genes), len(cts)))
    phi /= phi.sum(0, keepdims=True)
    from tissueresolve.results import ReferenceSignature

    ref = ReferenceSignature(gene_names=genes, cell_types=cts, phi=phi)
    bulk = pd.DataFrame(
        rng.poisson(40, (len(genes), 3)) + 1,
        index=genes, columns=[f"s{i}" for i in range(3)],
    )
    return bulk, ref


class TestDeconvBulk:
    def test_returns_bulk_result(self):
        bulk, ref = _toy_bulk_ref()
        result = tr.deconv_bulk(bulk, ref)
        assert isinstance(result.deconv, BulkDeconvResult)
        assert result.deconv.ESTIMATE_TYPE == "mRNA_proportion"


class TestDeconvSpatialAndDownstream:
    def _run(self):
        from tissueresolve.spatial.benchmark import simulate_visium
        from tissueresolve.config import TissueResolveConfig

        ds = simulate_visium(n_spots=40, n_types=3, n_genes=30, seed=0)
        cfg = TissueResolveConfig()
        cfg.spatial_solver.max_iter = 15
        result = tr.deconv_spatial(
            ds.dense_counts(), ds.to_reference(), ds.array_row, ds.array_col,
            ds.lib_sizes, list(ds.gene_names),
            marker_genes=list(ds.gene_names), config=cfg,
        )
        return ds, result

    def test_returns_spatial_result(self):
        _, result = self._run()
        assert isinstance(result.deconv, SpatialDeconvResult)

    def test_generate_report_dispatches_spatial(self, tmp_path):
        _, result = self._run()
        out = tr.generate_report(result, tmp_path / "report.html")
        assert out.exists()
        assert "spot-level" in out.read_text().lower() or "composition" in out.read_text().lower()

    def test_plot_results_spatial_requires_coords(self, tmp_path):
        _, result = self._run()
        with pytest.raises(ValueError, match="array_row"):
            tr.plot_results(result, tmp_path)

    def test_plot_results_spatial(self, tmp_path):
        pytest.importorskip("matplotlib")
        ds, result = self._run()
        figs = tr.plot_results(
            result, tmp_path, array_row=ds.array_row, array_col=ds.array_col
        )
        assert figs
        for f in figs:
            assert f.data_paths


class TestBuildReference:
    def test_unsupported_source_raises(self):
        with pytest.raises(TypeError):
            tr.build_reference(12345)
