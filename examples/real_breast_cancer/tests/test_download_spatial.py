"""
Offline tests for the spatial download path and skip-if-exists / --force logic
in 00_download_data.py.

No network: the Visium loader is injected, the Census module is faked, and the
tarfile compatibility shim is exercised against an injected fake module.
"""
from __future__ import annotations

import sys
import types

import numpy as np
import pytest


def _counter_clock(step=1.0):
    state = {"t": 0.0}

    def clock():
        state["t"] += step
        return state["t"]

    return clock


@pytest.fixture
def visium_adata():
    """Tiny Visium-like AnnData: count-like X, obsm['spatial'], no counts layer."""
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(0)
    n, g = 12, 10
    X = rng.poisson(4, (n, g)).astype("float32")
    obs = pd.DataFrame(index=[f"spot_{i}" for i in range(n)])
    var = pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(g)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.obsm["spatial"] = rng.integers(0, 100, (n, 2)).astype(float)
    return a


def _fake_census_module(*, adata=None):
    mod = types.ModuleType("cellxgene_census")

    class _Soma:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def open_soma(census_version=None, **kw):
        if adata is None:
            raise AssertionError("open_soma should not be called when skipping")
        return _Soma()

    def get_anndata(census, **kw):
        return adata

    mod.open_soma = open_soma
    mod.get_anndata = get_anndata
    return mod


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_parses_force(load_script):
    m00 = load_script("00_download_data.py")
    assert m00.build_parser().parse_args(["--force"]).force is True
    assert m00.build_parser().parse_args([]).force is False


def test_dry_run_accepts_force(load_script, harness, monkeypatch, tmp_path):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    monkeypatch.setattr(harness, "ALL_DIRS", [tmp_path / "d"])
    monkeypatch.setattr(harness, "MANIFEST_PATH", tmp_path / "m.json")
    m00 = load_script("00_download_data.py")
    assert m00.main(["--force"]) == 0


# ---------------------------------------------------------------------------
# tarfile compatibility shim
# ---------------------------------------------------------------------------


def test_tarfile_shim_patches_when_missing(load_script):
    m00 = load_script("00_download_data.py")
    fake = types.SimpleNamespace()  # no data_filter
    available = m00.ensure_tarfile_data_filter(fake)
    assert available is False
    assert hasattr(fake, "data_filter")
    assert callable(fake.data_filter)


def test_tarfile_shim_noop_when_present(load_script):
    m00 = load_script("00_download_data.py")
    sentinel = object()
    fake = types.SimpleNamespace(data_filter=sentinel)
    available = m00.ensure_tarfile_data_filter(fake)
    assert available is True
    assert fake.data_filter is sentinel


# ---------------------------------------------------------------------------
# Reference skip / force
# ---------------------------------------------------------------------------


def test_reference_skipped_when_exists(load_script, harness, monkeypatch, tmp_path):
    monkeypatch.setattr(harness, "REFERENCE_H5AD", tmp_path / "ref.h5ad")
    (tmp_path / "ref.h5ad").write_bytes(b"existing-reference")
    # Census that explodes if queried — proves we skip before querying.
    monkeypatch.setitem(sys.modules, "cellxgene_census", _fake_census_module())
    m00 = load_script("00_download_data.py")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_reference(manifest, force=False,
                                clock=_counter_clock(), now=lambda: "t")
    assert ok is True
    entry = manifest["datasets"]["reference"]
    assert entry["status"] == "already_exists"
    assert entry["method"] == "reuse_existing"
    assert entry["file_size_bytes"] > 0


def test_reference_force_redownloads(load_script, harness, monkeypatch, tmp_path,
                                     visium_adata):
    import anndata as ad
    import pandas as pd

    # build a tiny census-style adata with a cell_type column
    a = ad.AnnData(X=np.ones((20, 8), dtype="float32"),
                   obs=pd.DataFrame({"cell_type": ["A"] * 10 + ["B"] * 10}))
    monkeypatch.setattr(harness, "REFERENCE_H5AD", tmp_path / "ref.h5ad")
    (tmp_path / "ref.h5ad").write_bytes(b"old")
    monkeypatch.setitem(sys.modules, "cellxgene_census", _fake_census_module(adata=a))
    m00 = load_script("00_download_data.py")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_reference(manifest, force=True,
                                clock=_counter_clock(), now=lambda: "t")
    assert ok is True
    assert manifest["datasets"]["reference"]["status"] == "downloaded"


# ---------------------------------------------------------------------------
# Spatial download via injected loader
# ---------------------------------------------------------------------------


def test_spatial_download_success(load_script, harness, monkeypatch, tmp_path,
                                  visium_adata, capsys):
    monkeypatch.setattr(harness, "SPATIAL_H5AD", tmp_path / "sp.h5ad")
    m00 = load_script("00_download_data.py")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_spatial(
        manifest, loader=lambda sample: visium_adata,
        clock=_counter_clock(), now=lambda: "t",
    )
    assert ok is True
    assert (tmp_path / "sp.h5ad").exists()

    entry = manifest["datasets"]["spatial"]
    assert entry["spatial_status"] == "downloaded"
    assert entry["spatial_download_method"] == "scanpy.datasets.visium_sge"
    assert entry["spatial_has_counts_layer"] is True   # copied from count-like X
    assert entry["spatial_has_spatial_obsm"] is True
    assert entry["spatial_file_size_bytes"] > 0
    assert "python_version" in entry
    assert "tarfile_data_filter_available" in entry
    harness.validate_manifest(manifest)

    out = capsys.readouterr().out
    assert "[ 10%]" in out and "[100%]" in out


def test_spatial_non_count_X_not_copied(load_script, harness, monkeypatch, tmp_path):
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(1)
    a = ad.AnnData(X=rng.random((10, 6)).astype("float32"))  # floats, not counts
    a.obsm["spatial"] = rng.random((10, 2))
    monkeypatch.setattr(harness, "SPATIAL_H5AD", tmp_path / "sp.h5ad")
    m00 = load_script("00_download_data.py")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_spatial(manifest, loader=lambda s: a,
                              clock=_counter_clock(), now=lambda: "t")
    assert ok is True
    assert manifest["datasets"]["spatial"]["spatial_has_counts_layer"] is False


def test_spatial_skip_when_exists(load_script, harness, monkeypatch, tmp_path,
                                  visium_adata):
    spatial_path = tmp_path / "sp.h5ad"
    visium_adata.layers["counts"] = visium_adata.X.copy()
    visium_adata.write_h5ad(spatial_path)
    monkeypatch.setattr(harness, "SPATIAL_H5AD", spatial_path)
    m00 = load_script("00_download_data.py")

    def _boom(sample):
        raise AssertionError("loader must not run when skipping")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_spatial(manifest, force=False, loader=_boom,
                              clock=_counter_clock(), now=lambda: "t")
    assert ok is True
    entry = manifest["datasets"]["spatial"]
    assert entry["spatial_status"] == "already_exists"
    assert entry["spatial_has_counts_layer"] is True
    assert entry["spatial_has_spatial_obsm"] is True


def test_spatial_loader_failure_is_actionable(load_script, harness, monkeypatch,
                                              tmp_path, capsys):
    monkeypatch.setattr(harness, "SPATIAL_H5AD", tmp_path / "sp.h5ad")
    m00 = load_script("00_download_data.py")

    def _fail(sample):
        raise RuntimeError("module 'tarfile' has no attribute 'data_filter'")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_spatial(manifest, loader=_fail,
                              clock=_counter_clock(), now=lambda: "t")
    assert ok is False
    assert not (tmp_path / "sp.h5ad").exists()           # nothing written
    assert manifest["datasets"]["spatial"]["spatial_status"] == "failed"
    err = capsys.readouterr().err
    assert "spatial download failed" in err
    assert "Manual spatial download" in err
