"""
Offline tests for the CELLxGENE Census download path in 00_download_data.py.

The ``cellxgene_census`` module is replaced by a fake injected into
``sys.modules`` so no network access occurs.  Two scenarios are covered:
a successful (mocked) download and the SOMA-encoding incompatibility that the
fallback logic must turn into a clear, actionable message.
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
def tiny_census_adata():
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(0)
    n, g = 120, 30
    X = rng.poisson(3, (n, g)).astype("float32")
    obs = pd.DataFrame(
        {"cell_type": (["Tcell"] * 60 + ["Bcell"] * 60)},
        index=[f"c{i}" for i in range(n)],
    )
    var = pd.DataFrame(index=[f"GENE_{i:03d}" for i in range(g)])
    return ad.AnnData(X=X, obs=obs, var=var)


def _fake_census_module(*, adata=None, raise_soma=False):
    mod = types.ModuleType("cellxgene_census")

    class _Soma:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def open_soma(census_version=None, **kwargs):
        if raise_soma:
            raise ValueError("Unsupported SOMA object encoding version 1.1.0")
        return _Soma()

    def get_anndata(census, organism=None, obs_value_filter=None,
                    column_names=None, **kwargs):
        return adata

    mod.open_soma = open_soma
    mod.get_anndata = get_anndata
    mod.__version__ = "1.15.0-fake"
    return mod


# ---------------------------------------------------------------------------
# CLI argument
# ---------------------------------------------------------------------------


def test_census_version_argument_parsed(load_script):
    m00 = load_script("00_download_data.py")
    args = m00.build_parser().parse_args(["--census-version", "2023-07-25"])
    assert args.census_version == "2023-07-25"


def test_default_census_version(load_script):
    m00 = load_script("00_download_data.py")
    args = m00.build_parser().parse_args([])
    assert args.census_version == m00.DEFAULT_CENSUS_VERSION


def test_dry_run_accepts_census_version(load_script, harness, monkeypatch, tmp_path):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    monkeypatch.setattr(harness, "ALL_DIRS", [tmp_path / "d"])
    monkeypatch.setattr(harness, "MANIFEST_PATH", tmp_path / "m.json")
    m00 = load_script("00_download_data.py")
    assert m00.main(["--census-version", "2023-07-25"]) == 0


# ---------------------------------------------------------------------------
# Successful (mocked) download → manifest provenance
# ---------------------------------------------------------------------------


def test_download_reference_success_records_provenance(
    load_script, harness, monkeypatch, tmp_path, tiny_census_adata, capsys
):
    monkeypatch.setitem(sys.modules, "cellxgene_census",
                        _fake_census_module(adata=tiny_census_adata))
    monkeypatch.setattr(harness, "REFERENCE_H5AD", tmp_path / "ref.h5ad")
    m00 = load_script("00_download_data.py")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_reference(
        manifest, census_version="2024-07-01",
        clock=_counter_clock(), now=lambda: "2026-01-01T00:00:00Z",
    )
    assert ok is True
    assert (tmp_path / "ref.h5ad").exists()

    entry = manifest["datasets"]["reference"]
    assert entry["downloaded"] is True
    assert entry["census_version"] == "2024-07-01"
    assert entry["method"] == "cellxgene-census"
    assert entry["cellxgene_census_version"] != "not installed"
    assert entry["tiledbsoma_version"]  # recorded (real or 'not installed')
    assert entry["file_size_bytes"] > 0
    assert "elapsed_seconds" in entry
    assert entry["start_time"] and entry["end_time"]
    # Manifest still schema-valid.
    harness.validate_manifest(manifest)

    out = capsys.readouterr().out
    assert "[  5%]" in out and "[100%]" in out
    assert "cells" in out and "genes" in out


# ---------------------------------------------------------------------------
# SOMA-encoding incompatibility → actionable error, no silent continue
# ---------------------------------------------------------------------------


def test_download_reference_soma_error_actionable(
    load_script, harness, monkeypatch, tmp_path, capsys
):
    monkeypatch.setitem(sys.modules, "cellxgene_census",
                        _fake_census_module(raise_soma=True))
    monkeypatch.setattr(harness, "REFERENCE_H5AD", tmp_path / "ref.h5ad")
    m00 = load_script("00_download_data.py")

    manifest = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    ok = m00.download_reference(
        manifest, census_version="2024-07-01",
        clock=_counter_clock(), now=lambda: "2026-01-01T00:00:00Z",
    )
    assert ok is False
    assert not (tmp_path / "ref.h5ad").exists()  # never wrote a wrong dataset

    entry = manifest["datasets"]["reference"]
    # all candidate versions were attempted before giving up
    assert "2024-07-01" in entry["fallback_versions_attempted"]
    assert set(m00.CENSUS_FALLBACK_VERSIONS) <= set(entry["fallback_versions_attempted"])

    err = capsys.readouterr().err
    assert "cellxgene-census installed" in err
    assert "tiledbsoma installed" in err
    assert "census_version requested" in err
    assert "pip install -U cellxgene-census" in err
