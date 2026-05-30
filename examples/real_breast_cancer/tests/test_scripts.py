"""
Offline tests for the numbered harness scripts.

These confirm the scripts import, the dry-run path is network-free, missing
inputs give clear non-zero exits, and the summary works after a partial run.
No real data is downloaded.
"""
from __future__ import annotations

import socket

import pytest

SCRIPT_FILES = [
    "00_download_data.py",
    "01_prepare_reference.py",
    "02_make_pseudobulk.py",
    "03_run_bulk_validation.py",
    "04_run_spatial_validation.py",
    "05_summarize_results.py",
]


@pytest.mark.parametrize("name", SCRIPT_FILES)
def test_script_imports(load_script, name):
    mod = load_script(name)
    assert hasattr(mod, "main")


# ---------------------------------------------------------------------------
# 00 — download dry run
# ---------------------------------------------------------------------------


def test_download_dry_run_writes_manifest(load_script, harness, monkeypatch, tmp_path):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    monkeypatch.setattr(harness, "ALL_DIRS", [tmp_path / "d"])
    monkeypatch.setattr(harness, "MANIFEST_PATH", tmp_path / "manifest.json")
    m00 = load_script("00_download_data.py")

    rc = m00.main([])
    assert rc == 0
    assert (tmp_path / "manifest.json").exists()
    harness.validate_manifest(harness.read_manifest(tmp_path / "manifest.json"))


def test_download_dry_run_does_not_download(load_script, harness, monkeypatch, tmp_path):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    monkeypatch.setattr(harness, "ALL_DIRS", [tmp_path / "d"])
    monkeypatch.setattr(harness, "MANIFEST_PATH", tmp_path / "manifest.json")
    m00 = load_script("00_download_data.py")

    def _boom(*a, **k):
        raise AssertionError("download must not run in dry mode")

    monkeypatch.setattr(m00, "download_reference", _boom)
    monkeypatch.setattr(m00, "download_spatial", _boom)
    assert m00.main([]) == 0


def test_download_dry_run_is_offline(load_script, harness, monkeypatch, tmp_path):
    monkeypatch.delenv(harness.REAL_DATA_ENV, raising=False)
    monkeypatch.setattr(harness, "ALL_DIRS", [tmp_path / "d"])
    monkeypatch.setattr(harness, "MANIFEST_PATH", tmp_path / "manifest.json")

    def _blocked(*a, **k):
        raise OSError("network blocked in tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)
    m00 = load_script("00_download_data.py")
    assert m00.main([]) == 0


# ---------------------------------------------------------------------------
# 04 — spatial: clear non-zero exit when inputs missing
# ---------------------------------------------------------------------------


def test_spatial_missing_input_exits_clearly(load_script, tmp_path, capsys):
    m04 = load_script("04_run_spatial_validation.py")
    rc = m04.main(["--spatial-h5ad", str(tmp_path / "absent.h5ad")])
    assert rc != 0
    err = capsys.readouterr().err
    assert "not found" in err.lower()


# ---------------------------------------------------------------------------
# 05 — summary works after a partial / empty run
# ---------------------------------------------------------------------------


def test_summary_handles_missing_outputs(load_script, harness, monkeypatch, tmp_path):
    for attr in ("OUT_REFERENCE_DIR", "OUT_BULK_DIR", "OUT_SPATIAL_DIR",
                 "OUT_SUMMARY_DIR", "DERIVED_DIR"):
        d = tmp_path / attr
        d.mkdir()
        monkeypatch.setattr(harness, attr, d)
    monkeypatch.setattr(harness, "ALL_DIRS", [tmp_path / attr for attr in
                                              ("OUT_REFERENCE_DIR", "OUT_SUMMARY_DIR")])
    monkeypatch.setattr(harness, "PSEUDOBULK_TRUE_PROPS", tmp_path / "none.tsv")

    m05 = load_script("05_summarize_results.py")
    rc = m05.main([])
    assert rc == 0
    md = (tmp_path / "OUT_SUMMARY_DIR" / "validation_summary.md").read_text()
    assert "not available" in md
    assert (tmp_path / "OUT_SUMMARY_DIR" / "validation_summary.json").exists()
