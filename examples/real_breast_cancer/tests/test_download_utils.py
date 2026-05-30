"""
Offline tests for the download/progress utilities (_download_utils.py).

No network: download_file_with_progress is driven by an injected fake opener,
and the clock is injected for deterministic elapsed times.
"""
from __future__ import annotations

import io

import pytest


@pytest.fixture
def du():
    import _download_utils

    return _download_utils


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def test_format_bytes(du):
    assert du.format_bytes(0) == "0 B"
    assert du.format_bytes(512) == "512 B"
    assert du.format_bytes(1536) == "1.5 KB"
    assert du.format_bytes(int(317.5 * 1024 * 1024)) == "317.5 MB"
    assert du.format_bytes(2 * 1024 ** 3) == "2.0 GB"


def test_format_elapsed(du):
    assert du.format_elapsed(0) == "00:00:00"
    assert du.format_elapsed(201) == "00:03:21"
    assert du.format_elapsed(3 * 3600 + 21 * 60 + 5) == "03:21:05"
    assert du.format_elapsed(-5) == "00:00:00"


def test_get_package_version(du):
    assert du.get_package_version("numpy") != "not installed"
    assert du.get_package_version("definitely-not-a-real-package-xyz") == "not installed"


# ---------------------------------------------------------------------------
# ProgressLogger
# ---------------------------------------------------------------------------


def _counter_clock(step=1.0):
    state = {"t": 0.0}

    def clock():
        state["t"] += step
        return state["t"]

    return clock


def test_progress_logger_stage(du):
    buf = io.StringIO()
    log = du.ProgressLogger("reference", stream=buf, clock=_counter_clock())
    line = log.stage(35, "selecting breast-cancer cells")
    assert "[ 35%]" in line
    assert "reference:" in line
    assert "selecting breast-cancer cells" in line
    assert "elapsed" in line
    assert buf.getvalue().strip() == line


def test_progress_logger_done_is_100(du):
    buf = io.StringIO()
    log = du.ProgressLogger("reference", stream=buf, clock=_counter_clock())
    line = log.done("saved reference h5ad")
    assert "[100%]" in line


# ---------------------------------------------------------------------------
# progress_line
# ---------------------------------------------------------------------------


def test_progress_line_with_total(du):
    line = du.progress_line(int(317.5 * 1024 * 1024), int(750 * 1024 * 1024), 201)
    assert line.startswith("[42.3%]") or line.startswith("[42.")
    assert "317.5 MB" in line and "750.0 MB" in line
    assert "elapsed 00:03:21" in line


def test_progress_line_without_total(du):
    line = du.progress_line(100 * 1024 * 1024, None, 5)
    assert "downloaded" in line
    assert "100.0 MB" in line


# ---------------------------------------------------------------------------
# download_file_with_progress (mocked)
# ---------------------------------------------------------------------------


class _FakeResp:
    def __init__(self, data: bytes, content_length: bool = True):
        self._data = data
        self._pos = 0
        self.headers = {"Content-Length": str(len(data))} if content_length else {}

    def read(self, n: int) -> bytes:
        out = self._data[self._pos:self._pos + n]
        self._pos += len(out)
        return out

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _opener(data, content_length):
    def opener(url):
        return _FakeResp(data, content_length)

    return opener


def test_download_with_content_length(du, tmp_path):
    data = b"abc" * 1000  # 3000 bytes
    buf = io.StringIO()
    dest = tmp_path / "f.bin"
    result = du.download_file_with_progress(
        "http://x", dest, opener=_opener(data, True),
        chunk_size=1024, stream=buf, clock=_counter_clock(),
    )
    assert dest.read_bytes() == data
    assert result["bytes"] == 3000
    assert result["content_length_available"] is True
    assert result["total_bytes"] == 3000
    out = buf.getvalue()
    assert "%" in out and "/" in out


def test_download_without_content_length(du, tmp_path):
    data = b"z" * 2048
    buf = io.StringIO()
    dest = tmp_path / "f2.bin"
    result = du.download_file_with_progress(
        "http://x", dest, opener=_opener(data, False),
        chunk_size=512, stream=buf, clock=_counter_clock(),
    )
    assert dest.read_bytes() == data
    assert result["bytes"] == 2048
    assert result["content_length_available"] is False
    assert "downloaded" in buf.getvalue()


# ---------------------------------------------------------------------------
# SOMA error classification & manifest passthrough
# ---------------------------------------------------------------------------


def test_is_soma_encoding_error(du):
    assert du.is_soma_encoding_error(
        ValueError("Unsupported SOMA object encoding version 1.1.0"))
    assert du.is_soma_encoding_error(RuntimeError("Unsupported SOMA object encoding"))
    assert not du.is_soma_encoding_error(ValueError("some other error"))


def test_write_manifest_passthrough(du, harness, tmp_path):
    m = harness.default_manifest(access_date="2026-01-01T00:00:00Z")
    p = du.write_manifest(m, tmp_path / "manifest.json")
    assert p.exists()
    harness.validate_manifest(harness.read_manifest(p))
