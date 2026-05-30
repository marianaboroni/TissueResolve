"""
Lightweight, dependency-free download/progress utilities for the harness.

Nothing here touches the network at import time, and the streaming HTTP helper
accepts an injectable ``opener`` so it can be unit-tested fully offline.  No
tqdm dependency — progress is printed with plain stdlib.
"""
from __future__ import annotations

import sys
import time
import urllib.request
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Callable, Optional

__all__ = [
    "format_bytes",
    "format_elapsed",
    "get_package_version",
    "ProgressLogger",
    "progress_line",
    "download_file_with_progress",
    "is_soma_encoding_error",
    "write_manifest",
]


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def format_bytes(n: float) -> str:
    """Human-readable byte size, e.g. ``317.5 MB``."""
    n = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024.0 or unit == "TB":
            return f"{int(n)} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} TB"  # pragma: no cover


def format_elapsed(seconds: float) -> str:
    """Seconds → ``HH:MM:SS``."""
    seconds = int(max(0, seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def get_package_version(name: str) -> str:
    """Installed version string, or ``'not installed'``."""
    import importlib.metadata as md

    try:
        return md.version(name)
    except Exception:
        return "not installed"


# ---------------------------------------------------------------------------
# Stage-level progress logger
# ---------------------------------------------------------------------------


class ProgressLogger:
    """Stage-level progress reporter: ``[ 35%] dataset: message (elapsed ...)``.

    Parameters
    ----------
    dataset:
        Label printed on each line.
    stream:
        Output stream (default stdout).
    clock:
        Monotonic clock callable (injectable for deterministic tests).
    """

    def __init__(self, dataset: str, *, stream=None,
                 clock: Callable[[], float] = time.monotonic) -> None:
        self.dataset = dataset
        self.stream = stream if stream is not None else sys.stdout
        self._clock = clock
        self._t0 = clock()

    def elapsed(self) -> float:
        return self._clock() - self._t0

    def stage(self, pct: int, message: str) -> str:
        line = (f"[{pct:3d}%] {self.dataset}: {message} "
                f"(elapsed {format_elapsed(self.elapsed())})")
        print(line, file=self.stream)
        return line

    def info(self, message: str) -> str:
        line = f"       {self.dataset}: {message}"
        print(line, file=self.stream)
        return line

    def done(self, message: str) -> str:
        return self.stage(100, message)


# ---------------------------------------------------------------------------
# Streaming HTTP download with byte-level progress
# ---------------------------------------------------------------------------


def progress_line(downloaded: int, total: Optional[int], elapsed: float) -> str:
    """One progress line.  With *total*::

        [42.3%] 317.5 MB / 750.0 MB elapsed 00:03:21

    Without a known total (no Content-Length)::

        [ --- ] 317.5 MB downloaded elapsed 00:03:21
    """
    if total:
        pct = 100.0 * downloaded / total
        return (f"[{pct:.1f}%] {format_bytes(downloaded)} / {format_bytes(total)} "
                f"elapsed {format_elapsed(elapsed)}")
    return (f"[ --- ] {format_bytes(downloaded)} downloaded "
            f"elapsed {format_elapsed(elapsed)}")


def _content_length(resp: Any) -> Optional[int]:
    val = None
    headers = getattr(resp, "headers", None)
    if headers is not None:
        try:
            val = headers.get("Content-Length")
        except Exception:
            val = None
    if val is None and hasattr(resp, "getheader"):
        try:
            val = resp.getheader("Content-Length")
        except Exception:
            val = None
    try:
        return int(val) if val is not None else None
    except (TypeError, ValueError):
        return None


def download_file_with_progress(
    url: str,
    dest: Path | str,
    *,
    opener: Callable[[str], Any] = urllib.request.urlopen,
    chunk_size: int = 1 << 20,
    stream=None,
    clock: Callable[[], float] = time.monotonic,
) -> dict:
    """Stream *url* to *dest*, printing byte-level progress.

    Uses ``Content-Length`` for a percentage when the server provides it;
    otherwise reports downloaded MB and elapsed time.

    Returns
    -------
    dict
        ``{bytes, total_bytes, content_length_available, elapsed_seconds}``.
    """
    stream = stream if stream is not None else sys.stdout
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    t0 = clock()
    resp = opener(url)
    cm = resp if hasattr(resp, "__enter__") else nullcontext(resp)
    downloaded = 0
    total: Optional[int] = None
    with cm as r:
        total = _content_length(r)
        with open(dest, "wb") as fh:
            while True:
                chunk = r.read(chunk_size)
                if not chunk:
                    break
                fh.write(chunk)
                downloaded += len(chunk)
                print(progress_line(downloaded, total, clock() - t0),
                      end="\r", file=stream)
    # Final newline-terminated summary line.
    print(progress_line(downloaded, total, clock() - t0), file=stream)
    return {
        "bytes": downloaded,
        "total_bytes": total,
        "content_length_available": total is not None,
        "elapsed_seconds": clock() - t0,
    }


# ---------------------------------------------------------------------------
# Census compatibility helpers
# ---------------------------------------------------------------------------


def is_soma_encoding_error(exc: BaseException) -> bool:
    """True when *exc* is the SOMA-encoding-version incompatibility error."""
    msg = str(exc)
    return (
        "Unsupported SOMA object encoding" in msg
        or "SOMA object encoding version" in msg
        or "Unsupported SOMA" in msg
    )


def write_manifest(manifest: dict, path: Path | str) -> Path:
    """Validate-and-write the manifest (delegates to the harness writer)."""
    import _harness as H

    return H.write_manifest(manifest, Path(path))
