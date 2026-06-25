"""
Deprecated compatibility shim for report generation.

.. deprecated::
    This module's public functions delegate to the canonical
    :mod:`tissueresolve.report.orchestration` layer (which renders every report
    through the unified single-page shell).  New code should call
    :func:`tissueresolve.report.generate_report`.

The in-memory result section builders that used to live here now live in
:mod:`tissueresolve.report.result_sections`; they are re-exported here only for
backward compatibility.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Union

# Backward-compatibility re-exports (the real implementations moved out).
from tissueresolve.report.result_sections import (  # noqa: F401
    bulk_result_sections,
    spatial_result_sections,
    read_run_metadata as assets_read_metadata,
)

__all__ = [
    "generate_bulk_report", "generate_spatial_report", "generate_report",
    "bulk_result_sections", "spatial_result_sections", "assets_read_metadata",
]


def generate_bulk_report(source: Any, output_path: Union[str, Path, None] = None,
                         *, out: Union[str, Path, None] = None, **kwargs: Any) -> Path:
    """Deprecated shim — delegates to the canonical orchestration layer."""
    from tissueresolve.report.orchestration import generate_report as _gen
    return _gen("bulk", source, output_path, out=out, **kwargs)


def generate_spatial_report(source: Any, output_path: Union[str, Path, None] = None,
                            *, out: Union[str, Path, None] = None, **kwargs: Any) -> Path:
    """Deprecated shim — delegates to the canonical orchestration layer."""
    from tissueresolve.report.orchestration import generate_report as _gen
    return _gen("spatial", source, output_path, out=out, **kwargs)


def generate_report(modality: str, source: Any,
                    output_path: Union[str, Path, None] = None, *,
                    out: Union[str, Path, None] = None, **kwargs: Any) -> Path:
    """Deprecated shim — delegates to the canonical orchestration layer."""
    from tissueresolve.report.orchestration import generate_report as _gen
    return _gen(modality, source, output_path, out=out, **kwargs)
