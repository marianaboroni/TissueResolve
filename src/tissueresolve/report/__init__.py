"""
Report generation for TissueResolve.

Canonical path
--------------
``orchestration``
    The single entry point: builds ordered ``(title, body)`` sections from a
    results directory **or** an in-memory pipeline result and renders them
    through the unified single-page shell (``unified`` + ``components``).

``methods_text``
    Auto-generated methods text suitable for publication.

``sections`` / ``interpretation`` / ``assets``
    Results-directory section builders (content + data-driven prose).

``result_sections``
    In-memory result section builders (the counterpart of ``sections`` for a
    pipeline result), consumed by ``orchestration``.

``html``
    **Deprecated** compatibility shim — its public functions delegate to
    ``orchestration``; it re-exports the ``result_sections`` builders for
    backward compatibility only.
"""
from tissueresolve.report.orchestration import (
    generate_bulk_report,
    generate_report,
    generate_spatial_report,
)

__all__ = [
    "generate_bulk_report",
    "generate_spatial_report",
    "generate_report",
]
