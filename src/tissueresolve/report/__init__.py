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

``html``
    **Deprecated** compatibility shim — its public functions delegate to
    ``orchestration``.  Still hosts the in-memory result section builders
    (``bulk_result_sections`` / ``spatial_result_sections``) and shared HTML
    helpers consumed by ``orchestration``.
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
