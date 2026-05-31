"""
Report generation for TissueResolve.

Submodules (implemented in Stage 5):

``html``
    Unified HTMLReport with bulk and spatial variants.
    Self-contained HTML output with embedded CSS.

``methods_text``
    Auto-generated methods section text suitable for publication.
    Includes algorithm description, parameter values, and software versions.

``templates``, ``sections``, ``assets``
    HTML skeleton, section builders, and results-directory loaders for the
    publication-layer reports.
"""
from tissueresolve.report.html import (
    generate_bulk_report,
    generate_report,
    generate_spatial_report,
)

__all__ = [
    "generate_bulk_report",
    "generate_spatial_report",
    "generate_report",
]
