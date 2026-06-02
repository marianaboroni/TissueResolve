"""
Canonical report orchestration for TissueResolve.

This is the single entry point for HTML report generation.  It builds an
ordered list of ``(title, body_html)`` sections from **either** a results
directory **or** an in-memory pipeline result, and renders them through the
unified single-page shell (``report/unified.py`` + ``report/components.py``):
one ``report.html`` with a sticky sidebar, consistent design system, and
QC-first ordering.

Content builders live in the dedicated layers:

* results-directory sections → :mod:`tissueresolve.report.sections`
  (which uses :mod:`tissueresolve.report.interpretation` + ``assets``);
* in-memory result sections → :func:`tissueresolve.report.html.bulk_result_sections`
  / :func:`~tissueresolve.report.html.spatial_result_sections`.

``report/html.py``'s public functions are thin deprecation shims that delegate
here, so the project has a single canonical rendering path.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Optional, Sequence, Union

from tissueresolve.report.unified import Section, build_unified_report

__all__ = [
    "render_sections", "generate_report",
    "generate_bulk_report", "generate_spatial_report",
]

_BULK_TITLE = "TissueResolve — Bulk deconvolution report"
_SPATIAL_TITLE = "TissueResolve — Spatial deconvolution report"


def _anchor(title: str, index: int) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", str(title).lower()).strip("-") or "section"
    return f"{slug}-{index}"


def render_sections(out_path: Union[str, Path], page_title: str,
                    secs: Sequence[tuple[str, str]], *,
                    subtitle: str = "") -> Path:
    """Render an ordered ``(title, body_html)`` list to the unified single-page
    report and return the written path.  This is the one canonical renderer."""
    sections = [Section(_anchor(t, i), t, b) for i, (t, b) in enumerate(secs)]
    return build_unified_report(out_path, sections, title=page_title,
                                subtitle=subtitle)


def generate_report(
    modality: str,
    source: Any,
    output_path: Union[str, Path, None] = None,
    *,
    out: Union[str, Path, None] = None,
    run_metadata: Optional[dict] = None,
    warnings: Optional[list] = None,
    separability: Any = None,
    figures: Optional[Sequence[Any]] = None,
    output_files: Optional[Sequence[Any]] = None,
    methods: Optional[str] = None,
    title: Optional[str] = None,
) -> Path:
    """Generate a bulk or spatial HTML report (the canonical entry point).

    *source* is either a results-directory path/str or an in-memory pipeline
    result (``BulkPipelineResult`` / ``SpatialPipelineResult``).
    """
    if modality not in ("bulk", "spatial"):
        raise ValueError(f"modality must be 'bulk' or 'spatial', got {modality!r}.")
    page_title = title or (_BULK_TITLE if modality == "bulk" else _SPATIAL_TITLE)
    out_path = output_path or out

    if isinstance(source, (str, Path)):
        from tissueresolve.report import sections as S
        from tissueresolve.report.html import assets_read_metadata
        results_dir = Path(source)
        meta = run_metadata if run_metadata is not None else assets_read_metadata(results_dir)
        builder = S.bulk_sections if modality == "bulk" else S.spatial_sections
        secs = builder(results_dir, run_metadata=meta, warnings=warnings)
        out_path = Path(out_path) if out_path else (results_dir / "report.html")
    else:
        from tissueresolve.report import html as H
        if out_path is None:
            raise ValueError(
                "output_path/out is required when generating a report from an "
                "in-memory result.")
        out_path = Path(out_path)
        builder = (H.bulk_result_sections if modality == "bulk"
                   else H.spatial_result_sections)
        secs = builder(source, separability=separability, figures=figures,
                       output_files=output_files, methods=methods,
                       fig_dir=out_path.parent)

    return render_sections(out_path, page_title, secs)


def generate_bulk_report(source: Any, output_path: Union[str, Path, None] = None,
                         *, out: Union[str, Path, None] = None, **kwargs: Any) -> Path:
    return generate_report("bulk", source, output_path, out=out, **kwargs)


def generate_spatial_report(source: Any, output_path: Union[str, Path, None] = None,
                            *, out: Union[str, Path, None] = None, **kwargs: Any) -> Path:
    return generate_report("spatial", source, output_path, out=out, **kwargs)
