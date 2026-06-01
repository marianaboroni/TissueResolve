"""
Figure manifest for reproducible reports.

Accumulates one record per figure (title, paths, caption, methodology,
variables, status) and writes ``figures/figure_manifest.tsv`` so every figure in
a report is traceable to its source data and generation parameters.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

__all__ = ["FigureRecord", "FigureManifest"]

MANIFEST_COLUMNS = [
    "figure_id", "section", "title", "html_path", "png_path", "svg_path",
    "pdf_path", "source_data_path", "caption_path", "caption", "methodology",
    "variables_defined", "status", "reason_if_missing",
]

# status vocabulary used by the report harness
STATUS_GENERATED = "generated"
STATUS_MISSING_DATA = "missing_data"
STATUS_SKIPPED = "skipped"


@dataclass
class FigureRecord:
    figure_id: str
    section: str
    title: str
    html_path: str = ""
    png_path: str = ""
    svg_path: str = ""
    pdf_path: str = ""
    source_data_path: str = ""
    caption_path: str = ""
    caption: str = ""
    methodology: str = ""
    variables_defined: str = ""
    status: str = "generated"
    reason_if_missing: str = ""

    def row(self) -> dict:
        return {c: getattr(self, c) for c in MANIFEST_COLUMNS}


@dataclass
class FigureManifest:
    records: list = field(default_factory=list)

    def add(self, rec: FigureRecord) -> FigureRecord:
        self.records.append(rec)
        return rec

    def to_frame(self) -> pd.DataFrame:
        if not self.records:
            return pd.DataFrame(columns=MANIFEST_COLUMNS)
        return pd.DataFrame([r.row() for r in self.records], columns=MANIFEST_COLUMNS)

    def save(self, figures_dir: Path | str, name: str = "figure_manifest.tsv") -> Path:
        figures_dir = Path(figures_dir)
        figures_dir.mkdir(parents=True, exist_ok=True)
        path = figures_dir / name
        self.to_frame().to_csv(path, sep="\t", index=False)
        return path
