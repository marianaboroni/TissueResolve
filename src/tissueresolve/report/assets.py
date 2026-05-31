"""
Asset/data loading helpers for results-directory-driven reports.

Reads whatever tables/figures exist in a results directory (tolerant of
missing pieces — they render as "not available" rather than failing) and lists
figure HTML files for embedding.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

import pandas as pd

__all__ = [
    "read_tsv", "read_json", "find_table", "list_figures", "collect_warnings",
]


def read_tsv(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep="\t", comment="#", index_col=0)
    except Exception:
        return None


def read_json(path: Path) -> Optional[Any]:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None


def find_table(results_dir: Path, *names: str) -> Optional[pd.DataFrame]:
    """Return the first existing table among *names*, searching dir and tables/."""
    results_dir = Path(results_dir)
    for base in (results_dir / "tables", results_dir):
        for n in names:
            p = base / n
            if p.exists():
                df = read_tsv(p)
                if df is not None:
                    return df
    return None


def list_figures(results_dir: Path) -> dict[str, Path]:
    """Map ``figure_stem -> html Path`` for every ``figures/*.html``."""
    figdir = Path(results_dir) / "figures"
    out: dict[str, Path] = {}
    if figdir.is_dir():
        for p in sorted(figdir.glob("*.html")):
            out[p.stem] = p
    return out


def collect_warnings(results_dir: Path) -> list[str]:
    """Gather warnings from warnings.json / warnings.tsv if present."""
    results_dir = Path(results_dir)
    warns: list[str] = []
    for base in (results_dir / "tables", results_dir):
        wj = base / "warnings.json"
        if wj.exists():
            data = read_json(wj)
            if isinstance(data, dict):
                for k, v in data.items():
                    if "recommend" in k.lower() and isinstance(v, list):
                        warns.extend(str(x) for x in v)
            elif isinstance(data, list):
                warns.extend(str(x) for x in data)
    return warns
