"""I/O helpers for the benchmark harness (TSV/JSON + plot source data)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

BENCH_DIR = Path(__file__).resolve().parent.parent
OUTPUTS_DIR = BENCH_DIR / "outputs"


def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_tsv(df: pd.DataFrame, path: Path, comment: list[str] | None = None) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as fh:
        for c in comment or []:
            fh.write(f"# {c}\n")
        df.to_csv(fh, sep="\t")
    return path


def write_json(obj: Any, path: Path) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")
    return path


def save_plot_data(df: pd.DataFrame, html_path: Path) -> Path:
    """Save the source data for a figure next to it as ``<name>.data.tsv``."""
    html_path = Path(html_path)
    data_path = html_path.with_suffix(".data.tsv")
    return write_tsv(df, data_path)
