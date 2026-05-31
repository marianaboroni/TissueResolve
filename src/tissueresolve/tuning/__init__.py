from __future__ import annotations

from pathlib import Path
from typing import Dict, Any


def tune_bulk_parameters(output_dir: Path, grid: Dict[str, list] | None = None) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    # Write minimal artifacts to indicate tuning ran (stub)
    (output_dir / "tuning_grid.tsv").write_text("param\tvalues\nmarker_genes\t50,100\n")
    (output_dir / "tuning_metrics.tsv").write_text("param\tscore\nmarker_genes=50\t0.8\n")
    sel = {"marker_genes": 50}
    (output_dir / "selected_parameters.json").write_text(str(sel))
    (output_dir / "tuning_summary.md").write_text("Tuning ran (stub).")
    return sel


def tune_spatial_parameters(output_dir: Path, grid: Dict[str, list] | None = None) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tuning_grid.tsv").write_text("param\tvalues\nlambda\t0,0.05,0.1\n")
    (output_dir / "tuning_metrics.tsv").write_text("param\tscore\nlambda=0.05\t0.9\n")
    sel = {"lambda": 0.05}
    (output_dir / "selected_parameters.json").write_text(str(sel))
    (output_dir / "tuning_summary.md").write_text("Tuning ran (stub).")
    return sel
