#!/usr/bin/env python
"""
Run both benchmarks and build the single summary report that links to them.

Usage
-----
    python benchmarks/run_all.py --dry-run
    python benchmarks/run_all.py --toy
    python benchmarks/run_all.py --use-existing-real-data --max-spots 600
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.bulk import run_bulk_benchmark as B  # noqa: E402
from benchmarks.spatial import run_spatial_benchmark as S  # noqa: E402
from benchmarks.shared.io import OUTPUTS_DIR  # noqa: E402
from benchmarks.shared import report as R  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--toy", action="store_true")
    ap.add_argument("--use-existing-real-data", action="store_true")
    ap.add_argument("--no-external", action="store_true")
    ap.add_argument("--include-imported", action="store_true",
                    help="Include externally-run tool predictions imported via "
                         "benchmarks/import_external_results.py (counted as "
                         "executed/imported, not exported-only).")
    ap.add_argument("--max-spots", type=int, default=600)
    args = ap.parse_args(argv)

    common = []
    if args.dry_run:
        common.append("--dry-run")
    if args.use_existing_real_data:
        common.append("--use-existing-real-data")
    if args.no_external:
        common.append("--no-external")
    if args.include_imported:
        common.append("--include-imported")

    B.main(common + (["--toy"] if args.toy else []))
    S.main(common + (["--toy"] if args.toy else []) +
           (["--max-spots", str(args.max_spots)] if not args.dry_run else []))

    if args.dry_run:
        print("Dry run complete; summary report not written.")
        return 0

    import json
    bulk_html = OUTPUTS_DIR / "bulk" / "bulk_benchmark_report.html"
    spatial_html = OUTPUTS_DIR / "spatial" / "spatial_benchmark_report.html"

    highlights = {}
    for modality in ("bulk", "spatial"):
        meta_p = OUTPUTS_DIR / modality / "benchmark_metadata.json"
        if meta_p.exists():
            try:
                meta = json.loads(meta_p.read_text())
                best = meta.get("best_methods", {})
                if modality == "bulk":
                    highlights["best bulk fine-level (Pearson)"] = best.get("best_bulk_fine_pearson", "—")
                    highlights["best bulk family-level (Pearson)"] = best.get("best_bulk_family_pearson", "—")
                highlights[f"fastest executed ({modality})"] = best.get("best_runtime", "—")
                highlights[f"non-TissueResolve methods executed ({modality})"] = \
                    meta.get("n_nontissueresolve_executed", "—")
            except Exception:
                pass
    summary = R.build_summary_report(
        bulk_html if bulk_html.exists() else None,
        spatial_html if spatial_html.exists() else None,
        highlights=highlights or {"note": "See linked reports for details."},
        out_path=OUTPUTS_DIR / "benchmark_summary_report.html")

    # index.html + README so users know where results are
    (OUTPUTS_DIR / "index.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>Benchmark outputs</title>"
        "<h1>TissueResolve benchmark outputs</h1><ul>"
        "<li><a href='benchmark_summary_report.html'><b>benchmark_summary_report.html</b> "
        "— start here (links to bulk + spatial)</a></li>"
        "<li><a href='bulk/bulk_benchmark_report.html'>bulk/bulk_benchmark_report.html</a></li>"
        "<li><a href='spatial/spatial_benchmark_report.html'>spatial/spatial_benchmark_report.html</a></li>"
        "</ul>", encoding="utf-8")
    (OUTPUTS_DIR / "README.md").write_text(
        "# Benchmark outputs\n\n**Start here:** `benchmark_summary_report.html`.\n\n"
        "- `bulk/`, `spatial/` — per-modality reports, executive tables, fair "
        "metrics, capability matrix, figures (each figure has a `.data.tsv`).\n"
        "- `*/external/` — imported external-tool predictions (if any).\n"
        "- `*/exports/` — inputs exported for tools not run locally.\n",
        encoding="utf-8")
    print(f"Wrote summary report → {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
