#!/usr/bin/env python
"""
05 — Summarize the validation run into a Markdown + JSON report.

Reads whatever outputs exist (reference summary, pseudobulk, bulk metrics,
spatial outputs) and produces a concise summary.  Missing pieces are reported
as "not available" rather than causing a failure, so the summary works after a
partial run.

Usage
-----
    python scripts/05_summarize_results.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

import _harness as H


def _read_tsv(path: Path):
    return pd.read_csv(path, sep="\t", comment="#", index_col=0) if path.exists() else None


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def collect_summary() -> dict[str, Any]:
    s: dict[str, Any] = {"available": {}}

    ref_sum = _read_tsv(H.OUT_REFERENCE_DIR / "reference_summary.tsv")
    if ref_sum is not None:
        d = ref_sum["value"].to_dict() if "value" in ref_sum.columns else {}
        s["reference"] = {
            "n_cells": d.get("n_cells"),
            "n_genes": d.get("n_genes"),
            "n_cell_types": d.get("n_cell_types"),
            "cell_type_col": d.get("cell_type_col"),
        }
        s["available"]["reference"] = True
    else:
        s["available"]["reference"] = False

    true_props = _read_tsv(H.PSEUDOBULK_TRUE_PROPS)
    s["pseudobulk"] = {"n_samples": int(true_props.shape[0]) if true_props is not None else 0}
    s["available"]["pseudobulk"] = true_props is not None

    bulk_metrics = _read_tsv(H.OUT_BULK_DIR / "bulk_validation_metrics.tsv")
    bulk_warn = _read_json(H.OUT_BULK_DIR / "bulk_warnings.json")
    if bulk_metrics is not None:
        m = bulk_metrics["value"].to_dict() if "value" in bulk_metrics.columns else {}
        s["bulk"] = {k: m.get(k) for k in
                     ("pearson", "spearman", "rmse", "mae", "signed_bias")}
        if bulk_warn:
            s["bulk"]["gene_overlap"] = bulk_warn.get("gene_overlap")
            s["bulk"]["qc_recommendations"] = bulk_warn.get("qc_recommendations", [])
        s["available"]["bulk"] = True
    else:
        s["available"]["bulk"] = False

    sp_props = _read_tsv(H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv")
    sp_warn = _read_json(H.OUT_SPATIAL_DIR / "spatial_warnings.json")
    morans = _read_tsv(H.OUT_SPATIAL_DIR / "morans_i.tsv")
    if sp_props is not None:
        s["spatial"] = {
            "n_spots": int(sp_props.shape[0]),
            "converged": (sp_warn or {}).get("converged"),
            "lambda_spatial": (sp_warn or {}).get("lambda_spatial"),
            "gene_overlap": (sp_warn or {}).get("gene_overlap"),
            "morans_i_mean": float(morans["morans_i"].mean()) if morans is not None else None,
        }
        s["available"]["spatial"] = True
    else:
        s["available"]["spatial"] = False

    return s


def to_markdown(s: dict[str, Any]) -> str:
    lines = ["# TissueResolve real-data validation summary", ""]

    lines.append("## Reference")
    if s["available"].get("reference"):
        r = s["reference"]
        lines += [
            f"- cells: {r.get('n_cells')}",
            f"- genes: {r.get('n_genes')}",
            f"- cell types: {r.get('n_cell_types')}",
            f"- annotation column: {r.get('cell_type_col')}",
        ]
    else:
        lines.append("- not available (run 01_prepare_reference.py)")

    lines += ["", "## Bulk (pseudobulk) validation"]
    lines.append(f"- pseudobulk samples: {s['pseudobulk']['n_samples']}")
    if s["available"].get("bulk"):
        b = s["bulk"]
        lines += [
            f"- gene overlap (ref vs bulk): {b.get('gene_overlap')}",
            f"- Pearson: {b.get('pearson')}",
            f"- Spearman: {b.get('spearman')}",
            f"- RMSE: {b.get('rmse')}",
            f"- MAE: {b.get('mae')}",
            f"- signed bias: {b.get('signed_bias')}",
            f"- QC recommendations: {b.get('qc_recommendations')}",
        ]
        lines.append("- NOTE: bulk estimates are mRNA proportions, not cell fractions.")
    else:
        lines.append("- not available (run 03_run_bulk_validation.py)")

    lines += ["", "## Spatial validation"]
    if s["available"].get("spatial"):
        sp = s["spatial"]
        lines += [
            f"- spots: {sp.get('n_spots')}",
            f"- gene overlap (ref vs spatial): {sp.get('gene_overlap')}",
            f"- converged: {sp.get('converged')}",
            f"- lambda_spatial (smoothing, recorded): {sp.get('lambda_spatial')}",
            f"- mean Moran's I: {sp.get('morans_i_mean')}",
        ]
        lines.append("- NOTE: spatial estimates are spot-level RNA-derived "
                     "composition, not single-cell counts.")
    else:
        lines.append("- not available (run 04_run_spatial_validation.py)")

    lines += ["", "_Synthetic-free real-data validation; offline tests never "
              "download data._"]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    H.ensure_dirs()
    s = collect_summary()
    H.write_json(s, H.OUT_SUMMARY_DIR / "validation_summary.json")
    (H.OUT_SUMMARY_DIR / "validation_summary.md").write_text(
        to_markdown(s), encoding="utf-8")
    print(f"Saved summary -> {H.OUT_SUMMARY_DIR}")
    print(to_markdown(s))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
