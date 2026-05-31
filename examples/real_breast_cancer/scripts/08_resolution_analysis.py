#!/usr/bin/env python
"""
08 — Resolution analysis: recommend merge families + family-level predictions.

Reads existing outputs (no re-download, no re-fit):

    outputs/reference/breast_cancer_reference/   saved ReferenceSignature
    outputs/bulk/bulk_estimated_proportions.tsv  fine bulk estimates
    outputs/spatial/spatial_spot_proportions.tsv fine spatial estimates

computes the resolution-aware recommendation system, and writes to
``outputs/resolution/``:

    pairwise_separability.tsv      cell_type_families.tsv
    recommended_merges.tsv         unresolved_families.tsv
    bulk_family_proportions.tsv    spatial_family_proportions.tsv
    resolution_summary.md          resolution_mapping.json

Merging is **post-hoc** here: fine estimates are never overwritten — family
proportions are added as a safer interpretation, and the mapping + merge stage
are recorded.  ``--resolution-mode`` selects the behaviour (default: suggest).

Usage
-----
    python scripts/08_resolution_analysis.py
    python scripts/08_resolution_analysis.py --resolution-mode auto
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings as _warnings
from pathlib import Path

import pandas as pd

import _harness as H


def _read_tsv(p: Path):
    try:
        return pd.read_csv(p, sep="\t", comment="#", index_col=0)
    except Exception:
        return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--resolution-mode",
                        choices=["none", "suggest", "auto", "hierarchical"],
                        default="suggest",
                        help="How to handle non-separable types (default: suggest).")
    args = parser.parse_args(argv)

    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"ERROR: reference not found at {H.SAVED_REFERENCE_DIR}.  "
              "Run 01_prepare_reference.py first.", file=sys.stderr)
        return 1

    from tissueresolve.reference.hierarchy import aggregate_predictions_by_family
    from tissueresolve.reference.resolution import (
        ResolutionConfig, assign_resolution_families, build_resolution_report,
        recommend_cell_type_merges, summarize_resolution_report,
        write_recommended_merges)
    from tissueresolve.reference.separability import compute_separability
    from tissueresolve.results import ReferenceSignature

    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    cts = list(ref.cell_types)
    print(f"Loaded reference: {ref.n_genes} genes × {ref.n_cell_types} cell types")

    cfg = ResolutionConfig(resolution_mode=args.resolution_mode)
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")  # warning content is summarised below
        sep = compute_separability(ref, raise_on_critical=False)

    res = build_resolution_report(sep, cts, cfg)
    merges = recommend_cell_type_merges(sep, cts, cfg)
    mapping = assign_resolution_families(sep, cts, cfg)

    out = H.OUT_RESOLUTION_DIR
    out.mkdir(parents=True, exist_ok=True)

    # pairwise separability (raw metrics) + families + unresolved + merges
    res.pairwise.to_csv(out / "pairwise_separability.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "family": f.name, "resolvability": f.resolvability, "n_members": f.size,
        "mean_separability": round(f.mean_separability, 4),
        "members": "; ".join(f.members),
    } for f in res.families]).to_csv(out / "cell_type_families.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "family": f.name, "members": "; ".join(f.members),
        "resolvability": f.resolvability,
        "mean_separability": round(f.mean_separability, 4),
    } for f in res.unresolved_families]).to_csv(
        out / "unresolved_families.tsv", sep="\t", index=False)
    write_recommended_merges(merges, out / "recommended_merges.tsv")

    # Post-hoc family-level predictions (fine estimates untouched).
    n_family_outputs = 0
    bulk = _read_tsv(H.OUT_BULK_DIR / "bulk_estimated_proportions.tsv")
    if bulk is not None:
        fam = aggregate_predictions_by_family(bulk, mapping)
        fam.to_csv(out / "bulk_family_proportions.tsv", sep="\t")
        n_family_outputs += 1
        print(f"  bulk: {bulk.shape[1]} fine types → {fam.shape[1]} families")
    spatial = _read_tsv(H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv")
    if spatial is not None:
        fam = aggregate_predictions_by_family(spatial, mapping)
        fam.to_csv(out / "spatial_family_proportions.tsv", sep="\t")
        n_family_outputs += 1
        print(f"  spatial: {spatial.shape[1]} fine types → {fam.shape[1]} families")

    (out / "resolution_mapping.json").write_text(
        json.dumps({"resolution_mode": args.resolution_mode,
                    "merge_stage": "post_hoc_aggregation",
                    "fine_estimates_overwritten": False,
                    "mapping": mapping}, indent=2), encoding="utf-8")
    (out / "resolution_summary.md").write_text(
        summarize_resolution_report(sep, cts, cfg), encoding="utf-8")

    print(f"Recommended {len(merges)} merge family(ies); "
          f"{n_family_outputs} family-level prediction file(s) written.")
    print(f"Saved resolution outputs -> {out}")
    if len(merges):
        print("Recommended families:")
        for _, r in merges.iterrows():
            print(f"  {r['family_name']} ← {r['members']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
