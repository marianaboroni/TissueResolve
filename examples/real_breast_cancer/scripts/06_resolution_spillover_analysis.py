#!/usr/bin/env python
"""
06 — Resolution & spillover analysis (optional).

Turns the reference's separability diagnostics into an actionable,
resolution-aware report:

* group confusable cell types into families and classify resolvability,
* estimate a spillover matrix by deconvolving pure single-type pseudobulks
  (known truth = identity), and summarise per-type spillover risk + the main
  leaking partner,
* recommend merges and list unresolved families,
* annotate the bulk / spatial validation estimates (if present) with their
  resolvability and recommended interpretation level.

This is analysis tooling only; it does not modify the core algorithms or the
saved estimates.  Outputs are written under ``outputs/resolution/``.

Usage
-----
    python scripts/06_resolution_spillover_analysis.py
    python scripts/06_resolution_spillover_analysis.py --spillover-method expression
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

import _harness as H


def run_analysis(ref, *, spillover_method: str = "deconvolution"):
    """Compute the resolution + spillover artefacts for *ref*.  Returns a dict."""
    from tissueresolve.benchmark.spillover import (
        compute_spillover_risk,
        estimate_spillover_matrix,
        expression_spillover_proxy,
        summarize_spillover_partners,
    )
    from tissueresolve.reference.resolution import (
        build_resolution_report,
        suggest_merges,
    )
    from tissueresolve.reference.separability import compute_separability

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # separability warnings are summarised below
        sep = compute_separability(ref, warn_threshold=0.90, raise_on_critical=False)

    res = build_resolution_report(sep, list(ref.cell_types))
    merges = suggest_merges(sep, list(ref.cell_types))

    if spillover_method == "expression":
        spillover = expression_spillover_proxy(ref)
    else:
        spillover = estimate_spillover_matrix(ref)
    risk = compute_spillover_risk(spillover)
    partners = summarize_spillover_partners(spillover, min_fraction=0.05)

    return {
        "separability": sep, "resolution": res, "merges": merges,
        "spillover": spillover, "risk": risk, "partners": partners,
        "spillover_method": spillover_method,
    }


def _summary_md(ref, art) -> str:
    res = art["resolution"]
    counts = res.pairwise["resolvability"].value_counts().to_dict() if not res.pairwise.empty else {}
    unresolved = res.unresolved_families
    top_risk = art["risk"].sort_values("spillover_risk", ascending=False).head(10)
    lines = [
        "# Resolution & spillover analysis", "",
        f"- cell types: {ref.n_cell_types}",
        f"- families: {res.n_families} ({len(unresolved)} unresolved/poorly-resolved)",
        f"- spillover method: {art['spillover_method']}",
        "",
        "## Pair resolvability classes",
    ]
    for cls in ("resolved", "partially_resolved", "poorly_resolved", "unresolved"):
        lines.append(f"- {cls}: {counts.get(cls, 0)}")
    lines += ["", "## Highest spillover risk (top 10)"]
    for ct, row in top_risk.iterrows():
        lines.append(f"- {ct}: risk={row['spillover_risk']:.3f}, "
                     f"main partner={row['main_partner']} "
                     f"({row['main_partner_fraction']:.3f})")
    lines += ["", "## Unresolved / poorly-resolved families"]
    if unresolved:
        for f in unresolved:
            lines.append(f"- [{f.resolvability}] {f.name}: {', '.join(f.members)} "
                         f"(mean separability {f.mean_separability:.3f})")
        lines.append("")
        lines.append("Recommendation: report these families at the broad level, or "
                     "enable unresolved mode (allow_unresolved=True) so their mass "
                     "is reported as `unresolved_<family>` rather than confidently "
                     "split into subtypes.")
    else:
        lines.append("- none")
    lines += ["", "_Estimates are RNA-derived proportions (bulk: mRNA proportions; "
              "spatial: spot-level composition), never absolute cell counts._"]
    return "\n".join(lines) + "\n"


def _maybe_annotate_outputs(ref, art):
    """Annotate bulk/spatial estimates with resolvability + spillover, if present."""
    from tissueresolve.reference.resolution import annotate_resolution

    out = {}
    bulk_path = H.OUT_BULK_DIR / "bulk_estimated_proportions.tsv"
    if bulk_path.exists():
        est = pd.read_csv(bulk_path, sep="\t", comment="#", index_col=0)
        mean_est = est.mean(axis=0)
        out["bulk"] = annotate_resolution(mean_est, art["resolution"],
                                          spillover_risk=art["risk"], modality="bulk")
    sp_path = H.OUT_SPATIAL_DIR / "spatial_spot_proportions.tsv"
    if sp_path.exists():
        est = pd.read_csv(sp_path, sep="\t", comment="#", index_col=0)
        mean_est = est.mean(axis=0)
        out["spatial"] = annotate_resolution(mean_est, art["resolution"],
                                             spillover_risk=art["risk"], modality="spatial")
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spillover-method", choices=["deconvolution", "expression"],
                        default="deconvolution",
                        help="How to estimate spillover (default: deconvolution).")
    args = parser.parse_args(argv)

    if not H.SAVED_REFERENCE_DIR.exists():
        print(f"ERROR: reference not found at {H.SAVED_REFERENCE_DIR}.  "
              "Run 01_prepare_reference.py first.", file=sys.stderr)
        return 1

    from tissueresolve.results import ReferenceSignature

    ref = ReferenceSignature.load(H.SAVED_REFERENCE_DIR)
    print(f"Loaded reference: {ref.n_genes} genes × {ref.n_cell_types} cell types")
    print(f"Computing resolution & spillover (method={args.spillover_method}) …")
    art = run_analysis(ref, spillover_method=args.spillover_method)

    H.ensure_dirs()
    out = H.OUT_RESOLUTION_DIR
    art["resolution"].save(out)  # cell_type_families.tsv, pairwise_resolvability.tsv
    art["spillover"].to_csv(out / "spillover_matrix.tsv", sep="\t")
    art["risk"].to_csv(out / "spillover_risk_by_celltype.tsv", sep="\t")
    art["partners"].to_csv(out / "pairwise_spillover_report.tsv", sep="\t", index=False)
    art["merges"].to_csv(out / "recommended_merges.tsv", sep="\t", index=False)

    unresolved_rows = [{
        "family": f.name, "members": "; ".join(f.members),
        "resolvability": f.resolvability,
        "mean_separability": round(f.mean_separability, 4),
    } for f in art["resolution"].unresolved_families]
    pd.DataFrame(unresolved_rows, columns=[
        "family", "members", "resolvability", "mean_separability"
    ]).to_csv(out / "unresolved_families.tsv", sep="\t", index=False)

    annotations = _maybe_annotate_outputs(ref, art)
    for modality, ann in annotations.items():
        ann.to_csv(out / f"{modality}_resolution_annotation.tsv", sep="\t")

    (out / "resolution_summary.md").write_text(_summary_md(ref, art), encoding="utf-8")
    print(f"Saved resolution outputs -> {out}")
    print(art["resolution"].summary_str())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
