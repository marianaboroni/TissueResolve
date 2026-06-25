"""
Hierarchical reference construction and persistence.

Builds and saves the broad/fine reference bundle for hierarchical mode::

    reference/
    ├── fine_reference/             saved ReferenceSignature (subtypes)
    ├── family_reference/           saved ReferenceSignature (broad families)
    ├── cell_type_hierarchy.tsv     fine → broad mapping
    ├── family_reference_summary.tsv
    ├── broad_cell_type_counts.tsv
    ├── fine_cell_type_counts.tsv
    ├── hierarchy_summary.tsv       per-(broad, fine) cells + n_fine_in_family
    └── hierarchy_metadata.json     provenance + per-family subtype counts

The fine reference is preserved unchanged; the family reference is the
n_cells-weighted aggregate (gene order preserved, donor_cv dropped).  Families
with a single subtype are reported (and warned about) but never silently
merged.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from tissueresolve.reference.hierarchy import (
    aggregate_reference_by_family,
    build_cell_type_hierarchy,
    save_hierarchy_mapping,
)
from tissueresolve.results import ReferenceSignature

__all__ = ["build_hierarchical_reference_bundle", "save_hierarchical_reference_bundle"]


def build_hierarchical_reference_bundle(
    fine_ref: ReferenceSignature,
    mapping: Optional[dict[str, str]] = None,
) -> dict[str, Any]:
    """Build the broad/fine reference pair and hierarchy metadata.

    Returns a dict with ``fine_reference``, ``family_reference``, ``mapping``,
    ``groups`` (family → members), and ``summary`` (per-family subtype counts).
    """
    mapping = build_cell_type_hierarchy(list(fine_ref.cell_types), mapping)
    family_ref = aggregate_reference_by_family(fine_ref, mapping)

    groups: dict[str, list[str]] = {}
    for ct in fine_ref.cell_types:
        groups.setdefault(mapping[ct], []).append(ct)

    rows = []
    for fam in family_ref.cell_types:
        members = groups.get(fam, [])
        rows.append({
            "broad_cell_type": fam,
            "n_fine_subtypes": len(members),
            "n_cells": int(family_ref.n_cells_per_type.get(fam, 0)),
            "fine_subtypes": "; ".join(members),
            "single_subtype_only": len(members) <= 1,
        })
    summary = pd.DataFrame(rows).sort_values(
        "n_fine_subtypes", ascending=False).reset_index(drop=True)

    return {
        "fine_reference": fine_ref,
        "family_reference": family_ref,
        "mapping": mapping,
        "groups": groups,
        "summary": summary,
    }


def save_hierarchical_reference_bundle(
    fine_ref: ReferenceSignature,
    mapping: Optional[dict[str, str]],
    out_dir: Path | str,
    *,
    obs: Optional[pd.DataFrame] = None,
    broad_col: Optional[str] = None,
    fine_col: Optional[str] = None,
) -> dict[str, Path]:
    """Build and persist the full hierarchical reference bundle under *out_dir*.

    When *obs* (and the broad/fine column names) are supplied, the per-cell
    annotation summary (``hierarchy_summary.tsv`` plus broad/fine count tables)
    is derived from the real annotations; otherwise counts come from the
    reference's ``n_cells_per_type``.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    bundle = build_hierarchical_reference_bundle(fine_ref, mapping)
    written: dict[str, Path] = {}

    bundle["fine_reference"].save(out / "fine_reference")
    written["fine_reference"] = out / "fine_reference"
    bundle["family_reference"].save(out / "family_reference")
    written["family_reference"] = out / "family_reference"

    written["cell_type_hierarchy.tsv"] = save_hierarchy_mapping(
        bundle["mapping"], out / "cell_type_hierarchy.tsv")

    summary_path = out / "family_reference_summary.tsv"
    bundle["summary"].to_csv(summary_path, sep="\t", index=False)
    written["family_reference_summary.tsv"] = summary_path

    # broad / fine cell counts
    if obs is not None and broad_col and fine_col:
        from tissueresolve.io.validation import summarize_hierarchical_annotations
        hsum = summarize_hierarchical_annotations(obs, broad_col, fine_col)
    else:
        hsum = pd.DataFrame([
            {"broad_cell_type": bundle["mapping"][ct], "fine_cell_type": ct,
             "n_cells": int(fine_ref.n_cells_per_type.get(ct, 0))}
            for ct in fine_ref.cell_types
        ])
        nfine = hsum.groupby("broad_cell_type")["fine_cell_type"].transform("nunique")
        hsum["n_fine_in_family"] = nfine
        hsum = hsum.sort_values(["broad_cell_type", "n_cells"],
                                ascending=[True, False]).reset_index(drop=True)

    hsum_path = out / "hierarchy_summary.tsv"
    hsum.to_csv(hsum_path, sep="\t", index=False)
    written["hierarchy_summary.tsv"] = hsum_path

    broad_counts = (hsum.groupby("broad_cell_type")["n_cells"].sum()
                    .rename("n_cells").reset_index())
    bc_path = out / "broad_cell_type_counts.tsv"
    broad_counts.to_csv(bc_path, sep="\t", index=False)
    written["broad_cell_type_counts.tsv"] = bc_path

    fine_counts = hsum[["broad_cell_type", "fine_cell_type", "n_cells"]].copy()
    fc_path = out / "fine_cell_type_counts.tsv"
    fine_counts.to_csv(fc_path, sep="\t", index=False)
    written["fine_cell_type_counts.tsv"] = fc_path

    singletons = [f for f, m in bundle["groups"].items() if len(m) <= 1]
    if singletons:
        warnings.warn(
            f"{len(singletons)} broad family/ies have a single subtype: "
            f"{singletons[:8]}.  These are reported at the broad level only.",
            stacklevel=2,
        )
    meta = {
        "tissueresolve_object": "HierarchicalReference",
        "n_fine_cell_types": fine_ref.n_cell_types,
        "n_broad_families": bundle["family_reference"].n_cell_types,
        "n_genes": fine_ref.n_genes,
        "broad_cell_type_col": broad_col,
        "fine_cell_type_col": fine_col,
        "subtypes_per_family": {f: len(m) for f, m in bundle["groups"].items()},
        "single_subtype_families": singletons,
    }
    meta_path = out / "hierarchy_metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    written["hierarchy_metadata.json"] = meta_path
    return written
