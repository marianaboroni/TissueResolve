"""
Three-level annotation hierarchy: broad family → cell type → cell state.

Extends the existing two-level (fine → broad) model with an explicit middle
"cell type" level and a lower "cell state / subpopulation" level, inspired by
BayesPrism's state↔type distinction (without copying its Bayesian machinery).

Backward-compatible: a plain ``{fine → broad}`` mapping is lifted to a
three-level hierarchy with ``cell_type == fine`` and **no** state level
(``has_states=False``), so existing two-level callers are unaffected.

Validation never silently drops labels: every state must map to exactly one cell
type and every cell type to exactly one broad family; low-support states and
single-donor cell types are *warned about*, not removed.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Optional

import pandas as pd

__all__ = [
    "ThreeLevelHierarchy",
    "build_three_level_hierarchy",
    "build_three_level_from_two_level",
    "build_three_level_from_obs",
    "validate_three_level_hierarchy",
    "summarize_three_level_hierarchy",
    "write_three_level_outputs",
]


@dataclass
class ThreeLevelHierarchy:
    """A validated broad→cell_type→state hierarchy.

    ``state_to_celltype`` and ``celltype_to_broad`` are the two mapping dicts;
    ``has_states`` is False when the hierarchy was lifted from a two-level
    mapping (states are then identical to cell types).
    """
    state_to_celltype: dict[str, str]
    celltype_to_broad: dict[str, str]
    has_states: bool = True

    def states(self) -> list[str]:
        return sorted(self.state_to_celltype)

    def cell_types(self) -> list[str]:
        return sorted(set(self.celltype_to_broad))

    def broad_families(self) -> list[str]:
        return sorted(set(self.celltype_to_broad.values()))

    def state_to_broad(self) -> dict[str, str]:
        return {s: self.celltype_to_broad.get(ct, ct)
                for s, ct in self.state_to_celltype.items()}

    def states_of_celltype(self, cell_type: str) -> list[str]:
        return sorted(s for s, ct in self.state_to_celltype.items()
                      if ct == cell_type)

    def celltypes_of_broad(self, broad: str) -> list[str]:
        return sorted(ct for ct, b in self.celltype_to_broad.items() if b == broad)

    def frame(self) -> pd.DataFrame:
        rows = [{"state": s, "cell_type": ct,
                 "broad_family": self.celltype_to_broad.get(ct, ct)}
                for s, ct in sorted(self.state_to_celltype.items())]
        return pd.DataFrame(rows, columns=["state", "cell_type", "broad_family"])


def build_three_level_hierarchy(
    state_to_celltype: Mapping[str, str],
    celltype_to_broad: Mapping[str, str],
) -> ThreeLevelHierarchy:
    """Build + structurally validate a three-level hierarchy.

    Raises ``ValueError`` if a state's cell type is absent from
    ``celltype_to_broad`` (a missing mid-level label — never silently inferred).
    """
    s2c = {str(k): str(v) for k, v in state_to_celltype.items()}
    c2b = {str(k): str(v) for k, v in celltype_to_broad.items()}
    missing = sorted({ct for ct in s2c.values() if ct not in c2b})
    if missing:
        raise ValueError(
            "every cell type referenced by a state must map to a broad family; "
            f"missing cell_type→broad mappings for: {missing}")
    return ThreeLevelHierarchy(s2c, c2b, has_states=True)


def build_three_level_from_two_level(
    fine_to_broad: Mapping[str, str],
) -> ThreeLevelHierarchy:
    """Lift a two-level ``{fine → broad}`` mapping (cell_type = fine, no state)."""
    c2b = {str(k): str(v) for k, v in fine_to_broad.items()}
    s2c = {ct: ct for ct in c2b}            # states identical to cell types
    return ThreeLevelHierarchy(s2c, c2b, has_states=False)


def build_three_level_from_obs(
    obs: pd.DataFrame,
    broad_col: str,
    cell_type_col: str,
    state_col: Optional[str] = None,
) -> ThreeLevelHierarchy:
    """Derive the hierarchy from annotation columns in ``adata.obs``.

    When ``state_col`` is None the hierarchy has no state level (cell_type is the
    finest); the (broad, cell_type) pairing must be consistent (one broad per
    cell type) — a conflict raises ``ValueError``.
    """
    for col in (broad_col, cell_type_col):
        if col not in obs.columns:
            raise KeyError(f"annotation column '{col}' not in obs")
    c2b: dict[str, str] = {}
    for ct, grp in obs.groupby(obs[cell_type_col].astype(str)):
        broads = set(grp[broad_col].astype(str))
        if len(broads) > 1:
            raise ValueError(f"cell type '{ct}' maps to multiple broad families "
                             f"{sorted(broads)}")
        c2b[str(ct)] = next(iter(broads))
    if state_col and state_col in obs.columns:
        s2c: dict[str, str] = {}
        for st, grp in obs.groupby(obs[state_col].astype(str)):
            cts = set(grp[cell_type_col].astype(str))
            if len(cts) > 1:
                raise ValueError(f"state '{st}' maps to multiple cell types "
                                 f"{sorted(cts)}")
            s2c[str(st)] = next(iter(cts))
        return build_three_level_hierarchy(s2c, c2b)
    return build_three_level_from_two_level(c2b)


def validate_three_level_hierarchy(
    hierarchy: ThreeLevelHierarchy,
    *,
    state_counts: Optional[Mapping[str, int]] = None,
    celltype_donors: Optional[Mapping[str, object]] = None,
    min_cells_per_state: int = 20,
    emit_warnings: bool = True,
) -> pd.DataFrame:
    """Per-state validation table with low-support / single-donor flags.

    *state_counts*: ``{state: n_cells}``; *celltype_donors*:
    ``{cell_type: iterable of donor/batch ids}``.  Warnings are emitted (never
    hidden) and recorded as flag columns; nothing is dropped.
    """
    rows = []
    for st in hierarchy.states():
        ct = hierarchy.state_to_celltype[st]
        broad = hierarchy.celltype_to_broad.get(ct, ct)
        n = int(state_counts.get(st, -1)) if state_counts is not None else -1
        low = (state_counts is not None and 0 <= n < min_cells_per_state)
        donors = celltype_donors.get(ct) if celltype_donors is not None else None
        n_donors = len(set(donors)) if donors is not None else -1
        single_donor = (celltype_donors is not None and n_donors == 1)
        flags = []
        if low:
            flags.append(f"low_cells(<{min_cells_per_state})")
        if single_donor:
            flags.append("single_donor")
        rows.append({
            "state": st, "cell_type": ct, "broad_family": broad,
            "n_cells": n, "n_donors": n_donors,
            "low_support": bool(low), "single_donor": bool(single_donor),
            "flags": ";".join(flags),
        })
    df = pd.DataFrame(rows, columns=["state", "cell_type", "broad_family",
                                     "n_cells", "n_donors", "low_support",
                                     "single_donor", "flags"])
    if emit_warnings:
        n_low = int(df["low_support"].sum())
        n_sd = int(df["single_donor"].sum())
        if n_low:
            warnings.warn(f"{n_low} cell state(s) have < {min_cells_per_state} "
                          "cells; state-level resolution will be unreliable for "
                          "them.", stacklevel=2)
        if n_sd:
            warnings.warn(f"{n_sd} cell type(s) appear in a single donor/batch; "
                          "their signatures may be confounded.", stacklevel=2)
    return df


def summarize_three_level_hierarchy(
    hierarchy: ThreeLevelHierarchy,
    *,
    state_counts: Optional[Mapping[str, int]] = None,
) -> pd.DataFrame:
    """Per-broad-family summary: #cell types, #states, #cells."""
    rows = []
    for broad in hierarchy.broad_families():
        cts = hierarchy.celltypes_of_broad(broad)
        states = [s for ct in cts for s in hierarchy.states_of_celltype(ct)]
        n_cells = (sum(int(state_counts.get(s, 0)) for s in states)
                   if state_counts is not None else -1)
        rows.append({
            "broad_family": broad,
            "n_cell_types": len(cts),
            "n_states": len(states),
            "has_states": hierarchy.has_states,
            "n_reference_cells": n_cells,
        })
    return pd.DataFrame(rows, columns=["broad_family", "n_cell_types",
                                       "n_states", "has_states",
                                       "n_reference_cells"])


def write_three_level_outputs(
    hierarchy: ThreeLevelHierarchy,
    validation: pd.DataFrame,
    summary: pd.DataFrame,
    out_dir,
) -> dict[str, Path]:
    """Write the three Part-2 tables; returns ``{name: path}``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, df in (("hierarchy_three_level", hierarchy.frame()),
                     ("hierarchy_validation", validation),
                     ("hierarchy_summary", summary)):
        p = out / f"{name}.tsv"
        df.to_csv(p, sep="\t", index=False)
        paths[name] = p
    return paths
