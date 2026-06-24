"""Adaptive resolution reporting (P4): the best *supported* resolution per family.

This is **reporting / QC**, not a new deconvolution algorithm and not a change to
default estimation. It consumes existing diagnostics — pairwise Bhattacharyya
separability, the within-family state-similarity graph, discriminating-gene /
marker support, within-family signature conditioning (condition number / effective
rank), optional soft-gating subtype confidence, and optional unresolved mass — and
classifies every broad family into one of:

    resolved_fine            — subtypes are separable; report at fine level
    partially_resolved_group — some subtypes separable, others must be grouped
    broad_only               — all subtypes mutually collinear; report at family level
    unresolved_family        — most family mass is unresolved (gating abstained)
    diagnostic_only          — evidence insufficient to decide; flag for inspection

Rare/protected states are never auto-collapsed. The output is a tidy table:
``family | recommended_resolution | supported_states | grouped_states |
unresolved_mass | reason | confidence`` plus structured records.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd

from tissueresolve.results import ReferenceSignature

__all__ = ["FamilyResolution", "AdaptiveResolutionReport",
           "classify_family_resolution", "build_adaptive_resolution_report",
           "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"


@dataclass
class FamilyResolution:
    family: str
    recommended_resolution: str
    supported_states: list
    grouped_states: list           # list of lists (each an ambiguous group)
    unresolved_mass: float
    reason: str
    confidence: float
    n_members: int
    max_within_bc: float
    condition_number: float
    mean_subtype_confidence: float

    def as_row(self) -> dict:
        return {
            "family": self.family,
            "recommended_resolution": self.recommended_resolution,
            "supported_states": " | ".join(map(str, self.supported_states)),
            "grouped_states": " ; ".join("+".join(map(str, g)) for g in self.grouped_states),
            "unresolved_mass": round(self.unresolved_mass, 4),
            "n_members": self.n_members,
            "max_within_bc": round(self.max_within_bc, 4),
            "condition_number": round(self.condition_number, 2),
            "mean_subtype_confidence": round(self.mean_subtype_confidence, 4),
            "confidence": round(self.confidence, 4),
            "reason": self.reason,
        }


@dataclass
class AdaptiveResolutionReport:
    table: pd.DataFrame
    families: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    def write(self, out_dir) -> dict:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        p = out / "adaptive_resolution.tsv"
        self.table.to_csv(p, sep="\t", index=False)
        return {"adaptive_resolution": p}


def _within_family_condition_number(ref: ReferenceSignature, members: list) -> float:
    """Condition number of the within-family signature submatrix (log1p-CPM).

    Large κ ⇒ near-collinear columns ⇒ unstable fine split.
    """
    idx = {c: i for i, c in enumerate(ref.cell_types)}
    rows = [idx[m] for m in members if m in idx]
    if len(rows) < 2:
        return 1.0
    sub = ref.as_R_log()[rows, :].astype(np.float64)   # (m, G)
    sub = sub - sub.mean(axis=1, keepdims=True)
    s = np.linalg.svd(sub, compute_uv=False)
    s = s[s > 1e-12]
    if s.size < 2:
        return float("inf")
    return float(s[0] / s[-1])


def classify_family_resolution(
    family: str,
    members: list,
    bc_matrix: pd.DataFrame,
    *,
    ref: Optional[ReferenceSignature] = None,
    subtype_confidence: Optional[Mapping[str, float]] = None,
    unresolved_mass: float = 0.0,
    rare_protection: Optional[Iterable[str]] = None,
    bc_fine_threshold: float = 0.90,
    bc_collinear_threshold: float = 0.97,
    condition_number_threshold: float = 50.0,
    unresolved_family_threshold: float = 0.50,
) -> FamilyResolution:
    """Classify a single family's best supported resolution from diagnostics."""
    protected = {str(s) for s in (rare_protection or [])}
    members = [str(m) for m in members]
    n = len(members)
    conf_map = {str(k): float(v) for k, v in (subtype_confidence or {}).items()}
    mean_conf = float(np.mean([conf_map[m] for m in members if m in conf_map])) if any(
        m in conf_map for m in members) else float("nan")
    kappa = _within_family_condition_number(ref, members) if ref is not None else float("nan")

    # single-member family → trivially fine
    if n <= 1:
        return FamilyResolution(family, "resolved_fine", members, [], unresolved_mass,
                                "single cell type in family", 1.0, n, 0.0,
                                kappa if np.isfinite(kappa) else 1.0,
                                mean_conf if np.isfinite(mean_conf) else 1.0)

    # pairwise within-family BC
    pair_bc = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = members[i], members[j]
            if a in bc_matrix.index and b in bc_matrix.columns:
                pair_bc.append((a, b, float(bc_matrix.loc[a, b])))
    max_bc = max((w for _, _, w in pair_bc), default=0.0)

    # unresolved-dominated family
    if unresolved_mass >= unresolved_family_threshold:
        return FamilyResolution(family, "unresolved_family", [], [], unresolved_mass,
                                f"unresolved mass {unresolved_mass:.2f} ≥ "
                                f"{unresolved_family_threshold}", 1.0 - max_bc, n, max_bc,
                                kappa, mean_conf)

    # union-find grouping of collinear (BC > collinear_threshold) pairs, never
    # grouping a protected/rare state
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, w in pair_bc:
        if w > bc_collinear_threshold and a not in protected and b not in protected:
            parent[find(a)] = find(b)
    comp = {}
    for m in members:
        comp.setdefault(find(m), []).append(m)
    groups = [sorted(v) for v in comp.values() if len(v) > 1]
    singletons = [v[0] for v in comp.values() if len(v) == 1]

    all_collinear = (max_bc > bc_collinear_threshold) and (len(groups) == 1 and not singletons)
    well_conditioned = np.isfinite(kappa) and kappa <= condition_number_threshold
    separable = (max_bc < bc_fine_threshold) and well_conditioned
    conf = float(np.clip(1.0 - max_bc, 0.0, 1.0))
    if np.isfinite(mean_conf):
        conf = float(np.clip(0.5 * conf + 0.5 * mean_conf, 0.0, 1.0))

    if separable:
        return FamilyResolution(family, "resolved_fine", members, [], unresolved_mass,
                                f"all within-family BC < {bc_fine_threshold} "
                                f"(max {max_bc:.3f}), κ={kappa:.1f}", conf, n, max_bc,
                                kappa, mean_conf)
    if all_collinear:
        return FamilyResolution(family, "broad_only", [], [sorted(members)], unresolved_mass,
                                f"all subtypes mutually collinear (max BC {max_bc:.3f} > "
                                f"{bc_collinear_threshold}); report at family level",
                                conf, n, max_bc, kappa, mean_conf)
    if groups:
        return FamilyResolution(family, "partially_resolved_group", sorted(singletons), groups,
                                unresolved_mass,
                                f"{len(groups)} collinear group(s) merged, "
                                f"{len(singletons)} separable subtype(s) kept", conf, n, max_bc,
                                kappa, mean_conf)
    # 0.90 < max_bc ≤ 0.97 or ill-conditioned but no hard collinear pair → ambiguous
    return FamilyResolution(family, "diagnostic_only", members, [], unresolved_mass,
                            f"borderline separability (max BC {max_bc:.3f}, κ={kappa:.1f}); "
                            "inspect before trusting fine split", conf, n, max_bc, kappa, mean_conf)


def build_adaptive_resolution_report(
    ref: ReferenceSignature,
    mapping: Mapping[str, str],
    *,
    subtype_confidence: Optional[Mapping[str, float]] = None,
    unresolved_mass: Optional[Mapping[str, float]] = None,
    rare_protection: Optional[Iterable[str]] = None,
    bc_fine_threshold: float = 0.90,
    bc_collinear_threshold: float = 0.97,
    condition_number_threshold: float = 50.0,
) -> AdaptiveResolutionReport:
    """Build the per-family adaptive-resolution report from a reference + mapping.

    ``unresolved_mass`` maps family → mean ``unresolved_<family>`` mass (e.g. from a
    hierarchical run). ``subtype_confidence`` maps subtype → soft-gating confidence.
    Reporting only — does not alter any estimate.
    """
    from tissueresolve.reference.separability import compute_separability
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        rep = compute_separability(ref, warn_threshold=1.1)  # no warning spam here
    cts = [str(c) for c in ref.cell_types]
    bc = pd.DataFrame(np.eye(len(cts)), index=cts, columns=cts)
    for p in rep.pairs:
        if p.type_a in bc.index and p.type_b in bc.columns:
            bc.loc[p.type_a, p.type_b] = bc.loc[p.type_b, p.type_a] = p.bhattacharyya_coeff

    fam_members: dict[str, list] = {}
    for ct in cts:
        fam_members.setdefault(str(mapping.get(ct, ct)), []).append(ct)

    um = {str(k): float(v) for k, v in (unresolved_mass or {}).items()}
    families = []
    for fam, members in sorted(fam_members.items()):
        fr = classify_family_resolution(
            fam, members, bc, ref=ref, subtype_confidence=subtype_confidence,
            unresolved_mass=um.get(fam, um.get(f"unresolved_{fam}", 0.0)),
            rare_protection=rare_protection, bc_fine_threshold=bc_fine_threshold,
            bc_collinear_threshold=bc_collinear_threshold,
            condition_number_threshold=condition_number_threshold)
        families.append(fr)

    table = pd.DataFrame([f.as_row() for f in families])
    counts = table["recommended_resolution"].value_counts().to_dict() if not table.empty else {}
    return AdaptiveResolutionReport(
        table=table, families=families,
        metadata={"feature_status": FEATURE_STATUS, "n_families": len(families),
                  "resolution_counts": counts, "estimates_modified": False,
                  "used_for": "reporting_only"})
