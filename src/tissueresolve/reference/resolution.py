"""
Resolution-aware analysis of a cell-type reference.

This turns the existing pairwise *separability* diagnostics into an explicit,
actionable statement of **what can and cannot be resolved**:

* group confusable cell types into :class:`CellTypeFamily` objects,
* classify every pair into one of four resolvability classes,
* recommend merges, and
* (with :func:`apply_unresolved_mode`) collapse unresolvable subtype families
  into an honest ``unresolved_<family>`` mass rather than forcing confident
  subtype estimates.

It consumes the outputs of :func:`tissueresolve.reference.separability.
compute_separability` (Bhattacharyya coefficient, Jeffreys divergence,
Pearson r, number of discriminating genes) — it does **not** re-derive the
separability metrics or touch the core deconvolution algorithms.

Scientific-honesty rules enforced here:

* Non-separability is never hidden — unresolved families are labelled as such.
* Subtype accuracy is never overclaimed — a poorly/unresolved family is
  reported at the family (broad) level unless the user opts out.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from tissueresolve.results import PairSeparability, SeparabilityReport

__all__ = [
    "RESOLVABILITY_CLASSES",
    "ResolutionConfig",
    "CellTypeFamily",
    "ResolutionNode",
    "ResolutionReport",
    "classify_resolvability",
    "compute_pairwise_resolution_score",
    "build_similarity_graph",
    "infer_cell_type_families",
    "suggest_merges",
    "build_resolution_report",
    "apply_unresolved_mode",
    "annotate_resolution",
    # Resolution-aware recommendation system
    "RESOLUTION_MODES",
    "build_nonseparable_graph",
    "recommend_cell_type_merges",
    "assign_resolution_families",
    "family_label_for_members",
    "write_recommended_merges",
    "summarize_resolution_report",
]

RESOLUTION_MODES = ("none", "suggest", "auto", "hierarchical")

RESOLVABILITY_CLASSES = (
    "resolved", "partially_resolved", "poorly_resolved", "unresolved",
)


@dataclass
class ResolutionConfig:
    """Thresholds governing resolvability classification and the abstain mode.

    All thresholds are on the **separability score** ``s = 1 − Bhattacharyya``
    (1 = orthogonal/easy, 0 = identical), chosen to align with the existing
    :class:`~tissueresolve.results.PairSeparability` risk levels.
    """

    resolved_threshold: float = 0.20      # s ≥ 0.20  → resolved      (BC ≤ 0.80)
    partial_threshold: float = 0.10       # s ≥ 0.10  → partial        (BC ≤ 0.90)
    poor_threshold: float = 0.03          # s ≥ 0.03  → poorly         (BC ≤ 0.97)
    family_bc_threshold: float = 0.90     # pairs with BC ≥ this cluster into a family
    # Resolution mode — how non-separable types are handled:
    #   "none"         — report only (improved warning), estimates unchanged
    #   "suggest"      — compute recommended merges, estimates unchanged (default)
    #   "auto"         — apply recommended merges (aggregate to families)
    #   "hierarchical" — broad families first, subtypes only when resolvable
    resolution_mode: str = "suggest"
    # Abstain / unresolved mode
    allow_unresolved: bool = True
    unresolved_threshold: float = 0.10    # family mean s below → candidate to collapse
    spillover_threshold: float = 0.30     # spillover risk above → corroborates collapse
    uncertainty_threshold: float = 0.10   # bootstrap CI width above → corroborates collapse


@dataclass
class CellTypeFamily:
    """A group of cell types that may be confusable with one another."""

    name: str
    members: list[str]
    resolvability: str                    # worst-case class within the family
    mean_separability: float              # mean (1 − BC) over internal pairs
    representative: str                   # member used to label the family

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def is_resolvable(self) -> bool:
        return self.size == 1 or self.resolvability in ("resolved", "partially_resolved")


@dataclass
class ResolutionNode:
    """A node in the broad→subtype resolution hierarchy."""

    name: str
    members: list[str]
    is_family: bool = False
    resolvability: str = "resolved"
    children: list["ResolutionNode"] = field(default_factory=list)


@dataclass
class ResolutionReport:
    """Full resolvability picture for a reference."""

    cell_types: list[str]
    families: list[CellTypeFamily]
    pairwise: pd.DataFrame                # one row per cell-type pair
    config: ResolutionConfig

    @property
    def n_families(self) -> int:
        return len(self.families)

    @property
    def unresolved_families(self) -> list[CellTypeFamily]:
        return [f for f in self.families
                if f.size > 1 and f.resolvability in ("poorly_resolved", "unresolved")]

    def family_of(self, cell_type: str) -> Optional[CellTypeFamily]:
        for f in self.families:
            if cell_type in f.members:
                return f
        return None

    def summary_str(self) -> str:
        counts = self.pairwise["resolvability"].value_counts().to_dict()
        lines = [
            f"ResolutionReport: {len(self.cell_types)} cell types, "
            f"{self.n_families} families.",
            "  pair classes: " + ", ".join(
                f"{c}={counts.get(c, 0)}" for c in RESOLVABILITY_CLASSES),
            f"  unresolved families: {len(self.unresolved_families)}",
        ]
        for f in self.families:
            if f.size > 1:
                lines.append(
                    f"    [{f.resolvability:18s}] {f.name}: {', '.join(f.members)} "
                    f"(mean separability {f.mean_separability:.3f})")
        return "\n".join(lines)

    def save(self, path: Path | str) -> None:
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.pairwise.to_csv(path / "pairwise_resolvability.tsv", sep="\t", index=False)
        fam_rows = [{
            "family": f.name, "resolvability": f.resolvability,
            "n_members": f.size, "mean_separability": round(f.mean_separability, 4),
            "members": "; ".join(f.members),
        } for f in self.families]
        pd.DataFrame(fam_rows).to_csv(
            path / "cell_type_families.tsv", sep="\t", index=False)
        with (path / "resolution_config.json").open("w", encoding="utf-8") as fh:
            json.dump(asdict(self.config), fh, indent=2)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def classify_resolvability(separability_score: float,
                           cfg: Optional[ResolutionConfig] = None) -> str:
    """Map a separability score ``s = 1 − BC`` to a resolvability class."""
    cfg = cfg or ResolutionConfig()
    s = float(separability_score)
    if s >= cfg.resolved_threshold:
        return "resolved"
    if s >= cfg.partial_threshold:
        return "partially_resolved"
    if s >= cfg.poor_threshold:
        return "poorly_resolved"
    return "unresolved"


def compute_pairwise_resolution_score(
    pair: PairSeparability, *, marker_overlap: Optional[float] = None
) -> float:
    """A single ``[0, 1]`` resolvability score blending the separability metrics.

    Combines the separability score (``1 − BC``), profile decorrelation
    (``1 − max(Pearson, 0)``) and discriminating-gene density.  Optionally
    penalised by *marker_overlap* (fraction of shared markers).  Higher = more
    resolvable.  Classification still uses the raw ``1 − BC`` for transparency.
    """
    base = max(0.0, 1.0 - pair.bhattacharyya_coeff)
    decorr = 1.0 - max(0.0, pair.pearson_r)
    gene_density = min(1.0, pair.n_discriminating_genes / 100.0)
    score = 0.5 * base + 0.3 * decorr + 0.2 * gene_density
    if marker_overlap is not None:
        score *= max(0.0, 1.0 - float(marker_overlap))
    return float(np.clip(score, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Similarity graph & families
# ---------------------------------------------------------------------------


def build_similarity_graph(report: SeparabilityReport,
                           cell_types: list[str]) -> pd.DataFrame:
    """Return a ``K × K`` Bhattacharyya-coefficient similarity matrix (diag=1)."""
    idx = {ct: i for i, ct in enumerate(cell_types)}
    K = len(cell_types)
    M = np.eye(K, dtype=float)
    for p in report.pairs:
        if p.type_a in idx and p.type_b in idx:
            i, j = idx[p.type_a], idx[p.type_b]
            M[i, j] = M[j, i] = float(p.bhattacharyya_coeff)
    return pd.DataFrame(M, index=cell_types, columns=cell_types)


def infer_cell_type_families(
    report: SeparabilityReport,
    cell_types: list[str],
    cfg: Optional[ResolutionConfig] = None,
) -> list[CellTypeFamily]:
    """Cluster confusable cell types into families (connected components).

    Two types are linked when their Bhattacharyya coefficient is ≥
    ``cfg.family_bc_threshold`` (i.e. poorly separable).  Each connected
    component becomes a :class:`CellTypeFamily`; isolated types form
    singleton families classified ``resolved``.
    """
    cfg = cfg or ResolutionConfig()
    parent = {ct: ct for ct in cell_types}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    pair_lookup: dict[frozenset, PairSeparability] = {}
    for p in report.pairs:
        if p.type_a in parent and p.type_b in parent:
            pair_lookup[frozenset((p.type_a, p.type_b))] = p
            if p.bhattacharyya_coeff >= cfg.family_bc_threshold:
                union(p.type_a, p.type_b)

    groups: dict[str, list[str]] = {}
    for ct in cell_types:
        groups.setdefault(find(ct), []).append(ct)

    families: list[CellTypeFamily] = []
    for members in groups.values():
        members = sorted(members)
        rep = members[0]
        if len(members) == 1:
            families.append(CellTypeFamily(
                name=f"family:{rep}", members=members,
                resolvability="resolved", mean_separability=1.0, representative=rep))
            continue
        scores, classes = [], []
        for a, b in _pairs(members):
            p = pair_lookup.get(frozenset((a, b)))
            if p is None:
                continue
            s = max(0.0, 1.0 - p.bhattacharyya_coeff)
            scores.append(s)
            classes.append(classify_resolvability(s, cfg))
        worst = _worst_class(classes) if classes else "resolved"
        families.append(CellTypeFamily(
            name=f"family:{rep}+{len(members) - 1}", members=members,
            resolvability=worst,
            mean_separability=float(np.mean(scores)) if scores else 1.0,
            representative=rep))
    families.sort(key=lambda f: (f.size > 1, -f.size, f.name))
    return families


def suggest_merges(report: SeparabilityReport, cell_types: list[str],
                   cfg: Optional[ResolutionConfig] = None) -> pd.DataFrame:
    """Recommend cell-type groups to merge (poorly/unresolved multi-member families)."""
    cfg = cfg or ResolutionConfig()
    families = infer_cell_type_families(report, cell_types, cfg)
    rows = [{
        "family": f.name, "members": "; ".join(f.members), "n_members": f.size,
        "resolvability": f.resolvability,
        "mean_separability": round(f.mean_separability, 4),
        "recommended_merge": f.resolvability in ("poorly_resolved", "unresolved"),
    } for f in families if f.size > 1]
    return pd.DataFrame(rows, columns=[
        "family", "members", "n_members", "resolvability",
        "mean_separability", "recommended_merge"])


def build_resolution_report(
    report: SeparabilityReport,
    cell_types: list[str],
    cfg: Optional[ResolutionConfig] = None,
) -> ResolutionReport:
    """Assemble the full :class:`ResolutionReport` from a separability report."""
    cfg = cfg or ResolutionConfig()
    rows = []
    for p in report.pairs:
        s = max(0.0, 1.0 - p.bhattacharyya_coeff)
        rows.append({
            "type_a": p.type_a, "type_b": p.type_b,
            "bhattacharyya": round(p.bhattacharyya_coeff, 4),
            "separability_score": round(s, 4),
            "jeffreys": round(p.jeffreys_divergence, 4),
            "pearson": round(p.pearson_r, 4),
            "n_discriminating_genes": int(p.n_discriminating_genes),
            "resolution_score": round(compute_pairwise_resolution_score(p), 4),
            "resolvability": classify_resolvability(s, cfg),
        })
    pairwise = pd.DataFrame(rows)
    if not pairwise.empty:
        pairwise = pairwise.sort_values("separability_score").reset_index(drop=True)
    families = infer_cell_type_families(report, cell_types, cfg)
    return ResolutionReport(
        cell_types=list(cell_types), families=families,
        pairwise=pairwise, config=cfg)


# ---------------------------------------------------------------------------
# Abstain / unresolved mode
# ---------------------------------------------------------------------------


def apply_unresolved_mode(
    proportions: pd.DataFrame,
    resolution_report: ResolutionReport,
    *,
    spillover_risk: Optional[pd.Series] = None,
    bootstrap_uncertainty: Optional[pd.Series] = None,
    cfg: Optional[ResolutionConfig] = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Collapse unresolvable subtype families into honest ``unresolved_*`` mass.

    A multi-member family is collapsed when ``allow_unresolved`` is set and its
    mean separability is below ``unresolved_threshold``.  When *spillover_risk*
    and/or *bootstrap_uncertainty* are provided, they must *also* exceed their
    thresholds (defence in depth — collapse only when several signals agree).

    The collapsed column ``unresolved_<family>`` holds the summed mass of its
    members; **total mass per sample is preserved**.  Resolved families keep
    their subtype columns.  Returns ``(new_proportions, info)``.

    Never forces confident subtype estimates; never silently drops mass.
    """
    cfg = cfg or resolution_report.config or ResolutionConfig()
    cols = list(proportions.columns)
    collapsed: list[dict[str, Any]] = []
    kept: dict[str, list[str]] = {}

    if not cfg.allow_unresolved:
        return proportions.copy(), {"collapsed_families": [], "allow_unresolved": False}

    for fam in resolution_report.families:
        members = [m for m in fam.members if m in cols]
        if len(members) <= 1:
            continue
        below_sep = fam.mean_separability < cfg.unresolved_threshold
        risk_ok = True
        if spillover_risk is not None:
            fam_risk = float(np.nanmean([spillover_risk.get(m, np.nan)
                                         for m in members]))
            risk_ok = np.isnan(fam_risk) or fam_risk >= cfg.spillover_threshold
        unc_ok = True
        if bootstrap_uncertainty is not None:
            fam_unc = float(np.nanmean([bootstrap_uncertainty.get(m, np.nan)
                                        for m in members]))
            unc_ok = np.isnan(fam_unc) or fam_unc >= cfg.uncertainty_threshold
        if below_sep and risk_ok and unc_ok:
            kept[f"unresolved_{fam.representative}"] = members
            collapsed.append({"family": fam.name, "label": f"unresolved_{fam.representative}",
                              "members": members, "mean_separability": fam.mean_separability})

    if not collapsed:
        return proportions.copy(), {"collapsed_families": []}

    collapsed_members = {m for c in collapsed for m in c["members"]}
    out = pd.DataFrame(index=proportions.index)
    for col in cols:
        if col not in collapsed_members:
            out[col] = proportions[col]
    for label, members in kept.items():
        out[label] = proportions[members].sum(axis=1)

    return out, {"collapsed_families": collapsed,
                 "n_collapsed": len(collapsed),
                 "mass_preserved": True}


# ---------------------------------------------------------------------------
# Per-cell-type annotation (additive QC output for bulk & spatial)
# ---------------------------------------------------------------------------


def annotate_resolution(
    mean_estimate: pd.Series,
    resolution_report: ResolutionReport,
    *,
    spillover_risk: Optional[pd.DataFrame] = None,
    modality: str = "bulk",
) -> pd.DataFrame:
    """Per-cell-type resolution/spillover annotation (an *additional* output).

    Combines the (mean) estimate with each type's family, family resolvability,
    spillover risk + main leaking partner (from
    :func:`tissueresolve.benchmark.spillover.compute_spillover_risk`), and a
    recommended interpretation level (``subtype`` vs ``family/broad``).  It does
    **not** modify the core estimate — it annotates it.
    """
    cfg = resolution_report.config
    rows = []
    for ct in mean_estimate.index:
        fam = resolution_report.family_of(ct)
        fam_name = fam.name if fam else f"family:{ct}"
        fam_class = fam.resolvability if fam else "resolved"
        risk = partner = None
        partner_frac = np.nan
        if spillover_risk is not None and ct in spillover_risk.index:
            risk = float(spillover_risk.loc[ct, "spillover_risk"])
            partner = spillover_risk.loc[ct, "main_partner"]
            partner_frac = float(spillover_risk.loc[ct, "main_partner_fraction"])
        resolvable = (fam is None or fam.is_resolvable) and (
            risk is None or risk < cfg.spillover_threshold)
        rows.append({
            "cell_type": ct,
            "mean_estimate": round(float(mean_estimate[ct]), 4),
            "family": fam_name,
            "family_resolvability": fam_class,
            "spillover_risk": round(risk, 4) if risk is not None else np.nan,
            "main_partner": partner,
            "main_partner_fraction": round(partner_frac, 4)
            if not np.isnan(partner_frac) else np.nan,
            "recommended_level": "subtype" if resolvable else "family/broad",
            "modality": modality,
        })
    return pd.DataFrame(rows).set_index("cell_type")


# ---------------------------------------------------------------------------
# Resolution-aware recommendation system
# ---------------------------------------------------------------------------


def build_nonseparable_graph(
    report: SeparabilityReport,
    cell_types: list[str],
    cfg: Optional[ResolutionConfig] = None,
) -> dict[str, set[str]]:
    """Adjacency dict linking cell types whose pair is poorly separable.

    An edge ``a–b`` exists when ``Bhattacharyya ≥ cfg.family_bc_threshold``.
    Isolated (well-separated) types appear with an empty neighbour set.
    """
    cfg = cfg or ResolutionConfig()
    graph: dict[str, set[str]] = {ct: set() for ct in cell_types}
    present = set(cell_types)
    for p in report.pairs:
        if (p.type_a in present and p.type_b in present
                and p.bhattacharyya_coeff >= cfg.family_bc_threshold):
            graph[p.type_a].add(p.type_b)
            graph[p.type_b].add(p.type_a)
    return graph


# Friendly family names from the set of broad groups in a component.
_LYMPHOID_GROUPS = {"T cell", "NK cell", "Lymphoid", "B cell", "Plasma cell"}


def family_label_for_members(members: list[str]) -> str:
    """Heuristic family name for a group of confusable cell types.

    Examples: CD4/CD8/Treg/NK/lymphocyte → "T/NK lymphocytes";
    macrophage/monocyte/DC → "myeloid cells"; endothelial → "endothelial cells";
    luminal/basal/epithelial → "epithelial cells"; pericyte/smooth muscle →
    "mural cells".
    """
    from tissueresolve.reference.hierarchy import infer_broad_groups_from_labels

    if len(members) == 1:
        return members[0]
    groups = set(infer_broad_groups_from_labels(members).values())

    if groups <= _LYMPHOID_GROUPS and groups:
        if groups & {"T cell", "NK cell"}:
            return "T/NK lymphocytes"
        if groups == {"B cell"} or groups == {"Plasma cell"}:
            return "B/plasma cells"
        return "lymphocytes"
    if groups == {"Myeloid"}:
        return "myeloid cells"
    if groups == {"Endothelial"}:
        return "endothelial cells"
    if groups == {"Epithelial"}:
        return "epithelial cells"
    if groups == {"Perivascular"}:
        return "mural cells"
    if groups == {"Fibroblast"}:
        return "fibroblasts"
    if len(groups) == 1:
        return f"{next(iter(groups))} family"
    # Mixed lineages — name by the most common broad group.
    bg = infer_broad_groups_from_labels(members)
    counts: dict[str, int] = {}
    for g in bg.values():
        counts[g] = counts.get(g, 0) + 1
    dominant = max(counts, key=counts.get)
    return f"mixed family ({dominant})"


def _connected_components(graph: dict[str, set[str]]) -> list[list[str]]:
    seen: set[str] = set()
    comps: list[list[str]] = []
    for node in graph:
        if node in seen:
            continue
        stack, comp = [node], []
        seen.add(node)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in graph[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def _named_components(
    report: SeparabilityReport, cell_types: list[str], cfg: ResolutionConfig,
) -> list[tuple[str, list[str]]]:
    """Connected components with **unique** family labels.

    Two distinct components must never share a name, or post-hoc aggregation
    would merge unrelated families.  On collision the representative member is
    appended (and a counter if still colliding).
    """
    graph = build_nonseparable_graph(report, cell_types, cfg)
    used: set[str] = set()
    out: list[tuple[str, list[str]]] = []
    for comp in _connected_components(graph):
        if len(comp) == 1:
            label = comp[0]
        else:
            label = family_label_for_members(comp)
            if label in used:
                label = f"{label} ({comp[0]})"
            n = 2
            while label in used:
                label = f"{family_label_for_members(comp)} ({comp[0]} #{n})"
                n += 1
        used.add(label)
        out.append((label, comp))
    return out


def recommend_cell_type_merges(
    report: SeparabilityReport,
    cell_types: list[str],
    cfg: Optional[ResolutionConfig] = None,
) -> pd.DataFrame:
    """Recommend merge families from the non-separable graph.

    Returns one row per **multi-member** family with a unique suggested family
    name, its members, size, resolvability and mean separability.  Singletons
    (well separated types) are intentionally omitted — they should not be merged.
    """
    cfg = cfg or ResolutionConfig()
    res = build_resolution_report(report, cell_types, cfg)
    fam_lookup = {tuple(sorted(f.members)): f for f in res.families}
    rows = []
    for label, comp in _named_components(report, cell_types, cfg):
        if len(comp) < 2:
            continue
        fam = fam_lookup.get(tuple(comp))
        rows.append({
            "family_name": label,
            "members": "; ".join(comp),
            "n_members": len(comp),
            "resolvability": fam.resolvability if fam else "poorly_resolved",
            "mean_separability": round(fam.mean_separability, 4) if fam else float("nan"),
            "recommended_merge": True,
        })
    rows.sort(key=lambda r: -r["n_members"])
    return pd.DataFrame(rows, columns=[
        "family_name", "members", "n_members", "resolvability",
        "mean_separability", "recommended_merge"])


def assign_resolution_families(
    report: SeparabilityReport,
    cell_types: list[str],
    cfg: Optional[ResolutionConfig] = None,
) -> dict[str, str]:
    """Map every cell type → its (unique) family name (separated types map to self).

    This is the explicit, reproducible merge mapping used by
    :func:`tissueresolve.reference.hierarchy.merge_reference_cell_types` and
    :func:`aggregate_predictions_by_family`.
    """
    cfg = cfg or ResolutionConfig()
    mapping: dict[str, str] = {}
    for label, comp in _named_components(report, cell_types, cfg):
        for ct in comp:
            mapping[ct] = label
    for ct in cell_types:
        mapping.setdefault(ct, ct)
    return mapping


def write_recommended_merges(merges: pd.DataFrame, path: Path | str) -> Path:
    """Write the recommended-merges table (machine-readable TSV)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    merges.to_csv(path, sep="\t", index=False)
    return path


def summarize_resolution_report(
    report: SeparabilityReport,
    cell_types: list[str],
    cfg: Optional[ResolutionConfig] = None,
    *,
    spillover_risk: Optional[pd.Series] = None,
) -> str:
    """Markdown summary: high-risk pairs, recommended merge families, guidance."""
    cfg = cfg or ResolutionConfig()
    res = build_resolution_report(report, cell_types, cfg)
    merges = recommend_cell_type_merges(report, cell_types, cfg)
    counts = (res.pairwise["resolvability"].value_counts().to_dict()
              if not res.pairwise.empty else {})
    n_crit = sum(1 for p in report.pairs if p.bhattacharyya_coeff > 0.97)
    n_high = sum(1 for p in report.pairs
                 if cfg.family_bc_threshold < p.bhattacharyya_coeff <= 0.97)
    lines = [
        "# Resolution recommendation", "",
        f"- cell types: {len(cell_types)}",
        f"- HIGH separability pairs: {n_high}; CRITICAL: {n_crit}",
        f"- pair classes: " + ", ".join(
            f"{c}={counts.get(c, 0)}" for c in RESOLVABILITY_CLASSES),
        f"- recommended merge families: {len(merges)}",
        "",
        "## Recommended merge families",
    ]
    if len(merges):
        for _, r in merges.iterrows():
            lines.append(f"- **{r['family_name']}** ← {r['members']} "
                         f"({r['resolvability']}, mean separability "
                         f"{r['mean_separability']})")
    else:
        lines.append("- none (all pairs adequately separable)")
    lines += [
        "", "## Guidance",
        "- Subtype-level estimates for confusable pairs may be unreliable; "
        "interpret them at the family level.",
        "- Use `--resolution-mode suggest` (default) to get recommendations, "
        "`auto` to aggregate to families before deconvolution, or "
        "`hierarchical` for broad-then-subtype estimation.",
        "- Merging is always explicit and recorded — TissueResolve never "
        "merges cell types silently.",
        "- Bulk = mRNA-derived proportions (not cell fractions); spatial = "
        "spot-level composition (not cell counts).",
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pairs(items):
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            yield items[i], items[j]


def _worst_class(classes: list[str]) -> str:
    order = {c: i for i, c in enumerate(RESOLVABILITY_CLASSES)}
    return max(classes, key=lambda c: order.get(c, 0))
