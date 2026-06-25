"""Experimental state-similarity / Redeconve-*inspired* regularization (opt-in).

This module addresses a limitation that pure **spatial** smoothing could not:
diffuse mass spread across biologically *similar fine states* within a broad
family (inflated effective-N, weak conditional within-family recovery). The
spatial CAR penalty couples neighbouring *spots*; this module couples similar
*states*.

It is **inspired by the concept** of state-aware regularization (as discussed for
Redeconve and related methods). It is NOT a copy of Redeconve's code, API,
objective, or claims, and asserts no equivalence to Redeconve.

Design (see ``docs/STATE_SIMILARITY_REGULARIZATION_REPORT.md`` Phase 0):

* It does NOT modify the NB-CAR solver (that would be a major rewrite). Instead it
  is a **post-fit, constrained, within-family refinement** of the fitted spot×state
  proportions:
  - ``compute_state_similarity_graph`` — similarity graph over states (within-family
    by default, so unrelated families are never coupled);
  - ``apply_state_regularized_refinement`` — concentrates mass among *redundant*
    similar states (sharpening, not Laplacian smoothing) + sparsity shrinkage,
    **conserving each spot's broad-family mass exactly** and leaving
    ``unresolved_<family>`` columns untouched;
  - ``apply_sparsity_aware_refinement`` — sparsity shrinkage alone, with rare-state
    protection;
  - ``recommend_state_groups`` — diagnostic grouping of indistinguishable states.

All operations preserve non-negativity and valid compositions, conserve broad-family
mass, keep unresolved mass available, and never touch default behaviour
(everything here is opt-in and feature-flagged ``experimental``).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
import pandas as pd

__all__ = [
    "StateSimilarityGraph",
    "compute_state_similarity_graph",
    "apply_state_regularized_refinement",
    "apply_sparsity_aware_refinement",
    "recommend_state_groups",
    "FEATURE_STATUS",
]

FEATURE_STATUS = "experimental"


# ---------------------------------------------------------------------------
# State-similarity graph
# ---------------------------------------------------------------------------


@dataclass
class StateSimilarityGraph:
    """Undirected similarity graph over cell states/subtypes.

    ``state_i[e]``, ``state_j[e]``, ``weights[e]`` describe edge *e* between two
    states (indices into ``cell_types``). Edges are stored once (``i < j``), carry
    a non-negative weight in ``[0, 1]``, and contain no self-loops.
    """

    state_i: np.ndarray
    state_j: np.ndarray
    weights: np.ndarray
    cell_types: list[str]
    family_map: dict[str, str]
    metadata: dict = field(default_factory=dict)

    @property
    def n_states(self) -> int:
        return len(self.cell_types)

    @property
    def n_edges(self) -> int:
        return int(self.weights.size)

    def within_family_adjacency(self, states: list[str]) -> np.ndarray:
        """Dense row-normalised adjacency over the given ordered subset of states."""
        idx = {s: i for i, s in enumerate(states)}
        m = len(states)
        A = np.zeros((m, m), dtype=np.float64)
        for e in range(self.n_edges):
            a = self.cell_types[int(self.state_i[e])]
            b = self.cell_types[int(self.state_j[e])]
            if a in idx and b in idx:
                w = float(self.weights[e])
                A[idx[a], idx[b]] = w
                A[idx[b], idx[a]] = w
        return A


def _pairwise_similarity(profiles: np.ndarray, method: str) -> np.ndarray:
    """Symmetric (K, K) similarity in [-1, 1]. Correlation = centred cosine."""
    X = np.asarray(profiles, dtype=np.float64)
    if np.all(X >= 0):  # expression-like → log stabilises heavy tails (monotone)
        X = np.log1p(X)
    if method == "correlation":
        X = X - X.mean(axis=1, keepdims=True)
    elif method != "cosine":
        raise ValueError(f"Unknown similarity method {method!r}; use 'correlation' or 'cosine'.")
    norm = np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    U = X / norm
    S = U @ U.T
    np.clip(S, -1.0, 1.0, out=S)
    return S


def _connected_components(n: int, edges: list[tuple[int, int]]) -> int:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return len({find(i) for i in range(n)})


def compute_state_similarity_graph(
    reference_profiles: np.ndarray,
    cell_types: list[str],
    family_map: Optional[dict[str, str]] = None,
    method: str = "correlation",
    k_states: int = 5,
    min_similarity: float = 0.0,
    within_family_only: bool = True,
) -> StateSimilarityGraph:
    """Build a kNN similarity graph over cell states from reference profiles.

    Parameters
    ----------
    reference_profiles:
        ``(K, G)`` reference expression (e.g. CPM or log1p-CPM); rows are states.
    cell_types:
        ``K`` state names (row order of ``reference_profiles``).
    family_map:
        ``{state: broad_family}``. When ``within_family_only`` is True and this is
        ``None``, **no edges are created** (safe no-op) — unrelated families are
        never coupled without explicit family information.
    method:
        ``"correlation"`` (default) or ``"cosine"``.
    k_states:
        Keep up to this many strongest neighbours per state.
    min_similarity:
        Drop edges with similarity ``<=`` this value (negative similarities are
        always dropped — coupling weights are non-negative).
    within_family_only:
        Restrict edges to same-family state pairs (default, conservative).
    """
    cell_types = [str(c) for c in cell_types]
    K = len(cell_types)
    if reference_profiles.shape[0] != K:
        raise ValueError(
            f"reference_profiles has {reference_profiles.shape[0]} rows but "
            f"{K} cell_types were given.")

    if family_map is None:
        # No family information → each state is its own family. This makes
        # within-family operations (edges, refinement) a genuine no-op rather than
        # collapsing unrelated states into a single pseudo-family.
        fam = {c: c for c in cell_types}
        no_family_info = True
    else:
        fam = {c: str(family_map.get(c, c)) for c in cell_types}
        no_family_info = False

    S = _pairwise_similarity(reference_profiles, method)

    edges: dict[tuple[int, int], float] = {}
    # within-family-only with no family information → no edges (documented no-op)
    build = not (within_family_only and no_family_info)
    if build:
        for i in range(K):
            order = np.argsort(S[i])[::-1]  # descending similarity
            kept = 0
            for j in order:
                j = int(j)
                if j == i:
                    continue
                if within_family_only and fam[cell_types[i]] != fam[cell_types[j]]:
                    continue
                w = float(S[i, j])
                if w <= min_similarity or w <= 0.0:
                    continue
                a, b = (i, j) if i < j else (j, i)
                prev = edges.get((a, b))
                if prev is None or w > prev:
                    edges[(a, b)] = w
                kept += 1
                if kept >= k_states:
                    break

    if edges:
        keys = sorted(edges)
        si = np.array([k[0] for k in keys], dtype=np.int64)
        sj = np.array([k[1] for k in keys], dtype=np.int64)
        ww = np.array([edges[k] for k in keys], dtype=np.float64)
    else:
        si = np.zeros(0, dtype=np.int64)
        sj = np.zeros(0, dtype=np.int64)
        ww = np.zeros(0, dtype=np.float64)

    degree = np.zeros(K, dtype=np.int64)
    for a, b in zip(si, sj):
        degree[a] += 1
        degree[b] += 1
    isolated = [cell_types[i] for i in range(K) if degree[i] == 0]
    if degree.max(initial=0) > 0:
        order = np.argsort(degree)[::-1]
        highly_connected = [
            {"state": cell_types[int(i)], "degree": int(degree[int(i)])}
            for i in order[:5] if degree[int(i)] > 0
        ]
    else:
        highly_connected = []

    metadata = {
        "n_states": K,
        "n_edges": int(ww.size),
        "within_family_only": bool(within_family_only),
        "similarity_method": method,
        "k_states": int(k_states),
        "min_similarity": float(min_similarity),
        "weight_min": float(ww.min()) if ww.size else float("nan"),
        "weight_max": float(ww.max()) if ww.size else float("nan"),
        "weight_mean": float(ww.mean()) if ww.size else float("nan"),
        "weight_median": float(np.median(ww)) if ww.size else float("nan"),
        "n_connected_components": _connected_components(K, list(zip(si.tolist(), sj.tolist()))),
        "families": sorted(set(fam.values())),
        "n_families": len(set(fam.values())),
        "highly_connected_states": highly_connected,
        "isolated_states": isolated,
        "n_isolated_states": len(isolated),
        "no_family_info": no_family_info,
        "feature_status": FEATURE_STATUS,
    }
    return StateSimilarityGraph(
        state_i=si, state_j=sj, weights=ww,
        cell_types=cell_types, family_map=fam, metadata=metadata,
    )


# ---------------------------------------------------------------------------
# Within-family helpers (shared by refinement functions)
# ---------------------------------------------------------------------------


def _split_columns(columns: list[str], family_map: Optional[dict[str, str]],
                   known_states: Optional[set[str]]):
    """Return (fine_state_cols, unresolved_cols, family_of) for a proportions frame.

    ``unresolved_<family>`` columns are passthrough (never refined). A column is a
    refinable fine state only if it is not an ``unresolved_*`` column and (when
    ``known_states`` is provided) appears among the graph's states.
    """
    fine, unresolved, family_of = [], [], {}
    for c in columns:
        cs = str(c)
        if cs.startswith("unresolved_"):
            unresolved.append(c)
            continue
        if known_states is not None and cs not in known_states:
            # unknown column (not in the state graph) → leave untouched
            continue
        fine.append(c)
        if family_map is not None:
            family_of[c] = str(family_map.get(cs, cs))
        else:
            family_of[c] = cs  # each its own family → no within-family coupling
    return fine, unresolved, family_of


def _families_to_members(fine_cols: list[str], family_of: dict) -> dict[str, list[str]]:
    fam_members: dict[str, list[str]] = {}
    for c in fine_cols:
        fam_members.setdefault(family_of[c], []).append(c)
    return fam_members


# ---------------------------------------------------------------------------
# State-regularized refinement (Option B: within-family, mass-conserving)
# ---------------------------------------------------------------------------


def apply_state_regularized_refinement(
    proportions: pd.DataFrame,
    state_graph: Optional[StateSimilarityGraph] = None,
    family_map: Optional[dict[str, str]] = None,
    lambda_state: float = 0.01,
    lambda_sparse: float = 0.001,
    preserve_broad_mass: bool = True,
    rare_protection: Optional[Iterable[str]] = None,
) -> tuple[pd.DataFrame, dict]:
    """Constrained within-family refinement of fitted spot×state proportions.

    Two mass-conserving operations are applied **within each broad family**, per
    spot (never across families, never touching ``unresolved_*`` columns):

    1. **State-graph redundancy concentration** (``lambda_state``): for each
       within-family similar pair (i, j) with weight w, move
       ``lambda_state · w · (s_i − s_j)`` toward the *larger* of the two — i.e.
       concentrate mass among redundant similar states (sharpening, NOT Laplacian
       smoothing). This is mass-conserving within the family.
    2. **Sparsity shrinkage** (``lambda_sparse``): soft-threshold small within-family
       masses toward zero (protecting the dominant and any ``rare_protection``
       states) and renormalise to the original family mass.

    With ``lambda_state == 0`` and ``lambda_sparse == 0`` the input is returned
    unchanged (exact baseline). Broad-family mass per spot is conserved (deviation
    recorded). Rows remain valid compositions; outputs stay non-negative.

    Returns ``(refined_df, metadata)``.
    """
    df = proportions.copy()
    cols = list(df.columns)
    known = set(state_graph.cell_types) if state_graph is not None else None
    fmap = family_map if family_map is not None else (
        state_graph.family_map if state_graph is not None else None)
    fine_cols, unresolved_cols, family_of = _split_columns(cols, fmap, known)
    fam_members = _families_to_members(fine_cols, family_of)
    protected = {str(s) for s in (rare_protection or [])}

    P = df[fine_cols].to_numpy(dtype=np.float64) if fine_cols else np.zeros((len(df), 0))
    col_idx = {c: i for i, c in enumerate(fine_cols)}
    broad_before = {}
    mass_shifted_state = 0.0
    mass_shifted_sparse = 0.0

    if lambda_state > 0 and state_graph is not None and state_graph.n_edges:
        for fam, members in fam_members.items():
            if len(members) < 2:
                continue
            mi = [col_idx[m] for m in members]
            sub = P[:, mi]                      # (n_spots, |F|)
            fam_mass = sub.sum(axis=1, keepdims=True)
            A = state_graph.within_family_adjacency(members)  # (|F|,|F|) weights
            delta = np.zeros_like(sub)
            mm = len(members)
            for a in range(mm):
                for b in range(a + 1, mm):
                    w = A[a, b]
                    if w <= 0:
                        continue
                    move = lambda_state * w * (sub[:, a] - sub[:, b])  # toward larger
                    delta[:, a] += move
                    delta[:, b] -= move
            new = np.clip(sub + delta, 0.0, None)
            row = new.sum(axis=1, keepdims=True)
            scale = np.divide(fam_mass, np.maximum(row, 1e-12))
            new = new * scale                    # restore family mass exactly
            mass_shifted_state += float(np.abs(new - sub).sum())
            P[:, mi] = new

    if lambda_sparse > 0:
        for fam, members in fam_members.items():
            if len(members) < 2:
                continue
            mi = [col_idx[m] for m in members]
            sub = P[:, mi]
            fam_mass = sub.sum(axis=1, keepdims=True)
            thresh = lambda_sparse * fam_mass               # per-spot absolute floor
            keep_dom = np.zeros_like(sub, dtype=bool)
            dom = np.argmax(sub, axis=1)
            keep_dom[np.arange(sub.shape[0]), dom] = True    # never shrink dominant
            for a, m in enumerate(members):
                if m in protected:
                    keep_dom[:, a] = True
            shrunk = np.where(keep_dom, sub, np.clip(sub - thresh, 0.0, None))
            row = shrunk.sum(axis=1, keepdims=True)
            # avoid degenerate all-zero (only if a family somehow emptied)
            empty = row[:, 0] <= 1e-12
            scale = np.divide(fam_mass, np.maximum(row, 1e-12))
            new = shrunk * scale
            if empty.any():
                new[empty] = sub[empty]                      # fall back to baseline
            mass_shifted_sparse += float(np.abs(new - sub).sum())
            P[:, mi] = new

    # write back fine columns
    for c in fine_cols:
        df[c] = P[:, col_idx[c]]

    # broad-mass conservation check (fine + their unresolved sibling, per family)
    max_dev = 0.0
    if preserve_broad_mass:
        before = proportions
        all_fams = set(family_of.values())
        for fam in all_fams:
            members = [c for c in fine_cols if family_of[c] == fam]
            unres = [c for c in unresolved_cols if str(c) == f"unresolved_{fam}"]
            group = members + unres
            if not group:
                continue
            b = before[group].to_numpy(float).sum(axis=1)
            a = df[group].to_numpy(float).sum(axis=1)
            max_dev = max(max_dev, float(np.abs(a - b).max()))

    metadata = {
        "module": "state_regularized_refinement",
        "feature_status": FEATURE_STATUS,
        "lambda_state": float(lambda_state),
        "lambda_sparse": float(lambda_sparse),
        "preserve_broad_mass": bool(preserve_broad_mass),
        "within_family_only": True,
        "n_fine_states": len(fine_cols),
        "n_families_refined": sum(1 for m in fam_members.values() if len(m) >= 2),
        "n_unresolved_cols_preserved": len(unresolved_cols),
        "mass_shifted_state": mass_shifted_state,
        "mass_shifted_sparse": mass_shifted_sparse,
        "broad_mass_max_deviation": max_dev,
        "rare_protected": sorted(protected),
        "state_graph_edges": int(state_graph.n_edges) if state_graph is not None else 0,
    }
    return df, metadata


# ---------------------------------------------------------------------------
# Sparsity-aware refinement (sparsity only; rare-state protection)
# ---------------------------------------------------------------------------


def apply_sparsity_aware_refinement(
    proportions: pd.DataFrame,
    family_map: Optional[dict[str, str]] = None,
    state_graph: Optional[StateSimilarityGraph] = None,
    lambda_sparse: float = 0.001,
    rare_protection: Optional[Iterable[str]] = None,
    preserve_broad_mass: bool = True,
) -> tuple[pd.DataFrame, dict]:
    """Within-family sparsity shrinkage that reduces inflated effective-N.

    Soft-thresholds small within-family masses toward zero and renormalises to the
    original family mass. The per-spot dominant state and any ``rare_protection``
    states are never shrunk, so true supported rare states are preserved while
    unsupported diffuse mass is removed. Broad-family mass is conserved.

    Thin wrapper over :func:`apply_state_regularized_refinement` with
    ``lambda_state = 0`` (sparsity only).
    """
    return apply_state_regularized_refinement(
        proportions, state_graph=state_graph, family_map=family_map,
        lambda_state=0.0, lambda_sparse=lambda_sparse,
        preserve_broad_mass=preserve_broad_mass, rare_protection=rare_protection,
    )


# ---------------------------------------------------------------------------
# Adaptive resolution / state grouping (diagnostic)
# ---------------------------------------------------------------------------


def recommend_state_groups(
    state_similarity_graph: StateSimilarityGraph,
    separability_metrics: Optional[dict] = None,
    family_map: Optional[dict[str, str]] = None,
    threshold: float = 0.9,
    rare_protection: Optional[Iterable[str]] = None,
) -> dict:
    """Recommend grouping fine states that are not reliably distinguishable.

    Groups within-family state pairs whose similarity ``>= threshold`` (and, if
    ``separability_metrics`` is given as ``{(a, b): separability_score}``, pairs
    whose separability ``<= 1 - threshold``). Never groups across families; never
    auto-groups a ``rare_protection`` state. This is a **reporting/interpretation**
    recommendation — it does not modify estimates.
    """
    g = state_similarity_graph
    fam = family_map or g.family_map
    protected = {str(s) for s in (rare_protection or [])}
    ct = g.cell_types

    # union-find over high-similarity within-family edges
    parent = list(range(g.n_states))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    reasons: list[dict] = []
    idx = {c: i for i, c in enumerate(ct)}
    for e in range(g.n_edges):
        i, j = int(g.state_i[e]), int(g.state_j[e])
        a, b = ct[i], ct[j]
        if fam.get(a, a) != fam.get(b, b):
            continue
        if a in protected or b in protected:
            continue
        w = float(g.weights[e])
        sep_ok = False
        if separability_metrics is not None:
            sep = separability_metrics.get((a, b), separability_metrics.get((b, a)))
            sep_ok = sep is not None and sep <= (1.0 - threshold)
        if w >= threshold or sep_ok:
            union(i, j)
            reasons.append({"state_a": a, "state_b": b, "similarity": round(w, 4),
                            "family": fam.get(a, a),
                            "reason": "high_similarity" if w >= threshold else "low_separability"})

    comp: dict[int, list[str]] = {}
    for i in range(g.n_states):
        comp.setdefault(find(i), []).append(ct[i])
    groups = [sorted(members) for members in comp.values() if len(members) >= 2]
    groups.sort()

    grouped_states = {s for grp in groups for s in grp}
    families_unresolved = sorted({fam.get(grp[0], grp[0]) for grp in groups})

    return {
        "module": "adaptive_resolution",
        "feature_status": FEATURE_STATUS,
        "threshold": float(threshold),
        "recommended_groups": groups,
        "n_groups": len(groups),
        "n_grouped_states": len(grouped_states),
        "families_with_unresolved_fine_structure": families_unresolved,
        "protected_states": sorted(protected),
        "reasons": reasons,
        "used_for": "reporting_only",
        "used_for_estimation": False,
    }
