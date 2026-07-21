"""Identifiability certificate for reference-based deconvolution (Stage 2A, experimental).

Before deconvolving, ask what is *recoverable*. Deconvolution solves y ≈ R θ (R = reference
mean profiles restricted to query-detectable genes; θ = proportions). The recoverable
information is the projection of θ onto the well-conditioned, above-noise subspace of R; the
rest is structurally indistinguishable no matter the solver.

This module computes, from R restricted to query-detectable genes plus a query noise model:
  * singular-value spectrum, condition number, (stable + count) effective rank;
  * NOISE-ADJUSTED effective rank: directions recoverable above the query's counting-noise floor;
  * near-null right-singular directions and the cell types confounded along each;
  * per-type recoverability class + confounded groups + recommended merges;
  * a theoretical_detection_floor per direction/type (minimal detectable abundance shift).

It changes NO estimate and NO default. It reuses the reference means (as_R_cpm) — no new distance
metric. Outputs are RNA-proportion-space statements about what the (reference, query) pair can know.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["identifiability_certificate", "IdentifiabilityCertificate", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "identifiability_certificate-0.1.0"


@dataclass
class IdentifiabilityCertificate:
    cell_types: list
    singular_values: np.ndarray
    condition_number: float
    stable_rank: float
    effective_rank_count: int
    noise_adjusted_rank: int
    per_type: pd.DataFrame          # class, rec_fraction, detection_floor, closest_confounders, ...
    null_space: pd.DataFrame        # direction, singular_value, detection_floor, involved_types
    confounded_groups: list         # list of {group, recoverable_group, reason}
    query_gene_retention: dict
    metadata: dict = field(default_factory=dict)

    def save(self, out_dir):
        out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
        self.per_type.to_csv(out / "identifiability_certificate.tsv", sep="\t", index=False)
        self.per_type.to_csv(out / "recoverability_by_celltype.tsv", sep="\t", index=False)
        pd.DataFrame({"index": range(len(self.singular_values)),
                      "singular_value": self.singular_values}).to_csv(
            out / "singular_values.tsv", sep="\t", index=False)
        self.null_space.to_csv(out / "null_space_components.tsv", sep="\t", index=False)
        pd.DataFrame(self.confounded_groups or [{"group": None}]).to_csv(
            out / "confounded_celltype_groups.tsv", sep="\t", index=False)
        pd.DataFrame([self.query_gene_retention]).to_csv(
            out / "query_gene_retention.tsv", sep="\t", index=False)
        (out / "identifiability_manifest.json").write_text(json.dumps(
            {"condition_number": self.condition_number, "stable_rank": self.stable_rank,
             "effective_rank_count": self.effective_rank_count,
             "noise_adjusted_rank": self.noise_adjusted_rank,
             "query_gene_retention": self.query_gene_retention,
             "algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
             **self.metadata}, indent=2))


def identifiability_certificate(
    R_cpm: np.ndarray,                 # (K types, G genes) — reference means (as_R_cpm())
    cell_types: Sequence[str],
    gene_names: Sequence[str],
    *,
    query_detectable_genes: Optional[Sequence[str]] = None,
    library_size: float = 1e7,
    z: float = 2.0,                    # detection SNR
    detect_floor: float = 0.10,       # a direction is "recoverable" if detectable at ≤10% abundance shift
    null_sv_tol: float = 1e-3,        # relative singular-value tolerance for structural null
    involve_tol: float = 0.30,        # |v_i[k]| above this ⇒ type k involved in direction i
    rec_fraction_resolvable: float = 0.80,
    rec_fraction_weak: float = 0.50,
    n_donors_by_type: Optional[dict] = None,
    min_donors_testable: int = 2,
) -> IdentifiabilityCertificate:
    """Compute the identifiability certificate on query-detectable genes only."""
    R = np.asarray(R_cpm, dtype=float)                 # (K, G)
    K = R.shape[0]
    genes = [str(g) for g in gene_names]
    gset = set(genes)
    if query_detectable_genes is not None:
        keep = [g for g in query_detectable_genes if g in gset]
        idx = [genes.index(g) for g in keep]
        retention = {"n_reference_genes": len(genes), "n_query_detectable": len(set(query_detectable_genes)),
                     "n_retained": len(idx), "retention_frac": round(len(idx) / max(len(genes), 1), 4)}
    else:
        idx = list(range(len(genes)))
        retention = {"n_reference_genes": len(genes), "n_query_detectable": len(genes),
                     "n_retained": len(genes), "retention_frac": 1.0}
    Rq = R[:, idx].T                                   # (G_q, K) genes × types

    # --- SVD in θ-space: y = Rq θ, right singular vectors live in ℝ^K ---
    if Rq.shape[0] < 2 or K < 1:
        raise ValueError("too few query-detectable genes or cell types for a certificate")
    U, S, Vt = np.linalg.svd(Rq, full_matrices=False)  # Rq = U diag(S) Vt ; S len min(G_q,K)
    V = Vt.T                                            # (K, r) right singular vectors (columns)
    r = len(S)
    smax = float(S[0]) if r else 0.0
    cond = float(smax / S[-1]) if r and S[-1] > 0 else float("inf")
    stable_rank = float((S.sum() ** 2) / (np.square(S).sum())) if r else 0.0
    eff_rank_count = int((S > null_sv_tol * smax).sum())

    # --- noise-aware detection floor per singular direction (matched filter) ---
    pbar = R[:, idx].mean(axis=0)                      # mean bulk profile (CPM), G_q
    b = np.clip(library_size * pbar / 1e6, 1e-9, None) # baseline expected counts
    floors = np.empty(r)
    for i in range(r):
        s_i = library_size * (Rq @ V[:, i]) / 1e6      # count change per unit α along v_i
        snr1 = np.sqrt(np.sum(s_i ** 2 / b))           # SNR per unit abundance shift
        floors[i] = float(z / snr1) if snr1 > 0 else np.inf
    noise_adj_rank = int((floors < detect_floor).sum())
    recoverable = floors < detect_floor                # per-direction

    # --- recoverable-subspace projector ---
    Vrec = V[:, recoverable]
    P_rec = Vrec @ Vrec.T if Vrec.shape[1] else np.zeros((K, K))

    # --- near-null (confounded) directions + involved types ---
    null_rows, edges = [], []
    for i in range(r):
        if not recoverable[i]:
            involved = [cell_types[k] for k in range(K) if abs(V[k, i]) > involve_tol]
            if len(involved) >= 2:
                null_rows.append({"direction": i, "singular_value": round(float(S[i]), 4),
                                  "detection_floor": round(float(floors[i]), 4),
                                  "involved_types": ";".join(involved)})
                for a in range(len(involved)):
                    for bb in range(a + 1, len(involved)):
                        edges.append((involved[a], involved[bb]))
    # confounded groups = connected components over shared near-null directions
    groups = _connected_components(list(cell_types), edges)
    groups = [g for g in groups if len(g) >= 2]

    # --- per-type recoverability ---
    rows = []
    ct_index = {c: k for k, c in enumerate(cell_types)}
    group_of = {c: g for g in groups for c in g}
    for k, c in enumerate(cell_types):
        e = np.zeros(K); e[k] = 1.0
        rec_frac = float(e @ P_rec @ e)                # fraction of type-identity that is recoverable
        floor_k = _type_floor(e, V, floors)
        n_don = (n_donors_by_type or {}).get(c, None)
        if retention["n_retained"] < 2 or (n_don is not None and n_don < min_donors_testable):
            cls, reason = "NOT_TESTABLE", "too few query genes or reference donors"
            conf, merge, grp = "", "", ""
        elif rec_frac >= rec_fraction_resolvable and c not in group_of:
            cls, reason, conf, merge, grp = "RESOLVABLE", "well-conditioned & above noise", "", "", ""
        elif c in group_of:
            grp_members = group_of[c]
            grp_rec = _group_recoverable(grp_members, ct_index, P_rec)
            conf = ";".join([x for x in grp_members if x != c])
            grp = ";".join(sorted(grp_members))
            if grp_rec:
                cls, reason, merge = "GROUP_ONLY", "individual axis near-null; group sum recoverable", grp
            else:
                cls, reason, merge = "UNRESOLVABLE", "type and group both near-null", grp
        elif rec_frac >= rec_fraction_weak:
            cls, reason, conf, merge, grp = "WEAKLY_RESOLVABLE", "partially above noise", "", "", ""
        else:
            cls, reason, conf, merge, grp = "UNRESOLVABLE", "identity in near-null subspace", "", "", ""
        rows.append({"cell_type": c, "recoverability": cls, "rec_fraction": round(rec_frac, 4),
                     "theoretical_detection_floor": round(floor_k, 4),
                     "closest_confounders": conf, "recommended_merge": merge,
                     "recoverable_group": grp, "n_donors": n_don, "reason": reason})
    per_type = pd.DataFrame(rows)
    conf_groups = [{"group": ";".join(sorted(g)),
                    "recoverable_group": _group_recoverable(g, ct_index, P_rec),
                    "reason": "shared near-null direction(s)"} for g in groups]
    meta = {"library_size": library_size, "z": z, "detect_floor": detect_floor,
            "n_types": K, "algorithm_version": ALGORITHM_VERSION}
    return IdentifiabilityCertificate(
        cell_types=list(cell_types), singular_values=S, condition_number=cond,
        stable_rank=stable_rank, effective_rank_count=eff_rank_count,
        noise_adjusted_rank=noise_adj_rank, per_type=per_type,
        null_space=pd.DataFrame(null_rows, columns=["direction", "singular_value",
                                                    "detection_floor", "involved_types"]),
        confounded_groups=conf_groups,
        query_gene_retention=retention, metadata=meta)


def _type_floor(e, V, floors):
    """Detection floor for a type = the floor of its least-recoverable participating direction."""
    contrib = np.abs(V.T @ e)                           # |projection of e_k on each direction|
    signif = contrib > 0.1
    if not signif.any():
        return float("inf")
    return float(np.max(floors[signif]))


def _group_recoverable(group, ct_index, P_rec):
    """Is the group's SUM (aggregate proportion) recoverable? (uniform direction over the group)."""
    K = P_rec.shape[0]
    g = np.zeros(K)
    for c in group:
        g[ct_index[c]] = 1.0
    if g.sum() == 0:
        return False
    g = g / np.linalg.norm(g)
    return bool(float(g @ P_rec @ g) >= 0.5)


def _connected_components(nodes, edges):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        parent[find(a)] = find(b)
    comp = {}
    for n in nodes:
        comp.setdefault(find(n), []).append(n)
    return list(comp.values())
