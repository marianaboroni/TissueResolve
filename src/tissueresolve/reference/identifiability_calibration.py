"""Calibrated identifiability certificate (Stage 2B, experimental).

The linear certificate (``identifiability.py``) is an *unconstrained* SVD analysis: a conservative
lower bound that over-flags confounding the nonnegative solver actually resolves, and whose
recoverability outcome (absolute RMSE) is abundance-dominated. This module calibrates it along five
axes, using forward simulation from the reference itself (no query labels):

  1. NONNEGATIVE / SIMPLEX-AWARE — recoverability is measured by running the *nonnegative* production
     solver on simulated mixtures, so nonnegativity is honoured (geometry proposes candidate
     confusable clusters; the constrained solve disposes).
  2. SWAP-BASED OUTCOME — the recoverability outcome is the WITHIN-cluster conditional (swap) RMSE
     (member / cluster-sum), which isolates member confusion and is abundance-normalised — not the
     abundance-dominated absolute RMSE.
  3. DONOR-LEVEL UNCERTAINTY — the data-generating profile is resampled per draw using the reference
     cross-donor CV (``donor_cv``), so the predicted error SD reflects donor/biological variation
     (the dominant error source), not just counting noise.
  4. DEPTH / PANEL — ``library_size`` and ``query_detectable_genes`` parameterise the query regime.
  5. (tissue-agnostic: takes any ReferenceSignature.)

It is EXPERIMENTAL, changes no default, and adds no new solver — it reuses the existing production
solver in a forward simulation for analysis only.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from tissueresolve.reference.identifiability import _connected_components

__all__ = ["calibrated_identifiability_certificate", "CalibratedCertificate", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "calibrated_identifiability_certificate-0.1.0"

# Donor-shift inflation that gives ~nominal (90%) held-out coverage WITHOUT a per-run calibration
# split. Validated on breast + lung (Stage 2C, benchmarks/signatures/run_stage2c_default_shift.py);
# a per-run calibration split (Stage 2B) may refine it. Re-audit when a third tissue is available.
DEFAULT_SHIFT_SCALE = 2.0


@dataclass
class CalibratedCertificate:
    cell_types: list
    per_type: pd.DataFrame          # recoverability, sim_rmse, sim_bias, pred_sd, swap_rmse, cluster, ...
    clusters: pd.DataFrame          # candidate confusable clusters + swap RMSE + group-sum RMSE + verdict
    predicted: pd.DataFrame         # per-type predicted uncertainty (pred_sd) for held-out interval building
    metadata: dict = field(default_factory=dict)

    def save(self, out_dir):
        out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
        self.per_type.to_csv(out / "calibrated_recoverability_by_celltype.tsv", sep="\t", index=False)
        self.clusters.to_csv(out / "calibrated_confusable_clusters.tsv", sep="\t", index=False)
        self.predicted.to_csv(out / "predicted_uncertainty_by_celltype.tsv", sep="\t", index=False)
        (out / "calibrated_identifiability_manifest.json").write_text(json.dumps(
            {"algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
             **self.metadata}, indent=2, default=str))


def _candidate_clusters(Rq: np.ndarray, cell_types, tau: float):
    """Connected components of the graph {(a,b): cosine(R_a, R_b) >= tau} over query-detectable genes."""
    norm = Rq / (np.linalg.norm(Rq, axis=1, keepdims=True) + 1e-12)
    C = norm @ norm.T
    K = len(cell_types)
    edges = [(cell_types[i], cell_types[j]) for i in range(K) for j in range(i + 1, K) if C[i, j] >= tau]
    groups = [g for g in _connected_components(list(cell_types), edges) if len(g) >= 2]
    return groups, C


def _draws(K, clusters, cell_types, n_sim, rng):
    """Dirichlet draws over all types + cluster-stress draws that split mass within each cluster."""
    thetas = [rng.dirichlet(np.full(K, 0.5)) for _ in range(n_sim)]
    idx = {c: i for i, c in enumerate(cell_types)}
    per = max(8, n_sim // max(len(clusters), 1)) if clusters else 0
    for g in clusters:
        members = [idx[c] for c in g]
        for _ in range(per):
            th = rng.dirichlet(np.full(K, 0.05)) * 0.1
            a, b = rng.choice(members, size=2, replace=False)
            th[a] += rng.uniform(0.2, 0.6); th[b] += rng.uniform(0.2, 0.6)
            thetas.append(th / th.sum())
    return np.array(thetas)


def _swap_rmse(Theta, ThetaHat, member_idx, present_thresh=0.05):
    """Within-cluster conditional (swap) RMSE and group-sum RMSE for one cluster."""
    T = Theta[:, member_idx]; P = ThetaHat[:, member_idx]
    ts, ps = T.sum(1), P.sum(1)
    group_rmse = float(np.sqrt(np.mean((ps - ts) ** 2)))
    present = ts > present_thresh
    if present.sum() < 3:
        return np.nan, group_rmse, int(present.sum())
    tc = T[present] / ts[present, None]
    pc = P[present] / np.clip(ps[present, None], 1e-9, None)
    return float(np.sqrt(np.mean((tc - pc) ** 2))), group_rmse, int(present.sum())


def calibrated_identifiability_certificate(
    ref,
    *,
    query_detectable_genes: Optional[Sequence[str]] = None,
    library_size: float = 1e7,
    n_sim: int = 200,
    noise: str = "poisson",
    donor_uncertainty: bool = True,
    cv_clip: float = 2.0,
    shift_scale: float = 1.0,          # inflation of donor CV to match train→held-out shift
                                       # (calibrate on a donor-disjoint calibration split; see Stage 2B)
    tau_cluster: float = 0.90,
    swap_resolvable: float = 0.10,     # cluster swap RMSE below this ⇒ nonnegativity resolves it
    swap_confounded: float = 0.15,     # at/above this ⇒ within-cluster split is not recoverable
    rmse_resolvable: float = 0.10,     # standalone-type simulated RMSE thresholds
    rmse_weak: float = 0.20,
    group_sum_recoverable: float = 0.12,
    n_donors_by_type: Optional[dict] = None,
    min_donors_testable: int = 2,
    solver: str = "poisson",
    seed: int = 0,
) -> CalibratedCertificate:
    """Forward-simulation, nonnegative, swap-based, donor-uncertainty-aware identifiability."""
    from tissueresolve.solver import PoissonGLMSolver, NNLSSolver

    R = np.asarray(ref.as_R_cpm(), dtype=float)            # (K, G)
    K = R.shape[0]
    cell_types = [str(c) for c in ref.cell_types]
    genes = [str(g) for g in ref.gene_names]
    gset = set(genes)
    if query_detectable_genes is not None:
        detect = [g for g in query_detectable_genes if g in gset]
    else:
        detect = list(genes)
    didx = [genes.index(g) for g in detect]
    if len(didx) < 2:
        raise ValueError("need >=2 query-detectable genes for a calibrated certificate")
    Rq = R[:, didx]                                        # (K, Gq)

    donor_cv = getattr(ref, "donor_cv", None)
    have_cv = donor_uncertainty and donor_cv is not None
    cvq = None
    if have_cv:
        cvq = np.nan_to_num(np.asarray(donor_cv)[didx, :].T, nan=0.0)   # (K, Gq)
        cvq = np.clip(cvq, 0.0, cv_clip) * float(shift_scale)

    clusters, cos = _candidate_clusters(Rq, cell_types, tau_cluster)
    rng = np.random.default_rng(seed)
    Theta = _draws(K, clusters, cell_types, n_sim, rng)    # (S, K)
    S = Theta.shape[0]

    # generate simulated query counts (donor-perturbed generating profile, then counting noise).
    # Donor variation is modelled as multiplicative lognormal with the target CV: sigma from the CV
    # (sigma = sqrt(log1p(CV^2))), median-corrected so E[R_i] = R (no upward bias).
    sigma = np.sqrt(np.log1p(cvq ** 2)) if have_cv else None
    counts = np.empty((len(didx), S))
    for i in range(S):
        Ri = Rq
        if have_cv:
            z = rng.normal(0.0, 1.0, Rq.shape)
            Ri = Rq * np.exp(z * sigma - 0.5 * sigma ** 2)
        ycpm = Theta[i] @ Ri                                # (Gq,)
        mu = np.clip(library_size * ycpm / 1e6, 0.0, library_size)   # a type cannot exceed total depth
        if noise == "nb":
            phi = getattr(ref, "phi_g", None)
            r = 10.0 if phi is None else float(np.median(1.0 / np.clip(np.asarray(phi)[didx], 1e-3, None)))
            counts[:, i] = rng.negative_binomial(r, r / (r + mu))
        else:
            counts[:, i] = rng.poisson(mu)

    dfc = pd.DataFrame(counts, index=detect, columns=[f"s{i}" for i in range(S)])
    Solver = {"poisson": PoissonGLMSolver, "nnls": NNLSSolver}[solver]
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ThetaHat = (Solver(genes=detect).solve(dfc, ref).proportions
                    .reindex(columns=cell_types).fillna(0.0).to_numpy(float))

    # per-type simulated recoverability + donor-aware predicted uncertainty
    resid = ThetaHat - Theta
    sim_rmse = np.sqrt(np.mean(resid ** 2, axis=0))
    sim_bias = np.mean(resid, axis=0)
    pred_sd = np.std(resid, axis=0)

    # cluster swap analysis
    idx = {c: i for i, c in enumerate(cell_types)}
    cl_rows, member_cluster, cluster_swap, cluster_group_rmse = [], {}, {}, {}
    for gi, g in enumerate(clusters):
        members = [idx[c] for c in g]
        swap, grp, npres = _swap_rmse(Theta, ThetaHat, members)
        confounded = (swap == swap) and swap >= swap_confounded
        grp_ok = grp < group_sum_recoverable
        verdict = ("RESOLVED_BY_NONNEG" if (swap == swap and swap < swap_resolvable)
                   else "GROUP_ONLY" if (confounded and grp_ok)
                   else "UNRESOLVABLE" if confounded else "PARTIAL")
        cl_rows.append(dict(cluster_id=gi, group=";".join(g), n_members=len(g),
                            max_cosine=round(float(cos[np.ix_(members, members)][np.triu_indices(len(members), 1)].max()), 4)
                            if len(members) > 1 else np.nan,
                            swap_rmse=round(swap, 4) if swap == swap else np.nan,
                            group_sum_rmse=round(grp, 4), n_present=npres, verdict=verdict))
        for c in g:
            member_cluster[c] = ";".join(g); cluster_swap[c] = swap; cluster_group_rmse[c] = grp

    rows = []
    for k, c in enumerate(cell_types):
        n_don = (n_donors_by_type or {}).get(c)
        if n_don is not None and n_don < min_donors_testable:
            cls, reason = "NOT_TESTABLE", "too few reference donors"
        elif c in member_cluster:
            swap = cluster_swap[c]; grp = cluster_group_rmse[c]
            if swap == swap and swap < swap_resolvable:
                cls, reason = "RESOLVABLE", "in confusable cluster but nonnegativity resolves the split"
            elif swap == swap and swap >= swap_confounded:
                cls = "GROUP_ONLY" if grp < group_sum_recoverable else "UNRESOLVABLE"
                reason = "within-cluster split not recoverable (group sum recoverable)" if grp < group_sum_recoverable \
                    else "within-cluster split and group sum both unrecoverable"
            else:
                cls, reason = "WEAKLY_RESOLVABLE", "partial within-cluster recovery"
        else:
            cls = ("RESOLVABLE" if sim_rmse[k] < rmse_resolvable
                   else "WEAKLY_RESOLVABLE" if sim_rmse[k] < rmse_weak else "UNRESOLVABLE")
            reason = "standalone type; simulated recovery " + (
                "good" if cls == "RESOLVABLE" else "partial" if cls == "WEAKLY_RESOLVABLE" else "poor")
        # recommended_merge: for confounded (GROUP_ONLY/UNRESOLVABLE) types in a cluster, report the
        # cluster members whose individual split is unrecoverable — read out as an aggregate instead.
        merge = member_cluster.get(c, "") if cls in ("GROUP_ONLY", "UNRESOLVABLE") and c in member_cluster else ""
        rows.append(dict(cell_type=c, recoverability=cls,
                         sim_rmse=round(float(sim_rmse[k]), 4), sim_bias=round(float(sim_bias[k]), 4),
                         pred_sd=round(float(pred_sd[k]), 4),
                         swap_rmse=round(float(cluster_swap[c]), 4) if c in cluster_swap and cluster_swap[c] == cluster_swap[c] else np.nan,
                         cluster=member_cluster.get(c, ""), recommended_merge=merge,
                         n_donors=n_don, reason=reason))
    per_type = pd.DataFrame(rows)
    predicted = per_type[["cell_type", "pred_sd", "sim_bias"]].copy()
    meta = {"library_size": library_size, "n_sim_total": int(S), "n_sim_base": int(n_sim),
            "noise": noise, "donor_uncertainty": bool(have_cv), "shift_scale": float(shift_scale),
            "tau_cluster": tau_cluster,
            "n_detectable": len(didx), "n_clusters": len(clusters), "solver": solver,
            "swap_resolvable": swap_resolvable, "swap_confounded": swap_confounded,
            "algorithm_version": ALGORITHM_VERSION}
    return CalibratedCertificate(cell_types=cell_types, per_type=per_type,
                                 clusters=pd.DataFrame(cl_rows, columns=[
                                     "cluster_id", "group", "n_members", "max_cosine", "swap_rmse",
                                     "group_sum_rmse", "n_present", "verdict"]),
                                 predicted=predicted, metadata=meta)
