"""Shared two-resolution + identifiability-stratified scoring for the public benchmark.

Every method is scored on the SAME truth with the SAME functions. Broad = fine aggregated to
families via one shared fine→broad mapping. Conditional within-family and within-cluster *swap*
metrics isolate subpopulation confusion (a method can look good globally yet fail to separate
subtypes). Identifiability stratification groups per-type error by TissueResolve's recoverability
class (RESOLVABLE / WEAKLY_RESOLVABLE / GROUP_ONLY / UNRESOLVABLE / NOT_TESTABLE) — an evaluation
axis applied identically to every method, TR and external alike.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

RARE_BAND = (1e-9, 0.05)     # a present instance is "rare" if 0 < true prop <= 0.05
DETECT = 0.05                # present/detected threshold
RARE_DETECT = 0.01          # rare recall detection floor
OVERCONF = 0.15             # predicted mass above this on a wrong subtype = an overconfident call


def _align(truth, pred, cols):
    t = truth.reindex(columns=cols).fillna(0.0)
    p = pred.reindex(index=truth.index, columns=cols).fillna(0.0)
    return t, p


def _acc(t, p):
    tv, pv = t.to_numpy(float).ravel(), p.to_numpy(float).ravel()
    out = {"rmse": float(np.sqrt(np.mean((tv - pv) ** 2))),
           "mae": float(np.mean(np.abs(tv - pv))),
           "bias": float(np.mean(pv - tv))}
    out["pearson"] = float(np.corrcoef(tv, pv)[0, 1]) if tv.std() > 0 and pv.std() > 0 else np.nan
    from scipy.stats import spearmanr
    out["spearman"] = float(spearmanr(tv, pv).correlation) if tv.std() > 0 and pv.std() > 0 else np.nan
    return out


def _jsd(t, p):
    from benchmarks.shared import metrics as M
    vals = [M.jensen_shannon_divergence(t.iloc[i].to_numpy(float), p.iloc[i].to_numpy(float))
            for i in range(len(t))]
    return float(np.nanmean(vals))


def _dominant(t, p):
    return float((t.values.argmax(1) == p.values.argmax(1)).mean())


def unresolved_mass(pred):
    """Fraction of predicted mass a method leaves in unresolved_<family> columns (abstention).
    Zero for methods that always assign full fine mass; >0 only for resolution-aware TR-hierarchical."""
    uc = [c for c in pred.columns if str(c).startswith("unresolved_")]
    if not uc:
        return 0.0
    return float(pred[uc].to_numpy(float).sum(1).mean())


def broad_metrics(truth, pred, mapping):
    """Broad accuracy. Unresolved_<family> mass is CREDITED to its family (fair to resolution-aware
    abstention: hierarchical is confident at broad even when it abstains at fine)."""
    from benchmarks.shared import metrics as M
    pred = pred.reindex(index=truth.index)
    tb = M.aggregate_to_families(truth, mapping)
    fine_cols = [c for c in pred.columns if not str(c).startswith("unresolved_")]
    pb = M.aggregate_to_families(pred[fine_cols], mapping)
    for c in pred.columns:
        if str(c).startswith("unresolved_"):
            fam = str(c)[len("unresolved_"):]
            if fam in pb.columns:
                pb[fam] = pb[fam] + pred[c].to_numpy(float)
    cols = list(tb.columns)
    tb, pb = _align(tb, pb, cols)
    m = _acc(tb, pb)
    m.update({"jsd": _jsd(tb, pb), "dominant_acc": _dominant(tb, pb), "n_families": len(cols)})
    return {f"broad_{k}": v for k, v in m.items()}


def fine_metrics(truth, pred):
    cols = list(truth.columns)
    t, p = _align(truth, pred, cols)
    m = _acc(t, p)
    m.update({"jsd": _jsd(t, p), "dominant_acc": _dominant(t, p), "n_types": len(cols)})
    return {f"fine_{k}": v for k, v in m.items()}


def _conditional(truth, pred, groups, present_thresh=DETECT):
    """Within-group conditional composition (member/group-sum) RMSE & Pearson, pooled over groups."""
    T_all, P_all = [], []
    for members in groups:
        members = [m for m in members if m in truth.columns]
        if len(members) < 2:
            continue
        T = truth.reindex(columns=members).fillna(0).to_numpy(float)
        P = pred.reindex(index=truth.index, columns=members).fillna(0).to_numpy(float)
        ts, ps = T.sum(1), P.sum(1)
        present = ts > present_thresh
        if present.sum() < 1:
            continue
        T_all.append((T[present] / ts[present, None]).ravel())
        P_all.append((P[present] / np.clip(ps[present, None], 1e-9, None)).ravel())
    if not T_all:
        return {"rmse": np.nan, "pearson": np.nan, "n_groups": 0}
    tv, pv = np.concatenate(T_all), np.concatenate(P_all)
    return {"rmse": float(np.sqrt(np.mean((tv - pv) ** 2))),
            "pearson": float(np.corrcoef(tv, pv)[0, 1]) if tv.std() > 0 and pv.std() > 0 else np.nan,
            "n_groups": len(T_all)}


def conditional_within_family(truth, pred, mapping):
    fam = defaultdict(list)
    for c in truth.columns:
        fam[str(mapping.get(c, c))].append(c)
    r = _conditional(truth, pred, [m for m in fam.values() if len(m) > 1])
    return {f"condfam_{k}": v for k, v in r.items()}


def swap_within_clusters(truth, pred, clusters):
    r = _conditional(truth, pred, clusters)
    return {f"swap_{k}": v for k, v in r.items()}


def rare_and_absent(truth, pred):
    cols = list(truth.columns)
    t, p = _align(truth, pred, cols)
    T, P = t.to_numpy(float), p.to_numpy(float)
    rare = (T > RARE_BAND[0]) & (T <= RARE_BAND[1])
    present = T > DETECT
    absent = T <= 1e-9
    return {
        "rare_recall": float((P[rare] > RARE_DETECT).mean()) if rare.any() else np.nan,
        "rare_precision": float((T[P > RARE_DETECT] > RARE_BAND[0]).mean()) if (P > RARE_DETECT).any() else np.nan,
        "rare_fpr": float((P[absent] > RARE_DETECT).mean()) if absent.any() else np.nan,
        "detection_recall": float((P[present] > DETECT).mean()) if present.any() else np.nan,
        "absent_mass": float(P[absent].mean()) if absent.any() else np.nan,
        "false_positive_subtype_rate": float((P[absent] > DETECT).mean()) if absent.any() else np.nan,
        "eff_n_pred": float(np.mean([1.0 / np.sum(np.square(row / row.sum())) if row.sum() > 0 else 0
                                     for row in P])),
        "eff_n_true": float(np.mean([1.0 / np.sum(np.square(row / row.sum())) if row.sum() > 0 else 0
                                     for row in T])),
    }


def identifiability_stratified(truth, pred, class_by_type):
    """Per-type held-out RMSE grouped by TR recoverability class + overconfident-error / swap-risk."""
    cols = list(truth.columns)
    t, p = _align(truth, pred, cols)
    by_class = defaultdict(list)
    overconf_wrong = defaultdict(list)
    for j, c in enumerate(cols):
        cls = class_by_type.get(c, "NOT_TESTABLE")
        tv, pv = t[c].to_numpy(float), p[c].to_numpy(float)
        by_class[cls].append(float(np.sqrt(np.mean((pv - tv) ** 2))))
        if cls in ("GROUP_ONLY", "UNRESOLVABLE"):
            # overconfident: confidently predicts this individual subtype when it is not truly dominant
            wrong = (pv > OVERCONF) & (tv < pv - 0.05)
            overconf_wrong[cls].append(float(wrong.mean()))
    out = {}
    for cls in ("RESOLVABLE", "WEAKLY_RESOLVABLE", "GROUP_ONLY", "UNRESOLVABLE", "NOT_TESTABLE"):
        out[f"err_{cls}"] = float(np.mean(by_class[cls])) if by_class[cls] else np.nan
        out[f"n_{cls}"] = len(by_class[cls])
    oc = [x for v in overconf_wrong.values() for x in v]
    out["overconfident_error_rate"] = float(np.mean(oc)) if oc else np.nan
    return out


def score_all(truth, pred, mapping, clusters, class_by_type):
    """Full metric row for one (method, scenario) prediction."""
    row = {}
    row.update(broad_metrics(truth, pred, mapping))
    row.update(fine_metrics(truth, pred))
    row.update(conditional_within_family(truth, pred, mapping))
    row.update(swap_within_clusters(truth, pred, clusters))
    row.update(rare_and_absent(truth, pred))
    row.update(identifiability_stratified(truth, pred, class_by_type))
    row["unresolved_mass"] = unresolved_mass(pred)
    return row


def per_celltype_rows(truth, pred, class_by_type, mapping):
    cols = list(truth.columns)
    t, p = _align(truth, pred, cols)
    rows = []
    for c in cols:
        tv, pv = t[c].to_numpy(float), p[c].to_numpy(float)
        present = tv > DETECT; absent = tv <= 1e-9
        rows.append(dict(cell_type=c, family=str(mapping.get(c, c)),
                         recoverability=class_by_type.get(c, "NOT_TESTABLE"),
                         rmse=float(np.sqrt(np.mean((pv - tv) ** 2))),
                         mean_true=float(tv.mean()),
                         detection_recall=float((pv[present] > DETECT).mean()) if present.any() else np.nan,
                         absent_mass=float(pv[absent].mean()) if absent.any() else np.nan))
    return rows
