"""Compact high-granularity subtype strategy (experimental).

EXPERIMENTAL — not imported by default, not wired into any pipeline.  Tests
whether high-granularity fine-subpopulation prediction can be improved *without*
amplifying fragile subtype contrasts (which the FineGranularityRefiner showed
worsens spillover / false positives), by instead making the reference more
query-compatible and restricting the candidate subtype space, then taking a
consensus over robust panels, BEFORE the validated partial confidence-weighted
soft gate.

Strategy (soft gating applied exactly once, AFTER this procedure):

    π_f (broad mass, unchanged)
    → q^raw_{k|f}                              current fine solve
    → ReferenceCalibration (gene-level, conservative; preserves contrasts)
    → CandidateRestriction (plausible subtypes per family; never empties a family)
    → multi-panel fine deconvolution on the retained candidates (3–5 panels)
    → ConsensusStability (mean/median/weighted/conservative; stays on the simplex)
    → stability/candidate-aware confidence  c_k
    → θ_k = π_f q^consensus_{k|f};  θ^resolved_k = θ_k c_k;  u_f = π_f − Σ θ^resolved_k

Only the within-family conditional is changed; broad-family and total mass are
conserved exactly.  Unsupported families fall back to the raw conditional.

Reuses the pure helpers from ``fine_refiner`` (simplex projection, family
sub-reference, within-family WNNLS, donor residuals) to keep the file footprint
minimal.  No PyTorch/JAX/Pyro/NumPyro/OT/NB-CAR-VI, no deep learning, no L1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd

from .fine_refiner import (
    project_to_simplex, family_program_and_residual, contrast_weights,
    _subset_ref_celltypes, _wnnls_family_conditional,
)

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "high_granularity-0.1.0"

PRESENT = 0.01          # within-family conditional "present" threshold
RARE_RAW = 0.02         # protect subtypes whose raw conditional ≥ this (rare-but-real)


# ===========================================================================
# Config / result
# ===========================================================================
@dataclass
class HighGranularityConfig:
    """Settings for the high-granularity strategy (all components ablatable)."""
    # top-level mode
    mode: str = "none"   # none|candidate_consensus|reference_calibrated_consensus|auto
    # Component 1 — reference calibration
    calibration: str = "none"          # none|scale|logscale|residual_downweight
    scale_lo: float = 0.5
    scale_hi: float = 2.0
    # Component 2 — candidate restriction
    candidate_rule: str = "none"       # none|topk|evidence|hybrid
    top_k: int = 3
    min_candidates: int = 1
    max_candidates: int = 6
    min_raw_conditional: float = 0.05
    min_query_detectability: float = 0.10
    # Component 3 — panels
    panels: tuple = ("current", "donor_stable", "query_detectable", "marker")
    panel_size: int = 200
    # Component 4 — consensus
    consensus: str = "mean"            # mean|median|weighted|conservative
    min_panel_detection: int = 2       # conservative: subtype must appear in ≥ M panels
    # Component 5 — stability/candidate-aware confidence
    stability_cv_tol: float = 0.5      # CV above this → downweight subtype confidence
    # support gating (fallback)
    min_subtypes: int = 2
    min_family_donors: int = 3
    eps: float = 1e-9
    feature_status: str = FEATURE_STATUS
    version: str = ALGORITHM_VERSION


@dataclass
class HighGranularityResult:
    refined_conditional: pd.DataFrame                 # samples × subtypes; within-family Σ=1
    raw_conditional: pd.DataFrame
    stability_confidence: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    retained: dict = field(default_factory=dict)      # family -> [subtypes]
    dropped: dict = field(default_factory=dict)       # family -> [subtypes]
    refined_families: list = field(default_factory=list)
    skipped_families: dict = field(default_factory=dict)
    diagnostics: pd.DataFrame = field(default_factory=pd.DataFrame)
    calibration_status: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


# ===========================================================================
# Component 1 — reference calibration (conservative, contrast-preserving)
# ===========================================================================
def calibrate_gene_scale(sub_ref, query: pd.DataFrame, cfg: HighGranularityConfig,
                         eval_genes: Optional[list] = None):
    """Conservative per-gene scaling toward the query magnitude.

    ``S*_{gk} = a_g S_{gk}`` with ``a_g = clip(query_cpm_g / ref_mean_cpm_g)``.
    Scaling is identical across subtypes within a gene, so within-family subtype
    *contrasts* (ratios) are preserved exactly; only the cross-gene weighting
    shifts toward the query.  Returns (calibrated_sub_ref, diagnostics).
    """
    from tissueresolve.results import ReferenceSignature
    genes = [str(g) for g in sub_ref.gene_names]
    R = sub_ref.as_R_cpm().astype(float)            # (k, G)
    ref_mean = R.mean(axis=0) + cfg.eps             # family mean signature (G,)
    q = query.reindex(index=genes).fillna(0.0)
    q_cpm = (q.div(q.sum(axis=0).where(lambda s: s > 0, np.nan), axis=1).fillna(0.0)
             * 1e6).mean(axis=1).to_numpy()          # mean query CPM per gene (G,)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(ref_mean > cfg.eps, q_cpm / ref_mean, 1.0)
    a = np.clip(np.nan_to_num(a, nan=1.0, posinf=cfg.scale_hi, neginf=1.0),
                cfg.scale_lo, cfg.scale_hi)
    Rc = (R * a[None, :]).astype(np.float32)
    cal = ReferenceSignature(
        gene_names=genes, cell_types=list(sub_ref.cell_types),
        R_cpm=Rc, R_log=np.log1p(Rc).astype(np.float32),
        genome=sub_ref.genome,
        n_cells_per_type={m: sub_ref.n_cells_per_type.get(m, 0) for m in sub_ref.cell_types})
    # held-out reconstruction improvement (on eval_genes if given)
    ev = [g for g in (eval_genes or genes) if g in genes]
    gi = [genes.index(g) for g in ev]
    base_err = float(np.mean(np.abs(ref_mean[gi] - q_cpm[gi]))) if gi else float("nan")
    cal_err = float(np.mean(np.abs((R * a[None, :]).mean(axis=0)[gi] - q_cpm[gi]))) if gi else float("nan")
    diag = {"method": "scale", "n_genes_modified": int((np.abs(a - 1.0) > 1e-6).sum()),
            "heldout_recon_base": base_err, "heldout_recon_calibrated": cal_err,
            "median_scale": float(np.median(a))}
    return cal, diag


def calibrate_logscale(sub_ref, query: pd.DataFrame, cfg: HighGranularityConfig,
                       eval_genes: Optional[list] = None):
    """Affine log-scale calibration ``log S* = α + β log S`` (global α,β toward query)."""
    from tissueresolve.results import ReferenceSignature
    genes = [str(g) for g in sub_ref.gene_names]
    R = sub_ref.as_R_cpm().astype(float)
    ref_mean = R.mean(axis=0)
    q = query.reindex(index=genes).fillna(0.0)
    q_cpm = (q.div(q.sum(axis=0).where(lambda s: s > 0, np.nan), axis=1).fillna(0.0)
             * 1e6).mean(axis=1).to_numpy()
    x = np.log1p(ref_mean); y = np.log1p(q_cpm)
    mask = np.isfinite(x) & np.isfinite(y) & (ref_mean > cfg.eps)
    if mask.sum() >= 3:
        beta, alpha = np.polyfit(x[mask], y[mask], 1)
        beta = float(np.clip(beta, 0.5, 1.5))
    else:
        alpha, beta = 0.0, 1.0
    Rc = np.expm1(alpha + beta * np.log1p(R)).astype(np.float32)
    Rc = np.clip(Rc, 0.0, None)
    cal = ReferenceSignature(
        gene_names=genes, cell_types=list(sub_ref.cell_types),
        R_cpm=Rc, R_log=np.log1p(Rc).astype(np.float32), genome=sub_ref.genome,
        n_cells_per_type={m: sub_ref.n_cells_per_type.get(m, 0) for m in sub_ref.cell_types})
    diag = {"method": "logscale", "alpha": float(alpha), "beta": beta,
            "n_genes_modified": int(mask.sum())}
    return cal, diag


# ===========================================================================
# Component 2 — candidate restriction
# ===========================================================================
def select_candidates(members: list, raw_cond_fam: pd.DataFrame, cfg: HighGranularityConfig,
                      *, query_detect: Optional[Mapping[str, float]] = None,
                      donor_stable: Optional[Mapping[str, float]] = None):
    """Select plausible candidate subtypes within a family (never empties it).

    Returns (retained, dropped).  Rare-but-real subtypes (mean raw conditional ≥
    RARE_RAW) are always protected from exclusion to preserve rare sensitivity.
    """
    if cfg.candidate_rule == "none" or len(members) <= cfg.min_candidates:
        return list(members), []
    mean_raw = raw_cond_fam[members].mean(axis=0)
    order = list(mean_raw.sort_values(ascending=False).index)
    protected = [m for m in members if float(mean_raw[m]) >= RARE_RAW]

    if cfg.candidate_rule == "topk":
        keep = order[:max(cfg.top_k, cfg.min_candidates)]
    elif cfg.candidate_rule == "evidence":
        keep = []
        for m in members:
            ev_raw = float(mean_raw[m]) >= cfg.min_raw_conditional
            ev_q = (query_detect is None) or (float(query_detect.get(m, 1.0)) >= cfg.min_query_detectability)
            if ev_raw and ev_q:
                keep.append(m)
    elif cfg.candidate_rule == "hybrid":
        keep = list(order[:max(1, cfg.top_k)])
        for m in members:
            ev_raw = float(mean_raw[m]) >= cfg.min_raw_conditional
            ev_q = (query_detect is None) or (float(query_detect.get(m, 1.0)) >= cfg.min_query_detectability)
            if ev_raw and ev_q and m not in keep:
                keep.append(m)
    else:
        return list(members), []

    keep = list(dict.fromkeys(keep + protected))                # protect rare-but-real
    if not keep:                                                # never empty a family
        keep = order[:cfg.min_candidates]
    if len(keep) < cfg.min_candidates:
        keep += [m for m in order if m not in keep][:cfg.min_candidates - len(keep)]
    if len(keep) > cfg.max_candidates:
        keep = [m for m in order if m in keep][:cfg.max_candidates]
    retained = [m for m in members if m in keep]
    dropped = [m for m in members if m not in keep]
    return retained, dropped


# ===========================================================================
# Component 3 — panels
# ===========================================================================
def build_family_panels(sub_ref, query: pd.DataFrame, cfg: HighGranularityConfig,
                        *, donor_residuals: Optional[np.ndarray] = None,
                        omega: Optional[np.ndarray] = None) -> dict:
    """Return {panel_name: gene_list} for the requested robust panels."""
    genes = [str(g) for g in sub_ref.gene_names]
    G = len(genes)
    S = sub_ref.as_R_cpm().astype(float)
    n = min(cfg.panel_size, G)
    out: dict = {}
    qset = set(map(str, query.index))
    for name in cfg.panels:
        if name == "current":
            out["current"] = [g for g in genes if g in qset]
            continue
        if name == "marker":
            _, D = family_program_and_residual(S, omega)
            C = np.sqrt((D ** 2).mean(axis=0))
            idx = np.argsort(C)[::-1][:n]
        elif name == "donor_stable":
            if donor_residuals is not None and len(donor_residuals) >= 2:
                dr = np.asarray(donor_residuals, float)             # (n_d, k, G)
                per = np.sqrt((dr ** 2).mean(axis=1))               # (n_d, G)
                stab = 1.0 / (1.0 + per.var(axis=0))
                idx = np.argsort(stab)[::-1][:n]
            else:                                                   # fall back to contrast
                _, D = family_program_and_residual(S, omega)
                idx = np.argsort(np.sqrt((D ** 2).mean(axis=0)))[::-1][:n]
        elif name == "query_detectable":
            q = query.reindex(index=genes).fillna(0.0)
            det = (q > 0).mean(axis=1).to_numpy()
            idx = np.argsort(det)[::-1][:n]
        else:
            continue
        out[name] = [genes[i] for i in idx if genes[i] in qset]
    # never let a panel collapse to nothing
    out = {k: (v if len(v) >= 2 else [g for g in genes if g in qset]) for k, v in out.items()}
    return out


def multi_panel_conditional(query, sub_ref_candidates, panels: dict, members: list) -> dict:
    """Run within-family WNNLS per panel (restricted to candidate subtypes)."""
    res = {}
    for name, panel_genes in panels.items():
        ref_p = sub_ref_candidates.subset_genes(
            [g for g in panel_genes if g in set(map(str, sub_ref_candidates.gene_names))]) \
            if panel_genes else sub_ref_candidates
        try:
            cond = _wnnls_family_conditional(query.loc[query.index.intersection(panel_genes)]
                                             if panel_genes else query, ref_p, None)
        except Exception:
            continue
        res[name] = cond.reindex(columns=members).fillna(0.0)
    return res


# ===========================================================================
# Component 4 — consensus + stability
# ===========================================================================
def consensus_stability(panel_conds: dict, members: list, cfg: HighGranularityConfig,
                        *, weights: Optional[Mapping[str, float]] = None):
    """Combine panel conditionals into a stable consensus on the simplex.

    Returns (consensus_df, stability_by_subtype, diag).  ``stability_by_subtype``
    is the per-subtype detection frequency × (1 − clipped CV), in [0,1], used to
    modulate soft-gating confidence (low stability → routed to unresolved).
    """
    names = list(panel_conds)
    if not names:
        return None, pd.Series(1.0, index=members), {}
    idx = panel_conds[names[0]].index
    stack = np.stack([panel_conds[n].reindex(index=idx, columns=members).fillna(0.0).to_numpy()
                      for n in names], axis=0)        # (P, S, K)
    if cfg.consensus == "median":
        comb = np.median(stack, axis=0)
    elif cfg.consensus == "weighted" and weights:
        w = np.array([max(weights.get(n, 1.0), 0.0) for n in names], float)
        w = w / (w.sum() + cfg.eps)
        comb = (w[:, None, None] * stack).sum(axis=0)
    elif cfg.consensus == "conservative":
        det_freq = (stack > PRESENT).mean(axis=0)     # (S, K) fraction of panels detecting
        keep = det_freq >= (cfg.min_panel_detection / max(len(names), 1))
        comb = np.median(stack, axis=0) * keep        # drop subtypes not detected by ≥M panels
    else:  # mean
        comb = stack.mean(axis=0)
    comb = project_to_simplex(comb)
    consensus = pd.DataFrame(comb, index=idx, columns=members)

    mean_k = stack.mean(axis=0)                        # (S, K)
    std_k = stack.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(mean_k > cfg.eps, std_k / np.where(mean_k > cfg.eps, mean_k, 1.0), 0.0)
    det_freq = (stack > PRESENT).mean(axis=0)          # (S, K)
    stab = det_freq * np.clip(1.0 - cv / max(cfg.stability_cv_tol, cfg.eps), 0.0, 1.0)
    stability = pd.Series(np.clip(stab.mean(axis=0), 0.0, 1.0), index=members)
    diag = {"n_panels": len(names),
            "mean_cross_panel_cv": float(np.nanmean(cv)),
            "mean_detection_freq": float(det_freq.mean()),
            "consensus_entropy": float(_entropy(comb).mean()),
            "consensus_eff_n": float(_eff_n(comb).mean())}
    return consensus, stability, diag


def _entropy(P):
    P = np.clip(P, 1e-12, None); P = P / P.sum(axis=1, keepdims=True)
    return -(P * np.log(P)).sum(axis=1)


def _eff_n(P):
    P = np.clip(P, 1e-12, None); P = P / P.sum(axis=1, keepdims=True)
    return 1.0 / (P ** 2).sum(axis=1)


# ===========================================================================
# Orchestrator
# ===========================================================================
class HighGranularityStrategy:
    """ReferenceCalibration + CandidateRestriction + ConsensusStability (pre-gate)."""

    def __init__(self, config: Optional[HighGranularityConfig] = None):
        self.cfg = config or HighGranularityConfig()

    def refine(self, raw_conditional: pd.DataFrame, fine_ref, mapping: Mapping[str, str],
               query: pd.DataFrame, *,
               family_donor_counts: Optional[Mapping[str, int]] = None,
               donor_residuals: Optional[Mapping[str, np.ndarray]] = None) -> HighGranularityResult:
        cfg = self.cfg
        subtypes = [str(c) for c in raw_conditional.columns]
        group_of = {st: str(mapping.get(st, st)) for st in subtypes}
        families: dict[str, list] = {}
        for st in subtypes:
            families.setdefault(group_of[st], []).append(st)

        refined = raw_conditional.copy().astype(float)
        stability = pd.Series(1.0, index=subtypes)
        retained_d, dropped_d, refined_fams, skipped = {}, {}, [], {}
        diag_rows, warns, cal_status = [], [], {}

        if cfg.mode == "none":
            return HighGranularityResult(
                refined_conditional=refined, raw_conditional=raw_conditional.copy(),
                stability_confidence=stability,
                skipped_families={f: "mode=none" for f in families},
                metadata=self._meta(cfg, [], {f: "mode=none" for f in families}, 0.0))

        ct = [str(c) for c in fine_ref.cell_types]
        use_calib = (cfg.mode == "reference_calibrated_consensus"
                     or (cfg.mode == "auto" and cfg.calibration != "none"))

        for fam, members in families.items():
            members = [m for m in members if m in ct]
            if len(members) < cfg.min_subtypes:
                skipped[fam] = f"<{cfg.min_subtypes} subtypes"
                continue
            ndon = (family_donor_counts or {}).get(fam)
            if ndon is not None and ndon < cfg.min_family_donors:
                skipped[fam] = f"{ndon} donors (<{cfg.min_family_donors})"
                continue
            try:
                sub_ref = _subset_ref_celltypes(fine_ref, members)
                omega = np.array([max(fine_ref.n_cells_per_type.get(m, 1), 1) for m in members], float)
                # ---- Component 1: calibration ----
                if use_calib and cfg.calibration in ("scale", "logscale"):
                    cal_fn = calibrate_gene_scale if cfg.calibration == "scale" else calibrate_logscale
                    sub_ref, cdiag = cal_fn(sub_ref, query, cfg)
                    cal_status[fam] = cdiag
                # ---- Component 2: candidate restriction ----
                qdet = self._subtype_query_detect(fine_ref, members, mapping, query, omega)
                retained, dropped = select_candidates(
                    members, raw_conditional, cfg, query_detect=qdet)
                retained_d[fam], dropped_d[fam] = retained, dropped
                cand_ref = _subset_ref_celltypes(sub_ref, retained) \
                    if set(retained) != set(members) else sub_ref
                # ---- Component 3: multi-panel deconvolution on candidates ----
                dr = (donor_residuals or {}).get(fam)
                # align donor residuals to the candidate ordering
                dr_cand = self._align_residuals(dr, members, retained) if dr is not None else None
                panels = build_family_panels(cand_ref, query, cfg,
                                             donor_residuals=dr_cand,
                                             omega=np.array([omega[members.index(m)] for m in retained], float))
                panel_conds = multi_panel_conditional(query, cand_ref, panels, retained)
                if not panel_conds:
                    raise RuntimeError("no panel produced a solution")
                # ---- Component 4: consensus + stability ----
                cons, stab_fam, cdiag = consensus_stability(panel_conds, retained, cfg)
                if cons is None:
                    raise RuntimeError("empty consensus")
                # write consensus back over the full member set (dropped → 0)
                full = pd.DataFrame(0.0, index=raw_conditional.index, columns=members)
                full[retained] = cons.reindex(index=raw_conditional.index, columns=retained).fillna(0.0)
                full = pd.DataFrame(project_to_simplex(full.to_numpy(float)),
                                    index=raw_conditional.index, columns=members)
                refined[members] = full
                for m in members:
                    stability[m] = float(stab_fam.get(m, 0.0)) if m in retained else 0.0
                refined_fams.append(fam)
                diag_rows.append({"family": fam, "n_members": len(members),
                                  "n_retained": len(retained), "n_dropped": len(dropped),
                                  **cdiag})
            except Exception as exc:
                skipped[fam] = f"failed: {exc}"
                warns.append(f"family {fam}: {exc}; kept raw conditional")
                continue

        mass_err = self._mass_err(refined, families)
        return HighGranularityResult(
            refined_conditional=refined, raw_conditional=raw_conditional.copy(),
            stability_confidence=stability.clip(0.0, 1.0),
            retained=retained_d, dropped=dropped_d, refined_families=refined_fams,
            skipped_families=skipped, diagnostics=pd.DataFrame(diag_rows),
            calibration_status=cal_status, warnings=warns,
            metadata=self._meta(cfg, refined_fams, skipped, mass_err))

    # ---- helpers ----
    @staticmethod
    def _subtype_query_detect(fine_ref, members, mapping, query, omega):
        """Per-subtype query-detectable marker support (fraction of a subtype's top
        contrast genes detected in the query)."""
        sub = _subset_ref_celltypes(fine_ref, members)
        S = sub.as_R_cpm().astype(float)
        genes = [str(g) for g in sub.gene_names]
        _, D = family_program_and_residual(S, omega)
        q = query.reindex(index=genes).fillna(0.0)
        det = (q > 0).mean(axis=1).to_numpy()
        out = {}
        for j, m in enumerate(members):
            top = np.argsort(D[j])[::-1][:30]            # subtype's most over-expressed genes
            out[m] = float(det[top].mean()) if len(top) else 1.0
        return out

    @staticmethod
    def _align_residuals(dr, members, retained):
        if dr is None:
            return None
        idx = [members.index(m) for m in retained if m in members]
        return np.asarray(dr)[:, idx, :]

    @staticmethod
    def _mass_err(cond, families):
        err = 0.0
        for fam, members in families.items():
            members = [m for m in members if m in cond.columns]
            if members:
                err = max(err, float((cond[members].sum(axis=1) - 1.0).abs().max()))
        return err

    @staticmethod
    def _meta(cfg, refined_fams, skipped, mass_err):
        return {
            "feature_status": cfg.feature_status, "version": cfg.version,
            "high_granularity_mode": cfg.mode,
            "reference_calibration": cfg.calibration,
            "candidate_rule": cfg.candidate_rule,
            "panels": list(cfg.panels), "consensus": cfg.consensus,
            "refined_families": list(refined_fams),
            "skipped_families": dict(skipped),
            "within_family_mass_error": float(mass_err),
        }
