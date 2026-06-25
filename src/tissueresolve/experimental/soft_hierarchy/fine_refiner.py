"""Compact FineGranularityRefiner (experimental).

EXPERIMENTAL — not imported by default, not wired into any pipeline. Tests whether
a small fine-level refinement layer can improve *conditional within-family*
subtype prediction in families where subtypes are biologically distinguishable but
highly collinear in full-panel deconvolution (shared-lineage programs dominate the
loss while subtype contrasts are small, donor-fragile, and protocol-sensitive).

The refiner runs **before** the validated soft gate and changes only the
conditional fine proportions ``q_{k|f}`` within each family; it never touches the
broad-family mass ``π_f``.  The canonical order is:

    1. broad-family mass π_f                       (unchanged)
    2. raw conditional q^raw_{k|f}                 (current fine solve)
    3. q^refined_{k|f} = refine(q^raw, ...)         (THIS module)
    4. confidence c_k                              (original features + refiner diagnostics)
    5. θ^resolved_k = π_f q^refined_{k|f} c_k       (partial confidence-weighted soft gate)
    6. unresolved u_f = π_f − Σ_k θ^resolved_k

Mass is conserved exactly: ``Σ_k q^refined_{k|f} = 1`` within each family, so
``Σ_k π_f q^refined_{k|f} = π_f`` and the broad/total mass are untouched.  Soft
gating is applied exactly once, *after* the refiner — never before, because the
gate would remove the very mass the refiner is meant to recover.

Four optional, **independently ablatable** components (do not assume the combined
model is best):

  1. contrast-weighted WNNLS  — re-solve the within-family fine problem with genes
     weighted by donor-stable subtype contrast (down-weights shared-lineage genes).
  2. residual subtype signatures — decompose S = B + D (shared family program +
     centered subtype residual) and refine the conditional along the D contrasts.
  3. family-specific spillover calibration — learn a within-family confusion map on
     donor-held-out calibration pseudobulk and project the correction to the simplex.
  4. conservative query-adaptive gene reweighting — down-weight systematically
     mismatched genes (used only if 1–3 are insufficient).

Pure / deterministic.  No PyTorch/JAX/Pyro/NumPyro/OT/NB-CAR-VI, no L1 sparsity,
no reference adaptation, no expression reconstruction.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "fine_refiner-0.1.0"

# Technical/protocol-risk gene families (mitochondrial, ribosomal, haemoglobin,
# heat-shock/stress).  Matched case-insensitively against gene symbols.
_TECH_PATTERNS = re.compile(
    r"^(MT-|MTRNR|RPL|RPS|MRPL|MRPS|HB[ABDEGMQZ]\d?$|HSP|DNAJ|FOS|JUN|JUNB|EGR1|"
    r"MALAT1$|NEAT1$|XIST$)", re.IGNORECASE)


# ===========================================================================
# Config / result
# ===========================================================================
@dataclass
class FineRefinerConfig:
    """Settings for the FineGranularityRefiner (all components ablatable)."""
    mode: str = "none"                 # none|contrast_weighted|residual_contrast|spillover_calibrated|auto
    # support thresholds — families below these fall back to the raw conditional
    min_subtypes: int = 2
    min_family_donors: int = 3
    min_family_cells: int = 30
    # contrast-weight component toggles (Component 1)
    use_donor_stability: bool = True
    use_query_detectability: bool = True
    use_pairwise_support: bool = True
    use_shared_dominance: bool = True
    use_redundancy: bool = True
    use_tech_risk: bool = True
    weight_floor: float = 0.05         # fraction of mean weight kept as a floor (no hard zeros / no L1)
    # residual-contrast component (Component 2)
    residual_step: float = 1.0         # step size along the residual-contrast adjustment
    # spillover calibration (Component 3)
    calibration: str = "none"          # none|diagonal|full|ridge
    ridge_alpha: float = 1.0
    # query-adaptive reweighting (Component 4)
    query_adaptive: bool = False
    query_adaptive_quantile: float = 0.95
    eps: float = 1e-9
    feature_status: str = FEATURE_STATUS
    version: str = ALGORITHM_VERSION


@dataclass
class FineRefinerResult:
    """Output of :meth:`FineGranularityRefiner.refine`."""
    refined_conditional: pd.DataFrame          # samples × subtypes; within-family sums to 1
    raw_conditional: pd.DataFrame
    family_gene_weights: dict[str, pd.Series] = field(default_factory=dict)
    refined_families: list = field(default_factory=list)
    skipped_families: dict[str, str] = field(default_factory=dict)
    contrast_support: pd.DataFrame = field(default_factory=pd.DataFrame)
    calibration_status: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


# ===========================================================================
# Pure numerical helpers
# ===========================================================================
def project_to_simplex(V: np.ndarray) -> np.ndarray:
    """Euclidean projection of each row of *V* onto the probability simplex.

    Duchi et al. (2008).  Rows that are all ≤0 map to the uniform distribution.
    Always returns non-negative rows that sum to 1.
    """
    V = np.atleast_2d(np.asarray(V, float))
    n, d = V.shape
    if d == 0:
        return V
    U = np.sort(V, axis=1)[:, ::-1]
    css = np.cumsum(U, axis=1) - 1.0
    rng = np.arange(1, d + 1)
    cond = U - css / rng > 0
    rho = np.where(cond.any(axis=1), cond.cumsum(axis=1).argmax(axis=1), 0)
    theta = css[np.arange(n), rho] / (rho + 1.0)
    W = np.clip(V - theta[:, None], 0.0, None)
    s = W.sum(axis=1, keepdims=True)
    bad = (s <= 0).ravel()
    if bad.any():
        W[bad] = 1.0 / d
        s = W.sum(axis=1, keepdims=True)
    return W / s


def family_program_and_residual(
    S: np.ndarray, omega: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Decompose subtype signatures ``S = B + D`` within a family.

    Parameters
    ----------
    S : (members × genes) linear subtype signatures (e.g. CPM).
    omega : (members,) non-negative subtype weights (uniform or ∝ support).

    Returns
    -------
    (B, D)
        ``B`` (genes,) is the ω-weighted shared family program; ``D``
        (members × genes) is the centered subtype residual contrast, satisfying
        ``Σ_k ω_k D_{k} = 0`` (centering is exact by construction).  ``D`` may be
        negative and must NOT be read as a count-expression signature.
    """
    S = np.atleast_2d(np.asarray(S, float))
    k = S.shape[0]
    if omega is None:
        omega = np.ones(k)
    omega = np.asarray(omega, float)
    w = omega / (omega.sum() + 1e-12)
    B = (w[:, None] * S).sum(axis=0)        # ω-weighted mean over members
    D = S - B[None, :]
    return B, D


def contrast_weights(
    S: np.ndarray,
    omega: Optional[np.ndarray] = None,
    *,
    gene_names: Optional[list] = None,
    donor_residuals: Optional[np.ndarray] = None,
    query_detectability: Optional[np.ndarray] = None,
    config: Optional[FineRefinerConfig] = None,
) -> np.ndarray:
    """Per-gene donor-stable subtype-contrast weight ``w_{g,f}`` (Component 1).

        w_{g,f} = (C · D · Q · P) / (1 + S + V + R + T)

    All factors use simple robust definitions; the result is finite, non-negative,
    and floored at a fraction of its mean (no hard zeros, no L1).

    Parameters
    ----------
    S : (members × genes) linear subtype signatures.
    omega : (members,) subtype weights.
    donor_residuals : optional (n_donors × members × genes) per-donor centered
        residual contrasts ``D`` — used for donor stability/instability.
    query_detectability : optional (genes,) in [0,1], query detectability per gene.
    """
    cfg = config or FineRefinerConfig()
    S = np.atleast_2d(np.asarray(S, float))
    k, G = S.shape
    B, D = family_program_and_residual(S, omega)
    w = np.ones(k) if omega is None else np.asarray(omega, float)
    wn = w / (w.sum() + 1e-12)

    # C — within-family subtype contrast strength: ω-weighted variance of the
    #     residual across subtypes (between-subtype variance after removing mean).
    C = np.sqrt((wn[:, None] * D ** 2).sum(axis=0))                  # (G,)

    # D_stab — donor stability of the contrast: inverse donor variance of the
    #     per-donor contrast magnitude (1 when donor data absent).
    if cfg.use_donor_stability and donor_residuals is not None and len(donor_residuals) >= 2:
        dr = np.asarray(donor_residuals, float)                     # (n_d, k, G)
        per_donor_C = np.sqrt((wn[None, :, None] * dr ** 2).sum(axis=1))  # (n_d, G)
        donor_var = per_donor_C.var(axis=0)                          # (G,)
        D_stab = 1.0 / (1.0 + donor_var)
        V = donor_var / (donor_var + np.median(donor_var) + cfg.eps)  # instability ∈ [0,1)
    else:
        D_stab = np.ones(G)
        V = np.zeros(G)

    # Q — query detectability (1 when query info absent).
    if cfg.use_query_detectability and query_detectability is not None:
        Q = np.clip(np.asarray(query_detectability, float), 0.0, 1.0)
    else:
        Q = np.ones(G)

    # P — pairwise subtype support: fraction of subtype pairs the gene helps
    #     discriminate (|D_ag − D_bg| above the gene-wise median contrast).
    if cfg.use_pairwise_support and k >= 2:
        thr = np.median(C) + cfg.eps
        cnt = np.zeros(G)
        npairs = 0
        for a in range(k):
            for b in range(a + 1, k):
                cnt += (np.abs(D[a] - D[b]) > thr).astype(float)
                npairs += 1
        P = (cnt / max(npairs, 1)) if npairs else np.ones(G)
        P = 0.1 + 0.9 * P                                            # keep weakly-supported genes alive
    else:
        P = np.ones(G)

    # S_dom — shared-lineage dominance: |B| relative to subtype residual.
    if cfg.use_shared_dominance:
        S_dom = np.abs(B) / (np.abs(B) + C + cfg.eps)               # ∈ [0,1)
    else:
        S_dom = np.zeros(G)

    # R — redundancy with the single highest-contrast gene's subtype pattern.
    if cfg.use_redundancy and k >= 2 and G >= 2:
        anchor = int(np.argmax(C))
        ref_pat = D[:, anchor]
        rp = ref_pat - ref_pat.mean()
        denom = np.linalg.norm(rp) + cfg.eps
        Dc = D - D.mean(axis=0, keepdims=True)
        num = (Dc * rp[:, None]).sum(axis=0)
        R = np.abs(num) / (denom * (np.linalg.norm(Dc, axis=0) + cfg.eps))
        R = np.clip(np.nan_to_num(R, nan=0.0), 0.0, 1.0)
    else:
        R = np.zeros(G)

    # T — technical/protocol risk from gene symbols.
    if cfg.use_tech_risk and gene_names is not None:
        T = np.array([1.0 if _TECH_PATTERNS.match(str(g)) else 0.0 for g in gene_names])
    else:
        T = np.zeros(G)

    num = C * D_stab * Q * P
    den = 1.0 + S_dom + V + R + T
    w_g = np.nan_to_num(num / den, nan=0.0, posinf=0.0, neginf=0.0)
    w_g = np.clip(w_g, 0.0, None)
    mean_w = w_g.mean()
    if mean_w > 0:
        w_g = w_g + cfg.weight_floor * mean_w                       # floor: never fully zero a gene
    else:
        w_g = np.ones(G)
    return w_g


# ===========================================================================
# Spillover calibration (Component 3)
# ===========================================================================
@dataclass
class SpilloverCalibrator:
    """Per-family within-family confusion correction ``q_corr = project(C_f q_pred)``.

    Fitted ONLY on calibration data (caller guarantees donor/seed split
    separation).  Skips families with insufficient calibration support.
    """
    kind: str = "none"                  # none|diagonal|full|ridge
    ridge_alpha: float = 1.0
    maps_: dict = field(default_factory=dict)        # family -> (members, C matrix)
    skipped_: dict = field(default_factory=dict)
    min_rows: int = 6

    def fit(self, pred_by_family: Mapping[str, pd.DataFrame],
            true_by_family: Mapping[str, pd.DataFrame]) -> "SpilloverCalibrator":
        self.maps_ = {}
        self.skipped_ = {}
        if self.kind == "none":
            return self
        for fam, P in pred_by_family.items():
            T = true_by_family.get(fam)
            if T is None or P.shape[1] < 2:
                self.skipped_[fam] = "no truth / <2 subtypes"
                continue
            members = list(P.columns)
            T = T.reindex(columns=members).fillna(0.0)
            X = P.to_numpy(float)             # (n, k) predicted conditional
            Y = T.to_numpy(float)             # (n, k) true conditional
            if X.shape[0] < self.min_rows:
                self.skipped_[fam] = f"only {X.shape[0]} calibration rows (<{self.min_rows})"
                continue
            k = X.shape[1]
            if self.kind == "diagonal":
                # per-subtype shrinkage toward truth: scalar gain g_j minimising
                # ||g_j x_j − y_j||, then a stochastic column matrix is not formed;
                # store a diagonal gain applied then re-projected.
                g = np.ones(k)
                for j in range(k):
                    xj = X[:, j]
                    denom = float(xj @ xj) + self.ridge_alpha
                    g[j] = float(xj @ Y[:, j]) / denom if denom > 0 else 1.0
                C = np.diag(np.clip(g, 0.0, None))
            elif self.kind in ("full", "ridge"):
                # solve Y ≈ X C  (C: k×k), ridge-regularised; "full" uses a tiny ridge.
                alpha = self.ridge_alpha if self.kind == "ridge" else 1e-3
                XtX = X.T @ X + alpha * np.eye(k)
                C = np.linalg.solve(XtX, X.T @ Y)     # (k, k), maps pred→true
                C = np.clip(C, 0.0, None)
            else:
                self.skipped_[fam] = f"unknown calibration kind {self.kind!r}"
                continue
            self.maps_[fam] = (members, C)
        return self

    def transform_family(self, fam: str, q_pred: pd.DataFrame) -> pd.DataFrame:
        """Apply the calibration map for *fam* and project rows back to the simplex."""
        if fam not in self.maps_:
            return q_pred
        members, C = self.maps_[fam]
        q = q_pred.reindex(columns=members).fillna(0.0).to_numpy(float)
        corrected = q @ C                       # (n, k)
        proj = project_to_simplex(corrected)
        return pd.DataFrame(proj, index=q_pred.index, columns=members)


# ===========================================================================
# Reference helpers
# ===========================================================================
def _subset_ref_celltypes(ref, members: list):
    """Build a sub-:class:`ReferenceSignature` containing only *members*."""
    from tissueresolve.results import ReferenceSignature
    ct = [str(c) for c in ref.cell_types]
    idx = [ct.index(m) for m in members]
    R = ref.as_R_cpm()[idx].astype(np.float32)
    return ReferenceSignature(
        gene_names=list(ref.gene_names), cell_types=list(members),
        R_cpm=R, R_log=np.log1p(R).astype(np.float32),
        phi=(ref.as_phi()[:, idx]) if ref.phi is not None else None,
        phi_g=ref.phi_g, genome=ref.genome,
        n_cells_per_type={m: ref.n_cells_per_type.get(m, 0) for m in members})


def _wnnls_family_conditional(query: pd.DataFrame, sub_ref, gene_weights: Optional[pd.Series]):
    """Re-solve the within-family fine problem and return conditional proportions.

    *query* is genes × samples.  Returns samples × members with rows summing to 1.
    """
    from tissueresolve.bulk.solver import WNNLSSolver
    members = [str(c) for c in sub_ref.cell_types]
    panel = sorted(set(map(str, sub_ref.gene_names)) & set(map(str, query.index)))
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        res = WNNLSSolver().solve(query, sub_ref, panel, gene_weights=gene_weights)
    props = res.proportions.reindex(columns=members).fillna(0.0)
    totals = props.sum(axis=1)
    cond = props.div(totals.where(totals > 0, np.nan), axis=0)
    cond = cond.fillna(1.0 / len(members))            # uniform where no signal
    return cond


# ===========================================================================
# Orchestrator
# ===========================================================================
class FineGranularityRefiner:
    """Refine conditional within-family subtype proportions before soft gating."""

    def __init__(self, config: Optional[FineRefinerConfig] = None,
                 calibrator: Optional[SpilloverCalibrator] = None):
        self.cfg = config or FineRefinerConfig()
        self.calibrator = calibrator

    # ---- per-family contrast support diagnostics ----
    @staticmethod
    def _query_detectability(query: pd.DataFrame, genes: list) -> np.ndarray:
        q = query.reindex(index=genes).fillna(0.0)
        det = (q > 0).mean(axis=1).to_numpy(float)        # fraction of samples detecting the gene
        return det

    def refine(
        self,
        raw_conditional: pd.DataFrame,
        fine_ref,
        mapping: Mapping[str, str],
        query: Optional[pd.DataFrame] = None,
        *,
        family_donor_counts: Optional[Mapping[str, int]] = None,
        family_cell_counts: Optional[Mapping[str, int]] = None,
        donor_residuals: Optional[Mapping[str, np.ndarray]] = None,
    ) -> FineRefinerResult:
        """Return refined conditional proportions (within-family sums preserved).

        *raw_conditional* : samples × subtypes, within-family sums to 1 (the current
        fine solve, e.g. ``compute_conditional_subtype_proportions``).  Only the
        conditional is changed — broad-family mass is never touched here.
        """
        cfg = self.cfg
        subtypes = [str(c) for c in raw_conditional.columns]
        group_of = {st: str(mapping.get(st, st)) for st in subtypes}
        families: dict[str, list] = {}
        for st in subtypes:
            families.setdefault(group_of[st], []).append(st)

        refined = raw_conditional.copy().astype(float)
        weights_out: dict[str, pd.Series] = {}
        refined_fams, skipped = [], {}
        support_rows, warns = [], []

        if cfg.mode == "none":
            return FineRefinerResult(
                refined_conditional=refined, raw_conditional=raw_conditional.copy(),
                refined_families=[], skipped_families={f: "mode=none" for f in families},
                metadata=self._meta(cfg, [], {f: "mode=none" for f in families}, 0.0))

        ct = [str(c) for c in fine_ref.cell_types]
        R_cpm = fine_ref.as_R_cpm()                 # (K, G)
        gene_names = [str(g) for g in fine_ref.gene_names]

        for fam, members in families.items():
            members = [m for m in members if m in ct]
            # ---- support gating (fallback preserves the raw conditional) ----
            if len(members) < cfg.min_subtypes:
                skipped[fam] = f"<{cfg.min_subtypes} subtypes"
                continue
            ndon = (family_donor_counts or {}).get(fam, None)
            ncell = (family_cell_counts or {}).get(fam, None)
            if ndon is not None and ndon < cfg.min_family_donors:
                skipped[fam] = f"{ndon} donors (<{cfg.min_family_donors})"
                continue
            if ncell is not None and ncell < cfg.min_family_cells:
                skipped[fam] = f"{ncell} cells (<{cfg.min_family_cells})"
                continue

            idx = [ct.index(m) for m in members]
            S = R_cpm[idx]                          # (k, G)
            omega = np.array([max(fine_ref.n_cells_per_type.get(m, 1), 1)
                              for m in members], float)
            qd = (self._query_detectability(query, gene_names)
                  if (query is not None and cfg.use_query_detectability) else None)
            dr = (donor_residuals or {}).get(fam)

            w_g = contrast_weights(S, omega, gene_names=gene_names,
                                   donor_residuals=dr, query_detectability=qd, config=cfg)
            w_series = pd.Series(w_g, index=gene_names)
            weights_out[fam] = w_series
            B, Dres = family_program_and_residual(S, omega)
            support_rows.append({
                "family": fam, "n_subtypes": len(members),
                "n_contrast_genes": int((w_g > w_g.mean()).sum()),
                "mean_contrast_weight": float(w_g.mean()),
                "shared_lineage_fraction": float(
                    np.linalg.norm(B) / (np.linalg.norm(S, axis=1).mean() + cfg.eps)),
                "mean_query_detectability": float(np.mean(qd)) if qd is not None else float("nan"),
                "donor_stability_available": dr is not None,
            })

            q_fam = raw_conditional[members]
            try:
                if cfg.mode in ("contrast_weighted", "spillover_calibrated", "auto"):
                    if query is None:
                        raise ValueError("contrast_weighted requires the query matrix")
                    sub_ref = _subset_ref_celltypes(fine_ref, members)
                    q_ref = _wnnls_family_conditional(query, sub_ref, w_series)
                    if cfg.query_adaptive:
                        q_ref = self._query_adaptive(query, sub_ref, w_series, q_ref)
                elif cfg.mode == "residual_contrast":
                    q_ref = self._residual_refine(q_fam, S, omega, Dres, members, query, gene_names, cfg)
                else:
                    raise ValueError(f"unknown mode {cfg.mode!r}")
            except Exception as exc:                # fall back safely on any failure
                skipped[fam] = f"refine failed: {exc}"
                warns.append(f"family {fam}: refine failed ({exc}); kept raw conditional")
                continue

            # ---- optional spillover calibration on top ----
            if cfg.mode == "spillover_calibrated" and self.calibrator is not None:
                q_ref = self.calibrator.transform_family(fam, q_ref)

            q_ref = q_ref.reindex(index=raw_conditional.index, columns=members)
            # enforce simplex within family (non-negative, sums to 1)
            proj = project_to_simplex(q_ref.fillna(0.0).to_numpy(float))
            refined[members] = proj
            refined_fams.append(fam)

        support = pd.DataFrame(support_rows)
        mass_err = self._within_family_mass_error(refined, families)
        meta = self._meta(cfg, refined_fams, skipped, mass_err)
        if self.calibrator is not None:
            meta["calibration"] = self.calibrator.kind
            meta["calibration_families"] = sorted(self.calibrator.maps_.keys())
        return FineRefinerResult(
            refined_conditional=refined, raw_conditional=raw_conditional.copy(),
            family_gene_weights=weights_out, refined_families=refined_fams,
            skipped_families=skipped, contrast_support=support,
            calibration_status={"kind": getattr(self.calibrator, "kind", "none"),
                                "skipped": getattr(self.calibrator, "skipped_", {})},
            warnings=warns, metadata=meta)

    # ---- Component 2: residual-contrast refinement ----
    def _residual_refine(self, q_fam, S, omega, Dres, members, query, gene_names, cfg):
        """Adjust the conditional along donor-stable residual contrasts D.

        Uses the query's deviation from the family-mean reconstruction, projected
        onto the (contrast-weighted) residual directions, to nudge the raw
        conditional; result is projected back to the simplex.  When no query is
        available it returns the raw conditional unchanged (safe fallback).
        """
        if query is None:
            return q_fam
        w_g = contrast_weights(S, omega, gene_names=gene_names, config=cfg)
        genes = [g for g in gene_names]
        qd = query.reindex(index=genes).fillna(0.0)
        # L1-normalise each sample (CPM-like) to compare with linear signatures
        qn = qd.div(qd.sum(axis=0).where(lambda s: s > 0, np.nan), axis=1).fillna(0.0) * 1e6
        Sn = S / (S.sum(axis=1, keepdims=True) + cfg.eps) * 1e6      # (k, G) normalised
        Bn = (np.asarray(omega, float)[:, None] / (omega.sum() + cfg.eps) * Sn).sum(axis=0)
        Dn = Sn - Bn[None, :]                                        # centered residuals (k, G)
        Wg = np.sqrt(np.clip(w_g, 0, None))
        Dw = Dn * Wg[None, :]                                        # weighted residual basis
        out = q_fam.copy().astype(float)
        for samp in q_fam.index:
            y = qn[samp].to_numpy(float)
            resid = (y - Bn) * Wg
            # least-squares load of residual onto weighted contrasts → per-subtype score
            G = Dw @ Dw.T + cfg.eps * np.eye(Dw.shape[0])
            score = np.linalg.solve(G, Dw @ resid)
            score = score - score.mean()                            # centered adjustment
            base = q_fam.loc[samp, members].to_numpy(float)
            adj = base + cfg.residual_step * score * base.sum() / (np.abs(score).sum() + cfg.eps)
            out.loc[samp, members] = adj
        proj = project_to_simplex(out[members].to_numpy(float))
        return pd.DataFrame(proj, index=q_fam.index, columns=members)

    # ---- Component 4: conservative query-adaptive reweighting ----
    def _query_adaptive(self, query, sub_ref, w_series, q_ref):
        """One conservative reweight: down-weight genes with large reconstruction
        residuals, re-solve once.  Uses separate residual genes from the panel."""
        members = [str(c) for c in sub_ref.cell_types]
        panel = sorted(set(map(str, sub_ref.gene_names)) & set(map(str, query.index)))
        phi = pd.DataFrame(sub_ref.as_phi(), index=sub_ref.gene_names, columns=members).loc[panel]
        R = phi.to_numpy(float)
        Q = query.loc[panel].to_numpy(float)
        Qn = Q / (Q.sum(axis=0, keepdims=True) + self.cfg.eps)
        recon = R @ q_ref.reindex(columns=members).fillna(0.0).to_numpy(float).T   # (G, n)
        resid = np.abs(Qn - recon).mean(axis=1)                     # per-gene mean |residual|
        thr = np.quantile(resid, self.cfg.query_adaptive_quantile)
        down = pd.Series(1.0, index=panel)
        down[resid > thr] = 0.25                                    # conservative down-weight
        w2 = (w_series.reindex(panel).fillna(w_series.mean()) * down)
        return _wnnls_family_conditional(query, sub_ref, w2)

    # ---- diagnostics ----
    @staticmethod
    def _within_family_mass_error(cond: pd.DataFrame, families: Mapping[str, list]) -> float:
        err = 0.0
        for fam, members in families.items():
            members = [m for m in members if m in cond.columns]
            if not members:
                continue
            s = cond[members].sum(axis=1)
            err = max(err, float((s - 1.0).abs().max()))
        return err

    @staticmethod
    def _meta(cfg, refined_fams, skipped, mass_err):
        return {
            "feature_status": cfg.feature_status, "version": cfg.version,
            "fine_refinement": cfg.mode,
            "refined_families": list(refined_fams),
            "skipped_families": dict(skipped),
            "within_family_mass_error": float(mass_err),
            "contrast_components": {
                "donor_stability": cfg.use_donor_stability,
                "query_detectability": cfg.use_query_detectability,
                "pairwise_support": cfg.use_pairwise_support,
                "shared_dominance": cfg.use_shared_dominance,
                "redundancy": cfg.use_redundancy,
                "tech_risk": cfg.use_tech_risk,
            },
            "calibration": cfg.calibration,
            "query_adaptive": cfg.query_adaptive,
        }


def build_donor_residuals(adata, fine_ref, mapping, *, celltype_col: str,
                          donor_col: str, min_cells: int = 10) -> dict:
    """Per-family per-donor centered residual contrasts D (for donor stability).

    Returns ``{family: ndarray (n_donors × members × genes)}`` aligned to
    ``fine_ref.gene_names``/member order.  Donors lacking a member (or with too
    few cells) contribute the family-mean (zero residual) for that member, so the
    array is always well-formed.  Pure aggregation; no test-donor leakage is
    enforced by the caller passing only calibration/reference donors.
    """
    import scipy.sparse as sp
    ct = [str(c) for c in fine_ref.cell_types]
    genes = [str(g) for g in fine_ref.gene_names]
    gidx = {g: i for i, g in enumerate(genes)}
    var_names = [str(v) for v in adata.var_names]
    keep = [i for i, v in enumerate(var_names) if v in gidx]
    col_to_gene = [gidx[var_names[i]] for i in keep]
    fam_members: dict[str, list] = {}
    for c in ct:
        fam_members.setdefault(str(mapping.get(c, c)), []).append(c)

    labels = adata.obs[celltype_col].astype(str).to_numpy()
    donors = adata.obs[donor_col].astype(str).to_numpy()
    uniq_donors = sorted(set(donors))
    out: dict[str, np.ndarray] = {}
    for fam, members in fam_members.items():
        if len(members) < 2:
            continue
        omega = np.array([max(fine_ref.n_cells_per_type.get(m, 1), 1) for m in members], float)
        stacks = []
        for d in uniq_donors:
            S_d = np.zeros((len(members), len(genes)))
            ok = True
            for j, m in enumerate(members):
                rows = np.where((labels == m) & (donors == d))[0]
                if rows.size < min_cells:
                    ok = False
                    break
                X = adata.X[rows]
                X = X.toarray() if sp.issparse(X) else np.asarray(X, float)
                cpm = X / np.clip(X.sum(axis=1, keepdims=True), 1, None) * 1e6
                vec = cpm.mean(axis=0)
                S_d[j, col_to_gene] = vec[keep] if vec.shape[0] == len(var_names) else vec[:len(keep)]
            if not ok:
                continue
            _, D_d = family_program_and_residual(S_d, omega)
            stacks.append(D_d)
        if len(stacks) >= 2:
            out[fam] = np.stack(stacks, axis=0)
    return out
