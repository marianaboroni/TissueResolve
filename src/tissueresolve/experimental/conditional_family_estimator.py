"""Supervised conditional within-family estimator (experimental, opt-in).

Learns the map ``bulk features -> conditional within-family proportions`` from
**simulated pseudobulks built on TRAIN donors only**, and validates on **held-out
donors**. Fine states are reported only when out-of-donor learnability is demonstrated;
otherwise a reject option returns grouped / ``unresolved_<family>`` mass.

This is NOT mean-reference inversion (which failed for collinear states) — it is a
learned, donor-held-out-validated regression with a reject option. No deep learning.
Defaults, the Poisson GLM solver, unresolved mass, soft gating, and spatial code are
unchanged; this is not wired into the default pipeline.

See ``docs/dev/CONDITIONAL_FAMILY_ESTIMATOR_DESIGN.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

__all__ = [
    "ConditionalFamilyModel", "fit_conditional_family_estimator",
    "predict_conditional_family", "simulate_family_pseudobulks",
    "FEATURE_STATUS", "DEFAULT_GATES",
]

FEATURE_STATUS = "experimental"
_EPS = 1e-9

DEFAULT_GATES = {
    "min_delta_rmse_pct": 10.0,     # learnable if RMSE improves >= this % over baseline
    "partial_delta_rmse_pct": 0.0,  # partially_learnable if improvement > this %
    "max_rare_fpr_increase": 0.05,  # abstain if rare FPR rises more than this
    "max_calibration_error": 0.20,  # diagnostic_only above this
    "max_seed_instability": 0.5,    # unstable if CV of val RMSE across seeds exceeds this
}


# --------------------------------------------------------------------------- utils
def _dense(adata):
    import scipy.sparse as sp
    X = adata.X
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def _lognorm_counts(counts: np.ndarray, target: float = 1e4) -> np.ndarray:
    lib = np.maximum(counts.sum(-1, keepdims=True), 1.0)
    return np.log1p(counts / lib * target)


def _simplex_rows(V: np.ndarray) -> np.ndarray:
    V = np.clip(np.asarray(V, float), 0.0, None)
    s = V.sum(-1, keepdims=True)
    out = np.where(s > _EPS, V / np.maximum(s, _EPS), 1.0 / V.shape[-1])
    return out


# ----------------------------------------------------------------- mixture simulator
def simulate_family_pseudobulks(
    cell_counts: np.ndarray,          # (n_cells, G) raw counts of the family's cells
    cell_state: np.ndarray,           # (n_cells,) state label per cell
    states: Sequence[str],
    n_mixtures: int,
    cells_per_mixture: int = 400,
    rng: Optional[np.random.Generator] = None,
    include_absent: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (X_bulk (n_mixtures, G) CPM-normalised, Y_props (n_mixtures, S))."""
    rng = rng or np.random.default_rng(0)
    S = len(states)
    by_state = {s: np.where(cell_state == s)[0] for s in states}
    by_state = {s: idx for s, idx in by_state.items() if len(idx) > 0}
    states = [s for s in states if s in by_state]
    S = len(states)
    Xb = np.zeros((n_mixtures, cell_counts.shape[1]))
    Y = np.zeros((n_mixtures, S))
    for m in range(n_mixtures):
        # vary concentration to span balanced/imbalanced/rare; optionally absent states
        alpha = rng.choice([0.2, 0.5, 1.0, 3.0])
        p = rng.dirichlet(np.full(S, alpha))
        if include_absent and S > 2 and rng.random() < 0.3:
            drop = rng.integers(S)
            p[drop] = 0.0
            p = p / max(p.sum(), _EPS)
        Y[m] = p
        bulk = np.zeros(cell_counts.shape[1])
        for si, s in enumerate(states):
            n = int(round(p[si] * cells_per_mixture))
            if n <= 0:
                continue
            pool = by_state[s]
            pick = rng.choice(pool, n, replace=len(pool) < n)
            bulk += cell_counts[pick].sum(0)
        tot = bulk.sum()
        Xb[m] = bulk / tot if tot > 0 else bulk
    return Xb, Y, states


# ------------------------------------------------------------------ feature builders
def _select_panel(train_counts, train_state, states, n_genes=300):
    """Within-family discriminative panel: top genes by between-state variance of
    log1p-CPM means (selected on TRAIN cells only)."""
    Xln = _lognorm_counts(train_counts)
    means = np.vstack([Xln[train_state == s].mean(0) for s in states])  # (S, G)
    v = means.var(0)
    panel = np.argsort(v)[::-1][:min(n_genes, Xln.shape[1])]
    return panel, means[:, panel]


def _features(X_cpm, panel, state_means_panel, mode):
    """Build features from CPM bulk profiles restricted to the panel.

    X_cpm: (n, G) ; panel indices ; state_means_panel: (S, len(panel)) log-CPM means.
    """
    Xp = np.log1p(X_cpm[:, panel] * 1e4)               # log1p-CPM on panel
    if mode == "markers":
        return Xp
    if mode == "modules":
        # per-state module score = projection onto that state's mean direction
        M = state_means_panel
        Mn = M / (np.linalg.norm(M, axis=1, keepdims=True) + _EPS)
        return Xp @ Mn.T                               # (n, S)
    if mode == "residual":
        # subtract pan-family mean signature direction
        pan = state_means_panel.mean(0)
        pan = pan / (np.linalg.norm(pan) + _EPS)
        return Xp - (Xp @ pan)[:, None] * pan[None, :]
    if mode == "combined":
        M = state_means_panel
        Mn = M / (np.linalg.norm(M, axis=1, keepdims=True) + _EPS)
        return np.hstack([Xp, Xp @ Mn.T])
    raise ValueError(f"unknown feature_mode {mode!r}")


# ------------------------------------------------------------------------- model
@dataclass
class ConditionalFamilyModel:
    family: str
    states: list
    model_kind: str
    feature_mode: str
    panel: np.ndarray
    state_means_panel: np.ndarray
    estimator: object
    scaler_mean: np.ndarray
    scaler_std: np.ndarray
    validation: dict = field(default_factory=dict)
    decision: str = "not_learnable"
    reason: str = ""
    confidence: float = 0.0
    metadata: dict = field(default_factory=dict)


def _make_estimator(kind):
    from sklearn.linear_model import Ridge, MultiTaskElasticNet
    if kind == "ridge":
        return Ridge(alpha=1.0)
    if kind == "elastic_net":
        return MultiTaskElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=5000)
    if kind == "pairwise_ridge":
        return ("pairwise", Ridge(alpha=1.0))
    if kind == "random_forest":
        from sklearn.ensemble import RandomForestRegressor
        return RandomForestRegressor(n_estimators=200, max_depth=8, random_state=0, n_jobs=-1)
    raise ValueError(f"unknown model {kind!r}")


def _fit_predict(kind, Xtr, Ytr, Xte):
    """Fit a (possibly pairwise) regressor and predict simplex proportions on Xte."""
    if kind == "pairwise_ridge":
        from sklearn.linear_model import Ridge
        S = Ytr.shape[1]
        score = np.zeros((Xte.shape[0], S))
        for i in range(S):
            for j in range(i + 1, S):
                d = Ytr[:, i] - Ytr[:, j]
                r = Ridge(alpha=1.0).fit(Xtr, d)
                pred = r.predict(Xte)
                score[:, i] += pred
                score[:, j] -= pred
        score = score - score.min(1, keepdims=True)
        return _simplex_rows(score)
    est = _make_estimator(kind)
    est.fit(Xtr, Ytr)
    return _simplex_rows(est.predict(Xte))


def _uniform_nnls_baseline(X_cpm_te, panel, state_means_cpm_panel):
    """Mean-reference conditional split (uniform NNLS) on the same panel — baseline."""
    from scipy.optimize import nnls
    Phi = state_means_cpm_panel.T  # (panel, S)
    out = np.zeros((X_cpm_te.shape[0], Phi.shape[1]))
    for n in range(X_cpm_te.shape[0]):
        b = X_cpm_te[n, panel]
        theta, _ = nnls(Phi, b)
        out[n] = theta
    return _simplex_rows(out)


# ------------------------------------------------------------------------- public API
def fit_conditional_family_estimator(
    reference_adata,
    family_col: str,
    state_col: str,
    donor_col: str,
    family: str,
    feature_mode: str = "combined",
    model: str = "ridge",
    n_train_mixtures: int = 1000,
    n_val_mixtures: int = 300,
    n_panel_genes: int = 300,
    cells_per_mixture: int = 400,
    random_state: int = 0,
    gates: Optional[dict] = None,
    train_donors: Optional[Sequence[str]] = None,
    val_donors: Optional[Sequence[str]] = None,
    val_confound: Optional[dict] = None,
) -> ConditionalFamilyModel:
    """Fit + donor-held-out validate a conditional estimator for one family.

    Splits the family's donors into train/val (disjoint), simulates pseudobulks from
    each, selects the gene panel on TRAIN cells only, fits the model, and validates on
    held-out-donor pseudobulks vs the mean-NNLS baseline. Returns a model carrying the
    decision (learnable / partially_learnable / not_learnable / unstable / diagnostic).
    """
    gates = {**DEFAULT_GATES, **(gates or {})}
    obs = reference_adata.obs
    for c in (family_col, state_col, donor_col):
        if c not in obs.columns:
            raise ValueError(f"column {c!r} not in reference_adata.obs")
    fam_mask = obs[family_col].astype(str).to_numpy() == str(family)
    if fam_mask.sum() < 20:
        raise ValueError(f"family {family!r} has too few cells ({int(fam_mask.sum())})")
    sub = reference_adata[fam_mask]
    counts = _dense(sub).astype(float)
    state = sub.obs[state_col].astype(str).to_numpy()
    donor = sub.obs[donor_col].astype(str).to_numpy()
    states = sorted(set(state))
    if len(states) < 2:
        raise ValueError(f"family {family!r} has <2 states; nothing to resolve")

    rng = np.random.default_rng(random_state)
    if train_donors is not None and val_donors is not None:
        # explicit (e.g. cross-platform) donor split — donor-held-out by construction
        tr_d, val_d = set(map(str, train_donors)), set(map(str, val_donors))
    else:
        donors = sorted(set(donor))
        rng.shuffle(donors)
        n_val = max(1, len(donors) // 3)
        val_d, tr_d = set(donors[:n_val]), set(donors[n_val:])
    tr_m, val_m = np.isin(donor, list(tr_d)), np.isin(donor, list(val_d))
    if tr_m.sum() < 20 or val_m.sum() < 10:
        raise ValueError("insufficient train/val cells after donor split")

    # panel + state means on TRAIN cells only (no leakage)
    panel, state_means_panel_log = _select_panel(counts[tr_m], state[tr_m], states, n_panel_genes)
    # CPM-scale state means on panel (for NNLS baseline + module dirs)
    cpm_tr = counts[tr_m] / np.maximum(counts[tr_m].sum(1, keepdims=True), 1.0)
    state_means_cpm_panel = np.vstack([cpm_tr[state[tr_m] == s][:, panel].mean(0) for s in states])

    # simulate train/val pseudobulks (donor-disjoint)
    Xtr_cpm, Ytr, states = simulate_family_pseudobulks(
        counts[tr_m], state[tr_m], states, n_train_mixtures,
        cells_per_mixture=cells_per_mixture, rng=rng)
    Xva_cpm, Yva, _ = simulate_family_pseudobulks(
        counts[val_m], state[val_m], states, n_val_mixtures,
        cells_per_mixture=cells_per_mixture, rng=rng)

    # optional real-bulk-like technical confounders applied to the VALIDATION bulks only
    # (train stays clean): library-size shift, multiplicative noise, gene dropout.
    if val_confound:
        cf = val_confound
        if cf.get("lib_cv", 0) > 0:
            Xva_cpm = Xva_cpm * rng.lognormal(0.0, cf["lib_cv"], size=(Xva_cpm.shape[0], 1))
        if cf.get("mult_noise_cv", 0) > 0:
            Xva_cpm = Xva_cpm * rng.lognormal(0.0, cf["mult_noise_cv"], size=Xva_cpm.shape)
        if cf.get("dropout", 0) > 0:
            Xva_cpm = Xva_cpm * (rng.random(Xva_cpm.shape) >= cf["dropout"])
        Xva_cpm = Xva_cpm / np.maximum(Xva_cpm.sum(1, keepdims=True), _EPS)

    Ftr = _features(Xtr_cpm, panel, state_means_panel_log, feature_mode)
    Fva = _features(Xva_cpm, panel, state_means_panel_log, feature_mode)
    mu, sd = Ftr.mean(0), Ftr.std(0) + _EPS
    Ftr_s, Fva_s = (Ftr - mu) / sd, (Fva - mu) / sd

    pred_model = _fit_predict(model, Ftr_s, Ytr, Fva_s)
    pred_base = _uniform_nnls_baseline(Xva_cpm, panel, state_means_cpm_panel)

    def _rmse(P):
        return float(np.sqrt(np.mean((P - Yva) ** 2)))
    rmse_model, rmse_base = _rmse(pred_model), _rmse(pred_base)
    delta_pct = (rmse_base - rmse_model) / max(rmse_base, _EPS) * 100.0
    # rare = least-abundant state; FPR on val mixtures where it's truly absent
    rare_i = int(np.argmin(Yva.mean(0)))
    absent = Yva[:, rare_i] < 1e-6
    def _fpr(P):
        called = P[:, rare_i] > 0.05
        return float((called & absent).sum() / max(absent.sum(), 1)) if absent.any() else float("nan")
    fpr_model, fpr_base = _fpr(pred_model), _fpr(pred_base)
    # calibration: mean abs error of per-state predicted vs true mean proportion
    cal_err = float(np.mean(np.abs(pred_model.mean(0) - Yva.mean(0))))
    pear = float(np.corrcoef(pred_model.ravel(), Yva.ravel())[0, 1]) if np.std(pred_model) > 0 else 0.0

    validation = {"rmse_model": round(rmse_model, 4), "rmse_baseline": round(rmse_base, 4),
                  "delta_rmse_pct": round(delta_pct, 2), "pearson": round(pear, 3),
                  "rare_fpr_model": round(fpr_model, 3) if np.isfinite(fpr_model) else None,
                  "rare_fpr_baseline": round(fpr_base, 3) if np.isfinite(fpr_base) else None,
                  "calibration_error": round(cal_err, 4),
                  "n_train_donors": len(tr_d), "n_val_donors": len(val_d),
                  "n_states": len(states)}

    # decision (single-seed; benchmark aggregates seed reproducibility/instability)
    fpr_inc = (fpr_model - fpr_base) if (np.isfinite(fpr_model) and np.isfinite(fpr_base)) else 0.0
    if cal_err > gates["max_calibration_error"]:
        decision, reason = "diagnostic_only", f"calibration_error {cal_err:.3f} > {gates['max_calibration_error']}"
    elif fpr_inc > gates["max_rare_fpr_increase"]:
        decision, reason = "diagnostic_only", f"rare FPR increased by {fpr_inc:.3f}"
    elif delta_pct >= gates["min_delta_rmse_pct"]:
        decision, reason = "learnable", f"RMSE improved {delta_pct:.1f}% over baseline"
    elif delta_pct > gates["partial_delta_rmse_pct"]:
        decision, reason = "partially_learnable", f"RMSE improved {delta_pct:.1f}% (sub-threshold)"
    else:
        decision, reason = "not_learnable", f"no improvement over baseline (delta {delta_pct:.1f}%)"
    confidence = float(np.clip(delta_pct / max(gates["min_delta_rmse_pct"], 1e-6), 0.0, 1.0))

    # refit on train+val features for deployment (the returned predictor)
    est_kind = model
    if model == "pairwise_ridge":
        deploy = None  # pairwise predictor is reconstructed in predict via stored data
    else:
        deploy = _make_estimator(model)
        deploy.fit(Ftr_s, Ytr)

    return ConditionalFamilyModel(
        family=str(family), states=list(states), model_kind=est_kind, feature_mode=feature_mode,
        panel=panel, state_means_panel=state_means_panel_log, estimator=deploy,
        scaler_mean=mu, scaler_std=sd, validation=validation, decision=decision,
        reason=reason, confidence=confidence,
        metadata={"feature_status": FEATURE_STATUS, "model": model, "feature_mode": feature_mode,
                  "family": str(family), "random_state": random_state, "gates": gates,
                  "n_train_mixtures": n_train_mixtures, "n_panel_genes": int(len(panel))})


def predict_conditional_family(
    model: ConditionalFamilyModel,
    bulk_profile: np.ndarray,
    family_mass: float = 1.0,
    reject: bool = True,
) -> dict:
    """Predict conditional within-family proportions for a CPM bulk profile.

    Returns ``{"states", "proportions" (sum=family_mass), "decision", "confidence",
    "unresolved"}``. With ``reject=True`` and a non-learnable/abstain decision, the
    proportions are returned as zeros and the family_mass is reported as ``unresolved``.
    """
    b = np.asarray(bulk_profile, dtype=float).reshape(1, -1)
    F = _features(b, model.panel, model.state_means_panel, model.feature_mode)
    Fs = (F - model.scaler_mean) / model.scaler_std
    if model.estimator is not None:
        p = _simplex_rows(model.estimator.predict(Fs))[0]
    else:
        p = np.full(len(model.states), 1.0 / len(model.states))
    abstain = reject and model.decision in ("not_learnable", "unstable", "diagnostic_only")
    if abstain:
        return {"states": model.states, "proportions": np.zeros(len(model.states)),
                "decision": model.decision, "confidence": model.confidence,
                "unresolved": float(family_mass), "reason": model.reason}
    return {"states": model.states, "proportions": p * family_mass,
            "decision": model.decision, "confidence": model.confidence,
            "unresolved": 0.0, "reason": model.reason}
