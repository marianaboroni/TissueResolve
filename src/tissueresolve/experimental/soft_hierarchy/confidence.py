"""Confidence features + calibration for partial-confidence gating (Phase 2A).

EXPERIMENTAL. Computes per-subtype confidence *features* (reference- and
run-derived) and calibrates a monotonic map feature→confidence on held-out
labelled mixtures.  The calibrated confidence feeds
:func:`apply_partial_confidence_gating`.

Calibrators compared (pick the simplest well-calibrated one, rule 16):
- ``monotonic`` : a clipped linear ramp on the evidence score (no learning) — baseline;
- ``logistic``  : 1-D logistic regression score→P(correct);
- ``isotonic``  : isotonic regression score→P(correct).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Optional

import numpy as np
import pandas as pd


# --------------------------------------------------------------------- features
def compute_reference_confidence_features(
    fine_ref, mapping: Mapping[str, str], *,
    min_discriminating_genes: int = 10,
    family_gene_panels: Optional[dict] = None,
) -> pd.DataFrame:
    """Per-subtype reference-derived features: closest-sibling separability and
    discriminating-gene support (reuses the existing confidence computation)."""
    from tissueresolve.reference.hierarchy import compute_within_family_subtype_confidence
    conf = compute_within_family_subtype_confidence(
        fine_ref, dict(mapping), min_discriminating_genes=min_discriminating_genes,
        family_gene_panels=family_gene_panels)
    rows = []
    for st, d in conf.items():
        sep = float(d.get("confidence", 0.0))      # (1-BC) to closest sibling (×0.5 if few genes)
        disc = int(d.get("min_disc_genes", 0))
        rows.append({"subtype": st, "family": d.get("family", mapping.get(st, st)),
                     "separability": sep,
                     "min_disc_genes": disc,
                     "marker_support": float(np.log1p(max(disc, 0))),
                     # combined reference evidence score (specificity × marker support, scaled)
                     "ref_evidence": sep * float(np.tanh(max(disc, 0) / 20.0))})
    return pd.DataFrame(rows).set_index("subtype")


def add_run_confidence_features(
    feat: pd.DataFrame,
    conditional_by_solver: Mapping[str, pd.DataFrame],
    query_detectability: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Augment with run-derived features.

    conditional_by_solver : {solver_name: conditional P(subtype|family) samples×subtypes}
        cross-solver agreement = 1 − mean|condA − condB| per subtype (over samples).
    query_detectability : per-subtype fraction of marker genes detected in the query.
    """
    feat = feat.copy()
    names = list(conditional_by_solver)
    if len(names) >= 2:
        a, b = conditional_by_solver[names[0]], conditional_by_solver[names[1]]
        cols = [c for c in a.columns if c in b.columns]
        diff = (a[cols].reindex(index=a.index) - b[cols].reindex(index=a.index)).abs().mean(axis=0)
        feat["cross_solver_agreement"] = feat.index.map(
            lambda s: float(1.0 - diff.get(s, np.nan)))
    if query_detectability is not None:
        feat["query_detectability"] = feat.index.map(lambda s: float(query_detectability.get(s, np.nan)))
    return feat


# --------------------------------------------------------------------- calibration
@dataclass
class ConfidenceModel:
    """Calibrated map: evidence score(s) → confidence in [0,1]."""
    kind: str                                  # 'monotonic' | 'logistic' | 'isotonic'
    score_col: str
    params: dict = field(default_factory=dict)
    _model: object = None
    feature_status: str = "experimental"

    def predict(self, features: pd.DataFrame) -> pd.Series:
        x = features[self.score_col].to_numpy(float)
        x = np.where(np.isfinite(x), x, np.nan)
        if self.kind == "monotonic":
            lo, hi = self.params["lo"], self.params["hi"]
            c = np.clip((x - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
        elif self.kind == "logistic":
            c = self._model.predict_proba(np.nan_to_num(x, nan=self.params["xmean"]).reshape(-1, 1))[:, 1]
        elif self.kind == "isotonic":
            c = self._model.predict(np.nan_to_num(x, nan=self.params["xmean"]))
            c = np.clip(c, 0.0, 1.0)
        else:
            raise ValueError(f"unknown calibrator {self.kind!r}")
        out = pd.Series(c, index=features.index, name="confidence")
        # NaN score → 0 confidence (conservative); recorded by the gater
        out[~np.isfinite(x)] = np.nan
        return out


def fit_confidence_model(features: pd.DataFrame, labels: pd.Series, *,
                         kind: str = "isotonic", score_col: str = "ref_evidence") -> ConfidenceModel:
    """Fit a calibrator from *score_col* to a binary correct-resolution *labels*.

    Fit ONLY on calibration data (caller guarantees split separation)."""
    df = features.join(labels.rename("_label"), how="inner").dropna(subset=[score_col, "_label"])
    x = df[score_col].to_numpy(float)
    y = df["_label"].to_numpy(float)
    xmean = float(np.nanmean(x)) if x.size else 0.0
    if kind == "monotonic":
        # ramp between the score quantiles of negatives and positives
        lo = float(np.quantile(x[y < 0.5], 0.5)) if (y < 0.5).any() else float(np.min(x) if x.size else 0)
        hi = float(np.quantile(x[y >= 0.5], 0.5)) if (y >= 0.5).any() else float(np.max(x) if x.size else 1)
        if hi <= lo:
            hi = lo + 1e-6
        return ConfidenceModel("monotonic", score_col, {"lo": lo, "hi": hi, "xmean": xmean})
    if kind == "logistic":
        from sklearn.linear_model import LogisticRegression
        m = LogisticRegression(max_iter=1000)
        if len(np.unique(y)) < 2:
            # degenerate: constant label → constant confidence
            return ConfidenceModel("monotonic", score_col,
                                   {"lo": xmean - 1, "hi": xmean + 1, "xmean": xmean})
        m.fit(x.reshape(-1, 1), y)
        return ConfidenceModel("logistic", score_col, {"xmean": xmean}, _model=m)
    if kind == "isotonic":
        from sklearn.isotonic import IsotonicRegression
        m = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        if len(np.unique(y)) < 2 or x.size < 3:
            return ConfidenceModel("monotonic", score_col,
                                   {"lo": xmean - 1, "hi": xmean + 1, "xmean": xmean})
        m.fit(x, y)
        return ConfidenceModel("isotonic", score_col, {"xmean": xmean}, _model=m)
    raise ValueError(f"unknown calibrator kind {kind!r}")


def calibration_metrics(confidence: np.ndarray, labels: np.ndarray, n_bins: int = 10) -> dict:
    """Brier score + expected calibration error + reliability-curve source data."""
    c = np.asarray(confidence, float); y = np.asarray(labels, float)
    mask = np.isfinite(c) & np.isfinite(y)
    c, y = c[mask], y[mask]
    if c.size == 0:
        return {"brier": float("nan"), "ece": float("nan"), "n": 0, "reliability": []}
    brier = float(np.mean((c - y) ** 2))
    bins = np.linspace(0, 1, n_bins + 1)
    rel, ece = [], 0.0
    for i in range(n_bins):
        m = (c >= bins[i]) & (c < bins[i + 1] if i < n_bins - 1 else c <= bins[i + 1])
        if m.sum() == 0:
            continue
        conf_mean, acc = float(c[m].mean()), float(y[m].mean())
        rel.append({"bin_lo": float(bins[i]), "bin_hi": float(bins[i + 1]),
                    "mean_confidence": conf_mean, "observed_accuracy": acc, "n": int(m.sum())})
        ece += (m.sum() / c.size) * abs(conf_mean - acc)
    return {"brier": brier, "ece": float(ece), "n": int(c.size), "reliability": rel}
