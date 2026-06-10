"""Unit + invariant tests for experimental partial-confidence gating (Phase 2A).

Offline, deterministic. Covers mass conservation, confidence bounds, the
hard-gate limits (c=0 / c=1), partial behaviour, calibration train/test
separation, missing-feature handling, and backward-compatible defaults.
"""
import numpy as np
import pandas as pd
import pytest

from tissueresolve.experimental.soft_hierarchy import (
    apply_partial_confidence_gating, SoftGatingConfig,
    ConfidenceModel, fit_confidence_model, calibration_metrics,
    compute_reference_confidence_features,  # noqa: F401 (import smoke)
)

FAMILY_MAP = {"A1": "A", "A2": "A", "B1": "B", "B2": "B"}


def _toy(n=4, seed=0):
    """broad_mass (samples×families) and pre-gating fine = family×conditional."""
    rng = np.random.default_rng(seed)
    fams = ["A", "B"]
    broad = rng.dirichlet([1, 1], size=n)
    broad = pd.DataFrame(broad, columns=fams, index=[f"s{i}" for i in range(n)])
    # conditional within each family
    fine = pd.DataFrame(index=broad.index, columns=["A1", "A2", "B1", "B2"], dtype=float)
    for fam, members in {"A": ["A1", "A2"], "B": ["B1", "B2"]}.items():
        q = rng.dirichlet([1, 1], size=n)
        for j, m in enumerate(members):
            fine[m] = broad[fam].to_numpy() * q[:, j]
    return broad, fine


def test_confidence_clipped_to_unit_interval():
    broad, fine = _toy()
    conf = pd.Series({"A1": 1.5, "A2": -0.3, "B1": 0.5, "B2": 0.0})  # out of range
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=conf)
    v = r.confidence_by_subtype.to_numpy()
    assert v.min() >= 0.0 and v.max() <= 1.0


def test_total_and_family_mass_conserved():
    broad, fine = _toy(seed=1)
    conf = pd.Series({"A1": 0.7, "A2": 0.2, "B1": 0.9, "B2": 0.4})
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=conf)
    assert r.mass_conservation_error < 1e-9
    # total mass preserved (== broad row sum)
    tot = r.combined.sum(axis=1)
    assert np.allclose(tot, broad.sum(axis=1), atol=1e-9)
    # family mass preserved: resolved_f + unresolved_f == broad_f
    for fam, members in {"A": ["A1", "A2"], "B": ["B1", "B2"]}.items():
        recon = r.resolved_fine_estimates[members].sum(axis=1) + r.unresolved_by_family[f"unresolved_{fam}"]
        assert np.allclose(recon, broad[fam], atol=1e-9)


def test_no_negative_mass():
    broad, fine = _toy(seed=2)
    conf = pd.Series({"A1": 0.3, "A2": 0.3, "B1": 0.3, "B2": 0.3})
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=conf)
    assert (r.resolved_fine_estimates.to_numpy() >= -1e-12).all()
    assert (r.unresolved_by_family.to_numpy() >= -1e-12).all()


def test_confidence_zero_routes_all_mass_to_unresolved():
    broad, fine = _toy(seed=3)
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=0.0)
    assert np.allclose(r.resolved_fine_estimates.to_numpy(), 0.0, atol=1e-12)
    # all family mass becomes unresolved
    for fam in ["A", "B"]:
        assert np.allclose(r.unresolved_by_family[f"unresolved_{fam}"], broad[fam], atol=1e-9)


def test_confidence_one_preserves_ungated_estimate():
    broad, fine = _toy(seed=4)
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=1.0)
    # resolved == raw fine; unresolved ~ 0
    assert np.allclose(r.resolved_fine_estimates.to_numpy(),
                       fine.clip(lower=0).to_numpy(), atol=1e-9)
    assert np.allclose(r.unresolved_by_family.to_numpy(), 0.0, atol=1e-9)


def test_partial_confidence_gives_partial_unresolved():
    broad, fine = _toy(seed=5)
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=0.5)
    u = r.unresolved_by_family.to_numpy()
    assert (u > 1e-6).any() and (u < broad[["A", "B"]].to_numpy() - 1e-9).any()
    # resolved is exactly half the raw (c=0.5), within family mass
    assert np.allclose(r.resolved_fine_estimates.to_numpy(),
                       0.5 * fine.clip(lower=0).to_numpy(), atol=1e-9)


def test_missing_confidence_is_conservative_and_warns():
    broad, fine = _toy(seed=6)
    conf = pd.Series({"A1": 0.8, "A2": np.nan, "B1": 0.8, "B2": 0.8})
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=conf)
    assert any("missing confidence" in w for w in r.warnings)
    # the NaN subtype contributes 0 resolved mass
    assert np.allclose(r.resolved_fine_estimates["A2"].to_numpy(), 0.0, atol=1e-12)
    assert r.mass_conservation_error < 1e-9


def test_hard_gate_is_special_case_of_soft():
    """c in {0,1} reproduces a hard gate; soft generalises it."""
    broad, fine = _toy(seed=7)
    hard = pd.Series({"A1": 1.0, "A2": 0.0, "B1": 1.0, "B2": 0.0})
    r = apply_partial_confidence_gating(broad, fine, FAMILY_MAP, confidence=hard)
    assert np.allclose(r.resolved_fine_estimates["A2"].to_numpy(), 0.0)
    assert np.allclose(r.resolved_fine_estimates["A1"].to_numpy(),
                       fine["A1"].to_numpy(), atol=1e-9)


def test_calibration_train_test_separation_and_determinism():
    rng = np.random.default_rng(0)
    score = rng.random(200)
    label = (score + 0.2 * rng.standard_normal(200) > 0.5).astype(float)
    feat = pd.DataFrame({"ref_evidence": score}, index=[f"x{i}" for i in range(200)])
    lab = pd.Series(label, index=feat.index)
    cal_idx, test_idx = feat.index[:120], feat.index[120:]
    for kind in ["monotonic", "logistic", "isotonic"]:
        m = fit_confidence_model(feat.loc[cal_idx], lab.loc[cal_idx], kind=kind)
        p1 = m.predict(feat.loc[test_idx]); p2 = m.predict(feat.loc[test_idx])
        assert p1.equals(p2)                       # deterministic
        assert (p1.dropna() >= 0).all() and (p1.dropna() <= 1).all()
    cm = calibration_metrics(m.predict(feat.loc[test_idx]).to_numpy(), lab.loc[test_idx].to_numpy())
    assert 0 <= cm["brier"] <= 1 and cm["ece"] >= 0


def test_default_production_behaviour_unchanged():
    """A fresh `import tissueresolve` must NOT transitively import the experimental
    soft-gating module (it is opt-in only)."""
    import subprocess, sys
    code = ("import tissueresolve, tissueresolve.api, sys; "
            "print('exp' if 'tissueresolve.experimental.soft_hierarchy' in sys.modules "
            "else 'clean'); "
            "assert hasattr(tissueresolve.api, 'deconv_bulk')")
    out = subprocess.check_output([sys.executable, "-c", code], text=True).strip()
    assert out == "clean", f"experimental module imported by default: {out}"
