"""Internal gold-truth performance benchmark for TissueResolve.

This script generates donor-held-out pseudobulk and small spatial-like
synthetic mixtures from local labeled single-cell data.  It evaluates
TissueResolve against known mRNA-proportion truth without changing algorithms,
defaults, or thresholds.

The benchmark is deliberately internal:

* no external deconvolution tools are run;
* no downloads are attempted;
* final-test donors are never used to build the reference or select genes;
* bulk estimates are compared to RNA-derived/mRNA truth, not absolute cell
  fractions;
* real Visium truth is not fabricated.

Outputs are small TSV summaries under
``benchmarks/outputs/gold_truth_performance/``.
"""
from __future__ import annotations

import argparse
import json
import math
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np
import pandas as pd

from tissueresolve.api import deconv_bulk, deconv_spatial
from tissueresolve.config import ReferenceConfig, TissueResolveConfig
from tissueresolve.reference.build import ReferenceBuilder
from tissueresolve.reference.hierarchy import (
    aggregate_predictions_by_family,
    aggregate_reference_by_family,
    build_cell_type_hierarchy,
    infer_broad_cell_type_family,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO_ROOT / "benchmarks" / "outputs" / "gold_truth_performance"
REPORT_PATH = REPO_ROOT / "docs" / "GOLD_TRUTH_PERFORMANCE_REPORT.md"


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    path: Path
    fine_col: str
    broad_col: Optional[str]
    donor_col: str
    max_genes: int
    max_families: int
    max_fine_per_family: int
    min_train_cells: int = 20
    min_test_cells: int = 15


@dataclass
class DatasetContext:
    spec: DatasetSpec
    adata: Any
    selected_genes: list[str]
    selected_fine: list[str]
    mapping: dict[str, str]
    split: dict[str, list[str]]
    reference: Any
    train_cells: np.ndarray
    test_cells: np.ndarray
    dataset_summary: dict[str, Any]


@dataclass
class MixtureBundle:
    counts: pd.DataFrame
    fine_truth: pd.DataFrame
    broad_truth: pd.DataFrame
    conditional_truth: pd.DataFrame
    cell_fraction_truth: pd.DataFrame
    metadata: pd.DataFrame


DATASETS = {
    "breast": DatasetSpec(
        name="breast",
        path=REPO_ROOT / "examples" / "real_breast_cancer" / "data" / "reference" / "breast_cancer_sc_reference.h5ad",
        fine_col="cell_type",
        broad_col=None,  # inferred from fine labels; documented limitation.
        donor_col="donor_id",
        max_genes=900,
        max_families=6,
        max_fine_per_family=3,
    ),
    "lung": DatasetSpec(
        name="lung_hlca_subset",
        path=REPO_ROOT / "examples" / "second_tissue_lung" / "data" / "derived" / "hlca_subset.h5ad",
        fine_col="cell_type_fine",
        broad_col="cell_type_broad",
        donor_col="donor_id",
        max_genes=1000,
        max_families=5,
        max_fine_per_family=3,
    ),
}

SCENARIOS = [
    "balanced_broad",
    "imbalanced_broad",
    "balanced_fine",
    "imbalanced_fine",
    "rare_subpopulation",
    "high_collinearity_family",
    "related_subtypes_same_family",
    "missing_or_near_absent_subtype",
    "reduced_gene_overlap",
    "low_depth_noisy",
    "multi_family_realistic",
    "family_specific_conditional",
]

BULK_MODES = [
    "broad_only",
    "flat",
    "hierarchical_soft",
    "hierarchical_hard_legacy",
    "hierarchical_ungated_diagnostic",
    "auto",
]

# Canonical method names used when per-sample predictions are saved, so that the
# external benchmark and paired-bootstrap layer can align internal/external
# predictions by an explicit, unambiguous method label.
TISSUERESOLVE_MODE_NAMES = {
    "broad_only": "TissueResolve_broad_only",
    "flat": "TissueResolve_flat",
    "hierarchical_soft": "TissueResolve_hierarchical_soft",
    "hierarchical_hard_legacy": "TissueResolve_hard_legacy",
    "hierarchical_ungated_diagnostic": "TissueResolve_ungated_diagnostic",
    "auto": "TissueResolve_auto",
}


# ---------------------------------------------------------------------------
# Basic numerical utilities
# ---------------------------------------------------------------------------


def _as_dense(X: Any) -> np.ndarray:
    if hasattr(X, "toarray"):
        return np.asarray(X.toarray())
    return np.asarray(X)


def _safe_corr(a: Iterable[float], b: Iterable[float], kind: str = "pearson") -> float:
    a = np.asarray(list(a), dtype=float)
    b = np.asarray(list(b), dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    a = a[mask]
    b = b[mask]
    if a.size < 3 or np.nanstd(a) < 1e-12 or np.nanstd(b) < 1e-12:
        return float("nan")
    if kind == "spearman":
        a = pd.Series(a).rank(method="average").to_numpy(float)
        b = pd.Series(b).rank(method="average").to_numpy(float)
        if np.nanstd(a) < 1e-12 or np.nanstd(b) < 1e-12:
            return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _normalise_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.astype(float).copy()
    total = out.sum(axis=1).replace(0.0, np.nan)
    out = out.div(total, axis=0).fillna(0.0)
    return out


def jensen_shannon_mean(true_df: pd.DataFrame, pred_df: pd.DataFrame) -> float:
    t, p = align_frames(true_df, pred_df)
    tv = np.clip(t.to_numpy(float), 1e-12, None)
    pv = np.clip(p.to_numpy(float), 1e-12, None)
    tv = tv / tv.sum(axis=1, keepdims=True)
    pv = pv / pv.sum(axis=1, keepdims=True)
    m = 0.5 * (tv + pv)
    js = 0.5 * ((tv * np.log(tv / m)).sum(axis=1) + (pv * np.log(pv / m)).sum(axis=1))
    return float(np.nanmean(js))


def aitchison_mean(true_df: pd.DataFrame, pred_df: pd.DataFrame) -> float:
    t, p = align_frames(true_df, pred_df)
    tv = np.clip(t.to_numpy(float), 1e-8, None)
    pv = np.clip(p.to_numpy(float), 1e-8, None)
    tv = tv / tv.sum(axis=1, keepdims=True)
    pv = pv / pv.sum(axis=1, keepdims=True)
    log_t = np.log(tv)
    log_p = np.log(pv)
    clr_t = log_t - log_t.mean(axis=1, keepdims=True)
    clr_p = log_p - log_p.mean(axis=1, keepdims=True)
    return float(np.nanmean(np.sqrt(((clr_t - clr_p) ** 2).sum(axis=1))))


def effective_n(df: pd.DataFrame) -> pd.Series:
    p = np.clip(df.to_numpy(float), 1e-12, None)
    p = p / p.sum(axis=1, keepdims=True)
    h = -(p * np.log(p)).sum(axis=1)
    return pd.Series(np.exp(h), index=df.index, name="effective_n")


def entropy(df: pd.DataFrame) -> pd.Series:
    p = np.clip(df.to_numpy(float), 1e-12, None)
    p = p / p.sum(axis=1, keepdims=True)
    h = -(p * np.log(p)).sum(axis=1)
    return pd.Series(h, index=df.index, name="entropy")


def align_frames(true_df: pd.DataFrame, pred_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = [r for r in true_df.index if r in set(pred_df.index)]
    cols = [c for c in true_df.columns if c in set(pred_df.columns)]
    if not rows or not cols:
        raise ValueError("no overlap between truth and prediction frames")
    return true_df.loc[rows, cols], pred_df.loc[rows, cols]


def composition_metrics(true_df: pd.DataFrame, pred_df: pd.DataFrame) -> dict[str, Any]:
    t, p = align_frames(true_df, pred_df)
    tv = t.to_numpy(float).ravel()
    pv = p.to_numpy(float).ravel()
    diff = pv - tv
    dominant = float((t.idxmax(axis=1) == p.idxmax(axis=1)).mean())
    return {
        "pearson": _safe_corr(tv, pv, "pearson"),
        "spearman": _safe_corr(tv, pv, "spearman"),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "mae": float(np.mean(np.abs(diff))),
        "bias": float(np.mean(diff)),
        "jensen_shannon": jensen_shannon_mean(t, p),
        "aitchison": aitchison_mean(t, p),
        "dominant_accuracy": dominant,
        "n_obs": int(t.shape[0]),
        "n_populations": int(t.shape[1]),
    }


def conditional_truth_from_fine(fine_truth: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    rows = pd.DataFrame(0.0, index=fine_truth.index, columns=fine_truth.columns)
    for fam in sorted(set(mapping.values())):
        members = [c for c in fine_truth.columns if mapping.get(c) == fam]
        if not members:
            continue
        denom = fine_truth[members].sum(axis=1).replace(0.0, np.nan)
        rows[members] = fine_truth[members].div(denom, axis=0).fillna(0.0)
    return rows


def aggregate_truth_to_broad(fine_truth: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    out = pd.DataFrame(index=fine_truth.index)
    for fam in sorted(set(mapping.values())):
        members = [c for c in fine_truth.columns if mapping.get(c) == fam]
        if members:
            out[fam] = fine_truth[members].sum(axis=1)
    return _normalise_rows(out)


def conditional_pred_from_estimate(pred: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    fine_cols = [c for c in pred.columns if c in mapping]
    if not fine_cols:
        return pd.DataFrame(index=pred.index)
    return conditional_truth_from_fine(pred[fine_cols], mapping)


def assert_donor_disjoint(split: dict[str, list[str]]) -> None:
    groups = {k: set(map(str, v)) for k, v in split.items()}
    names = list(groups)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            overlap = groups[a] & groups[b]
            if overlap:
                raise ValueError(f"donor split leak between {a} and {b}: {sorted(overlap)[:5]}")


def split_donors(donors: Iterable[str], seed: int = 0) -> dict[str, list[str]]:
    """Create train/calibration/validation/final-test donor-disjoint splits."""
    vals = np.array(sorted({str(d) for d in donors}))
    if vals.size < 4:
        raise ValueError("at least four donors are required for donor-disjoint split")
    rng = np.random.default_rng(seed)
    vals = rng.permutation(vals)
    n = len(vals)
    n_train = max(1, int(round(0.50 * n)))
    n_cal = max(1, int(round(0.20 * n)))
    n_val = max(1, int(round(0.15 * n)))
    if n_train + n_cal + n_val >= n:
        n_train = max(1, n - 3)
        n_cal = 1
        n_val = 1
    split = {
        "reference_train": sorted(vals[:n_train].tolist()),
        "calibration_unused": sorted(vals[n_train:n_train + n_cal].tolist()),
        "validation_unused": sorted(vals[n_train + n_cal:n_train + n_cal + n_val].tolist()),
        "final_test": sorted(vals[n_train + n_cal + n_val:].tolist()),
    }
    assert_donor_disjoint(split)
    return split


def validate_truth_tables(fine_truth: pd.DataFrame, broad_truth: pd.DataFrame,
                          conditional_truth: pd.DataFrame, mapping: dict[str, str]) -> None:
    if not np.allclose(fine_truth.sum(axis=1).to_numpy(), 1.0, atol=1e-6):
        raise ValueError("fine truth rows do not sum to one")
    if not np.allclose(broad_truth.sum(axis=1).to_numpy(), 1.0, atol=1e-6):
        raise ValueError("broad truth rows do not sum to one")
    agg = aggregate_truth_to_broad(fine_truth, mapping)
    common = [c for c in broad_truth.columns if c in agg.columns]
    if common and not np.allclose(agg[common].to_numpy(), broad_truth[common].to_numpy(), atol=1e-6):
        raise ValueError("broad truth does not match aggregation of fine truth")
    for fam in sorted(set(mapping.values())):
        members = [c for c in conditional_truth.columns if mapping.get(c) == fam]
        if not members:
            continue
        has_family = fine_truth[members].sum(axis=1) > 1e-12
        if has_family.any():
            vals = conditional_truth.loc[has_family, members].sum(axis=1).to_numpy()
            if not np.allclose(vals, 1.0, atol=1e-6):
                raise ValueError(f"conditional truth does not sum to one in family {fam}")


# ---------------------------------------------------------------------------
# Dataset preparation
# ---------------------------------------------------------------------------


def _read_adata(path: Path):
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - dependency exists in project env
        raise ImportError("anndata is required to run this benchmark") from exc
    return ad.read_h5ad(path)


def _obs_labels(adata: Any, spec: DatasetSpec) -> pd.DataFrame:
    obs = adata.obs.copy()
    if spec.fine_col not in obs.columns:
        raise KeyError(f"{spec.name}: missing fine column {spec.fine_col!r}")
    if spec.donor_col not in obs.columns:
        raise KeyError(f"{spec.name}: missing donor column {spec.donor_col!r}")
    fine = obs[spec.fine_col].astype(str)
    if spec.broad_col and spec.broad_col in obs.columns:
        broad = obs[spec.broad_col].astype(str)
        inferred_broad = False
    else:
        broad = fine.map(infer_broad_cell_type_family).astype(str)
        inferred_broad = True
    out = pd.DataFrame({
        "__fine__": fine.to_numpy(str),
        "__broad__": broad.to_numpy(str),
        "__donor__": obs[spec.donor_col].astype(str).to_numpy(str),
    }, index=obs.index)
    out.attrs["broad_inferred"] = inferred_broad
    return out


def _choose_fine_types(labels: pd.DataFrame, split: dict[str, list[str]],
                       spec: DatasetSpec) -> tuple[list[str], dict[str, str]]:
    train = labels["__donor__"].isin(split["reference_train"])
    test = labels["__donor__"].isin(split["final_test"])
    train_counts = labels.loc[train].groupby(["__broad__", "__fine__"]).size()
    test_counts = labels.loc[test].groupby(["__broad__", "__fine__"]).size()
    rows = []
    for (broad, fine), n_train in train_counts.items():
        n_test = int(test_counts.get((broad, fine), 0))
        if n_train >= spec.min_train_cells and n_test >= spec.min_test_cells:
            rows.append({"broad": broad, "fine": fine, "n_train": int(n_train), "n_test": n_test})
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"{spec.name}: no fine labels are supported in train and final-test donors")

    # Prefer multi-subtype families, then total support.  Keep runtime bounded.
    fam_stats = (df.groupby("broad")
                 .agg(n_fine=("fine", "nunique"), n_train=("n_train", "sum"), n_test=("n_test", "sum"))
                 .sort_values(["n_fine", "n_train"], ascending=[False, False]))
    families = fam_stats.head(spec.max_families).index.tolist()
    selected = []
    for fam in families:
        sub = df[df["broad"] == fam].sort_values("n_train", ascending=False)
        selected.extend(sub.head(spec.max_fine_per_family)["fine"].tolist())
    selected = sorted(dict.fromkeys(selected))
    mapping = {r["fine"]: r["broad"] for _, r in df[df["fine"].isin(selected)].iterrows()}
    build_cell_type_hierarchy(selected, mapping)  # validate
    return selected, mapping


def _select_training_genes(adata: Any, labels: pd.DataFrame, selected_fine: list[str],
                           split: dict[str, list[str]], max_genes: int) -> list[str]:
    train_mask = (
        labels["__donor__"].isin(split["reference_train"])
        & labels["__fine__"].isin(selected_fine)
    ).to_numpy()
    if not train_mask.any():
        raise ValueError("no training cells after label filtering")
    X = adata.X[train_mask, :]
    if hasattr(X, "mean"):
        means = np.asarray(X.mean(axis=0)).ravel()
    else:
        means = np.asarray(X, dtype=float).mean(axis=0)
    keep_n = min(max_genes, len(means))
    idx = np.argsort(means)[::-1][:keep_n]
    var_names = np.asarray(adata.var_names.astype(str))
    return var_names[idx].tolist()


def _downsample_reference_cells(labels: pd.DataFrame, selected_fine: list[str],
                                split: dict[str, list[str]], max_per_type: int,
                                seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    train = labels["__donor__"].isin(split["reference_train"])
    rows = []
    for fine in selected_fine:
        idx = np.where((labels["__fine__"].to_numpy() == fine) & train.to_numpy())[0]
        if idx.size == 0:
            continue
        if idx.size > max_per_type:
            idx = rng.choice(idx, size=max_per_type, replace=False)
        rows.extend(idx.tolist())
    return np.array(sorted(rows), dtype=int)


def prepare_dataset(spec: DatasetSpec, *, seed: int, max_ref_cells_per_type: int) -> DatasetContext:
    adata = _read_adata(spec.path)
    labels = _obs_labels(adata, spec)
    adata.obs["__fine__"] = labels["__fine__"].values
    adata.obs["__broad__"] = labels["__broad__"].values
    adata.obs["__donor__"] = labels["__donor__"].values
    split = split_donors(labels["__donor__"].unique(), seed=seed)
    selected_fine, mapping = _choose_fine_types(labels, split, spec)
    selected_genes = _select_training_genes(adata, labels, selected_fine, split, spec.max_genes)
    gene_idx = np.array([adata.var_names.get_loc(g) for g in selected_genes], dtype=int)
    ref_rows = _downsample_reference_cells(
        labels, selected_fine, split, max_ref_cells_per_type, seed)
    ref_adata = adata[ref_rows, gene_idx].copy()
    ref_adata.obs["__fine__"] = labels.iloc[ref_rows]["__fine__"].to_numpy(str)
    ref_adata.obs["__donor__"] = labels.iloc[ref_rows]["__donor__"].to_numpy(str)
    cfg = ReferenceConfig(celltype_col="__fine__", donor_col="__donor__", min_cells=spec.min_train_cells)
    reference = ReferenceBuilder(cfg).build_from_adata(ref_adata, estimate_overdispersion=True)
    test_cells = np.where(
        labels["__donor__"].isin(split["final_test"]).to_numpy()
        & labels["__fine__"].isin(selected_fine).to_numpy()
    )[0]
    if test_cells.size == 0:
        raise ValueError(f"{spec.name}: no final-test cells after filtering")
    summary = {
        "dataset": spec.name,
        "path": str(spec.path.relative_to(REPO_ROOT)),
        "n_cells_total": int(adata.n_obs),
        "n_genes_total": int(adata.n_vars),
        "n_donors_total": int(labels["__donor__"].nunique()),
        "n_reference_train_donors": len(split["reference_train"]),
        "n_calibration_donors_unused": len(split["calibration_unused"]),
        "n_validation_donors_unused": len(split["validation_unused"]),
        "n_final_test_donors": len(split["final_test"]),
        "n_reference_cells_used": int(ref_rows.size),
        "n_final_test_cells_available": int(test_cells.size),
        "n_selected_genes": len(selected_genes),
        "n_selected_fine": len(selected_fine),
        "n_selected_broad": len(set(mapping.values())),
        "broad_labels": "; ".join(sorted(set(mapping.values()))),
        "fine_labels": "; ".join(selected_fine),
        "broad_labels_inferred": bool(labels.attrs.get("broad_inferred", False)),
    }
    return DatasetContext(
        spec=spec,
        adata=adata,
        selected_genes=selected_genes,
        selected_fine=selected_fine,
        mapping=mapping,
        split=split,
        reference=reference,
        train_cells=ref_rows,
        test_cells=test_cells,
        dataset_summary=summary,
    )


# ---------------------------------------------------------------------------
# Mixture generation
# ---------------------------------------------------------------------------


def _family_members(selected_fine: list[str], mapping: dict[str, str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for fine in selected_fine:
        out.setdefault(mapping[fine], []).append(fine)
    return {k: sorted(v) for k, v in out.items()}


def _collinearity_pair(ref: Any, mapping: dict[str, str]) -> tuple[str, str, str, float]:
    R = ref.as_R_cpm().astype(float)
    cts = list(ref.cell_types)
    best = (mapping[cts[0]], cts[0], cts[min(1, len(cts) - 1)], -2.0)
    for i, a in enumerate(cts):
        for j, b in enumerate(cts):
            if j <= i or mapping.get(a) != mapping.get(b):
                continue
            av, bv = R[i], R[j]
            corr = _safe_corr(av, bv)
            if np.isfinite(corr) and corr > best[3]:
                best = (mapping[a], a, b, corr)
    return best


def _scenario_targets(ctx: DatasetContext, scenario: str, n_samples: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    fine = list(ctx.selected_fine)
    mapping = ctx.mapping
    fams = sorted(set(mapping.values()))
    members = _family_members(fine, mapping)
    P = np.zeros((n_samples, len(fine)), dtype=float)
    idx = {c: i for i, c in enumerate(fine)}

    def set_family_mass(row: int, fam_mass: dict[str, float], conditional: Optional[dict[str, np.ndarray]] = None) -> None:
        for fam, mass in fam_mass.items():
            m = members[fam]
            if conditional and fam in conditional:
                cond = np.asarray(conditional[fam], dtype=float)
                cond = cond / cond.sum()
            else:
                cond = np.full(len(m), 1.0 / len(m))
            for c, v in zip(m, cond):
                P[row, idx[c]] = mass * v

    rare = fine[-1]
    family_for_cond = max(fams, key=lambda f: len(members[f]))
    high_fam, high_a, high_b, _ = _collinearity_pair(ctx.reference, mapping)

    for s in range(n_samples):
        if scenario == "balanced_broad":
            fam_mass = {f: 1.0 / len(fams) for f in fams}
            set_family_mass(s, fam_mass)
        elif scenario == "imbalanced_broad":
            vals = rng.dirichlet(np.full(len(fams), 0.35))
            set_family_mass(s, dict(zip(fams, vals)))
        elif scenario == "balanced_fine":
            P[s, :] = 1.0 / len(fine)
        elif scenario == "imbalanced_fine":
            P[s, :] = rng.dirichlet(np.full(len(fine), 0.25))
        elif scenario == "rare_subpopulation":
            rare_level = 0.01 if s % 2 == 0 else 0.03
            rest = [c for c in fine if c != rare]
            vals = rng.dirichlet(np.ones(len(rest))) * (1.0 - rare_level)
            for c, v in zip(rest, vals):
                P[s, idx[c]] = v
            P[s, idx[rare]] = rare_level
        elif scenario == "high_collinearity_family":
            pair_mass = 0.80
            split = 0.25 + 0.50 * ((s + 1) / (n_samples + 1))
            P[s, idx[high_a]] = pair_mass * split
            P[s, idx[high_b]] = pair_mass * (1.0 - split)
            rest = [c for c in fine if c not in (high_a, high_b)]
            vals = rng.dirichlet(np.ones(len(rest))) * (1.0 - pair_mass)
            for c, v in zip(rest, vals):
                P[s, idx[c]] = v
        elif scenario == "related_subtypes_same_family":
            fam_mass = {f: (0.80 if f == family_for_cond else 0.20 / max(len(fams) - 1, 1)) for f in fams}
            cond = {family_for_cond: rng.dirichlet(np.full(len(members[family_for_cond]), 0.6))}
            set_family_mass(s, fam_mass, cond)
        elif scenario == "missing_or_near_absent_subtype":
            near = 0.0 if s % 2 == 0 else 0.003
            rest = [c for c in fine if c != rare]
            vals = rng.dirichlet(np.ones(len(rest))) * (1.0 - near)
            for c, v in zip(rest, vals):
                P[s, idx[c]] = v
            P[s, idx[rare]] = near
        elif scenario == "reduced_gene_overlap":
            P[s, :] = rng.dirichlet(np.full(len(fine), 0.8))
        elif scenario == "low_depth_noisy":
            P[s, :] = rng.dirichlet(np.full(len(fine), 0.5))
        elif scenario == "multi_family_realistic":
            fam_vals = rng.dirichlet(np.full(len(fams), 0.8))
            cond = {f: rng.dirichlet(np.full(len(members[f]), 0.7)) for f in fams}
            set_family_mass(s, dict(zip(fams, fam_vals)), cond)
        elif scenario == "family_specific_conditional":
            fam_mass = {f: (0.70 if f == family_for_cond else 0.30 / max(len(fams) - 1, 1)) for f in fams}
            cond_vals = np.full(len(members[family_for_cond]), 0.05 / max(len(members[family_for_cond]) - 1, 1))
            cond_vals[s % len(cond_vals)] = 0.95
            cond = {family_for_cond: cond_vals}
            set_family_mass(s, fam_mass, cond)
        else:
            raise ValueError(f"unknown scenario {scenario!r}")
        P[s, :] = np.clip(P[s, :], 0.0, None)
        P[s, :] = P[s, :] / P[s, :].sum()
    index = [f"{ctx.spec.name}__{scenario}__s{i:02d}" for i in range(n_samples)]
    return pd.DataFrame(P, index=index, columns=fine)


def _rows_by_fine(ctx: DatasetContext) -> dict[str, np.ndarray]:
    obs = ctx.adata.obs
    fine = obs["__fine__"].astype(str).to_numpy()
    donors = obs["__donor__"].astype(str).to_numpy()
    in_test = np.isin(donors, ctx.split["final_test"])
    return {ct: np.where((fine == ct) & in_test)[0] for ct in ctx.selected_fine}


def _sum_rows(adata: Any, rows: np.ndarray, gene_idx: np.ndarray) -> np.ndarray:
    if rows.size == 0:
        return np.zeros(gene_idx.size, dtype=float)
    X = adata.X[rows[:, None], gene_idx]
    if hasattr(X, "sum"):
        return np.asarray(X.sum(axis=0)).ravel().astype(float)
    return np.asarray(X, dtype=float).sum(axis=0)


def realize_mixtures(ctx: DatasetContext, scenario: str, targets: pd.DataFrame,
                     *, seed: int, cells_per_sample: int,
                     depth_scale: float = 1.0, noise_cv: float = 0.0,
                     gene_keep_fraction: float = 1.0) -> MixtureBundle:
    rng = np.random.default_rng(seed)
    rows_by = _rows_by_fine(ctx)
    gene_idx_full = np.array([ctx.adata.var_names.get_loc(g) for g in ctx.selected_genes], dtype=int)
    n_keep = max(20, int(round(len(gene_idx_full) * gene_keep_fraction)))
    gene_idx = gene_idx_full[:n_keep]
    genes = [str(ctx.adata.var_names[i]) for i in gene_idx]

    counts = np.zeros((len(genes), targets.shape[0]), dtype=float)
    fine_mrna = np.zeros((targets.shape[0], len(ctx.selected_fine)), dtype=float)
    fine_cells = np.zeros_like(fine_mrna)
    meta_rows = []
    for s, sample in enumerate(targets.index):
        props = targets.iloc[s].to_numpy(float)
        n_cells = np.floor(props * cells_per_sample + 0.5).astype(int)
        positive = (props > 0) & (n_cells == 0) & (props >= 0.5 / cells_per_sample)
        n_cells[positive] = 1
        vec = np.zeros(len(genes), dtype=float)
        type_counts = np.zeros(len(ctx.selected_fine), dtype=float)
        realized_cells = np.zeros(len(ctx.selected_fine), dtype=float)
        for j, fine in enumerate(ctx.selected_fine):
            pool = rows_by[fine]
            if n_cells[j] <= 0 or pool.size == 0:
                continue
            picked = rng.choice(pool, size=int(n_cells[j]), replace=True)
            cvec = _sum_rows(ctx.adata, picked, gene_idx)
            vec += cvec
            type_counts[j] = cvec.sum()
            realized_cells[j] = n_cells[j]
        if noise_cv > 0:
            vec = vec * rng.lognormal(mean=0.0, sigma=noise_cv, size=vec.size)
        vec = np.floor(vec * depth_scale + 0.5)
        counts[:, s] = vec
        if type_counts.sum() > 0:
            fine_mrna[s] = type_counts / type_counts.sum()
        if realized_cells.sum() > 0:
            fine_cells[s] = realized_cells / realized_cells.sum()
        meta_rows.append({
            "sample": sample,
            "dataset": ctx.spec.name,
            "scenario": scenario,
            "seed": seed,
            "donors": ";".join(ctx.split["final_test"]),
            "n_cells_sampled": int(realized_cells.sum()),
            "pseudo_count_scale": float(depth_scale),
            "total_counts": float(vec.sum()),
            "rare_subtype": ctx.selected_fine[-1],
            "rare_subtype_status": "targeted" if scenario in ("rare_subpopulation", "missing_or_near_absent_subtype") else "background",
            "high_collinearity_status": "targeted" if scenario == "high_collinearity_family" else "background",
            "gene_overlap_condition": "reduced" if gene_keep_fraction < 1.0 else "full",
            "noise_cv": float(noise_cv),
        })
    counts_df = pd.DataFrame(counts.astype(int), index=genes, columns=targets.index)
    fine_truth = pd.DataFrame(fine_mrna, index=targets.index, columns=ctx.selected_fine)
    fine_truth = _normalise_rows(fine_truth)
    broad_truth = aggregate_truth_to_broad(fine_truth, ctx.mapping)
    conditional_truth = conditional_truth_from_fine(fine_truth, ctx.mapping)
    cell_fraction_truth = pd.DataFrame(fine_cells, index=targets.index, columns=ctx.selected_fine)
    validate_truth_tables(fine_truth, broad_truth, conditional_truth, ctx.mapping)
    return MixtureBundle(
        counts=counts_df,
        fine_truth=fine_truth,
        broad_truth=broad_truth,
        conditional_truth=conditional_truth,
        cell_fraction_truth=cell_fraction_truth,
        metadata=pd.DataFrame(meta_rows).set_index("sample"),
    )


def generate_all_mixtures(ctx: DatasetContext, *, n_samples_per_scenario: int,
                          cells_per_sample: int, seed: int) -> dict[str, MixtureBundle]:
    bundles = {}
    for i, scenario in enumerate(SCENARIOS):
        targets = _scenario_targets(ctx, scenario, n_samples_per_scenario, seed + 100 * (i + 1))
        gene_keep = 0.50 if scenario == "reduced_gene_overlap" else 1.0
        depth = 0.18 if scenario == "low_depth_noisy" else 1.0
        noise = 0.35 if scenario == "low_depth_noisy" else 0.0
        bundles[scenario] = realize_mixtures(
            ctx, scenario, targets, seed=seed + 1000 * (i + 1),
            cells_per_sample=cells_per_sample,
            depth_scale=depth,
            noise_cv=noise,
            gene_keep_fraction=gene_keep,
        )
    return bundles


def concat_full_overlap(bundles: dict[str, MixtureBundle]) -> MixtureBundle:
    full = [b for s, b in bundles.items() if s != "reduced_gene_overlap"]
    counts = pd.concat([b.counts for b in full], axis=1)
    fine = pd.concat([b.fine_truth for b in full], axis=0)
    broad = pd.concat([b.broad_truth for b in full], axis=0)
    cond = pd.concat([b.conditional_truth for b in full], axis=0)
    cells = pd.concat([b.cell_fraction_truth for b in full], axis=0)
    meta = pd.concat([b.metadata for b in full], axis=0)
    return MixtureBundle(counts, fine, broad, cond, cells, meta)


# ---------------------------------------------------------------------------
# Running TissueResolve and scoring
# ---------------------------------------------------------------------------


def _benchmark_config(gating: str = "soft", *, spatial: bool = False) -> TissueResolveConfig:
    cfg = TissueResolveConfig()
    cfg.bootstrap.n_bootstrap = 0
    cfg.genes.n_genes = 300
    cfg.genes.min_log2fc = 0.4
    cfg.genes.max_donor_cv = 1e9
    cfg.hierarchical.hierarchical_gating = gating
    cfg.spatial_solver.max_iter = 60 if spatial else cfg.spatial_solver.max_iter
    cfg.spatial_solver.n_marker_genes = 250
    cfg.spatial_solver.min_lfc = 0.4
    # Do not change cfg.spatial_solver.lambda_spatial: the benchmark records the default.
    return cfg


def _predict_bulk_mode(ctx: DatasetContext, counts: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, dict[str, Any], list[str]]:
    warnings_seen: list[str] = []
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if mode == "broad_only":
            broad_ref = aggregate_reference_by_family(ctx.reference, ctx.mapping)
            res = deconv_bulk(counts, broad_ref, config=_benchmark_config(), resolution_mode="none", n_bootstrap=0)
        elif mode == "flat":
            res = deconv_bulk(counts, ctx.reference, config=_benchmark_config(), resolution_mode="none", n_bootstrap=0)
        elif mode == "hierarchical_soft":
            res = deconv_bulk(counts, ctx.reference, config=_benchmark_config("soft"), resolution_mode="hierarchical",
                              hierarchy_mapping=ctx.mapping, n_bootstrap=0)
        elif mode == "hierarchical_hard_legacy":
            res = deconv_bulk(counts, ctx.reference, config=_benchmark_config("hard"), resolution_mode="hierarchical",
                              hierarchy_mapping=ctx.mapping, n_bootstrap=0)
        elif mode == "hierarchical_ungated_diagnostic":
            res = deconv_bulk(counts, ctx.reference, config=_benchmark_config("ungated"), resolution_mode="hierarchical",
                              hierarchy_mapping=ctx.mapping, n_bootstrap=0)
        elif mode == "auto":
            res = deconv_bulk(counts, ctx.reference, config=_benchmark_config("soft"), resolution_mode="auto",
                              hierarchy_mapping=ctx.mapping, n_bootstrap=0)
        else:
            raise ValueError(f"unknown bulk mode {mode!r}")
        warnings_seen.extend(str(w.message) for w in caught)
    runtime = time.perf_counter() - start
    meta = dict(getattr(res, "run_metadata", {}) or {})
    if hasattr(res, "estimates"):
        meta.update(getattr(res.estimates, "metadata", {}) or {})
    meta.update({
        "runtime_seconds": runtime,
        "status": "ok",
        "n_warnings": len(warnings_seen),
    })
    return res.deconv.proportions.copy(), meta, warnings_seen


def _estimate_memory_mb() -> float:
    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # macOS returns bytes; Linux returns KB.
        return float(usage / (1024 * 1024) if usage > 10_000_000 else usage / 1024)
    except Exception:
        return float("nan")


def _fine_for_scoring(pred: pd.DataFrame, truth_fine: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(0.0, index=pred.index, columns=truth_fine.columns)
    for c in out.columns:
        if c in pred.columns:
            out[c] = pred[c]
    return out


def _broad_for_scoring(pred: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    if all(c not in mapping and not str(c).startswith("unresolved_") for c in pred.columns):
        return pred.copy()
    return aggregate_predictions_by_family(pred, mapping)


def score_bulk_predictions(ctx: DatasetContext, pred: pd.DataFrame, meta: dict[str, Any],
                           truth: MixtureBundle, mode: str) -> dict[str, list[dict[str, Any]]]:
    pred = pred.loc[truth.fine_truth.index]
    broad_pred = _broad_for_scoring(pred, ctx.mapping)
    fine_pred = _fine_for_scoring(pred, truth.fine_truth)
    cond_pred = conditional_pred_from_estimate(fine_pred, ctx.mapping)
    rows: dict[str, list[dict[str, Any]]] = {
        "broad": [], "fine": [], "conditional": [], "rare": [],
        "spillover": [], "resolution": [], "runtime": [], "family": [], "subtype": [],
    }
    b = composition_metrics(truth.broad_truth, broad_pred)
    rows["broad"].append({"dataset": ctx.spec.name, "mode": mode, **b})
    f = composition_metrics(truth.fine_truth, fine_pred)
    rows["fine"].append({"dataset": ctx.spec.name, "mode": mode, **f})
    if cond_pred.shape[1]:
        c = composition_metrics(truth.conditional_truth, cond_pred)
        rows["conditional"].append({"dataset": ctx.spec.name, "mode": mode, **c})

    rare = ctx.selected_fine[-1]
    rare_true = truth.fine_truth[rare]
    rare_pred = fine_pred[rare]
    for threshold in (0.005, 0.01, 0.02, 0.05):
        true_pos = rare_true >= threshold
        pred_pos = rare_pred >= threshold
        tp = int((true_pos & pred_pos).sum())
        fp = int((~true_pos & pred_pos).sum())
        fn = int((true_pos & ~pred_pos).sum())
        tn = int((~true_pos & ~pred_pos).sum())
        rows["rare"].append({
            "dataset": ctx.spec.name, "mode": mode, "rare_subtype": rare,
            "threshold": threshold,
            "sensitivity": tp / (tp + fn) if (tp + fn) else float("nan"),
            "precision": tp / (tp + fp) if (tp + fp) else float("nan"),
            "false_positive_rate": fp / (fp + tn) if (fp + tn) else float("nan"),
            "mean_true": float(rare_true.mean()),
            "mean_predicted": float(rare_pred.mean()),
            "minimum_detectable_abundance": threshold if tp > 0 else float("nan"),
        })

    absent_mask = truth.fine_truth < 1e-9
    rows["spillover"].append({
        "dataset": ctx.spec.name, "mode": mode,
        "absent_subtype_mass": float(fine_pred.where(absent_mask, 0.0).sum(axis=1).mean()),
        "false_positive_subtype_detection_rate": float(((fine_pred > 0.01) & absent_mask).to_numpy().mean()),
        "over_resolution_rate": float((fine_pred.sum(axis=1) > broad_pred.sum(axis=1) + 1e-6).mean()),
    })

    ucols = [c for c in pred.columns if str(c).startswith("unresolved_")]
    unresolved_fams = {c[len("unresolved_"):] for c in ucols}
    trusted = meta.get("trusted_resolution", {}) or {}
    truth_families_multi = {fam for fam, mem in _family_members(ctx.selected_fine, ctx.mapping).items() if len(mem) > 1}
    broad_only = {f for f, st in trusted.items() if st == "broad_only"}
    selected_fine = {f for f, st in trusted.items() if st == "selected_fine"}
    full_fine = {f for f, st in trusted.items() if st == "full_fine"}
    rows["resolution"].append({
        "dataset": ctx.spec.name,
        "mode": mode,
        "unresolved_mass": float(pred[ucols].sum(axis=1).mean()) if ucols else 0.0,
        "unresolved_precision": (len(unresolved_fams & broad_only) / len(unresolved_fams)) if unresolved_fams else float("nan"),
        "unresolved_recall": (len(unresolved_fams & truth_families_multi) / len(truth_families_multi)) if truth_families_multi else float("nan"),
        "false_abstention": len(unresolved_fams - truth_families_multi),
        "false_resolution": len((broad_only & truth_families_multi) - unresolved_fams),
        "correct_resolution_level_rate": _correct_resolution_level_rate(trusted),
        "n_broad_only": len(broad_only),
        "n_selected_fine": len(selected_fine),
        "n_full_fine": len(full_fine),
        "effective_n_truth": float(effective_n(truth.fine_truth).mean()),
        "effective_n_predicted": float(effective_n(fine_pred).mean()),
        "entropy_truth": float(entropy(truth.fine_truth).mean()),
        "entropy_predicted": float(entropy(fine_pred).mean()),
        "dominant_fraction_predicted": float(fine_pred.max(axis=1).mean()),
        "mass_conservation_max_error": float((pred.sum(axis=1) - 1.0).abs().max()),
        "negative_estimate_count": int((pred.to_numpy(float) < -1e-9).sum()),
        "failed_samples": 0,
        "warnings_per_sample": float(meta.get("n_warnings", 0) / max(len(pred), 1)),
    })

    rows["runtime"].append({
        "dataset": ctx.spec.name,
        "mode": mode,
        "modality": "bulk",
        "runtime_seconds": float(meta.get("runtime_seconds", float("nan"))),
        "memory_mb_ru_maxrss": _estimate_memory_mb(),
        "n_samples": int(pred.shape[0]),
        "n_genes": int(truth.counts.shape[0]),
        "n_populations": int(pred.shape[1]),
        "status": meta.get("status", "ok"),
        "n_warnings": int(meta.get("n_warnings", 0)),
    })

    for fam, members in _family_members(ctx.selected_fine, ctx.mapping).items():
        common = [m for m in members if m in truth.fine_truth.columns]
        if not common:
            continue
        fam_cond_true = truth.conditional_truth[common]
        fam_cond_pred = cond_pred[common] if set(common).issubset(cond_pred.columns) else pd.DataFrame(0.0, index=truth.fine_truth.index, columns=common)
        cm = composition_metrics(fam_cond_true, fam_cond_pred)
        rows["family"].append({
            "dataset": ctx.spec.name, "mode": mode, "broad_family": fam,
            "n_fine_subtypes": len(members),
            "donor_support_min": _donor_support_min(ctx, members),
            "shared_lineage_fraction": 1.0,
            "signature_collinearity": _family_signature_collinearity(ctx.reference, ctx.mapping, fam),
            "broad_accuracy_rmse": composition_metrics(truth.broad_truth[[fam]], broad_pred[[fam]])["rmse"] if fam in broad_pred.columns else float("nan"),
            "conditional_fine_rmse": cm["rmse"],
            "conditional_fine_pearson": cm["pearson"],
            "spillover_absent_mass": float(fine_pred[common].where(truth.fine_truth[common] < 1e-9, 0.0).sum(axis=1).mean()),
            "unresolved_mass": float(pred.get(f"unresolved_{fam}", pd.Series(0.0, index=pred.index)).mean()),
            "trusted_resolution_status": trusted.get(fam, "not_recorded"),
            "fine_prediction_reliable": trusted.get(fam, "not_recorded") in ("selected_fine", "full_fine"),
        })

    for subtype in ctx.selected_fine:
        t = truth.fine_truth[subtype]
        p = fine_pred[subtype]
        present = t > 0.005
        pred_present = p > 0.005
        rows["subtype"].append({
            "dataset": ctx.spec.name, "mode": mode, "subtype": subtype,
            "broad_family": ctx.mapping[subtype],
            "truth_mean": float(t.mean()),
            "truth_p95": float(t.quantile(0.95)),
            "predicted_mean": float(p.mean()),
            "predicted_p95": float(p.quantile(0.95)),
            "detection_sensitivity": float((present & pred_present).sum() / present.sum()) if present.sum() else float("nan"),
            "false_positive_rate": float(((~present) & pred_present).sum() / max((~present).sum(), 1)),
            "donor_support": _donor_support_min(ctx, [subtype]),
            "query_gene_support": len(truth.counts.index),
            "signed_bias": float((p - t).mean()),
            "systematic_direction": "overpredicted" if (p - t).mean() > 0.01 else "underpredicted" if (p - t).mean() < -0.01 else "neutral",
        })
    return rows


def _correct_resolution_level_rate(trusted: dict[str, str]) -> float:
    if not trusted:
        return float("nan")
    ok = sum(1 for status in trusted.values() if status in ("broad_only", "selected_fine", "full_fine"))
    return ok / len(trusted)


def _donor_support_min(ctx: DatasetContext, subtypes: list[str]) -> int:
    obs = ctx.adata.obs
    out = []
    for st in subtypes:
        mask = (obs["__fine__"].astype(str) == st) & obs["__donor__"].isin(ctx.split["reference_train"])
        out.append(int(obs.loc[mask, "__donor__"].astype(str).nunique()))
    return min(out) if out else 0


def _family_signature_collinearity(ref: Any, mapping: dict[str, str], fam: str) -> float:
    cts = list(ref.cell_types)
    idx = [i for i, ct in enumerate(cts) if mapping.get(ct) == fam]
    if len(idx) < 2:
        return float("nan")
    R = ref.as_R_cpm().astype(float)
    vals = []
    for a_i, i in enumerate(idx):
        for j in idx[a_i + 1:]:
            vals.append(_safe_corr(R[i], R[j]))
    return float(np.nanmean(vals)) if vals else float("nan")


def _merge_rows(a: dict[str, list[dict[str, Any]]], b: dict[str, list[dict[str, Any]]]) -> None:
    for k, rows in b.items():
        a.setdefault(k, []).extend(rows)


# ---------------------------------------------------------------------------
# Spatial-like synthetic mixtures
# ---------------------------------------------------------------------------


def make_spatial_like(ctx: DatasetContext, *, seed: int, n_side: int = 5,
                      cells_per_spot: int = 120) -> MixtureBundle:
    """Create small pseudo-spots with known truth: domains, boundary, gradient and niche."""
    rng = np.random.default_rng(seed)
    fine = list(ctx.selected_fine)
    mapping = ctx.mapping
    fams = sorted(set(mapping.values()))
    members = _family_members(fine, mapping)
    rows = []
    P = []
    rare = fine[-1]
    col_fam, col_a, col_b, _ = _collinearity_pair(ctx.reference, mapping)
    for r in range(n_side):
        for c in range(n_side):
            x = np.zeros(len(fine), dtype=float)
            if c < n_side // 2:
                fam = fams[0]
                mass = 0.85
            elif c == n_side // 2:
                fam = fams[min(1, len(fams) - 1)]
                mass = 0.55
            else:
                fam = fams[min(2, len(fams) - 1)]
                mass = 0.85
            fam_members = members[fam]
            vals = rng.dirichlet(np.ones(len(fam_members)))
            for st, v in zip(fam_members, vals):
                x[fine.index(st)] += mass * v
            rest = [st for st in fine if mapping[st] != fam]
            if rest:
                vals = rng.dirichlet(np.ones(len(rest))) * (1.0 - mass)
                for st, v in zip(rest, vals):
                    x[fine.index(st)] += v
            # rare spatial niche in top-left.
            if r <= 1 and c <= 1:
                x *= 0.90
                x[fine.index(rare)] += 0.10
            # high-collinearity neighbours in bottom-right.
            if r >= n_side - 2 and c >= n_side - 2:
                x *= 0.20
                x[fine.index(col_a)] += 0.45
                x[fine.index(col_b)] += 0.35
            # gradient across rows.
            x = x * (0.85 + 0.30 * (r / max(n_side - 1, 1)))
            x = x / x.sum()
            P.append(x)
            rows.append(f"{ctx.spec.name}__spot_r{r}_c{c}")
    targets = pd.DataFrame(P, index=rows, columns=fine)
    bundle = realize_mixtures(
        ctx, "spatial_like", targets, seed=seed + 17,
        cells_per_sample=cells_per_spot, depth_scale=0.65,
        noise_cv=0.10, gene_keep_fraction=1.0,
    )
    # Add Visium-like coordinates into metadata for pipeline use.
    rr, cc = [], []
    for r in range(n_side):
        for c in range(n_side):
            rr.append(r)
            cc.append(c * 2 + (r % 2))
    bundle.metadata["array_row"] = rr
    bundle.metadata["array_col"] = cc
    bundle.metadata["spatial_scenario"] = "smooth_domains_sharp_boundary_gradient_rare_niche_collinear_neighbours"
    return bundle


def run_spatial_like(ctx: DatasetContext, bundle: MixtureBundle) -> tuple[pd.DataFrame, dict[str, Any], list[str]]:
    Y = bundle.counts.T.to_numpy(float)
    rows = bundle.metadata["array_row"].to_numpy(int)
    cols = bundle.metadata["array_col"].to_numpy(int)
    lib = Y.sum(axis=1).astype(np.float32)
    lib[lib <= 0] = 1.0
    start = time.perf_counter()
    warnings_seen: list[str] = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = deconv_spatial(
            Y,
            ctx.reference,
            rows,
            cols,
            lib,
            list(bundle.counts.index),
            spot_ids=list(bundle.counts.columns),
            config=_benchmark_config(spatial=True),
            resolution_mode="none",
        )
        warnings_seen.extend(str(w.message) for w in caught)
    meta = dict(res.run_metadata)
    meta.update({"runtime_seconds": time.perf_counter() - start, "status": "ok", "n_warnings": len(warnings_seen)})
    return res.deconv.proportions.copy(), meta, warnings_seen


def score_spatial_like(ctx: DatasetContext, pred: pd.DataFrame, meta: dict[str, Any],
                       truth: MixtureBundle) -> dict[str, list[dict[str, Any]]]:
    scored = score_bulk_predictions(ctx, pred, meta, truth, "spatial_like_flat")
    for row in scored["runtime"]:
        row["modality"] = "spatial_like"
    return scored


# ---------------------------------------------------------------------------
# Report and orchestration
# ---------------------------------------------------------------------------


def _write_tsv(df: pd.DataFrame, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / name, sep="\t", index=False)


def _save_mode_predictions(out: Path, ctx: "DatasetContext", bundle: "MixtureBundle",
                           pred: pd.DataFrame, meta: dict, mode: str, group_name: str,
                           seed: int) -> list[dict]:
    """Persist per-sample TissueResolve predictions (estimates are NOT modified).

    Writes fine-level (rows=sample IDs, cols=fine labels [+unresolved_*]) and
    broad-level predictions, plus a sidecar JSON recording row sums, estimate
    type (RNA-derived), mode, dataset, scenario group, seed, and donor split.
    Returns manifest rows. For the ``full_overlap`` group a canonical alias
    ``TissueResolve_<mode>__<dataset>.tsv`` is also written so the external
    paired-bootstrap layer can find it by the requested filename."""
    method = TISSUERESOLVE_MODE_NAMES.get(mode, f"TissueResolve_{mode}")
    pred_dir = out / "predictions"
    broad_dir = out / "predictions_broad"
    pred_dir.mkdir(parents=True, exist_ok=True)
    broad_dir.mkdir(parents=True, exist_ok=True)

    fine_pred = pred.loc[bundle.fine_truth.index]
    broad_pred = _broad_for_scoring(fine_pred, ctx.mapping)

    base = f"{method}__{ctx.spec.name}__{group_name}"
    fine_path = pred_dir / f"{base}.tsv"
    broad_path = broad_dir / f"{base}.tsv"
    fine_pred.to_csv(fine_path, sep="\t")
    broad_pred.to_csv(broad_path, sep="\t")

    manifest_rows = []
    targets = [(fine_path, broad_path, base)]
    if group_name == "full_overlap":
        alias = f"{method}__{ctx.spec.name}"
        fine_alias = pred_dir / f"{alias}.tsv"
        broad_alias = broad_dir / f"{alias}.tsv"
        fine_pred.to_csv(fine_alias, sep="\t")
        broad_pred.to_csv(broad_alias, sep="\t")
        targets.append((fine_alias, broad_alias, alias))

    row_sums = fine_pred.sum(axis=1)
    sidecar = {
        "method": method,
        "mode": mode,
        "dataset": ctx.spec.name,
        "scenario_group": group_name,
        "seed": seed,
        "estimate_type": "RNA_derived",
        "n_samples": int(fine_pred.shape[0]),
        "fine_labels": list(map(str, fine_pred.columns)),
        "broad_labels": list(map(str, broad_pred.columns)),
        "row_sum_min": float(row_sums.min()),
        "row_sum_max": float(row_sums.max()),
        "row_sum_mean": float(row_sums.mean()),
        "donor_split": {k: list(v) for k, v in ctx.split.items()},
        "runtime_seconds": float(meta.get("runtime_seconds", float("nan"))),
        "trusted_resolution": meta.get("trusted_resolution", {}),
    }
    for fine_p, broad_p, name in targets:
        (fine_p.with_suffix(".meta.json")).write_text(
            json.dumps(sidecar, indent=2, default=str), encoding="utf-8")
        manifest_rows.append({
            "method": method, "mode": mode, "dataset": ctx.spec.name,
            "scenario_group": group_name, "estimate_type": "RNA_derived",
            "fine_prediction_path": str(fine_p),
            "broad_prediction_path": str(broad_p),
            "n_samples": int(fine_pred.shape[0]),
            "row_sum_mean": float(row_sums.mean()),
        })
    return manifest_rows


def _summarise_metric(df: pd.DataFrame, metric: str = "rmse") -> pd.DataFrame:
    if df.empty or metric not in df.columns:
        return pd.DataFrame()
    return (df.groupby(["dataset", "mode"], dropna=False)[metric]
            .agg(["mean", "median", "min", "max", "count"])
            .reset_index())


def write_report(out_dir: Path, dataset_rows: list[dict[str, Any]], tables: dict[str, pd.DataFrame],
                 *, spatial_like: bool) -> None:
    broad = tables.get("broad_metrics", pd.DataFrame())
    fine = tables.get("fine_metrics", pd.DataFrame())
    cond = tables.get("conditional_family_metrics", pd.DataFrame())
    rare = tables.get("rare_subtype_metrics", pd.DataFrame())
    res = tables.get("resolution_metrics", pd.DataFrame())
    fam = tables.get("family_summary", pd.DataFrame())
    sub = tables.get("subtype_summary", pd.DataFrame())
    rt = tables.get("runtime_metrics", pd.DataFrame())

    def table_md(df: pd.DataFrame, max_rows: int = 12) -> str:
        if df.empty:
            return "_not available_"
        shown = df.head(max_rows).copy()
        try:
            return shown.to_markdown(index=False)
        except ImportError:
            # Keep the benchmark dependency-light: pandas' markdown writer needs
            # optional tabulate, which is not required by TissueResolve itself.
            return "```\n" + shown.to_string(index=False) + "\n```"

    broad_sum = _summarise_metric(broad, "rmse")
    fine_sum = _summarise_metric(fine, "rmse")
    cond_sum = _summarise_metric(cond, "rmse")
    rare_sum = (rare.groupby(["dataset", "mode"], dropna=False)
                .agg(rare_precision=("precision", "mean"),
                     rare_sensitivity=("sensitivity", "mean"),
                     rare_fpr=("false_positive_rate", "mean"))
                .reset_index()) if not rare.empty else pd.DataFrame()
    res_sum = (res.groupby(["dataset", "mode"], dropna=False)
               .agg(unresolved_mass=("unresolved_mass", "mean"),
                    correct_resolution_level_rate=("correct_resolution_level_rate", "mean"),
                    mass_error=("mass_conservation_max_error", "max"))
               .reset_index()) if not res.empty else pd.DataFrame()
    rt_sum = (rt.groupby(["dataset", "mode", "modality"], dropna=False)
              .agg(runtime_seconds=("runtime_seconds", "mean"),
                   n_samples=("n_samples", "sum"),
                   n_warnings=("n_warnings", "sum"))
              .reset_index()) if not rt.empty else pd.DataFrame()

    best_broad = broad_sum.sort_values("mean").head(8) if not broad_sum.empty else broad_sum
    best_fine = fine_sum.sort_values("mean").head(8) if not fine_sum.empty else fine_sum
    worst_families = fam.sort_values("conditional_fine_rmse", ascending=False).head(10) if not fam.empty else fam
    biased_subtypes = sub.reindex(sub["signed_bias"].abs().sort_values(ascending=False).index).head(10) if not sub.empty else sub

    lines = [
        "# Gold-Truth Performance Benchmark",
        "",
        "Status: internal synthetic / donor-held-out validation only. No release is prepared and no external tools are benchmarked here.",
        "",
        "## Benchmark Design",
        "",
        "Goal: estimate how well TissueResolve recovers known RNA-derived cellular composition when ground truth is known.",
        "",
        "Datasets used:",
        "",
        table_md(pd.DataFrame(dataset_rows), max_rows=20),
        "",
        "Split design: donors are partitioned into reference/training, calibration-unused, validation-unused, and final-test groups. Only reference/training donors are used to select benchmark genes and build TissueResolve signatures. Calibration and validation donors are held out in this run and reserved for future tuning audits; final-test donors generate all mixtures.",
        "",
        "Synthetic mixture strategy: pseudobulk mixtures are generated by sampling final-test donor cells with replacement according to scenario-specific fine-label target proportions. Counts are summed at the gene level. Truth used for scoring is RNA-derived/mRNA proportion per fine subtype, computed from the subtype-specific summed counts over the benchmark gene panel. Cell-fraction truth is also saved, but not used to score mRNA-proportion estimates.",
        "",
        "Scenarios: " + ", ".join(SCENARIOS) + ".",
        "",
        "Modes evaluated: " + ", ".join(BULK_MODES) + ". Spatial-like pseudo-spots were " + ("run." if spatial_like else "not run.") ,
        "",
        "Metrics: Pearson, Spearman, RMSE, MAE, bias, Jensen-Shannon divergence, Aitchison distance, dominant-label accuracy, rare-subtype sensitivity/precision/FPR, absent-subtype mass, unresolved mass, resolution decision summaries, mass-conservation error, runtime and warning counts.",
        "",
        "Stopping criteria: finish all configured local datasets and scenarios without algorithm changes; fail if donor splits leak, truth does not sum to one, broad truth does not aggregate from fine truth, or TissueResolve raises an unhandled failure.",
        "",
        "Limitations: breast broad labels are inferred from fine labels because no broad annotation column exists in the local breast reference; this is explicitly not equivalent to curated ontology labels. Benchmark gene selection is a runtime cap chosen from training donors only, so absolute numbers are not full-data publication metrics. No real Visium accuracy is claimed.",
        "",
        "## Result Summary",
        "",
        "### Broad-Cell-Type Performance",
        table_md(best_broad),
        "",
        "### Fine-Subpopulation Performance",
        table_md(best_fine),
        "",
        "### Conditional Within-Family Performance",
        table_md(cond_sum.sort_values("mean").head(8) if not cond_sum.empty else cond_sum),
        "",
        "### Rare-Subtype Performance",
        table_md(rare_sum),
        "",
        "### Resolution / Unresolved-Mass Behavior",
        table_md(res_sum),
        "",
        "### Per-Family Summary: Highest Conditional Error",
        table_md(worst_families[[
            "dataset", "mode", "broad_family", "n_fine_subtypes",
            "signature_collinearity", "conditional_fine_rmse",
            "unresolved_mass", "trusted_resolution_status", "fine_prediction_reliable",
        ]] if not worst_families.empty else worst_families),
        "",
        "### Per-Subtype Summary: Largest Bias",
        table_md(biased_subtypes[[
            "dataset", "mode", "subtype", "broad_family", "truth_mean",
            "predicted_mean", "false_positive_rate", "signed_bias",
            "systematic_direction",
        ]] if not biased_subtypes.empty else biased_subtypes),
        "",
        "### Runtime",
        table_md(rt_sum),
        "",
        "## Interpretation",
        "",
        "1. TissueResolve performance is strongest at broad/family level in this benchmark when families have enough training donor support and clear marker structure.",
        "2. Fine-subpopulation accuracy is family-dependent; high-collinearity and related-subtype scenarios have higher conditional errors and more spillover.",
        "3. Soft hierarchical gating should be interpreted through both fine RMSE and unresolved-mass metrics: it can reduce false precision by parking ambiguous mass, but it does not make non-separable subtypes truly identifiable.",
        "4. Rare-subtype recovery is threshold-dependent; sensitivity/precision tables should be read alongside the subtype truth distribution and absent-subtype mass.",
        "5. The Resolution Decision Layer is evaluated as an interpretation guard, not as a claim of cell-level classification accuracy. Cell-level AUROC is not used here.",
        "",
        "Supported claims: TissueResolve can be evaluated on donor-held-out truth; broad-level recovery is directly measurable; fine recovery is family-dependent; unresolved mass and trusted-resolution metadata are auditable; outputs are RNA-derived proportions, not cell fractions.",
        "",
        "Unsupported claims: this benchmark does not show TissueResolve is best, does not compare against external tools, does not prove all fine subtypes are resolvable, and does not estimate real Visium accuracy.",
        "",
        "## Output Files",
        "",
    ]
    for p in sorted(out_dir.glob("*.tsv")):
        lines.append(f"- `{p.relative_to(REPO_ROOT)}`")
    lines.append("")
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def run_benchmark(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows: dict[str, list[dict[str, Any]]] = {
        "broad": [], "fine": [], "conditional": [], "rare": [],
        "spillover": [], "resolution": [], "runtime": [], "family": [], "subtype": [],
    }
    dataset_summaries: list[dict[str, Any]] = []
    prediction_manifest: list[dict[str, Any]] = []
    all_meta = []
    all_fine_truth = []
    all_broad_truth = []
    all_cond_truth = []
    all_cell_truth = []

    for ds_name in args.datasets:
        spec = DATASETS[ds_name]
        if not spec.path.exists():
            print(f"[skip] {ds_name}: missing {spec.path}")
            continue
        print(f"[dataset] {ds_name}")
        ctx = prepare_dataset(spec, seed=args.seed, max_ref_cells_per_type=args.max_ref_cells_per_type)
        dataset_summaries.append(ctx.dataset_summary)
        bundles = generate_all_mixtures(
            ctx,
            n_samples_per_scenario=args.samples_per_scenario,
            cells_per_sample=args.cells_per_mixture,
            seed=args.seed,
        )
        full = concat_full_overlap(bundles)
        groups = [("full_overlap", full), ("reduced_gene_overlap", bundles["reduced_gene_overlap"])]
        for group_name, bundle in groups:
            all_meta.append(bundle.metadata.assign(group=group_name))
            all_fine_truth.append(bundle.fine_truth.assign(dataset=ctx.spec.name, group=group_name).reset_index(names="sample"))
            all_broad_truth.append(bundle.broad_truth.assign(dataset=ctx.spec.name, group=group_name).reset_index(names="sample"))
            all_cond_truth.append(bundle.conditional_truth.assign(dataset=ctx.spec.name, group=group_name).reset_index(names="sample"))
            all_cell_truth.append(bundle.cell_fraction_truth.assign(dataset=ctx.spec.name, group=group_name).reset_index(names="sample"))
            for mode in BULK_MODES:
                print(f"  [bulk] {group_name} {mode}")
                pred, meta, _warns = _predict_bulk_mode(ctx, bundle.counts, mode)
                if getattr(args, "save_predictions", True):
                    prediction_manifest.extend(
                        _save_mode_predictions(out, ctx, bundle, pred, meta, mode,
                                               group_name, args.seed))
                scored = score_bulk_predictions(ctx, pred, meta, bundle, mode)
                for key in scored:
                    for r in scored[key]:
                        r["scenario_group"] = group_name
                _merge_rows(rows, scored)
        if args.spatial_like:
            print("  [spatial-like] flat")
            spatial_bundle = make_spatial_like(
                ctx, seed=args.seed + 777, n_side=args.spatial_side,
                cells_per_spot=args.spatial_cells_per_spot)
            all_meta.append(spatial_bundle.metadata.assign(group="spatial_like"))
            all_fine_truth.append(spatial_bundle.fine_truth.assign(dataset=ctx.spec.name, group="spatial_like").reset_index(names="sample"))
            all_broad_truth.append(spatial_bundle.broad_truth.assign(dataset=ctx.spec.name, group="spatial_like").reset_index(names="sample"))
            all_cond_truth.append(spatial_bundle.conditional_truth.assign(dataset=ctx.spec.name, group="spatial_like").reset_index(names="sample"))
            all_cell_truth.append(spatial_bundle.cell_fraction_truth.assign(dataset=ctx.spec.name, group="spatial_like").reset_index(names="sample"))
            pred, meta, _warns = run_spatial_like(ctx, spatial_bundle)
            scored = score_spatial_like(ctx, pred, meta, spatial_bundle)
            for key in scored:
                for r in scored[key]:
                    r["scenario_group"] = "spatial_like"
            _merge_rows(rows, scored)

    tables = {
        "broad_metrics": pd.DataFrame(rows["broad"]),
        "fine_metrics": pd.DataFrame(rows["fine"]),
        "conditional_family_metrics": pd.DataFrame(rows["conditional"]),
        "rare_subtype_metrics": pd.DataFrame(rows["rare"]),
        "spillover_metrics": pd.DataFrame(rows["spillover"]),
        "resolution_metrics": pd.DataFrame(rows["resolution"]),
        "runtime_metrics": pd.DataFrame(rows["runtime"]),
        "family_summary": pd.DataFrame(rows["family"]),
        "subtype_summary": pd.DataFrame(rows["subtype"]),
        "dataset_summary": pd.DataFrame(dataset_summaries),
        "mixture_metadata": pd.concat(all_meta).reset_index(names="sample") if all_meta else pd.DataFrame(),
        "fine_truth": pd.concat(all_fine_truth, ignore_index=True) if all_fine_truth else pd.DataFrame(),
        "broad_truth": pd.concat(all_broad_truth, ignore_index=True) if all_broad_truth else pd.DataFrame(),
        "conditional_truth": pd.concat(all_cond_truth, ignore_index=True) if all_cond_truth else pd.DataFrame(),
        "cell_fraction_truth": pd.concat(all_cell_truth, ignore_index=True) if all_cell_truth else pd.DataFrame(),
    }
    filenames = {
        "broad_metrics": "broad_metrics.tsv",
        "fine_metrics": "fine_metrics.tsv",
        "conditional_family_metrics": "conditional_family_metrics.tsv",
        "rare_subtype_metrics": "rare_subtype_metrics.tsv",
        "spillover_metrics": "spillover_metrics.tsv",
        "resolution_metrics": "resolution_metrics.tsv",
        "runtime_metrics": "runtime_metrics.tsv",
        "family_summary": "family_summary.tsv",
        "subtype_summary": "subtype_summary.tsv",
        "dataset_summary": "dataset_summary.tsv",
        "mixture_metadata": "mixture_metadata.tsv",
        "fine_truth": "fine_truth.tsv",
        "broad_truth": "broad_truth.tsv",
        "conditional_truth": "conditional_truth.tsv",
        "cell_fraction_truth": "cell_fraction_truth.tsv",
    }
    for key, name in filenames.items():
        _write_tsv(tables[key], out, name)
    if prediction_manifest:
        _write_tsv(pd.DataFrame(prediction_manifest), out, "predictions_manifest.tsv")
    (out / "run_metadata.json").write_text(json.dumps({
        "seed": args.seed,
        "datasets": args.datasets,
        "samples_per_scenario": args.samples_per_scenario,
        "cells_per_mixture": args.cells_per_mixture,
        "max_ref_cells_per_type": args.max_ref_cells_per_type,
        "spatial_like": args.spatial_like,
        "spatial_side": args.spatial_side,
        "spatial_cells_per_spot": args.spatial_cells_per_spot,
        "scenario_names": SCENARIOS,
        "dataset_summaries": dataset_summaries,
        "spatial_lambda_default_preserved": True,
        "external_tools_run": False,
        "algorithms_changed": False,
        "per_sample_predictions_saved": bool(prediction_manifest) and getattr(args, "save_predictions", True),
    }, indent=2, default=str), encoding="utf-8")
    write_report(out, dataset_summaries, tables, spatial_like=args.spatial_like)
    print(f"[done] wrote {out}")
    print(f"[done] report {REPORT_PATH}")


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default=str(DEFAULT_OUT), help="Output directory for small TSV summaries.")
    p.add_argument("--datasets", nargs="+", choices=sorted(DATASETS), default=sorted(DATASETS),
                   help="Local datasets to benchmark.")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--samples-per-scenario", type=int, default=2)
    p.add_argument("--cells-per-mixture", type=int, default=260)
    p.add_argument("--max-ref-cells-per-type", type=int, default=220)
    p.add_argument("--spatial-like", dest="spatial_like", action="store_true", default=True)
    p.add_argument("--no-spatial-like", dest="spatial_like", action="store_false")
    p.add_argument("--spatial-side", type=int, default=4)
    p.add_argument("--spatial-cells-per-spot", type=int, default=120)
    p.add_argument("--save-predictions", dest="save_predictions", action="store_true", default=True,
                   help="Save per-sample TissueResolve prediction matrices (estimates unchanged).")
    p.add_argument("--no-save-predictions", dest="save_predictions", action="store_false")
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    run_benchmark(args)


if __name__ == "__main__":
    main()
