"""
Spatial auto-parameter selection.

Chooses the CAR spatial-smoothing strength ``lambda_spatial`` by trading off
reconstruction quality against **over-smoothing** (loss of genuine spot-level
structure).  No ground truth is required: each candidate is run, the observed
marker expression is reconstructed from the estimated proportions, and an
over-smoothing penalty discourages collapsing real spatial variance.

Spatial accuracy is **not** claimed — without ground truth this selects a
parameter that reconstructs well without destroying spatial structure, and
records the comparison transparently.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

__all__ = [
    "graph_parameter_diagnostics", "over_smoothing_penalty",
    "reconstruction_score", "select_lambda_spatial", "save_spatial_auto_outputs",
]

DEFAULT_LAMBDAS = (0.0, 0.1, 0.5, 1.0)


def graph_parameter_diagnostics(array_row, array_col) -> dict:
    """Basic spatial-graph diagnostics (spot count, neighbour degree)."""
    info = {"n_spots": int(len(array_row))}
    try:
        from tissueresolve.spatial.graph import build_hex_graph_from_arrays
        G = build_hex_graph_from_arrays(np.asarray(array_row), np.asarray(array_col))
        deg = np.asarray(G.degree).ravel()
        info.update({"mean_degree": round(float(deg.mean()), 3),
                     "min_degree": int(deg.min()), "max_degree": int(deg.max()),
                     "isolated_spots": int((deg == 0).sum())})
    except Exception:
        info["mean_degree"] = None
    return info


def reconstruction_score(deconv, ref, Y, gene_names) -> float:
    """Per-spot-normalised Pearson between observed and reconstructed marker expr."""
    import scipy.sparse as sp
    from scipy.stats import pearsonr
    markers = [str(g) for g in deconv.marker_genes if str(g) in set(map(str, ref.gene_names))]
    if not markers:
        return float("nan")
    sub = ref.subset_genes(markers)
    cols = [c for c in deconv.proportions.columns if c in set(sub.cell_types)]
    idx = [list(sub.cell_types).index(c) for c in cols]
    recon = deconv.proportions[cols].to_numpy(float) @ sub.as_R_cpm()[idx]  # spots × markers
    gi = {str(g): i for i, g in enumerate(gene_names)}
    keep = [g for g in markers if g in gi]
    obs_idx = [gi[g] for g in keep]
    Yd = (Y.tocsr()[:, obs_idx].toarray() if sp.issparse(Y)
          else np.asarray(Y)[:, obs_idx]).astype(float)
    rec = recon[:, [markers.index(g) for g in keep]]
    Yd = Yd / (Yd.sum(axis=1, keepdims=True) + 1e-9)
    rec = rec / (rec.sum(axis=1, keepdims=True) + 1e-9)
    o, r = Yd.ravel(), rec.ravel()
    if o.std() < 1e-12 or r.std() < 1e-12:
        return float("nan")
    return float(pearsonr(o, r)[0])


def over_smoothing_penalty(proportions: pd.DataFrame,
                           baseline: pd.DataFrame) -> float:
    """Fractional loss of per-cell-type spatial variance vs the *baseline* (λ=0).

    1.0 = all spatial variance removed (maximal over-smoothing); 0.0 = none."""
    cols = [c for c in proportions.columns if c in baseline.columns]
    v = proportions[cols].var(axis=0).to_numpy()
    v0 = baseline[cols].var(axis=0).to_numpy()
    denom = float(v0.sum())
    if denom <= 0:
        return 0.0
    return float(np.clip(1.0 - v.sum() / denom, 0.0, 1.0))


def select_lambda_spatial(Y, ref, array_row, array_col, lib_sizes, gene_names,
                          spot_ids=None, *, candidates=DEFAULT_LAMBDAS,
                          config=None, max_iter: int = 15,
                          oversmooth_weight: float = 0.3):
    """Sweep ``lambda_spatial`` candidates and pick the best by
    ``reconstruction − oversmooth_weight × over-smoothing penalty``.

    Returns ``(best_lambda, comparison_df, selected_dict)``.  Each candidate is
    run with a reduced ``max_iter`` for tractability (recorded)."""
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.spatial.pipeline import SpatialPipeline

    cfg = config or TissueResolveConfig()
    base_props = None
    rows, props_by_lambda = [], {}
    for lam in candidates:
        c = TissueResolveConfig()
        c.spatial_solver = cfg.spatial_solver
        c.spatial_solver.lambda_spatial = float(lam)
        c.spatial_solver.max_iter = int(max_iter)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = SpatialPipeline(c).run(Y, ref, array_row, array_col, lib_sizes,
                                         gene_names, spot_ids, run_neighbourhood=False)
        props_by_lambda[lam] = res.deconv.proportions
        if lam == 0.0 or base_props is None:
            base_props = res.deconv.proportions
        recon = reconstruction_score(res.deconv, ref, Y, gene_names)
        rows.append({"lambda_spatial": float(lam), "reconstruction_pearson": recon,
                     "converged": bool(res.deconv.converged),
                     "n_iter": int(res.deconv.n_iter)})
    # over-smoothing penalty relative to the smallest-lambda baseline
    base_props = props_by_lambda[min(candidates)]
    for row in rows:
        pen = over_smoothing_penalty(props_by_lambda[row["lambda_spatial"]], base_props)
        row["over_smoothing_penalty"] = round(pen, 4)
        rec = row["reconstruction_pearson"]
        row["objective"] = (float("nan") if (rec is None or np.isnan(rec))
                            else round(rec - oversmooth_weight * pen, 4))
    comp = pd.DataFrame(rows).set_index("lambda_spatial")
    valid = comp["objective"].dropna()
    best_lambda = float(valid.idxmax()) if not valid.empty else float(min(candidates))
    selected = {
        "selected_lambda_spatial": best_lambda,
        "candidates": list(map(float, candidates)),
        "max_iter_used": int(max_iter),
        "oversmooth_weight": oversmooth_weight,
        "reason": (f"selected lambda_spatial={best_lambda} maximising "
                   "(reconstruction − over-smoothing penalty); spatial accuracy "
                   "is not claimed without ground truth."),
        "graph_diagnostics": graph_parameter_diagnostics(array_row, array_col),
    }
    return best_lambda, comp, selected


def save_spatial_auto_outputs(comparison: pd.DataFrame, selected: dict,
                              out_dir: Path | str) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    p1 = out / "spatial_parameter_comparison.tsv"
    comparison.to_csv(p1, sep="\t")
    p2 = out / "selected_spatial_parameters.json"
    p2.write_text(json.dumps(selected, indent=2), encoding="utf-8")
    return {"comparison": p1, "selected": p2}
