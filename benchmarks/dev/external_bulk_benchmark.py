#!/usr/bin/env python
"""Bulk benchmark (Category A — donor-held-out pseudobulk with known truth).

Runs a FRESH internal head-to-head on breast + HLCA/lung pseudobulk:
  NNLS_baseline (plain scipy NNLS on the CPM signature),
  WNNLS_baseline (weighted-NNLS solver, flat),
  TissueResolve_flat (pipeline flat fine),
  TissueResolve_hierarchical_soft (broad->fine + soft gating + resolution decision).

External R/Python tools (MuSiC, BisqueRNA, BayesPrism, cell2location, RCTD, CARD) are
INSTALLED (see docs/EXTERNAL_TOOL_INSTALLATION_NOTES.md + tool_registry.tsv) and were
executed in a prior harness run (benchmarks/outputs/real_external_method_status.tsv).
They are NOT re-executed here (runtime/stability on this machine); their status and
prior runtimes are folded into the run manifest, clearly labelled, never hidden.

Donor- and seed-disjoint; query/test donors never used for gene selection. No model
defaults changed. Outputs (git-ignored): external_bulk_metrics.tsv,
external_bulk_per_family_metrics.tsv, external_bulk_runtime.tsv,
external_tools/run_manifest.tsv.

Usage:  python external_bulk_benchmark.py --run-real-data --tissue both
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
from benchmarks.shared import synthetic_holdout as SH  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import stats as ST  # noqa: E402

OUT = REPO / "benchmarks" / "outputs"
EXT = REPO / "benchmarks" / "external_tools"
SEEDS = list(range(12))        # expanded for statistical power (12 donor-disjoint splits/tissue)
N_SAMPLES, CELLS, MIN_CELLS = 12, 500, 20
PRESENT, ABSENT = 0.01, 0.005
TR_HIER = "TissueResolve_hierarchical_soft"
LOWER = {"broad_rmse", "broad_mae", "broad_jsd", "broad_aitchison",
         "fine_rmse", "cond_rmse", "pairwise_spillover", "false_positive_rate"}
SUPERIORITY_METRICS = ["cond_rmse", "pairwise_spillover", "false_positive_rate",
                       "rare_precision", "fine_pearson", "broad_pearson"]

_BREAST_SUPPLEMENT = {
    "IgG plasma cell": "B/Plasma", "naive B cell": "B/Plasma",
    "unswitched memory B cell": "B/Plasma",
    "activated CD8-positive, alpha-beta T cell": "T/NK", "gamma-delta T cell": "T/NK",
    "mast cell": "Myeloid", "myeloid dendritic cell": "Myeloid",
    "neutrophil": "Myeloid", "plasmacytoid dendritic cell": "Myeloid",
}
TISSUES = {
    "breast": {"h5ad": REPO / "examples/real_breast_cancer/data/reference/breast_cancer_sc_reference.h5ad",
               "hmap": REPO / "examples/real_breast_cancer/config/breast_cancer_cell_type_hierarchy.tsv",
               "fine_col": "cell_type", "donor_col": "donor_id", "supplement": _BREAST_SUPPLEMENT},
    "lung": {"h5ad": REPO / "examples/second_tissue_lung/data/derived/hlca_subset.h5ad",
             "hmap": REPO / "examples/second_tissue_lung/data/derived/hlca_hierarchy.tsv",
             "fine_col": "cell_type_fine", "donor_col": "donor_id"},
}


def _load(cfg):
    import anndata as ad
    adata = ad.read_h5ad(cfg["h5ad"])
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).to_numpy()
        adata.var_names_make_unique()
    hmap = pd.read_csv(cfg["hmap"], sep="\t")
    mapping = dict(zip(hmap["fine_cell_type"].astype(str), hmap["broad_cell_type"].astype(str)))
    mapping.update(cfg.get("supplement", {}))
    return adata, mapping


def _conditional(df, mp, cols):
    fam_of = {c: mp.get(c, c) for c in cols}
    out = df[cols] * 0.0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        s = df[mem].sum(axis=1)
        for c in mem:
            out[c] = np.where(s > 0, df[c] / s.replace(0, np.nan), 0.0)
    return out.fillna(0.0)


def _spillover(tcond, pcond, mp, cols):
    fam_of = {c: mp.get(c, c) for c in cols}
    leak, n = 0.0, 0
    for fam in set(fam_of.values()):
        mem = [c for c in cols if fam_of[c] == fam]
        if len(mem) < 2:
            continue
        absent = ~(tcond[mem] > PRESENT)
        leaked = (pcond[mem].to_numpy() * absent.to_numpy()).sum(axis=1)
        active = (tcond[mem] > PRESENT).any(axis=1).to_numpy()
        if active.any():
            leak += float(leaked[active].sum()); n += int(active.sum())
    return leak / n if n else 0.0


def _jsd(t, p):
    t = t.to_numpy(float); p = p.to_numpy(float)
    t = t / np.clip(t.sum(1, keepdims=True), 1e-12, None)
    p = p / np.clip(p.sum(1, keepdims=True), 1e-12, None)
    m = 0.5 * (t + p)

    def _kl(a, b):
        a = np.clip(a, 1e-12, None); b = np.clip(b, 1e-12, None)
        return (a * np.log(a / b)).sum(1)
    return float(np.mean(0.5 * _kl(t, m) + 0.5 * _kl(p, m)))


def _aitchison(t, p):
    def clr(x):
        x = np.clip(x.to_numpy(float), 1e-6, None)
        x = x / x.sum(1, keepdims=True)
        g = np.exp(np.log(x).mean(1, keepdims=True))
        return np.log(x / g)
    return float(np.mean(np.sqrt(((clr(t) - clr(p)) ** 2).sum(1))))


def _metrics(truth, pred, mp, ref_types):
    cols = [c for c in truth.columns if c in ref_types]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    acc = M.accuracy_metrics(t, p)
    tfam = M.aggregate_to_families(t, mp)
    pfam = M.aggregate_to_families(pred.reindex(index=t.index), mp)
    bacc = M.accuracy_metrics(tfam, pfam)
    tcond, pcond = _conditional(t, mp, cols), _conditional(p, mp, cols)
    cond = M.accuracy_metrics(tcond, pcond)
    tv, pv = t.to_numpy(float), p.to_numpy(float)
    absent = tv <= ABSENT
    rare = (tv >= 0.005) & (tv <= 0.05)
    cx, cxt = M.composition_complexity(p), M.composition_complexity(t)
    dom_acc = float((tfam.values.argmax(1) == pfam.reindex(columns=tfam.columns).values.argmax(1)).mean())
    return {
        "broad_pearson": bacc["pearson"], "broad_spearman": bacc["spearman"],
        "broad_rmse": bacc["rmse"], "broad_mae": float(np.abs(tfam.values - pfam.reindex(
            index=tfam.index, columns=tfam.columns).values).mean()),
        "broad_jsd": _jsd(tfam, pfam.reindex(index=tfam.index, columns=tfam.columns).fillna(0)),
        "broad_aitchison": _aitchison(tfam, pfam.reindex(index=tfam.index, columns=tfam.columns).fillna(0)),
        "dominant_family_accuracy": dom_acc,
        "fine_rmse": acc["rmse"], "fine_pearson": acc["pearson"], "fine_spearman": acc["spearman"],
        "cond_rmse": cond["rmse"], "cond_pearson": cond["pearson"],
        "pairwise_spillover": _spillover(tcond, pcond, mp, cols),
        "false_positive_rate": float((pv[absent] > PRESENT).mean()) if absent.any() else 0.0,
        "rare_sensitivity": float((pv[rare] > ABSENT).mean()) if rare.any() else float("nan"),
        "rare_precision": float(((pv > ABSENT) & rare).sum() / max((pv[rare] > ABSENT).sum() +
                          (pv[absent] > PRESENT).sum(), 1)) if rare.any() else float("nan"),
        "eff_n_pred": float(cx["effective_n_populations"].mean()),
        "eff_n_truth": float(cxt["effective_n_populations"].mean()),
        "entropy_pred": float(cx["entropy"].mean()) if "entropy" in cx.columns else float("nan"),
    }


def _perfamily(truth, pred, mp, ref_types):
    cols = [c for c in truth.columns if c in ref_types]
    t = truth[cols]
    p = pred.reindex(index=t.index, columns=cols).fillna(0.0)
    tcond, pcond = _conditional(t, mp, cols), _conditional(p, mp, cols)
    fam_of = {c: mp.get(c, c) for c in cols}
    rows = []
    for fam in sorted(set(fam_of.values())):
        mem = [c for c in cols if fam_of[c] == fam]
        if len(mem) < 2:
            continue
        present = (t[mem].sum(1) > PRESENT).to_numpy()
        if not present.any():
            continue
        err = tcond[mem].to_numpy()[present] - pcond[mem].to_numpy()[present]
        rows.append({"family": fam, "n_subtypes": len(mem),
                     "cond_rmse": float(np.sqrt((err ** 2).mean()))})
    return rows


def run_tissue(tissue):
    from tissueresolve.api import deconv_bulk
    from tissueresolve.reference.hierarchy import build_cell_type_hierarchy
    cfg = TISSUES[tissue]
    fine_col, donor_col = cfg["fine_col"], cfg["donor_col"]
    adata, mapping = _load(cfg)
    rows, perfam, runtime = [], [], []
    for seed in SEEDS:
        ref_d, qry = SH.split_donors(adata, donor_col, ref_frac=0.5, seed=seed)
        rmask = adata.obs[donor_col].astype(str).isin(set(ref_d)).to_numpy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = H.prepare_reference(adata[rmask].copy(), cell_type_col=fine_col,
                                      min_cells=MIN_CELLS, estimate_overdispersion=False).reference
        ref_types = [str(c) for c in ref.cell_types]
        mp = dict(build_cell_type_hierarchy(ref_types, mapping))
        tgt = SH.build_target_proportions(ref_types, N_SAMPLES, "imbalanced", seed=1000 + seed)
        ds = SH.realize_pseudobulk(adata, tgt, celltype_col=fine_col, donor_col=donor_col,
                                   query_donors=qry, seed=seed, cells_per_sample=CELLS)
        truth = ds.true_mrna_proportions

        # plain NNLS baseline on the CPM signature (external-style reference baseline)
        def nnls_baseline():
            R = ref.as_R_cpm().T  # (G, K)
            genes = [str(g) for g in ref.gene_names]
            common = [g for g in genes if g in set(ds.counts.index)]
            gi = [genes.index(g) for g in common]
            Rm = R[gi]; Rm = Rm / np.clip(Rm.sum(0, keepdims=True), 1e-9, None)
            B = ds.counts.loc[common].to_numpy(float); B = B / np.clip(B.sum(0, keepdims=True), 1e-9, None)
            P = np.zeros((B.shape[1], Rm.shape[1]))
            for j in range(B.shape[1]):
                x, _ = nnls(Rm, B[:, j]); s = x.sum(); P[j] = x / s if s > 0 else x
            return pd.DataFrame(P, index=ds.counts.columns, columns=ref_types)

        methods = {}
        t0 = time.time(); methods["NNLS_baseline"] = (nnls_baseline(), time.time() - t0)
        for name, kw in [("WNNLS_baseline", dict(solver="weighted_nnls", resolution_mode="none")),
                         ("TissueResolve_flat", dict(solver="nnls", resolution_mode="none")),
                         ("TissueResolve_hierarchical_soft",
                          dict(solver="nnls", resolution_mode="hierarchical", hierarchy_mapping=mp))]:
            t0 = time.time()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = deconv_bulk(ds.counts, ref, **kw)
            pred = res.deconv.proportions if name != "TissueResolve_hierarchical_soft" \
                else res.estimates.combined_fine
            methods[name] = (pred, time.time() - t0)

        for name, (pred, rt) in methods.items():
            rows.append({"tissue": tissue, "seed": seed, "method": name,
                         **_metrics(truth, pred, mp, ref_types)})
            runtime.append({"tissue": tissue, "seed": seed, "method": name,
                            "runtime_s": round(rt, 3), "status": "executed", "source": "fresh"})
            for r in _perfamily(truth, pred, mp, ref_types):
                perfam.append({"tissue": tissue, "seed": seed, "method": name, **r})
        print(f"[{tissue}] seed {seed} done", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(perfam), pd.DataFrame(runtime)


def _external_manifest():
    """Fold external-tool status + prior-run evidence into the run manifest (honest)."""
    rows = []
    prior = OUT / "real_external_method_status.tsv"
    prior_df = pd.read_csv(prior, sep="\t") if prior.exists() else pd.DataFrame()
    ext_tools = {"MuSiC": "bulk", "BisqueRNA": "bulk", "BayesPrism": "bulk",
                 "cell2location": "spatial", "CARD": "spatial", "spacexr_RCTD": "spatial",
                 "DWLS": "bulk", "SCDC": "bulk", "CIBERSORTx": "bulk",
                 "SPOTlight": "spatial", "Tangram": "spatial"}
    for tool, modality in ext_tools.items():
        key = "spacexr" if tool == "spacexr_RCTD" else tool
        m = prior_df[prior_df["method"].astype(str).str.contains(key, case=False)] if not prior_df.empty else prior_df
        if len(m):
            r = m.iloc[0]
            rows.append({"tool": tool, "modality": modality, "executed_this_session": False,
                         "prior_status": r.get("status", ""), "prior_runtime_s": r.get("runtime_seconds", ""),
                         "reason": "installed; not re-executed this session (runtime/stability); prior-run evidence used"})
        else:
            rows.append({"tool": tool, "modality": modality, "executed_this_session": False,
                         "prior_status": "no_prior_record",
                         "prior_runtime_s": "",
                         "reason": "not installed / not runnable (see EXTERNAL_TOOL_INSTALLATION_NOTES.md)"
                         if tool in ("DWLS", "SCDC", "CIBERSORTx", "SPOTlight", "Tangram")
                         else "installed; not re-executed this session"})
    return pd.DataFrame(rows)


def _bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment for a list of p-values (NaN-safe)."""
    p = np.asarray([np.nan if v is None else v for v in pvals], float)
    out = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx])]
    m = idx.size
    prev = 1.0
    for rank, j in enumerate(order[::-1]):
        k = m - rank
        prev = min(prev, p[j] * m / k)
        out[j] = min(prev, 1.0)
    return out


def robust_statistics(md: pd.DataFrame) -> tuple:
    """Paired robustness stats over the (tissue × seed) replicates.

    For each superiority-claim metric, compares ``TissueResolve_hierarchical_soft``
    against every other method with: paired mean difference, percentile bootstrap
    95% CI of the difference, paired Wilcoxon p, matched-pairs effect size, and a
    BH-adjusted p across all (metric × comparison) tests.  Also: best method per
    metric + TR rank, and rank-stability (mean/std rank) across the replicates.
    """
    methods = sorted(md["method"].unique())
    others = [m for m in methods if m != TR_HIER]
    keys = ["tissue", "seed"]
    pivots = {met: md.pivot_table(index=keys, columns="method", values=met) for met in SUPERIORITY_METRICS}

    rows = []
    for met in SUPERIORITY_METRICS:
        piv = pivots[met].dropna()
        asc = met in LOWER                       # lower better → TR "wins" if diff<0
        for other in others:
            if TR_HIER not in piv.columns or other not in piv.columns:
                continue
            a = piv[TR_HIER].to_numpy(float); b = piv[other].to_numpy(float)
            mask = np.isfinite(a) & np.isfinite(b); a, b = a[mask], b[mask]
            if a.size < 3:
                continue
            diff = a - b
            ci = ST.bootstrap_ci(diff, statistic=np.mean, n_boot=4000, seed=0)
            wil = ST.paired_wilcoxon(a, b)
            tr_better = (np.median(diff) < 0) if asc else (np.median(diff) > 0)
            # CI excludes 0 in the favourable direction?
            ci_sig = (ci["hi"] < 0) if asc else (ci["lo"] > 0)
            rows.append({
                "metric": met, "direction": "lower_better" if asc else "higher_better",
                "comparison": f"{TR_HIER} vs {other}", "n_pairs": int(a.size),
                "tr_mean": round(float(a.mean()), 4), "other_mean": round(float(b.mean()), 4),
                "mean_diff": round(float(diff.mean()), 4),
                "diff_ci_lo": round(ci["lo"], 4), "diff_ci_hi": round(ci["hi"], 4),
                "wilcoxon_p": round(wil["p_value"], 5) if np.isfinite(wil["p_value"]) else None,
                "effect_size": round(wil["effect_size"], 3),
                "tr_better_direction": bool(tr_better),
                "ci_excludes_zero_favourably": bool(ci_sig),
            })
    stats = pd.DataFrame(rows)
    if not stats.empty:
        stats["wilcoxon_p_bh"] = np.round(_bh_adjust(stats["wilcoxon_p"].tolist()), 5)
        stats["significant_bh_0.05"] = stats["wilcoxon_p_bh"].lt(0.05) & stats["ci_excludes_zero_favourably"]

    # best method per metric + TR rank (over all metrics)
    g = md.groupby("method").mean(numeric_only=True)
    best_rows = []
    all_metrics = [c for c in g.columns]
    for met in all_metrics:
        s = g[met].dropna()
        if s.empty:
            continue
        asc = met in LOWER
        order = s.sort_values(ascending=asc)
        tr_rank = (list(order.index).index(TR_HIER) + 1) if TR_HIER in order.index else None
        best_rows.append({"metric": met, "direction": "lower_better" if asc else "higher_better",
                          "best_method": order.index[0], "best_value": round(float(order.iloc[0]), 4),
                          "TR_hier_value": round(float(g.loc[TR_HIER, met]), 4) if TR_HIER in g.index else None,
                          "TR_hier_rank": f"{tr_rank}/{len(order)}" if tr_rank else None})
    best = pd.DataFrame(best_rows)

    # rank stability across the (tissue,seed) replicates for the key metrics
    rk_rows = []
    for met in SUPERIORITY_METRICS:
        asc = met in LOWER
        ranks = pivots[met].rank(axis=1, ascending=asc)
        for meth in methods:
            if meth in ranks.columns:
                rk_rows.append({"metric": met, "method": meth,
                                "mean_rank": round(float(ranks[meth].mean()), 2),
                                "rank_std": round(float(ranks[meth].std()), 2),
                                "n_replicates": int(ranks[meth].notna().sum())})
    rank_stab = pd.DataFrame(rk_rows)
    return stats, best, rank_stab


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--tissue", choices=["breast", "lung", "both"], default="both")
    args = ap.parse_args(argv)
    if not args.run_real_data:
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    OUT.mkdir(parents=True, exist_ok=True); EXT.mkdir(parents=True, exist_ok=True)
    tissues = ["breast", "lung"] if args.tissue == "both" else [args.tissue]
    M_all, P_all, R_all = [], [], []
    for tissue in tissues:
        if not Path(TISSUES[tissue]["h5ad"]).exists():
            print(f"Missing {TISSUES[tissue]['h5ad']}; skipping {tissue}.", file=sys.stderr); continue
        m, p, r = run_tissue(tissue)
        M_all.append(m); P_all.append(p); R_all.append(r)
    md = pd.concat(M_all, ignore_index=True)
    pf = pd.concat(P_all, ignore_index=True)
    rt = pd.concat(R_all, ignore_index=True)
    md.to_csv(OUT / "external_bulk_metrics.tsv", sep="\t", index=False)
    pf.to_csv(OUT / "external_bulk_per_family_metrics.tsv", sep="\t", index=False)
    rt.to_csv(OUT / "external_bulk_runtime.tsv", sep="\t", index=False)
    # run manifest: fresh internal runs + external-tool status
    fresh = rt.assign(executed_this_session=True, prior_status="", prior_runtime_s="",
                      reason="fresh internal run").rename(columns={"method": "tool"})[
        ["tool", "executed_this_session", "prior_status", "prior_runtime_s", "reason"]].drop_duplicates("tool")
    manifest = pd.concat([fresh, _external_manifest()[["tool", "executed_this_session",
                          "prior_status", "prior_runtime_s", "reason"]]], ignore_index=True)
    manifest.to_csv(EXT / "run_manifest.tsv", sep="\t", index=False)

    # robust paired statistics (expanded seeds → power)
    stats, best, rank_stab = robust_statistics(md)
    stats.to_csv(OUT / "external_benchmark_statistical_tests.tsv", sep="\t", index=False)
    best.to_csv(OUT / "external_benchmark_best_per_metric.tsv", sep="\t", index=False)
    rank_stab.to_csv(OUT / "external_benchmark_rank_stability.tsv", sep="\t", index=False)

    g = md.groupby(["tissue", "method"]).mean(numeric_only=True)
    n_rep = md.groupby("method").size().max()
    print(f"\n=== bulk benchmark — method means ({len(SEEDS)} seeds × {len(tissues)} tissues "
          f"= up to {n_rep} replicates/method) ===")
    print(md.groupby("method").mean(numeric_only=True)[
        ["broad_pearson", "broad_rmse", "fine_pearson", "fine_rmse", "cond_rmse",
         "pairwise_spillover", "rare_sensitivity", "rare_precision", "false_positive_rate"]].round(4).to_string())
    print("\n=== robust paired superiority stats (TR_hier vs each baseline) ===")
    if not stats.empty:
        print(stats[["metric", "comparison", "n_pairs", "mean_diff", "diff_ci_lo",
                     "diff_ci_hi", "wilcoxon_p_bh", "significant_bh_0.05"]].to_string(index=False))
    print(f"\nwrote metrics/per_family/runtime/statistical_tests/best_per_metric/rank_stability; "
          f"run_manifest.tsv ({len(manifest)} tools)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
