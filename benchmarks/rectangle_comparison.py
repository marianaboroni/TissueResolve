#!/usr/bin/env python
"""Rectangle vs TissueResolve — single, reproducible, config-driven benchmark.

One command produces one long-format result table, one summary, one manifest, one
HTML report, and a handful of figures under a single output directory:

  python benchmarks/rectangle_comparison.py \
      --config benchmarks/rectangle_config.yaml \
      --output benchmarks/results/rectangle_benchmark --run-real-data

Design goals (see the continuation brief):
  * fair, scaled, reproducible, easy-to-read; NO proliferation of folders/reports;
  * everything in ONE long table (benchmark_results.parquet) + a readable summary;
  * reproducibility manifest (git, versions, seeds, config, timings, failures);
  * --dry-run (show the run matrix + estimates) and --resume (skip completed runs);
  * changes NO production default and NO core algorithm — bulk scenarios are
    query-side confounders applied only inside this harness.

Rectangle runs out-of-process in benchmarks/envs/rectangle_py311 (Python 3.11) via
benchmarks/rectangle_env_runner.py; it is skipped honestly if not provisioned.

Ground truth is mRNA proportions (what bulk deconvolvers estimate); Rectangle is
scored on its `raw` (bias-uncorrected) output and its `mrnacorrected` default is
recorded for transparency.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.shared import environment as ENV  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

ENV_RUNNER = REPO / "benchmarks" / "rectangle_env_runner.py"
BREAST_HIERARCHY = (REPO / "examples" / "real_breast_cancer" / "config"
                    / "breast_cancer_cell_type_hierarchy.tsv")
MIN_CELLS = 30
RARE_MAX = 0.05          # truth in (0, RARE_MAX] counts as a rare presence


# ===========================================================================
# config / manifest
# ===========================================================================
def load_config(path: Path) -> dict:
    return yaml.safe_load(Path(path).read_text())


def _git(*args) -> str:
    try:
        return subprocess.run(["git", *args], cwd=str(REPO), capture_output=True,
                              text=True, timeout=15).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def build_manifest(config: dict, argv: list) -> dict:
    def _ver(mod):
        try:
            import importlib.metadata as m
            return m.version(mod)
        except Exception:
            return "not installed"
    rect_py = ENV.rectangle_env_python()
    rect_ver = ""
    if rect_py is not None:
        try:
            rect_ver = subprocess.run(
                [str(rect_py), "-c", "import rectanglepy as r; print(r.version('rectanglepy'))"],
                capture_output=True, text=True, timeout=60).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            rect_ver = "unknown"
    return {
        "git_commit": _git("rev-parse", "HEAD"),
        "git_status_dirty": bool(_git("status", "--porcelain")),
        "python": platform.python_version(),
        "os": platform.platform(),
        "cpu": platform.processor() or platform.machine(),
        "versions": {m: _ver(m) for m in ("numpy", "scipy", "pandas", "anndata",
                                          "scikit-learn", "tissueresolve")},
        "rectanglepy": rect_ver,
        "rectangle_env": str(rect_py) if rect_py else None,
        "config": config,
        "command": "python " + " ".join(argv),
    }


def run_id(dataset, scenario, seed, method) -> str:
    return f"{dataset}__{scenario}__seed{seed}__{method}"


# ===========================================================================
# small numeric utilities
# ===========================================================================
def _tpm(counts_genes_by_samples: pd.DataFrame) -> pd.DataFrame:
    c = counts_genes_by_samples
    return c.T.div(c.sum(axis=0).replace(0, np.nan), axis=0).fillna(0.0) * 1e6


def _mapping_for(tr_ref, D, dataset):
    types = [str(c) for c in tr_ref.cell_types]
    if dataset == "breast":
        from tissueresolve.reference.hierarchy import (
            load_hierarchy_mapping, build_cell_type_hierarchy)
        raw = load_hierarchy_mapping(BREAST_HIERARCHY)
        return build_cell_type_hierarchy(types, {t: raw.get(t, t) for t in types})
    return {t: D["mapping"].get(t, t) for t in types}


def _build_ref(adata_sub, celltype_col):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return H.prepare_reference(adata_sub.copy(), min_cells=MIN_CELLS,
                                   cell_type_col=celltype_col,
                                   estimate_overdispersion=True).reference


def _specificity_panel(tr_ref, size):
    """Top-`size` genes by row-max specificity from the reference CPM matrix."""
    R = tr_ref.as_R_cpm()                       # (K, G)
    genes = list(tr_ref.gene_names)
    col_max = R.max(axis=0)
    col_sum = R.sum(axis=0) + 1e-9
    spec = col_max / col_sum                    # 1 = perfectly specific
    order = np.argsort(-spec)[:size]
    return [genes[i] for i in order]


# ===========================================================================
# scenario reference/query construction (+ query-side confounders)
# ===========================================================================
def _perturb(counts, scenario, cfg, seed):
    """Apply a query-side technical confounder to genes×samples counts."""
    rng = np.random.default_rng(1000 + seed)
    p = cfg["perturbations"]
    X = counts.to_numpy(dtype=np.float64)
    if scenario == "low_depth":
        X = rng.binomial(np.rint(X).astype(np.int64),
                         p["low_depth_fraction"]).astype(np.float64)
    elif scenario == "overdispersion":
        cv = p["overdispersion_cv"]; k = 1.0 / (cv * cv)
        g = rng.gamma(k, 1.0 / k, size=(X.shape[0], 1))     # per-gene multiplier
        X = np.rint(X * g)
    elif scenario == "technical_noise":
        cv = p["technical_noise_cv"]
        f = rng.lognormal(0.0, cv, size=X.shape)
        X = np.rint(X * f)
    elif scenario == "gene_dropout":
        mask = rng.random(X.shape) < p["gene_dropout_fraction"]
        X[mask] = 0.0
    out = pd.DataFrame(X, index=counts.index, columns=counts.columns)
    return out


def _scenario_refquery(scenario, D, adata, celltype_col, donor_col):
    """Return (ref_adata_sub, tr_reference, query_adata_sub, note) or None if N/A."""
    obs = adata.obs
    ref_donors, query_donors = set(D["ref_donors"]), set(D["query_donors"])
    ref_mask = obs[donor_col].astype(str).isin(ref_donors).to_numpy()
    q_mask = obs[donor_col].astype(str).isin(query_donors).to_numpy()

    # query-side confounder scenarios reuse the complete-reference split
    if scenario in ("complete_reference", "high_collinearity", "low_depth",
                    "overdispersion", "technical_noise", "gene_dropout"):
        return (adata[ref_mask].copy(), D["ref"], adata[q_mask].copy(),
                scenario.replace("_", " "))
    if scenario == "cross_platform":
        if "assay" not in obs.columns or obs["assay"].nunique() < 2:
            return None
        assay = obs["assay"].astype(str)
        pa, pb = assay.value_counts().index.tolist()[:2]
        ref_sub = adata[(assay == pa).to_numpy()].copy()
        return (ref_sub, _build_ref(ref_sub, celltype_col),
                adata[(assay == pb).to_numpy()].copy(), f"train {pa} -> test {pb}")
    if scenario in ("missing_rare", "missing_abundant", "missing_activated_state"):
        types = obs.loc[ref_mask, celltype_col].astype(str)
        if scenario == "missing_rare":
            drop = D["rare_type"]
        elif scenario == "missing_abundant":
            drop = types.value_counts().index[0]
        else:
            cand = [t for t in types.unique()
                    if any(k in t.lower() for k in ("activat", "exhaust", "prolifer"))]
            if not cand:
                return None
            drop = cand[0]
        keep = ref_mask & (obs[celltype_col].astype(str) != drop).to_numpy()
        ref_sub = adata[keep].copy()
        return (ref_sub, _build_ref(ref_sub, celltype_col), adata[q_mask].copy(),
                f"removed from reference: {drop}")
    return None


# ===========================================================================
# method runners  ->  (pred_raw_df, resolution, config_label, audit)
# ===========================================================================
def _tr_solver(counts, reference, name):
    from tissueresolve.solver import (NNLSSolver, WeightedNNLSSolver, AutoSolver)
    if name == "nnls":
        return NNLSSolver().solve(counts, reference), {}
    if name == "wnnls":
        return WeightedNNLSSolver().solve(counts, reference), {}
    auto = AutoSolver(n_splits=2)
    best, comp, reason = auto.select(counts, reference)
    res = best.solve(counts, reference)
    top = comp.reset_index().iloc[0].to_dict()
    audit = {"selected_solver": best.name,
             "auto_objective": float(top.get("objective", np.nan)),
             "condition_number": float(top.get("condition_number", np.nan)),
             "masked_gene_cv": float(top.get("cv_masked_gene_pearson", np.nan))}
    return res, audit


def run_method(method, counts, tr_ref, ref_sub, scen_mapping, celltype_col,
               cfg, rect_cache):
    """Return (pred_raw, resolution, config_label, audit). pred_raw may carry
    unresolved_/Unknown columns. rect_cache holds Rectangle variants per run."""
    import tissueresolve as tr
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if method == "TissueResolve_flat_default":
            r = tr.deconv_bulk(counts, tr_ref, resolution_mode="flat", n_bootstrap=0)
            return r.deconv.proportions, "flat", "solver=default;flat", {}
        if method == "TissueResolve_hierarchical":
            r = tr.deconv_bulk(counts, tr_ref, resolution_mode="hierarchical",
                               hierarchy_mapping=scen_mapping, n_bootstrap=0)
            return r.deconv.proportions, "hierarchical", "soft_gating", {}
        if method == "TissueResolve_auto":
            res, audit = _tr_solver(counts, tr_ref, "auto")
            return res.proportions, "flat", f"solver=auto({audit.get('selected_solver')})", audit
        if method == "NNLS_baseline":
            res, _ = _tr_solver(counts, tr_ref, "nnls")
            return res.proportions, "flat", "nnls", {}
        if method == "Weighted_NNLS_baseline":
            res, _ = _tr_solver(counts, tr_ref, "wnnls")
            return res.proportions, "flat", "weighted_nnls", {}
        if method in ("Poisson_GLM", "NB_GLM"):
            from tissueresolve.experimental.nb_bulk_solver import NBGLMBulkSolver
            loss = "poisson" if method == "Poisson_GLM" else "nb"
            panel = _specificity_panel(tr_ref, cfg["glm"]["panel_size"])
            res = NBGLMBulkSolver(loss=loss).solve(counts, tr_ref, gene_panel=panel)
            return res.proportions, "flat", f"glm:{loss};panel={len(panel)}", {}
        if method in ("Rectangle_raw", "Rectangle_mrnacorrected"):
            variant = "raw" if method.endswith("raw") else "mrnacorrected"
            r = rect_cache
            if variant in r:
                return r[variant], "flat", f"rectangle:{variant}", {"runtime_s": r.get("runtime_s", 0.0)}
            return None, "flat", f"rectangle:{variant}", {"status": r.get("status", "failed")}
    raise ValueError(f"unknown method {method}")


def run_rectangle(ref_adata_sub, counts, celltype_col, cfg, seed):
    """Run Rectangle once (both variants) out-of-process; return a cache dict."""
    py = ENV.rectangle_env_python()
    if py is None:
        return {"status": "skipped"}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        ref_adata_sub.write_h5ad(tmp / "ref.h5ad")
        _tpm(counts).to_csv(tmp / "bulk_tpm.tsv", sep="\t")
        out_dir = tmp / "out"
        cmd = [str(py), str(ENV_RUNNER), "--ref-h5ad", str(tmp / "ref.h5ad"),
               "--cell-type-col", celltype_col, "--bulk-tpm", str(tmp / "bulk_tpm.tsv"),
               "--out-dir", str(out_dir), "--cells-per-type",
               str(cfg["rectangle"]["ref_cells_per_type"]), "--seed", str(seed),
               "--bootstraps", str(cfg["rectangle"]["bootstraps"])]
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))
        rt = time.perf_counter() - t0
        meta_f = out_dir / "rectangle_metadata.json"
        meta = json.loads(meta_f.read_text()) if meta_f.exists() else {"status": "failed"}
        cache = {"status": meta.get("status", "failed"), "runtime_s": rt}
        for v in ("raw", "mrnacorrected"):
            f = out_dir / f"predictions_{v}.tsv"
            if f.exists():
                cache[v] = pd.read_csv(f, sep="\t", index_col=0)
        return cache


# ===========================================================================
# scoring -> long rows
# ===========================================================================
def _align(truth, pred):
    keep = [c for c in pred.columns
            if not str(c).startswith("unresolved_") and str(c) != "Unknown"]
    return pred[keep].reindex(index=truth.index, columns=truth.columns).fillna(0.0)


def _conditional(df, cols, mapping):
    fam_members = {}
    for c in cols:
        fam_members.setdefault(str(mapping.get(c, c)), []).append(c)
    out = pd.DataFrame(0.0, index=df.index, columns=cols)
    multi = []
    for mem in fam_members.values():
        mem = [m for m in mem if m in df.columns]
        if len(mem) < 2:
            continue
        multi += mem
        s = df[mem].sum(axis=1).replace(0, np.nan)
        for m in mem:
            out[m] = (df[m] / s).fillna(0.0)
    return out, multi, fam_members


def _coverage_risk(truth, pred_raw, mapping):
    cols = list(truth.columns)
    _, _, fam_members = _conditional(truth, cols, mapping)
    multi = {f: [m for m in mem if m in cols] for f, mem in fam_members.items()
             if len([m for m in mem if m in cols]) > 1}
    if not multi:
        return float("nan"), float("nan")
    accepted = {}
    for fam, mem in multi.items():
        unres = f"unresolved_{fam}"
        fine = pred_raw[[m for m in mem if m in pred_raw.columns]].sum(axis=1).mean() \
            if any(m in pred_raw.columns for m in mem) else 0.0
        ures = pred_raw[unres].mean() if unres in pred_raw.columns else 0.0
        accepted[fam] = fine >= ures
    coverage = float(np.mean(list(accepted.values())))
    acc_members = [m for fam, mem in multi.items() if accepted[fam] for m in mem]
    if not acc_members:
        return coverage, float("nan")
    p = _align(truth, pred_raw)
    tc, _, _ = _conditional(truth, acc_members, mapping)
    pc, _, _ = _conditional(p, acc_members, mapping)
    risk = float(np.sqrt(np.mean((tc[acc_members].to_numpy(float)
                                  - pc[acc_members].to_numpy(float)) ** 2)))
    return coverage, risk


def score_rows(truth, pred_raw, mapping, thresholds, base):
    """Emit long-format rows for one method run."""
    pred = _align(truth, pred_raw)
    rows = []

    def add(metric, value, cell_type="", family="", threshold=np.nan):
        rows.append({**base, "metric": metric, "cell_type": cell_type,
                     "family": family, "threshold": threshold, "value": float(value)})

    m = H.compute_bulk_metrics(truth, pred)
    for k in ("rmse", "mae", "pearson", "spearman", "signed_bias"):
        add(k, m[k])
    pct = H.per_celltype_metrics(truth, pred)
    add("mean_per_type_rmse", float(pct["rmse"].mean()))
    for ct in pct.index:
        add("per_type_rmse", pct.loc[ct, "rmse"], cell_type=str(ct))

    tc, multi, _ = _conditional(truth, list(truth.columns), mapping)
    pc, _, _ = _conditional(pred, list(truth.columns), mapping)
    if multi:
        add("conditional_rmse", float(np.sqrt(np.mean(
            (tc[multi].to_numpy(float) - pc[multi].to_numpy(float)) ** 2))))
    cov, risk = _coverage_risk(truth, pred_raw, mapping)
    add("coverage", cov)
    add("selective_risk", risk)

    t = truth.to_numpy(float)
    p = pred.to_numpy(float)
    absent = t < 1e-9
    present = ~absent
    rare_present = present & (t < RARE_MAX)
    add("absent_mass_mean", float(np.where(absent, p, 0.0).sum(axis=1).mean()))
    if "Unknown" in pred_raw.columns:
        add("unknown_mass_mean", float(pred_raw["Unknown"].mean()))
    for thr in thresholds:
        pos = p > thr
        fp = int((pos & absent).sum()); tn = int((~pos & absent).sum())
        tp = int((pos & present).sum()); fn = int((~pos & present).sum())
        add("rare_fpr", fp / (fp + tn) if (fp + tn) else np.nan, threshold=thr)
        add("presence_recall", tp / (tp + fn) if (tp + fn) else np.nan, threshold=thr)
        add("presence_precision", tp / (tp + fp) if (tp + fp) else np.nan, threshold=thr)
        add("rare_recall", (int((pos & rare_present).sum()) / int(rare_present.sum()))
            if rare_present.sum() else np.nan, threshold=thr)
    return rows


# ===========================================================================
# matrix / driver
# ===========================================================================
def build_matrix(config):
    runs = []
    for ds in config["datasets"]:
        for scen in config["bulk_scenarios"]:
            for seed in config["seeds"]:
                for method in config["methods"]:
                    runs.append(dict(modality="bulk", dataset=ds, scenario=scen,
                                     seed=seed, method=method,
                                     run_id=run_id(ds, scen, seed, method)))
    return runs


def _estimate(config, matrix):
    n = len(matrix)
    rect = sum(1 for r in matrix if r["method"].startswith("Rectangle"))
    tr = n - rect
    # rough per-run seconds (breast fast, lung slower); Rectangle dominates
    est_s = rect * 120 + tr * 8
    return {"n_runs": n, "rectangle_runs": rect, "other_runs": tr,
            "rough_runtime_min": round(est_s / 60, 1),
            "parquet_est_mb": round(n * 60 * 200 / 1e6, 2)}  # ~60 rows/run, ~200 B/row


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "benchmarks" / "rectangle_config.yaml"))
    ap.add_argument("--output", default=str(REPO / "benchmarks" / "results" / "rectangle_benchmark"))
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=None, help="override config datasets")
    ap.add_argument("--scenarios", nargs="+", default=None)
    ap.add_argument("--seeds", type=int, nargs="+", default=None)
    ap.add_argument("--methods", nargs="+", default=None)
    args = ap.parse_args(argv)

    config = load_config(args.config)
    for key, ov in (("datasets", args.datasets), ("bulk_scenarios", args.scenarios),
                    ("seeds", args.seeds), ("methods", args.methods)):
        if ov is not None:
            config[key] = ov
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    matrix = build_matrix(config)
    est = _estimate(config, matrix)

    if args.dry_run:
        print(f"RUN MATRIX: {est['n_runs']} runs "
              f"({est['rectangle_runs']} Rectangle, {est['other_runs']} other)")
        print(f"Rough runtime ~{est['rough_runtime_min']} min; "
              f"parquet ~{est['parquet_est_mb']} MB")
        print("datasets:", config["datasets"], "| scenarios:", config["bulk_scenarios"])
        print("seeds:", config["seeds"], "| methods:", config["methods"])
        return 0

    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    results_path = out / "benchmark_results.parquet"
    done = set()
    prior_rows = []
    if args.resume and results_path.exists():
        prior = pd.read_parquet(results_path)
        done = set(prior["run_id"].unique())
        prior_rows = prior.to_dict("records")
        print(f"--resume: {len(done)} runs already complete; skipping them.")

    rows = list(prior_rows)
    failures = []
    t_start = time.perf_counter()
    thresholds = config["presence_thresholds"]

    for ds in config["datasets"]:
        D = _load_dataset(ds)
        adata, celltype_col, donor_col = D["adata"], D["celltype_col"], D["donor_col"]
        for scen in config["bulk_scenarios"]:
            rq = _scenario_refquery(scen, D, adata, celltype_col, donor_col)
            if rq is None:
                print(f"[{ds}/{scen}] N/A (skipped)")
                continue
            ref_sub, tr_ref, q_sub, note = rq
            scen_mapping = _mapping_for(tr_ref, D, ds)
            for seed in config["seeds"]:
                counts, truth, _ = H.generate_pseudobulk(
                    q_sub, celltype_col, n_per_regime=config["pseudobulk"]["n_per_regime"],
                    n_cells=config["pseudobulk"]["n_cells"], seed=seed)
                counts = _perturb(counts, scen, config, seed)
                rect_cache = None
                need_rect = any(m.startswith("Rectangle") for m in config["methods"]) \
                    and any(run_id(ds, scen, seed, m) not in done
                            for m in config["methods"] if m.startswith("Rectangle"))
                if need_rect:
                    rect_cache = run_rectangle(ref_sub, counts, celltype_col, config, seed)
                for method in config["methods"]:
                    rid = run_id(ds, scen, seed, method)
                    if rid in done:
                        continue
                    base = dict(modality="bulk", dataset=ds, scenario=scen, seed=seed,
                                method=method, run_id=rid, note=note)
                    t0 = time.perf_counter()
                    try:
                        pred, resolution, cfg_label, audit = run_method(
                            method, counts, tr_ref, ref_sub, scen_mapping,
                            celltype_col, config, rect_cache or {})
                        if pred is None:
                            status = audit.get("status", "failed")
                            failures.append((rid, status))
                            rows.append({**base, "configuration": "", "resolution": resolution,
                                         "metric": "status", "cell_type": "", "family": "",
                                         "threshold": np.nan, "value": np.nan,
                                         "runtime_seconds": 0.0, "status": status})
                            continue
                        rt = audit.get("runtime_s", time.perf_counter() - t0)
                        mrows = score_rows(truth, pred, scen_mapping, thresholds, base)
                        for r in mrows:
                            r.update(configuration=cfg_label, resolution=resolution,
                                     runtime_seconds=round(rt, 3), status="ok")
                        rows += mrows
                        for k, v in audit.items():
                            if k in ("selected_solver",):
                                continue
                            if isinstance(v, (int, float)) and k != "runtime_s":
                                rows.append({**base, "configuration": cfg_label,
                                             "resolution": resolution, "metric": f"audit_{k}",
                                             "cell_type": "", "family": "", "threshold": np.nan,
                                             "value": float(v), "runtime_seconds": round(rt, 3),
                                             "status": "ok"})
                        if "selected_solver" in audit:
                            rows.append({**base, "configuration": cfg_label, "resolution": resolution,
                                         "metric": "audit_selected_solver_is_nnls",
                                         "cell_type": "", "family": "", "threshold": np.nan,
                                         "value": float(audit["selected_solver"] == "nnls"),
                                         "runtime_seconds": round(rt, 3), "status": "ok"})
                        print(f"[{ds}/{scen}/seed{seed}] {method}: ok ({round(rt,1)}s)")
                    except Exception as exc:  # noqa: BLE001
                        failures.append((rid, f"{type(exc).__name__}: {exc}"))
                        rows.append({**base, "configuration": "", "resolution": "",
                                     "metric": "status", "cell_type": "", "family": "",
                                     "threshold": np.nan, "value": np.nan,
                                     "runtime_seconds": round(time.perf_counter() - t0, 3),
                                     "status": f"failed:{type(exc).__name__}"})
                        print(f"[{ds}/{scen}/seed{seed}] {method}: FAILED {exc}")

    df = pd.DataFrame(rows)
    for c in ("cell_type", "family", "note", "configuration", "resolution", "status"):
        if c not in df.columns:
            df[c] = ""
    df.to_parquet(results_path, index=False)

    # readable summary: scalar metrics, seed means
    scal = df[(df["status"] == "ok") & (df["cell_type"] == "") & (df["threshold"].isna())]
    if not scal.empty:
        summary = (scal.pivot_table(index=["dataset", "scenario", "method"],
                                    columns="metric", values="value", aggfunc="mean")
                   .round(4).reset_index())
        summary.to_csv(out / "benchmark_summary.tsv", sep="\t", index=False)

    manifest = build_manifest(config, ["benchmarks/rectangle_comparison.py"] + argv)
    manifest.update(total_runtime_s=round(time.perf_counter() - t_start, 1),
                    n_result_rows=len(df), failures=failures, estimate=est)
    (out / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2))

    try:
        from benchmarks.rectangle_report import build_report  # optional
        build_report(out)
    except Exception as exc:  # noqa: BLE001
        print(f"(report/figures step skipped: {exc})")

    print(f"\nWrote {results_path} ({len(df)} rows); {len(failures)} failures; "
          f"{round(manifest['total_runtime_s']/60,1)} min total.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
