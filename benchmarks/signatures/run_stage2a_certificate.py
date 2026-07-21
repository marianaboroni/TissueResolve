#!/usr/bin/env python
"""Stage 2A — Experiment A validation: does the identifiability certificate PREDICT, before any
deconvolution, which cell types are recoverable on held-out donor-disjoint pseudobulk?

Donor-disjoint TRAIN / CALIBRATION / TEST split:
  * TRAIN donors  -> reference (R). The certificate is computed from R + query-detectable genes
    + query depth + per-type reference donor counts. It uses NO held-out labels.
  * CALIBRATION donors -> the class->P(recoverable) map is fit here (empirical recoverable rate per
    class). Nothing is fit on TEST.
  * TEST donors   -> every reported metric (held-out RMSE, conditional RMSE, spillover, absent-type
    false positives, detection & rare recall, uncertainty coverage) is measured here.

Identifiability is depth/panel dependent, so we sweep query regimes (bulk_full / bulk_shallow /
spot_like / panel_1000) to force a spread of recoverability classes. Deconvolver = production
per-sample Poisson GLM. Uncertainty = gene-subsampling bootstrap around that same solver (validation
tooling; no new solver, no default change). Thresholds are pre-specified and never tuned on TEST.

Usage:
  PYTHONPATH=src:. python benchmarks/signatures/run_stage2a_certificate.py --run-real-data \
      --datasets breast lung --seeds 0 1 2 3 4
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.reference.identifiability import identifiability_certificate  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "signatures" / "stage2a" / "certificate"

# ---- pre-specified thresholds (NOT tuned on test) ----
REC_RMSE_THRESH = 0.10        # a type is empirically "recoverable" if held-out RMSE < this
DETECT_THRESH = 0.05          # a type is "present"/"detected" in a sample if proportion > this
RARE_BAND = (1e-9, 0.05)      # a present instance is "rare" if 0 < true proportion <= 0.05
RARE_DETECT = 0.01            # rare instance counted as recalled if predicted proportion > this
NOMINAL_COVERAGE = 0.90       # target for uncertainty intervals
BOOT_B = 20                   # gene-subsampling bootstrap replicates
BOOT_FRAC = 0.8               # fraction of panel genes per bootstrap replicate
BOOT_MAX_SAMPLES = 24         # cap bootstrap solves for runtime
CLASS_PROB_PRIOR = {"RESOLVABLE": 0.90, "WEAKLY_RESOLVABLE": 0.60,
                    "GROUP_ONLY": 0.25, "UNRESOLVABLE": 0.10}   # fallback if calibration is thin

CONDITIONS = [
    ("bulk_full", 1.0, None),        # ~1e6 counts, all genes
    ("bulk_shallow", 0.02, None),    # ~2% depth
    ("spot_like", 0.005, None),      # Visium-spot-like shallow depth
    ("panel_1000", 1.0, 1000),       # small arbitrary gene panel, full depth
]
STRESSED = {"breast": "spot_like", "lung": "panel_1000"}   # condition with a real class spread


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _thin(counts, factor, rng):
    if factor >= 1.0:
        return counts
    vals = rng.binomial(counts.to_numpy().astype(np.int64), factor)
    return pd.DataFrame(vals, index=counts.index, columns=counts.columns)


def _cond_rmse_by_type(truth, pred, mapping, types):
    """Within-family conditional-composition RMSE per type (member / family-sum, family present)."""
    fam = defaultdict(list)
    for c in types:
        fam[str(mapping.get(c, c))].append(c)
    out = {c: np.nan for c in types}
    for members in fam.values():
        if len(members) < 2:
            continue
        T = truth.reindex(columns=members).fillna(0).to_numpy(float)
        P = pred.reindex(columns=members).fillna(0).to_numpy(float)
        ts, ps = T.sum(1), P.sum(1)
        present = ts > DETECT_THRESH
        if present.sum() < 3:
            continue
        tc = T[present] / ts[present, None]
        pc = P[present] / np.clip(ps[present, None], 1e-9, None)
        for i, m in enumerate(members):
            out[m] = float(np.sqrt(np.mean((tc[:, i] - pc[:, i]) ** 2)))
    return out


def _per_type_metrics(truth, pred, mapping, types):
    cond = _cond_rmse_by_type(truth, pred, mapping, types)
    rows = {}
    for k in types:
        tv = truth[k].to_numpy(float) if k in truth else np.zeros(len(truth))
        pv = pred[k].to_numpy(float) if k in pred else np.zeros(len(truth))
        present = tv > DETECT_THRESH
        absent = tv <= 1e-9
        rare = (tv > RARE_BAND[0]) & (tv <= RARE_BAND[1])
        rows[k] = dict(
            held_out_rmse=float(np.sqrt(np.mean((pv - tv) ** 2))),
            conditional_rmse=cond[k],
            spillover_when_absent=float(pv[absent].mean()) if absent.any() else np.nan,
            detection_recall=float((pv[present] > DETECT_THRESH).mean()) if present.any() else np.nan,
            rare_recall=float((pv[rare] > RARE_DETECT).mean()) if rare.any() else np.nan,
            false_positive_rate=float((pv[absent] > DETECT_THRESH).mean()) if absent.any() else np.nan,
            mean_true=float(tv.mean()))
    return rows


def _fit_class_prob(cal_metrics):
    """Empirical P(recoverable | class) on CALIBRATION donors; prior fallback if a class is thin."""
    out = dict(CLASS_PROB_PRIOR)
    d = cal_metrics[cal_metrics.recoverability.isin(CLASS_PROB_PRIOR)].copy()
    d["y"] = (d.held_out_rmse < REC_RMSE_THRESH).astype(float)
    for cls, grp in d.groupby("recoverability"):
        if len(grp) >= 5:
            out[cls] = float(np.clip(grp.y.mean(), 0.02, 0.98))
    return out


def _calibration(test_metrics, class_prob):
    d = test_metrics[test_metrics.recoverability.isin(class_prob)].copy()
    if d.empty:
        return {}
    d["p"] = d.recoverability.map(class_prob)
    d["y"] = (d.held_out_rmse < REC_RMSE_THRESH).astype(float)
    brier = float(np.mean((d.p - d.y) ** 2))
    ece = float(sum((len(g) / len(d)) * abs(p - g.y.mean()) for p, g in d.groupby("p")))
    rec = d[d.y == 1.0]
    coverage = float(rec.recoverability.isin(["RESOLVABLE", "WEAKLY_RESOLVABLE"]).mean()) if len(rec) else np.nan
    return {"brier": round(brier, 4), "ece": round(ece, 4),
            "resolved_coverage": round(coverage, 4), "n_testable": int(len(d))}


def _false_rates(test_metrics):
    d = test_metrics[test_metrics.recoverability.isin(CLASS_PROB_PRIOR)].copy()
    d["y"] = d.held_out_rmse < REC_RMSE_THRESH
    resolved = d[d.recoverability == "RESOLVABLE"]
    merged = d[d.recoverability.isin(["UNRESOLVABLE", "GROUP_ONLY"])]
    fr = float((~resolved.y).mean()) if len(resolved) else np.nan
    fm = float(merged.y.mean()) if len(merged) else np.nan
    return {"false_resolution_rate": round(fr, 4) if fr == fr else np.nan,
            "false_merge_rate": round(fm, 4) if fm == fm else np.nan,
            "n_resolvable": int(len(resolved)), "n_merged": int(len(merged))}


def _group_aggregate_recovery(cert, truth, pred, types):
    """GROUP_ONLY claims the group SUM is recoverable but the WITHIN-group split is not.
    Within-group conditional (swap) RMSE isolates member swaps; abundance-normalised."""
    out = []
    for g in cert.confounded_groups:
        members = [m for m in g["group"].split(";") if m in types]
        if len(members) < 2:
            continue
        T = truth.reindex(columns=members).fillna(0).to_numpy(float)
        P = pred.reindex(columns=members).fillna(0).to_numpy(float)
        ts, ps = T.sum(1), P.sum(1)
        group_rmse = float(np.sqrt(np.mean((ps - ts) ** 2)))
        present = ts > 0.05
        within = np.nan
        if present.sum() >= 3:
            tc = T[present] / ts[present, None]
            pc = P[present] / np.clip(ps[present, None], 1e-9, None)
            within = float(np.sqrt(np.mean((tc - pc) ** 2)))
        indiv = float(np.mean([np.sqrt(np.mean((P[:, i] - T[:, i]) ** 2)) for i in range(len(members))]))
        out.append(dict(group=";".join(members), n_members=len(members),
                        group_sum_rmse=round(group_rmse, 4),
                        within_group_conditional_rmse=round(within, 4) if within == within else np.nan,
                        n_samples_group_present=int(present.sum()),
                        mean_individual_rmse=round(indiv, 4),
                        group_recoverable_claim=g["recoverable_group"]))
    return out


def _bootstrap_coverage(counts, ref, truth, cert, types, panel, rng, PoissonGLMSolver):
    """Gene-subsampling bootstrap around the production solver; empirical coverage of nominal CIs,
    by recoverability class. Subsampling (m-out-of-n, no replacement) avoids duplicate-gene issues."""
    avail = [g for g in (panel or list(ref.gene_names)) if g in counts.index]
    if len(avail) < 10:
        return {}, np.nan
    cols = list(counts.columns)[:BOOT_MAX_SAMPLES]
    c_sub = counts[cols]
    tr = truth.loc[cols]
    preds = []
    m = max(5, int(len(avail) * BOOT_FRAC))
    for _ in range(BOOT_B):
        gsel = list(rng.choice(avail, size=m, replace=False))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = PoissonGLMSolver(genes=gsel).solve(c_sub, ref).proportions.reindex(columns=types).fillna(0.0)
        preds.append(p.reindex(index=cols).to_numpy(float))
    stack = np.stack(preds, axis=0)                       # (B, S, K)
    lo = np.percentile(stack, (1 - NOMINAL_COVERAGE) / 2 * 100, axis=0)
    hi = np.percentile(stack, (1 + NOMINAL_COVERAGE) / 2 * 100, axis=0)
    ct = cert.per_type.set_index("cell_type")
    by_class = defaultdict(list)
    Tv = tr.reindex(columns=types).fillna(0).to_numpy(float)
    for j, k in enumerate(types):
        cls = ct.loc[k, "recoverability"]
        inside = (Tv[:, j] >= lo[:, j] - 1e-9) & (Tv[:, j] <= hi[:, j] + 1e-9)
        by_class[cls].extend(inside.tolist())
    cov = {c: round(float(np.mean(v)), 4) for c, v in by_class.items() if v}
    overall = round(float(np.mean([x for v in by_class.values() for x in v])), 4)
    return cov, overall


def _pool(query, ctc, seeds, n_per_regime, H):
    counts_l, truth_l = [], []
    for seed in seeds:
        c, t, _ = H.generate_pseudobulk(query, ctc, n_per_regime=n_per_regime, n_cells=200, seed=seed)
        counts_l.append(c); truth_l.append(t)
    counts = pd.concat(counts_l, axis=1)
    counts.columns = [f"m{i}" for i in range(counts.shape[1])]
    truth = pd.concat(truth_l, axis=0).reset_index(drop=True)
    truth.index = counts.columns
    return counts, truth


def run(dataset, seeds, outdir, hashes):
    import _harness as H
    from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset, LUNG_H5AD
    from benchmarks.signatures import validation as V
    from tissueresolve.solver import PoissonGLMSolver
    h5ad_path = {"breast": Path(H.REFERENCE_H5AD), "lung": Path(LUNG_H5AD)}.get(dataset)

    D = _load_dataset(dataset)
    ad, ctc, dc, ref, mapping = D["adata"], D["celltype_col"], D["donor_col"], D["ref"], D["mapping"]
    types = [str(c) for c in ref.cell_types]
    train = set(map(str, D["ref_donors"]))

    # ---- donor-disjoint 3-way split: TRAIN (reference) / CALIBRATION / TEST ----
    q = sorted(map(str, D["query_donors"]))
    perm = np.random.default_rng(2024).permutation(q)
    half = max(1, len(perm) // 2)
    cal_donors, test_donors = sorted(perm[:half]), sorted(perm[half:])
    assert train.isdisjoint(cal_donors) and train.isdisjoint(test_donors)
    assert set(cal_donors).isdisjoint(test_donors)
    V.validate_donor_disjoint(train, set(cal_donors) | set(test_donors))
    V.validate_selection_donors(train, set(cal_donors) | set(test_donors))

    R_cpm = ref.as_R_cpm()
    genes = list(ref.gene_names)
    q_cal = ad[ad.obs[dc].astype(str).isin(set(cal_donors)).to_numpy()].copy()
    q_test = ad[ad.obs[dc].astype(str).isin(set(test_donors)).to_numpy()].copy()
    n_don = {t_: int(q_test.obs.loc[q_test.obs[ctc].astype(str) == t_, dc].nunique()) for t_ in types}

    counts_test_full, truth_test = _pool(q_test, ctc, seeds, 3, H)
    counts_cal_full, truth_cal = _pool(q_cal, ctc, seeds[:3], 2, H)
    truth_test = truth_test.reindex(columns=types).fillna(0.0)
    truth_cal = truth_cal.reindex(columns=types).fillna(0.0)

    panel_1000 = list(np.random.default_rng(12345).choice(genes, size=min(1000, len(genes)), replace=False))
    stressed = STRESSED[dataset]
    cal_conditions = {"bulk_full", stressed}

    test_rows, cal_rows, cond_rows, group_rows = [], [], [], []
    cert_saved = False
    cert_stressed = None
    for cname, depth, n_panel in CONDITIONS:
        rng = np.random.default_rng(abs(hash((dataset, cname))) % (2**32))
        panel = panel_1000 if n_panel else genes
        counts_t = _thin(counts_test_full, depth, rng)
        det_frac = (counts_t.loc[panel] > 0).mean(axis=1)
        detectable = [g for g in panel if det_frac.get(g, 0) >= 0.10]
        if len(detectable) < 2:
            continue
        lib = float(counts_t.loc[detectable].sum(axis=0).median())
        cert = identifiability_certificate(R_cpm, types, genes, query_detectable_genes=detectable,
                                           library_size=lib, n_donors_by_type=n_don)
        if cname == "bulk_full" and not cert_saved:
            cert.save(outdir); cert_saved = True
        if cname == stressed:
            cert_stressed = (cert, counts_t, panel)
        ct = cert.per_type.set_index("cell_type")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pred_t = PoissonGLMSolver(genes=(panel if n_panel else None)).solve(
                counts_t, ref).proportions.reindex(columns=types).fillna(0.0)
        V.validate_truth_pred_alignment(truth_test, pred_t)
        per_t = _per_type_metrics(truth_test, pred_t, mapping, types)
        for k in types:
            test_rows.append(dict(dataset=dataset, condition=cname, cell_type=k,
                                  recoverability=ct.loc[k, "recoverability"],
                                  rec_fraction=ct.loc[k, "rec_fraction"],
                                  theoretical_detection_floor=ct.loc[k, "theoretical_detection_floor"],
                                  n_donors=ct.loc[k, "n_donors"], **per_t[k]))
        cc = pd.Series([ct.loc[k, "recoverability"] for k in types]).value_counts()
        cond_rows.append(dict(dataset=dataset, condition=cname, depth_factor=depth,
                              n_detectable=len(detectable), library_size=round(lib),
                              condition_number=round(cert.condition_number, 1),
                              noise_adjusted_rank=cert.noise_adjusted_rank,
                              n_confounded_groups=len(cert.confounded_groups),
                              class_counts=cc.to_dict()))
        for g in _group_aggregate_recovery(cert, truth_test, pred_t, types):
            group_rows.append(dict(dataset=dataset, condition=cname, **g))

        if cname in cal_conditions:
            counts_c = _thin(counts_cal_full, depth, np.random.default_rng(7))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pred_c = PoissonGLMSolver(genes=(panel if n_panel else None)).solve(
                    counts_c, ref).proportions.reindex(columns=types).fillna(0.0)
            per_c = _per_type_metrics(truth_cal, pred_c, mapping, types)
            for k in types:
                cal_rows.append(dict(cell_type=k, recoverability=ct.loc[k, "recoverability"],
                                     held_out_rmse=per_c[k]["held_out_rmse"]))

    test_metrics = pd.DataFrame(test_rows)
    cal_metrics = pd.DataFrame(cal_rows)
    test_metrics.to_csv(outdir / "metrics_by_celltype.tsv", sep="\t", index=False)
    pd.DataFrame(cond_rows).to_csv(outdir / "condition_summary.tsv", sep="\t", index=False)
    pd.DataFrame(group_rows or [{"group": None}]).to_csv(
        outdir / "group_aggregate_recovery.tsv", sep="\t", index=False)

    by_class = test_metrics[test_metrics.recoverability != "NOT_TESTABLE"].groupby("recoverability")[
        ["held_out_rmse", "conditional_rmse", "spillover_when_absent", "detection_recall",
         "rare_recall", "false_positive_rate", "rec_fraction"]].mean().round(4)
    by_class["n_types"] = test_metrics.groupby("recoverability").size()
    by_class.to_csv(outdir / "metrics_by_class.tsv", sep="\t")

    class_prob = _fit_class_prob(cal_metrics)             # fit on CALIBRATION
    calib = _calibration(test_metrics, class_prob)        # evaluate on TEST
    frates = _false_rates(test_metrics)
    pd.DataFrame([{**class_prob}]).to_csv(outdir / "calibrated_class_probabilities.tsv", sep="\t", index=False)
    pd.DataFrame([{**calib}]).to_csv(outdir / "calibration_metrics.tsv", sep="\t", index=False)
    pd.DataFrame([{**frates}]).to_csv(outdir / "false_merge_false_resolution.tsv", sep="\t", index=False)

    # uncertainty coverage (gene-subsampling bootstrap) on the stressed condition, TEST samples
    cov_by_class, cov_overall = {}, np.nan
    cov_cal = np.nan
    if cert_stressed is not None:
        cert_s, counts_s, panel_s = cert_stressed
        cov_by_class, cov_overall = _bootstrap_coverage(
            counts_s, ref, truth_test, cert_s, types, panel_s if panel_s is not genes else None,
            np.random.default_rng(99), PoissonGLMSolver)
        # nominal-level check on CALIBRATION donors (verification, not tuning)
        depth_s = dict((c, d) for c, d, _ in CONDITIONS)[stressed]
        counts_cs = _thin(counts_cal_full, depth_s, np.random.default_rng(7))
        _, cov_cal = _bootstrap_coverage(counts_cs, ref, truth_cal, cert_s, types,
                                         panel_s if panel_s is not genes else None,
                                         np.random.default_rng(98), PoissonGLMSolver)
    pd.DataFrame([{"class": c, "coverage": v} for c, v in cov_by_class.items()]
                 + [{"class": "OVERALL_test", "coverage": cov_overall},
                    {"class": "OVERALL_calibration", "coverage": cov_cal}]).to_csv(
        outdir / "uncertainty_coverage.tsv", sep="\t", index=False)

    pd.DataFrame([{"split": "train", "n_donors": len(train), "donors": ";".join(sorted(train))},
                  {"split": "calibration", "n_donors": len(cal_donors), "donors": ";".join(cal_donors)},
                  {"split": "test", "n_donors": len(test_donors), "donors": ";".join(test_donors)}]).to_csv(
        outdir / "train_calibration_test_donors.tsv", sep="\t", index=False)
    hashes.append({"dataset": dataset, "h5ad": str(h5ad_path),
                   "md5": _md5(h5ad_path) if h5ad_path and h5ad_path.exists() else "n/a"})

    print(f"\n=== {dataset}: class distribution per condition ===")
    for r in cond_rows:
        print(f"  {r['condition']:13s} depth×{r['depth_factor']:<6} cond#={r['condition_number']:>7} "
              f"nrank={r['noise_adjusted_rank']:>3} groups={r['n_confounded_groups']:>2}  {r['class_counts']}")
    print(f"--- {dataset}: held-out metrics by class (TEST, pooled over conditions) ---\n", by_class.to_string())
    print(f"class_prob (fit on CAL): { {k: round(v,3) for k,v in class_prob.items()} }")
    print(f"calibration (TEST): {calib}\nfalse rates (TEST): {frates}")
    print(f"uncertainty coverage @ nominal {NOMINAL_COVERAGE:.0%} (stressed={stressed}): "
          f"by_class={cov_by_class} overall_test={cov_overall} overall_cal={cov_cal}")
    return dict(dataset=dataset, **calib, **frates, uncertainty_overall_test=cov_overall,
                class_counts=test_metrics.recoverability.value_counts().to_dict())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    import _harness as H
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    versions = {m: __import__(m).__version__ for m in ("numpy", "pandas", "scipy", "anndata")}
    versions["python"] = platform.python_version()
    t0 = time.perf_counter()
    summary, hashes = [], []
    for ds in args.datasets:
        od = OUT / ds; od.mkdir(parents=True, exist_ok=True)
        summary.append(run(ds, args.seeds, od, hashes))
    pd.DataFrame([{"package": k, "version": v} for k, v in versions.items()]).to_csv(
        OUT / "software_versions.tsv", sep="\t", index=False)
    pd.DataFrame(hashes).to_csv(OUT / "input_hashes.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(OUT / "stage2a_certificate_summary.tsv", sep="\t", index=False)
    (OUT / "benchmark_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "rec_rmse_thresh": REC_RMSE_THRESH,
         "detect_thresh": DETECT_THRESH, "rare_band": RARE_BAND, "nominal_coverage": NOMINAL_COVERAGE,
         "boot_B": BOOT_B, "class_prob_prior": CLASS_PROB_PRIOR, "conditions": [c[0] for c in CONDITIONS],
         "runtime_s": round(time.perf_counter() - t0, 1)}, indent=2, default=str))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
