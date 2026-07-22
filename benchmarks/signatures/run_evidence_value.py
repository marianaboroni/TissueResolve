#!/usr/bin/env python
"""Strategy E core (Stage 1.5): do sibling / rare_confirmation panels add information
BEYOND the fine_global abundance estimate? Leakage-safe, donor-disjoint, real-data.

For each (mixture, subtype in a multi-member family) we collect:
  abundance      = fine_global Poisson estimate
  sibling_score  = sibling_evidence (within-family target-vs-siblings signal)
  rare_score     = rare_confirmation_evidence (fraction of confirmation genes observed)
  truth, is_false_positive (pred>thr & true~0), abs_error

Then, splitting MIXTURES into dev/test (never using test truth to fit anything), we ask:
  (1) FP detection: does [abundance + evidence] beat [abundance] at flagging false positives
      (test ROC-AUC)?
  (2) error prediction: does [abundance + evidence] beat [abundance] at predicting |error|
      (test R2)?
A positive, replicated Δ on held-out mixtures ⇒ the extra panels add value (Strategy E works).

Usage:
  PYTHONPATH=src:. python benchmarks/signatures/run_evidence_value.py \
      --run-real-data --datasets breast lung --seeds 0 1 2 3 4
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.spatial.run_weak_smoothing_grid import _load_dataset  # noqa: E402

OUT = REPO / "benchmarks" / "results" / "signatures"
PRESENCE_THR = 0.01           # pre-specified (not tuned on test)
DEV_SEEDS = 3                 # first N seeds = dev, rest = test


def _collect(dataset, seeds):
    from tissueresolve.solver import PoissonGLMSolver
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    from tissueresolve.reference.signature_evidence import (
        rare_confirmation_evidence, sibling_evidence)
    D = _load_dataset(dataset)
    ad, ctc, dc, ref, types = (D["adata"], D["celltype_col"], D["donor_col"],
                               D["ref"], D["ref_types"])
    mp = {t: D["mapping"].get(t, t) for t in types}
    train = ad[ad.obs[dc].astype(str).isin(set(D["ref_donors"])).to_numpy()].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ReferenceSignatureOptimizer(ctc, dc, mapping=mp, seed=0).optimize(train)
    subs_with_evidence = set(m.sibling_genes_by_subtype)   # subtypes that HAVE evidence panels
    rows = []
    for seed in seeds:
        counts, truth, _ = H.generate_pseudobulk(
            ad[ad.obs[dc].astype(str).isin(set(D["query_donors"])).to_numpy()].copy(),
            ctc, n_per_regime=2, n_cells=200, seed=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            props = PoissonGLMSolver(genes=m.fine_global_genes).solve(counts, ref).proportions
            sib = sibling_evidence(counts, m.sibling_genes_by_subtype, mp)
            rare = rare_confirmation_evidence(counts, m.rare_confirmation)
        for s in truth.index:
            for sub in subs_with_evidence:
                if sub not in truth.columns:
                    continue
                ab = float(props.loc[s, sub]) if sub in props.columns else 0.0
                tv = float(truth.loc[s, sub])
                rows.append(dict(
                    dataset=dataset, seed=seed, mixture=f"{seed}:{s}", subtype=sub,
                    abundance=ab, sibling=float(sib.loc[s, sub]) if sub in sib.columns else np.nan,
                    rare=float(rare.loc[s, sub]) if sub in rare.columns else np.nan,
                    truth=tv, present=int(tv > 1e-9),
                    is_fp=int(ab > PRESENCE_THR and tv < 1e-9),
                    abs_error=abs(ab - tv)))
    return pd.DataFrame(rows)


def _fit_eval(df):
    """Leakage-safe: fit on dev seeds, evaluate on test seeds. Returns metric dict."""
    from sklearn.linear_model import LogisticRegression, LinearRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score, r2_score
    d = df.dropna(subset=["sibling", "rare"]).copy()
    dev = d[d.seed < DEV_SEEDS]; test = d[d.seed >= DEV_SEEDS]
    res = {"n_dev": len(dev), "n_test": len(test)}
    if len(dev) < 30 or len(test) < 30:
        return {**res, "note": "insufficient rows"}

    # (1) FP detection among predicted-present entries
    def _fp_auc(cols):
        dv = dev[dev.abundance > PRESENCE_THR]; te = test[test.abundance > PRESENCE_THR]
        if dv.is_fp.nunique() < 2 or te.is_fp.nunique() < 2 or len(dv) < 20:
            return np.nan
        sc = StandardScaler().fit(dv[cols])
        clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(
            sc.transform(dv[cols]), dv.is_fp)
        return float(roc_auc_score(te.is_fp, clf.predict_proba(sc.transform(te[cols]))[:, 1]))
    res["fp_auc_abundance"] = _fp_auc(["abundance"])
    res["fp_auc_abund+evid"] = _fp_auc(["abundance", "sibling", "rare"])

    # (2) |error| prediction (all entries)
    def _err_r2(cols):
        sc = StandardScaler().fit(dev[cols])
        reg = LinearRegression().fit(sc.transform(dev[cols]), dev.abs_error)
        return float(r2_score(test.abs_error, reg.predict(sc.transform(test[cols]))))
    res["err_r2_abundance"] = _err_r2(["abundance"])
    res["err_r2_abund+evid"] = _err_r2(["abundance", "sibling", "rare"])
    res["delta_fp_auc"] = res["fp_auc_abund+evid"] - res["fp_auc_abundance"]
    res["delta_err_r2"] = res["err_r2_abund+evid"] - res["err_r2_abundance"]
    # partial association of each evidence with error, controlling abundance (test set)
    for name in ("sibling", "rare"):
        r_raw = float(np.corrcoef(test[name], test.abs_error)[0, 1])
        res[f"corr_{name}_abserr"] = r_raw
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--datasets", nargs="+", default=["breast", "lung"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows, summary = [], []
    for ds in args.datasets:
        df = _collect(ds, args.seeds)
        all_rows.append(df)
        r = _fit_eval(df); r["dataset"] = ds; summary.append(r)
        print(f"\n[{ds}] {r}", flush=True)
    pd.concat(all_rows, ignore_index=True).to_csv(OUT / "evidence_value_rows.tsv", sep="\t", index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(OUT / "evidence_value_summary.tsv", sep="\t", index=False)
    (OUT / "evidence_value_manifest.json").write_text(json.dumps(
        {"datasets": args.datasets, "seeds": args.seeds, "dev_seeds": DEV_SEEDS,
         "presence_thr": PRESENCE_THR,
         "question": "does sibling/rare evidence add info beyond fine_global abundance"}, indent=2))
    print("\n=== EVIDENCE VALUE SUMMARY (test-set; Δ>0 ⇒ panels add info) ===")
    cols = [c for c in ["dataset", "fp_auc_abundance", "fp_auc_abund+evid", "delta_fp_auc",
                        "err_r2_abundance", "err_r2_abund+evid", "delta_err_r2",
                        "corr_sibling_abserr", "corr_rare_abserr"] if c in sm.columns]
    print(sm[cols].round(4).to_string(index=False))
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
