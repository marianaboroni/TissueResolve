#!/usr/bin/env python
"""Pilot: does DISTRIBUTION-AWARE reference modeling (low-rank covariance /
donor×state metacells) help separate COLLINEAR T/NK states, where mean profiles
(BC>0.98) say they are non-identifiable?

Two probes on the breast reference (donor-held-out):

Probe A — headroom (cell level): for each collinear T/NK pair, donor-held-out
classification AUROC using only the MEAN (NearestCentroid) vs distribution/
covariance-aware models (LDA shared-cov, QDA in low-rank PC space = "low-rank
covariance", Logistic). If covariance/distribution-aware >> mean-only on HELD-OUT
donors, the mean discards real separating signal → a distribution-aware deconvolution
has headroom.

Probe B — metacell deconvolution (bulk level): T/NK-restricted pseudobulk from
held-out query donors; deconvolve (NNLS) against (i) one MEAN profile per state vs
(ii) a donor×state METACELL dictionary (aggregated back to states). Does the metacell
dictionary improve the conditional within-T/NK split vs the mean profile?

Experimental / opt-in; no defaults changed; outputs gitignored.

Usage:
  PYTHONPATH=src:. python benchmarks/dev/distribution_aware_tnk_pilot.py --run-real-data
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

OUT = REPO / "benchmarks" / "outputs" / "distribution_aware_tnk"
N_HVG = 800
N_PC = 30
N_TOP_PAIRS = 5
MIN_CELLS_PER_CLASS = 60
MIN_CELLS_PER_METACELL = 15


def _lognorm(X):
    lib = np.maximum(X.sum(1, keepdims=True), 1.0)
    return np.log1p(X / lib * 1e4)


def _dense(adata):
    import scipy.sparse as sp
    X = adata.X
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2

    from sklearn.decomposition import PCA
    from sklearn.discriminant_analysis import (
        LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA)
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import NearestCentroid
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import roc_auc_score
    from scipy.optimize import nnls

    D = _load_dataset("breast")
    adata, mapping, ctc, dc = D["adata"], D["mapping"], D["celltype_col"], D["donor_col"]
    tnk = [t for t, f in mapping.items() if f == "T/NK"]
    obs = adata.obs
    ct = obs[ctc].astype(str).to_numpy()
    dn = obs[dc].astype(str).to_numpy()
    tnk_mask = np.isin(ct, tnk)

    # HVG within T/NK (variance of log1p-CPM)
    Xt = _dense(adata[tnk_mask])
    Xln = _lognorm(Xt)
    hvg = np.argsort(Xln.var(0))[::-1][:N_HVG]
    genes = np.array(adata.var_names)[hvg]

    # collinear T/NK pairs by mean BC
    ref = D["ref"]; cts = [str(c) for c in ref.cell_types]
    R = ref.as_R_cpm(); Pn = R / (R.sum(1, keepdims=True) + 1e-10)
    gi = {c: i for i, c in enumerate(cts)}
    present = [t for t in tnk if t in gi]
    counts = {t: int((ct == t).sum()) for t in present}
    pairs = []
    for i in range(len(present)):
        for j in range(i + 1, len(present)):
            a, b = present[i], present[j]
            if counts[a] < MIN_CELLS_PER_CLASS or counts[b] < MIN_CELLS_PER_CLASS:
                continue
            bc = float(np.sum(np.sqrt(Pn[gi[a]] * Pn[gi[b]])))
            pairs.append((bc, a, b))
    pairs.sort(reverse=True)
    pairs = pairs[:N_TOP_PAIRS]
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------- Probe A: headroom (donor-held-out cell classification) ----------
    rng = np.random.default_rng(args.seed)
    rowsA = []
    for bc, a, b in pairs:
        ma = (ct == a) & tnk_mask if False else (ct == a)
        mb = (ct == b)
        idx = np.where(ma | mb)[0]
        X = _lognorm(_dense(adata[idx]))[:, hvg]
        y = (ct[idx] == b).astype(int)
        groups = dn[idx]
        ndon = len(np.unique(groups))
        nsplits = min(5, ndon)
        if nsplits < 2:
            continue
        gkf = GroupKFold(n_splits=nsplits)
        aucs = {"mean_centroid": [], "lda_sharedcov": [], "qda_lowrank_cov": [], "logistic": []}
        for tr_i, te_i in gkf.split(X, y, groups):
            if len(np.unique(y[tr_i])) < 2 or len(np.unique(y[te_i])) < 2:
                continue
            pca = PCA(n_components=min(N_PC, len(tr_i) - 1, X.shape[1])).fit(X[tr_i])
            Ztr, Zte = pca.transform(X[tr_i]), pca.transform(X[te_i])
            # mean-only: nearest centroid → score = signed distance margin
            nc = NearestCentroid().fit(Ztr, y[tr_i])
            d0 = np.linalg.norm(Zte - nc.centroids_[0], axis=1)
            d1 = np.linalg.norm(Zte - nc.centroids_[1], axis=1)
            aucs["mean_centroid"].append(roc_auc_score(y[te_i], d0 - d1))
            try:
                lda = LDA().fit(Ztr, y[tr_i])
                aucs["lda_sharedcov"].append(roc_auc_score(y[te_i], lda.predict_proba(Zte)[:, 1]))
            except Exception: pass
            try:
                qda = QDA(reg_param=0.1).fit(Ztr, y[tr_i])
                aucs["qda_lowrank_cov"].append(roc_auc_score(y[te_i], qda.predict_proba(Zte)[:, 1]))
            except Exception: pass
            try:
                lr = LogisticRegression(max_iter=1000, C=1.0).fit(Ztr, y[tr_i])
                aucs["logistic"].append(roc_auc_score(y[te_i], lr.predict_proba(Zte)[:, 1]))
            except Exception: pass
        rowsA.append({"pair": f"{a} || {b}", "BC_mean": round(bc, 4),
                      "n_a": counts[a], "n_b": counts[b], "n_donors": int(ndon),
                      **{k: round(float(np.mean(v)), 3) if v else float("nan")
                         for k, v in aucs.items()}})
        r = rowsA[-1]
        print(f"  BC={bc:.3f} {a[:22]:<22}|{b[:22]:<22} centroid={r['mean_centroid']:.3f} "
              f"lda={r['lda_sharedcov']:.3f} qda_lowrank={r['qda_lowrank_cov']:.3f} logit={r['logistic']:.3f}")
    dfA = pd.DataFrame(rowsA)
    dfA.to_csv(OUT / "probeA_headroom_auroc.tsv", sep="\t", index=False)

    # ---------- Probe B: metacell vs mean NNLS on T/NK-restricted pseudobulk ----------
    ref_d = set(D["ref_donors"]); qry_d = set(D["query_donors"])
    usable = [t for t in present if counts[t] >= MIN_CELLS_PER_CLASS]
    # build reference cell pools (ref donors) and query pools (query donors), T/NK genes
    refmask = np.isin(dn, list(ref_d))
    qmask = np.isin(dn, list(qry_d))
    G = hvg
    def state_cells(mask, state):
        idx = np.where(mask & (ct == state))[0]
        return idx
    # mean profile per state (ref donors), CPM on HVG
    def cpm_rows(idx):
        X = _dense(adata[idx])[:, G].astype(float)
        return X  # raw counts on HVG
    mean_ref = {}
    metacell_ref = {}   # state -> list of (donor mean profile)
    for s in usable:
        ridx = state_cells(refmask, s)
        if len(ridx) < MIN_CELLS_PER_CLASS:
            continue
        Xs = cpm_rows(ridx)
        prof = Xs.sum(0); prof = prof / (prof.sum() + 1e-9)
        mean_ref[s] = prof
        # donor metacells
        dons = dn[ridx]
        mcs = []
        for d in np.unique(dons):
            sub = Xs[dons == d]
            if sub.shape[0] >= MIN_CELLS_PER_METACELL:
                p = sub.sum(0); p = p / (p.sum() + 1e-9); mcs.append(p)
        if len(mcs) >= 2:
            metacell_ref[s] = np.array(mcs)
    states = [s for s in usable if s in mean_ref and s in metacell_ref]
    rowsB = []
    if len(states) >= 3:
        Phi_mean = np.array([mean_ref[s] for s in states]).T           # (G, S)
        mc_cols, mc_state = [], []
        for s in states:
            for p in metacell_ref[s]:
                mc_cols.append(p); mc_state.append(s)
        Phi_mc = np.array(mc_cols).T                                    # (G, M)
        mc_state = np.array(mc_state)
        # generate T/NK-only pseudobulk from query donors
        n_samples = 16
        for k in range(n_samples):
            true = rng.dirichlet(np.ones(len(states)))
            # realize counts: sample cells per state from query pool
            bulk = np.zeros(len(G))
            for si, s in enumerate(states):
                qidx = state_cells(qmask, s)
                if len(qidx) == 0:
                    continue
                ncell = max(1, int(round(true[si] * 400)))
                pick = rng.choice(qidx, ncell, replace=len(qidx) < ncell)
                bulk += _dense(adata[pick])[:, G].astype(float).sum(0)
            bn = bulk / (bulk.sum() + 1e-9)
            # mean NNLS
            wm, _ = nnls(Phi_mean, bn); wm = wm / (wm.sum() + 1e-9)
            # metacell NNLS → aggregate to states
            wc, _ = nnls(Phi_mc, bn)
            agg = np.array([wc[mc_state == s].sum() for s in states]); agg = agg / (agg.sum() + 1e-9)
            rowsB.append({"sample": k, "metric": "mean", **{s: wm[i] for i, s in enumerate(states)},
                          "_true": json.dumps({s: float(true[i]) for i, s in enumerate(states)})})
            rowsB.append({"sample": k, "metric": "metacell", **{s: agg[i] for i, s in enumerate(states)},
                          "_true": json.dumps({s: float(true[i]) for i, s in enumerate(states)})})
        dfB = pd.DataFrame(rowsB)
        dfB.to_csv(OUT / "probeB_metacell_vs_mean_preds.tsv", sep="\t", index=False)
        # score conditional RMSE (within T/NK) per method
        def rmse(method):
            sub = dfB[dfB.metric == method]
            errs = []
            for _, r in sub.iterrows():
                tr = json.loads(r["_true"])
                p = np.array([r[s] for s in states]); t = np.array([tr[s] for s in states])
                errs.append((p - t) ** 2)
            return float(np.sqrt(np.mean(np.concatenate([e[None, :] for e in errs]))))
        mean_rmse, mc_rmse = rmse("mean"), rmse("metacell")
        # correlation to truth
        def corr(method):
            sub = dfB[dfB.metric == method]; P, T = [], []
            for _, r in sub.iterrows():
                tr = json.loads(r["_true"]); P += [r[s] for s in states]; T += [tr[s] for s in states]
            return float(np.corrcoef(P, T)[0, 1])
        probeB = {"n_states": len(states), "n_metacells": int(Phi_mc.shape[1]),
                  "mean_rmse": round(mean_rmse, 4), "metacell_rmse": round(mc_rmse, 4),
                  "mean_pearson": round(corr("mean"), 3), "metacell_pearson": round(corr("metacell"), 3)}
    else:
        probeB = {"status": "insufficient states/metacells for Probe B"}
    (OUT / "probeB_summary.json").write_text(json.dumps(probeB, indent=2), encoding="utf-8")

    (OUT / "manifest.json").write_text(json.dumps({
        "dataset": "breast", "family": "T/NK", "n_hvg": N_HVG, "n_pc": N_PC,
        "pairs": [f"{a} || {b}" for _, a, b in pairs], "seed": args.seed,
        "note": "pilot; experimental; no defaults changed",
    }, indent=2), encoding="utf-8")

    print("\n=== Probe A (donor-held-out AUROC; mean vs covariance/distribution-aware) ===")
    if not dfA.empty:
        print(dfA[["pair", "BC_mean", "mean_centroid", "lda_sharedcov",
                   "qda_lowrank_cov", "logistic"]].to_string(index=False))
    print("\n=== Probe B (metacell vs mean NNLS on T/NK pseudobulk) ===")
    print(json.dumps(probeB, indent=2))
    print(f"\nWrote -> {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
