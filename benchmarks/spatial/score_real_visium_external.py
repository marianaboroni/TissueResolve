#!/usr/bin/env python
"""Add RCTD + cell2location to the real-Visium no-GT evidence + cross-method concordance.

Reloads the same deterministic 500-spot scenario, scores the external methods'
evidence panel (marker recovery, Moran's I, etc.), and computes FULL pairwise
cross-method concordance (fine + family) across all executed methods.  Concordance
is method agreement, NOT accuracy.

Usage:  python benchmarks/spatial/score_real_visium_external.py --run-real-data
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "benchmarks" / "spatial"))
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
import _harness as H  # noqa: E402
from benchmarks.shared import metrics as M  # noqa: E402
import run_real_visium_evidence as RV  # noqa: E402

OUT = REPO / "benchmarks" / "outputs" / "spatial_real_visium"
EI = OUT / "external_inputs"
RAW = OUT / "raw"
INTERNAL = ["NNLS_per_spot", "TissueResolve_spatial_flat", "TissueResolve_spatial_hierarchical"]
EXTERNAL = ["RCTD", "cell2location"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data.", file=sys.stderr); return 2

    sc = RV._load_scenario()
    ref, mapping = sc["ref"], sc["mapping"]
    ref_types = list(ref.cell_types)
    marker_sets = RV._marker_sets(ref, sc["gene_names"])
    all_markers = sorted({g for v in marker_sets.values() for g in v})
    gidx = {g: i for i, g in enumerate(sc["gene_names"])}
    Ym = sc["Y"][:, [gidx[g] for g in all_markers]]
    Ym_cpm = np.log1p(Ym / np.clip(Ym.sum(1, keepdims=True), 1, None) * 1e6)
    marker_expr = pd.DataFrame(Ym_cpm, index=sc["spot_ids"], columns=all_markers)

    preds, ev_rows, status = {}, [], []
    # internal (re-score for a single combined table)
    for m in INTERNAL:
        pf = RAW / f"{m}_proportions.tsv"
        if not pf.exists():
            continue
        p = pd.read_csv(pf, sep="\t", index_col=0).reindex(index=sc["spot_ids"])
        preds[m] = p
        ev_rows.append({"method": m, "family": "TissueResolve" if "TissueResolve" in m else "baseline",
                        **RV._evidence(p, sc, marker_expr, marker_sets)})
    # external
    for m in EXTERNAL:
        sfile = EI / f"{m}_status.tsv"
        st = pd.read_csv(sfile, sep="\t").iloc[0].to_dict() if sfile.exists() else {"status": "not_run"}
        status.append({"method": m, **{k: st.get(k) for k in ("status", "version", "runtime_seconds")}})
        pf = EI / f"{m}_pred.tsv"
        if st.get("status") != "executed" or not pf.exists():
            print(f"{m}: {st.get('status')} — not scored"); continue
        p = pd.read_csv(pf, sep="\t", index_col=0)
        p.index = [str(i) for i in p.index]
        p = p.reindex(index=sc["spot_ids"]).fillna(0.0)
        preds[m] = p
        ev_rows.append({"method": m, "family": "external",
                        **RV._evidence(p, sc, marker_expr, marker_sets)})
        print(f"{m}: scored marker_rec={ev_rows[-1]['marker_recovery_mean']:.3f}")

    # full pairwise concordance (fine + family)
    fine = {k: v.reindex(columns=ref_types).fillna(0.0) for k, v in preds.items()}
    fam = {k: M.aggregate_to_families(v, mapping) for k, v in preds.items()}
    conc = M.pairwise_method_correlation(fine)
    conc_fam = M.pairwise_method_correlation(fam)

    ev = pd.DataFrame(ev_rows)
    with open(OUT / "spatial_real_metrics_all.tsv", "w") as fh:
        fh.write("# Real Visium evidence (NO GROUND TRUTH), internal + external methods. "
                 "marker recovery=proxy; concordance!=accuracy.\n")
        ev.to_csv(fh, sep="\t", index=False)
    pd.DataFrame([{"pair": k, "fine_pearson": v} for k, v in conc.items()]).to_csv(
        OUT / "method_concordance_fine_all.tsv", sep="\t", index=False)
    pd.DataFrame([{"pair": k, "family_pearson": v} for k, v in conc_fam.items()]).to_csv(
        OUT / "method_concordance_family_all.tsv", sep="\t", index=False)
    pd.DataFrame(status).to_csv(OUT / "external_real_status.tsv", sep="\t", index=False)

    # concordance heatmap (fine) across all methods
    names = list(preds)
    Cm = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for k, v in conc.items():
        a, b = k.split("__vs__"); Cm.loc[a, b] = v; Cm.loc[b, a] = v
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(Cm.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{Cm.iloc[i,j]:.2f}", ha="center", va="center", fontsize=7,
                    color="black" if abs(Cm.iloc[i,j]) < 0.6 else "white")
    fig.colorbar(im, ax=ax, label="fine Pearson (method agreement)")
    ax.set_title("Real Visium cross-method concordance (NOT accuracy)")
    fig.tight_layout()
    figd = OUT / "figures"; figd.mkdir(exist_ok=True)
    fig.savefig(figd / "figM_real_concordance.png", dpi=150)
    fig.savefig(figd / "figM_real_concordance.svg"); plt.close(fig)
    Cm.to_csv(figd / "figM_real_concordance.data.tsv", sep="\t")
    (figd / "figM_real_concordance.caption.txt").write_text(
        "Pairwise fine-level Pearson between methods on the real Visium section "
        "(500 spots). Method AGREEMENT, not accuracy. Includes established tools "
        "RCTD and cell2location.\n", encoding="utf-8")

    print("\n=== evidence (all methods, NO GROUND TRUTH) ===")
    print(ev[["method", "marker_recovery_mean", "mean_morans_i", "mean_entropy",
              "near_zero_fraction", "recon_pearson"]].round(3).to_string(index=False))
    print("\n=== fine concordance (agreement, NOT accuracy) ===")
    for k, v in sorted(conc.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:.3f}")
    print(f"\nWrote -> {OUT}/ (figM, *_all.tsv)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
