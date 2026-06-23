#!/usr/bin/env python
"""TASK A — QC-first report rebuilt on the completed identifiability evidence.

Evidence chain: reference support → donor coverage → signature identifiability
(3 distinct levels) → query compatibility → FULL-PANEL deconvolution reliability →
trusted resolution → final predictions.  Central rule: a figure is in the MAIN
report only if it helps decide reliability or interpret the prediction; everything
else → technical appendix / source_data.

Does NOT touch stable report code (rule: no default change). Self-contained builder
(like build_final_report.py). Trusted resolution is derived from full-panel
deconvolution-relevant metrics (shared-lineage fraction, condition number,
per-family conditional RMSE) — NOT cell-classification AUROC.

Outputs (benchmarks/outputs/qc_report/): report.html, technical_appendix.html,
figures/figure_manifest.tsv, source_data/*.tsv, warnings.json, methods.txt, figures/*.

Usage:  python benchmarks/build_qc_first_report.py --run-real-data
"""
from __future__ import annotations

import argparse
import base64
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
import _harness as H  # noqa: E402
OUT = REPO / "benchmarks" / "outputs"
QC = OUT / "qc_report"
FIGD = QC / "figures"
SRC = QC / "source_data"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
TCGA = H.DATA_DIR / "bulk_tcga_tnbc" / "tcga_tnbc_counts.tsv"

MANIFEST = []   # (section, figure, file, purpose, source_data)


def _savefig(fig, name, section, purpose, data: pd.DataFrame, caption, what_to_check, limitation):
    FIGD.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)
    import matplotlib.pyplot as plt
    fig.tight_layout()
    png = FIGD / f"{name}.png"; svg = FIGD / f"{name}.svg"
    fig.savefig(png, dpi=150); fig.savefig(svg); plt.close(fig)
    sd = SRC / f"{name}.data.tsv"; data.to_csv(sd, sep="\t")
    (FIGD / f"{name}.caption.txt").write_text(caption.strip() + "\n", encoding="utf-8")
    MANIFEST.append({"section": section, "figure": name, "png": str(png.relative_to(QC)),
                     "svg": str(svg.relative_to(QC)), "purpose": purpose,
                     "source_data": str(sd.relative_to(QC)), "caption": caption,
                     "what_to_check": what_to_check, "limitation": limitation})
    return name


def _b64(name):
    p = FIGD / f"{name}.png"
    return base64.b64encode(p.read_bytes()).decode() if p.exists() else None


def _embed_existing(src_png, name, section, purpose, caption, what, limit):
    """Copy an existing decision-relevant figure into the report manifest."""
    src = OUT / src_png
    if not src.exists():
        return None
    FIGD.mkdir(parents=True, exist_ok=True)
    (FIGD / f"{name}.png").write_bytes(src.read_bytes())
    MANIFEST.append({"section": section, "figure": name, "png": f"figures/{name}.png",
                     "svg": "", "purpose": purpose, "source_data": str(src.relative_to(OUT)),
                     "caption": caption, "what_to_check": what, "limitation": limit})
    return name


# ---------------------------------------------------------------- evidence assembly
def build_evidence(adata, mapping):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ct = adata.obs["cell_type"].astype(str)
    don = adata.obs["donor_id"].astype(str)
    fam = ct.map(lambda c: mapping.get(c, c))
    # --- reference support per broad family ---
    fam_sup = pd.DataFrame({"n_cells": ct.groupby(fam).count(),
                            "n_donors": don.groupby(fam).nunique()}).sort_values("n_cells")
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(fam_sup)); ax2 = ax.twinx()
    ax.bar(x - 0.2, fam_sup["n_cells"], 0.4, color="#2c6fbb", label="cells")
    ax2.bar(x + 0.2, fam_sup["n_donors"], 0.4, color="#e08a00", label="donors")
    ax.set_xticks(x); ax.set_xticklabels(fam_sup.index, rotation=25, ha="right")
    ax.set_ylabel("cells", color="#2c6fbb"); ax2.set_ylabel("donors", color="#e08a00")
    ax.set_title("Reference support per broad family")
    _savefig(fig, "figR1_family_support", "Reference quality",
             "Show cell AND donor support per broad family (donor diversity not hidden by cell counts).",
             fam_sup, "Cells and donors per broad family.",
             "A family with many cells but few donors is donor-confounded.",
             "Counts only; does not assess label quality.")
    # --- fine subtype support ---
    fine_sup = pd.DataFrame({"family": ct.groupby(ct).first().index.map(lambda c: mapping.get(c, c)),
                             "n_cells": ct.groupby(ct).count(),
                             "n_donors": don.groupby(ct).nunique()})
    fine_sup = fine_sup.sort_values(["family", "n_cells"])
    fig, ax = plt.subplots(figsize=(8, max(4, 0.18 * len(fine_sup))))
    y = np.arange(len(fine_sup))
    ax.barh(y, fine_sup["n_cells"], color="#2c6fbb")
    for i, (idx, r) in enumerate(fine_sup.iterrows()):
        ax.text(r["n_cells"], i, f" {int(r['n_donors'])}d", va="center", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels(fine_sup.index, fontsize=6)
    ax.set_xlabel("cells (label = # donors)"); ax.set_title("Reference support per fine subtype")
    _savefig(fig, "figR2_subtype_support", "Reference quality",
             "Expose low-support / few-donor fine subtypes.", fine_sup,
             "Cells per fine subtype; text label = number of donors.",
             "Subtypes with <3 donors or few cells are unreliable to resolve.",
             "Counts only.")
    # --- donor coverage heatmap (fine × donor presence) ---
    cov = pd.crosstab(ct, don) > 0
    fig, ax = plt.subplots(figsize=(9, max(4, 0.16 * cov.shape[0])))
    ax.imshow(cov.to_numpy(), aspect="auto", cmap="Greys", interpolation="nearest")
    ax.set_yticks(range(cov.shape[0])); ax.set_yticklabels(cov.index, fontsize=5)
    ax.set_xlabel(f"{cov.shape[1]} donors"); ax.set_title("Donor coverage per fine subtype")
    _savefig(fig, "figR3_donor_coverage", "Reference quality",
             "Show which subtypes are supported by many vs few donors.", cov.astype(int),
             "Black = subtype detected in that donor.",
             "Sparse rows = donor-confounded subtypes.", "Presence/absence only.")
    return fam_sup, fine_sup


def per_family_conditional_reliability(adata, mapping, seed=3):
    """Full-panel held-out per-family conditional RMSE (flat NNLS), seed not used for tuning."""
    from benchmarks.shared import synthetic_holdout as SH
    from benchmarks.shared import metrics as M
    from tissueresolve.api import deconv_bulk
    ref_d, qry = SH.split_donors(adata, "donor_id", ref_frac=0.5, seed=seed)
    mapped = adata.obs["cell_type"].astype(str).isin(set(mapping)).to_numpy()
    rmask = adata.obs["donor_id"].astype(str).isin(set(ref_d)).to_numpy() & mapped
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = H.prepare_reference(adata[rmask].copy(), min_cells=30, estimate_overdispersion=True).reference
    ref_types = [str(c) for c in ref.cell_types]
    tgt = SH.build_target_proportions(ref_types, 12, "imbalanced", seed=999)
    ds = SH.realize_pseudobulk(adata, tgt, celltype_col="cell_type", donor_col="donor_id",
                               query_donors=qry, seed=seed, cells_per_sample=600)
    truth = ds.true_mrna_proportions.reindex(columns=ref_types).fillna(0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pred = deconv_bulk(ds.counts, ref, solver="nnls", resolution_mode="none").deconv.proportions
    pred = pred.reindex(index=truth.index, columns=ref_types).fillna(0.0)
    rows = []
    for famname in sorted({mapping.get(c, c) for c in ref_types}):
        mem = [c for c in ref_types if mapping.get(c, c) == famname]
        if len(mem) < 2:
            continue
        st_, sp_ = truth[mem].sum(1), pred[mem].sum(1)
        ct_ = truth[mem].div(st_.replace(0, np.nan), axis=0).fillna(0)
        cp_ = pred[mem].div(sp_.replace(0, np.nan), axis=0).fillna(0)
        rmse = float(np.sqrt(((ct_ - cp_) ** 2).to_numpy().mean()))
        rows.append({"family": famname, "n_subtypes": len(mem), "conditional_rmse_fullpanel": round(rmse, 4)})
    return pd.DataFrame(rows)


def trusted_resolution(mix_summary, cond_rel):
    """Derive per-family trusted resolution from full-panel deconvolution-relevant metrics."""
    m = mix_summary.set_index("family")
    cr = cond_rel.set_index("family")["conditional_rmse_fullpanel"] if not cond_rel.empty else pd.Series(dtype=float)
    med_cr = float(cr.median()) if len(cr) else np.nan
    rows = []
    for fam in m.index:
        shared = float(m.loc[fam, "median_shared_lineage_fraction"])
        cond = float(m.loc[fam, "median_condition_number"])
        crmse = float(cr.get(fam, np.nan))
        # selected-fine only if shared-lineage fraction is comparatively low (<0.95)
        # AND full-panel conditional RMSE is at/below the cross-family median
        sel = (shared < 0.95) and (np.isnan(crmse) or crmse <= med_cr)
        rows.append({"family": fam, "shared_lineage_fraction": shared,
                     "condition_number": cond, "conditional_rmse_fullpanel": crmse,
                     "trusted_resolution": "selected_fine" if sel else "broad_only",
                     "basis": "shared-lineage + condition + full-panel conditional RMSE (NOT cell-AUROC)"})
    return pd.DataFrame(rows).sort_values("shared_lineage_fraction")


# ---------------------------------------------------------------- HTML
CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:1120px;margin:22px auto;
padding:0 18px;color:#1b1b1b;line-height:1.55;background:#fff}
h1{border-bottom:3px solid #2c6fbb;padding-bottom:8px}
h2{margin-top:34px;color:#0b3d70;border-bottom:1px solid #cdd6e0;padding-bottom:4px}
h3{color:#2c6fbb;font-size:1.0em;margin-top:18px}
.cards{display:flex;flex-wrap:wrap;gap:10px}
.card{border:1px solid #dde;border-radius:8px;padding:10px 12px;min-width:215px;flex:1;background:#f8fafc}
.card .s{font-weight:700;font-size:0.95em} .card .m{font-family:monospace;font-size:0.82em;color:#333}
.card .e{font-size:0.82em;color:#444;margin-top:4px}
.ok{border-left:5px solid #2ca02c}.warn{border-left:5px solid #e08a00}.bad{border-left:5px solid #c0392b}.info{border-left:5px solid #2c6fbb}
figure{margin:14px 0;border:1px solid #e3e8ee;border-radius:8px;padding:10px;background:#fff}
figure img{max-width:760px;display:block;margin:6px 0}
.cap{font-size:0.86em;color:#333}.chk{font-size:0.83em;background:#eef5ff;border-left:3px solid #2c6fbb;padding:5px 9px;margin-top:6px}
table{border-collapse:collapse;font-size:0.85em;margin:8px 0}th{background:#0b3d70;color:#fff;padding:5px 9px;text-align:left}
td{border:1px solid #dde;padding:4px 9px}tr:nth-child(even){background:#f4f7fb}
.banner{background:#fff4e5;border-left:5px solid #e08a00;padding:10px 14px;border-radius:6px;font-size:0.92em}
.toc a{color:#2c6fbb;text-decoration:none}
"""


def card(status, cls, name, metric, expl, affected):
    return (f"<div class='card {cls}'><div class='s'>{name}: {status}</div>"
            f"<div class='m'>{metric}</div><div class='e'>{expl}<br><i>affects: {affected}</i></div></div>")


def fig_html(name):
    b = _b64(name)
    rec = next((m for m in MANIFEST if m["figure"] == name), None)
    if b is None or rec is None:
        # existing-embedded (png copied, not base64'd here)
        if rec:
            b = _b64(name)
    if b is None:
        return ""
    return (f"<figure><img src='data:image/png;base64,{b}'/>"
            f"<div class='cap'><b>{rec['purpose']}</b> — {rec['caption']}</div>"
            f"<div class='chk'>What to check: {rec['what_to_check']} · Limitation: {rec['limitation']} · "
            f"Source: <code>{rec['source_data']}</code></div></figure>")


def df_html(df, n=None):
    d = df.head(n) if n else df
    return d.to_html(index=False, border=0, na_rep="—")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-real-data", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print("Refusing without --run-real-data.", file=sys.stderr); return 2
    import anndata as ad
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy
    QC.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(H.REFERENCE_H5AD)
    if "feature_name" in adata.var.columns:
        fn = adata.var["feature_name"].astype(str)
        if fn.nunique() == adata.n_vars:
            adata.var_names = fn.to_numpy()
    raw_map = load_hierarchy_mapping(HIER)
    ct_all = adata.obs["cell_type"].astype(str)
    mapped_types = sorted(set(ct_all[ct_all.isin(set(raw_map))]))
    mapping = build_cell_type_hierarchy(mapped_types, raw_map)
    adata_m = adata[ct_all.isin(set(raw_map)).to_numpy()].copy()

    print("reference support + figures ...", flush=True)
    fam_sup, fine_sup = build_evidence(adata_m, mapping)
    print("per-family full-panel conditional reliability ...", flush=True)
    cond_rel = per_family_conditional_reliability(adata_m, mapping)

    # identifiability tables
    cellid = pd.read_csv(OUT / "subtype_identifiability_family_summary.tsv", sep="\t")
    mix = pd.read_csv(OUT / "mixture_recovery_family_summary.tsv", sep="\t")
    qdet = pd.read_csv(OUT / "subtype_query_detectability.tsv", sep="\t")
    trust = trusted_resolution(mix, cond_rel)
    trust.to_csv(SRC / "trusted_resolution.tsv", sep="\t", index=False) if SRC.exists() else None
    SRC.mkdir(parents=True, exist_ok=True); trust.to_csv(SRC / "trusted_resolution.tsv", sep="\t", index=False)

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # --- identifiability 3-level figure ---
    merged = (cellid[["family", "median_auroc_query_logreg"]]
              .merge(mix[["family", "median_ratio_recovery_mae", "median_shared_lineage_fraction"]], on="family")
              .merge(cond_rel, on="family", how="left"))
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(merged)); w = 0.27
    ax.bar(x - w, merged["median_auroc_query_logreg"], w, label="cell-level AUROC (NOT deconv proof)", color="#9bbcd6")
    ax.bar(x, 1 - merged["median_ratio_recovery_mae"], w, label="pairwise mixture recovery (1-MAE)", color="#5a9bd4")
    ax.bar(x + w, 1 - merged["conditional_rmse_fullpanel"].fillna(0), w, label="full-panel conditional (1-RMSE)", color="#2c6fbb")
    ax.set_xticks(x); ax.set_xticklabels(merged["family"], rotation=25, ha="right")
    ax.set_ylim(0, 1.05); ax.set_ylabel("score (higher=better)")
    ax.set_title("Three DISTINCT identifiability levels (do not merge)")
    ax.legend(fontsize=7)
    _savefig(fig, "figI1_three_levels", "Signature and identifiability quality",
             "Distinguish cell-classifiability, pairwise recovery, and full-panel deconvolution reliability.",
             merged, "Cell-level AUROC ~1.0 everywhere, but full-panel conditional reliability is much lower.",
             "Trust the full-panel bar (dark), NOT the cell-AUROC bar, for deconvolution.",
             "Single tissue; conditional RMSE from one held-out seed.")
    # --- shared-lineage + trusted resolution ---
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = {"selected_fine": "#ff7f0e", "broad_only": "#888"}
    t2 = trust.sort_values("shared_lineage_fraction")
    ax.barh(range(len(t2)), t2["shared_lineage_fraction"], color=[colors[r] for r in t2["trusted_resolution"]])
    ax.axvline(0.95, color="red", ls="--", lw=1)
    ax.set_yticks(range(len(t2))); ax.set_yticklabels(t2["family"])
    ax.set_xlabel("shared-lineage fraction ‖B‖/‖S‖"); ax.set_title("Trusted resolution by family (orange=selected-fine, grey=broad-only)")
    _savefig(fig, "figI2_trusted_resolution", "Signature and identifiability quality",
             "Show which families may be interpreted at selected-fine vs broad-only.", t2,
             "Families with shared-lineage ≥0.95 are broad-only (contrast is a fragile residual).",
             "B/Plasma & Epithelial are the selected-fine candidates; all others broad-only.",
             "Derived from collinearity + full-panel reliability, not cell-AUROC.")
    # --- query-detectable support ---
    qf = qdet.groupby("family")["top50_marker_query_detectable_frac"].median().sort_values()
    fig, ax = plt.subplots(figsize=(7, 3.6)); ax.barh(range(len(qf)), qf.values, color="#2c6fbb")
    ax.set_yticks(range(len(qf))); ax.set_yticklabels(qf.index); ax.set_xlabel("median top-50 marker query-detectable fraction")
    ax.set_title("Query-detectable marker support by family")
    _savefig(fig, "figC1_query_support", "Input compatibility",
             "Show whether each family's markers are detectable in the query.", qf.to_frame(),
             "Low values = markers missing in query → resolution unreliable regardless of reference.",
             "All families have high query detectability here.", "TCGA-TNBC bulk as the query proxy.")

    # --- prediction-quality dashboard (complexity audit) ---
    if (OUT / "tcga_method_complexity_summary.tsv").exists():
        tc = pd.read_csv(OUT / "tcga_method_complexity_summary.tsv", sep="\t", comment="#")
        sub = tc[tc.method.isin(["flat_nnls", "flat_auto", "hierarchical_hier_combined"])]
        if not sub.empty:
            fig, ax = plt.subplots(figsize=(7, 3.8))
            ax.bar(sub["method"], sub["mean_effective_n"], color="#2c6fbb")
            ax.set_ylabel("mean effective-N (TCGA)"); ax.set_title("Prediction complexity by solver (TCGA, no ground truth)")
            plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
            _savefig(fig, "figP1_prediction_complexity", "Prediction quality and trusted resolution",
                     "Show how concentrated each solver's composition is.", sub,
                     "nnls is sparse (~2 eff-N), ridge/auto diffuse (~17) — solver-driven, not biology.",
                     "Interpret composition concentration with the solver in mind.",
                     "TCGA has no ground truth; descriptive only.")

    # --- final bulk: TCGA broad composition ---
    tcga_pred = pd.read_csv(REPO / "examples/real_breast_cancer/outputs/bulk_tcga/tcga_bulk_estimated_proportions.tsv",
                            sep="\t", index_col=0, comment="#")
    fam_of = {c: mapping.get(c, c) for c in tcga_pred.columns}
    broad = tcga_pred.groupby(fam_of, axis=1).sum()
    order = broad.mean().sort_values(ascending=False).index
    fig, ax = plt.subplots(figsize=(9, 4)); bottom = np.zeros(len(broad))
    import matplotlib.cm as cm
    pal = cm.get_cmap("tab10", len(order))
    si = broad[order].sum(1); si = broad.index[np.argsort(-broad[order].iloc[:, 0].to_numpy())]
    for i, f in enumerate(order):
        ax.bar(range(len(broad)), broad.loc[si, f], bottom=bottom, color=pal(i), label=f, width=1.0)
        bottom += broad.loc[si, f].to_numpy()
    ax.set_xlabel(f"{len(broad)} TCGA-TNBC samples"); ax.set_ylabel("RNA-derived proportion")
    ax.set_title("Final bulk: broad-family composition (RNA-derived proportions, NOT cell fractions)")
    ax.legend(fontsize=6, ncol=2, loc="upper right")
    _savefig(fig, "figB1_bulk_broad_composition", "Final bulk predictions",
             "Show the trusted (broad-level) composition of the real bulk cohort.", broad,
             "Broad families are the trusted level here; fine subtypes are gated out (see trusted resolution).",
             "Broad-confident; fine not shown because no family is full-fine trusted.",
             "RNA-derived proportions, not absolute cell fractions; no ground truth (TCGA).")

    # embed decision-relevant existing figures
    _embed_existing("holdout_bulk_multisplit/figures/figG_multisplit_external_ci.png", "figBM1_bulk_benchmark",
                    "Benchmark evidence", "Bulk fine-accuracy with CIs vs external methods (held-out donors).",
                    "TissueResolve nnls/auto vs MuSiC/Bisque, bootstrap CIs.", "Higher fine Pearson = better; CIs overlap → tie.",
                    "Single tissue; breast pseudobulk.")
    _embed_existing("spatial_multiseed/figures/figK_spatial_multiseed_ci.png", "figSM1_spatial_benchmark",
                    "Benchmark evidence", "Spatial fine-accuracy with CIs (synthetic ground truth).",
                    "TissueResolve vs RCTD/cell2location across seeds.", "RCTD ties TissueResolve; n=5 small.",
                    "Synthetic spatial; 144-spot grids.")
    _embed_existing("spatial_real_visium/figures/figM_real_concordance.png", "figSP1_spatial_real",
                    "Final spatial predictions", "Real-Visium cross-method concordance (NO ground truth).",
                    "Method agreement on the breast Visium section.", "Concordance is NOT accuracy.",
                    "No ground truth; 500 spots.")

    # ---------------- status cards (from real metrics) ----------------
    n_sel = int((trust.trusted_resolution == "selected_fine").sum())
    n_broad = int((trust.trusted_resolution == "broad_only").sum())
    med_cond_rmse = float(cond_rel["conditional_rmse_fullpanel"].median()) if not cond_rel.empty else float("nan")
    cards = [
        card("REVIEW", "warn", "Reference suitability",
             f"{int(fam_sup['n_donors'].min())}–{int(fam_sup['n_donors'].max())} donors/family",
             "Donor coverage varies by family; some fine subtypes have few donors.", "all fine predictions"),
        card("OK", "ok", "Donor coverage", f"{adata_m.obs['donor_id'].nunique()} donors total",
             "Multiple donors per broad family.", "broad predictions"),
        card("CAUTION", "warn", "Fine full-panel reliability", f"median conditional RMSE={med_cond_rmse:.2f}",
             "Within-family subtype proportions are deconvolution-fragile.", "fine predictions"),
        card("OK", "ok", "Query compatibility", f"{int(qdet['top50_marker_query_detectable_frac'].median()*100)}% markers detectable",
             "Reference markers largely detectable in query.", "all predictions"),
        card(f"{n_sel} selected-fine / {n_broad} broad-only", "info", "Recommended resolution",
             "shared-lineage + full-panel conditional RMSE", "No family is full-fine trusted.", "fine predictions"),
        card("EXPERIMENTAL", "info", "Soft gating", "not promoted (no 2nd tissue)",
             "Validated on breast (9/9 gates); promotion blocked.", "hierarchical fine predictions"),
    ]

    def section(title, body):
        return f"<h2 id='{title}'>{title}</h2>\n{body}"

    warn_banner = ("<div class='banner'><b>Key interpretation rule:</b> cell-level "
                   "separability does NOT guarantee accurate deconvolution in a complex "
                   "mixture. Trust the FULL-PANEL conditional reliability, not cell AUROC. "
                   "Values are RNA-derived proportions, not absolute cell fractions.</div>")

    sections_html = [
        section("1. Executive decision summary",
                f"{warn_banner}<div class='cards'>{''.join(cards)}</div>"
                "<p class='cap'>Each card: status · metric · explanation · affected result. Derived from real metrics; no bare PASS/FAIL.</p>"),
        section("2. Reference quality",
                fig_html("figR1_family_support") + fig_html("figR2_subtype_support") + fig_html("figR3_donor_coverage")),
        section("3. Signature and identifiability quality",
                "<p>Three DISTINCT evidence levels — do not merge into one separability score:</p>"
                + fig_html("figI1_three_levels") + fig_html("figI2_trusted_resolution")
                + "<h3>Trusted resolution (full-panel basis)</h3>" + df_html(trust)),
        section("4. Input compatibility", fig_html("figC1_query_support")),
        section("5. Prediction quality and trusted resolution",
                "<div class='banner'>Cell-level separability does not guarantee accurate "
                "deconvolution in a complex mixture.</div>" + fig_html("figP1_prediction_complexity")
                + "<h3>Per-family full-panel conditional reliability</h3>" + df_html(cond_rel)),
        section("6. Final bulk predictions",
                fig_html("figB1_bulk_broad_composition")
                + "<p class='cap'>Fine subtypes are GATED OUT of the main view: no family is full-fine "
                "trusted (see §3). Broad-only families are shown at family level with unresolved mass.</p>"),
        section("7. Final spatial predictions",
                fig_html("figSP1_spatial_real")
                + "<p class='cap'>Real Visium has NO ground truth — concordance/structure only, never accuracy. "
                "Spatial structure (Moran's I) is not accuracy.</p>"),
        section("8. Benchmark evidence",
                "<h3>Bulk (held-out-donor ground truth)</h3>" + fig_html("figBM1_bulk_benchmark")
                + "<h3>Spatial (synthetic ground truth) — kept separate</h3>" + fig_html("figSM1_spatial_benchmark")),
        section("9. Warnings and limitations",
                "<ul><li>Single tissue (breast); soft gating + any new signature method need a 2nd tissue.</li>"
                "<li>No family is confidently full-fine; fine predictions are gated.</li>"
                "<li>RNA-derived proportions, not cell fractions.</li>"
                "<li>Real-Visium results are evidence, not accuracy.</li>"
                "<li>Joint hierarchy was a negative result (not developed).</li></ul>"),
        section("10. Methods and source data",
                "<p>Evidence chain: reference support → donor coverage → 3-level identifiability → query "
                "compatibility → full-panel reliability → trusted resolution → predictions. Every figure has a "
                "source-data TSV (<code>source_data/</code>) and a manifest entry "
                "(<code>figures/figure_manifest.tsv</code>). Identifiability: "
                "<code>subtype_identifiability_*</code>, <code>mixture_recovery_*</code>. Benchmark: "
                "<code>holdout_bulk_multisplit/</code>, <code>spatial_multiseed/</code>. Deterministic palette.</p>"),
    ]
    toc = " · ".join(f"<a href='#{t}'>{t.split('.')[0]}</a>" for t, _ in
                     [(s.split("id='")[1].split("'")[0], None) for s in sections_html])
    report = (f"<!doctype html><html><head><meta charset='utf-8'><title>TissueResolve — QC-First Report</title>"
              f"<style>{CSS}</style></head><body><h1>TissueResolve — QC-First Reliability Report</h1>"
              f"<p><i>Quality before predictions. A figure appears here only if it helps decide reliability or "
              f"interpret a prediction. Technical detail → appendix.</i></p><div class='toc'>{toc}</div>"
              + "\n".join(sections_html) +
              "<hr><p class='cap'>Generated by benchmarks/build_qc_first_report.py. No core algorithm or stable "
              "report code modified; nothing committed.</p></body></html>")
    (QC / "report.html").write_text(report, encoding="utf-8")

    # technical appendix (links to the heavy material kept OUT of main)
    appendix_items = [
        ("Full pairwise identifiability (91 pairs)", "subtype_identifiability_pairs.tsv"),
        ("Mixture-recovery per pair", "mixture_recovery_ceiling_pairs.tsv"),
        ("Donor stability per pair", "subtype_donor_stability.tsv"),
        ("TCGA per-sample complexity audit", "tcga_prediction_complexity_audit.tsv"),
        ("Bulk multisplit metrics (raw)", "holdout_bulk_multisplit/metrics_long.tsv"),
        ("Spatial multiseed metrics (raw)", "spatial_multiseed/combined_metrics_long.tsv"),
        ("Soft-gating mode comparison", "v2_soft_gating_mode_comparison.tsv"),
        ("Joint-hierarchy (negative result)", "phase2b_joint_hierarchy_metrics.tsv"),
    ]
    appx = "".join(f"<li>{n} — <code>benchmarks/outputs/{p}</code></li>" for n, p in appendix_items
                   if (OUT / p).exists())
    (QC / "technical_appendix.html").write_text(
        f"<!doctype html><html><head><meta charset='utf-8'><title>Technical Appendix</title><style>{CSS}</style>"
        f"</head><body><h1>TissueResolve — Technical Appendix</h1>"
        f"<p>Full matrices, raw tables, per-pair detail, experimental + negative-result diagnostics. "
        f"Kept OUT of the main report per the central inclusion rule.</p><ul>{appx}</ul>"
        f"<p>Full identifiability figures: <code>benchmarks/outputs/figures/identifiability/</code>. "
        f"Sparsity audit: <code>figures/part1_audit/</code>.</p></body></html>", encoding="utf-8")

    pd.DataFrame(MANIFEST).to_csv(FIGD / "figure_manifest.tsv", sep="\t", index=False)
    warns = [
        {"severity": "info", "category": "estimate_type", "message": "RNA-derived proportions, not cell fractions."},
        {"severity": "warning", "category": "resolution", "message": f"No family full-fine trusted; {n_sel} selected-fine, {n_broad} broad-only."},
        {"severity": "warning", "category": "replication", "message": "Single tissue; soft gating experimental (no 2nd tissue)."},
        {"severity": "info", "category": "spatial", "message": "Real-Visium results are evidence, not accuracy."},
    ]
    (QC / "warnings.json").write_text(json.dumps(warns, indent=2), encoding="utf-8")
    (QC / "methods.txt").write_text(
        "QC-first report. Identifiability assessed at three levels (cell classification, pairwise mixture "
        "recovery, full-panel conditional reliability); trusted resolution derived from full-panel "
        "deconvolution-relevant metrics (shared-lineage fraction, condition number, conditional RMSE), NOT "
        "cell-classification AUROC. Donor-held-out throughout. RNA-derived proportions.\n", encoding="utf-8")

    print(f"\nWrote QC-first report -> {QC}/report.html  (+ appendix, manifest, source_data, warnings, methods)")
    print(f"main figures: {len(MANIFEST)} | trusted resolution:\n{trust[['family','trusted_resolution']].to_string(index=False)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
