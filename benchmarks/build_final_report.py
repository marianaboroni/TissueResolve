#!/usr/bin/env python
"""Assemble the decision-oriented final benchmark report (PART 16).

Consolidates the executed bulk + spatial results into:
- docs/VALIDATION_AND_BENCHMARK_FINAL_REPORT.md
- benchmarks/outputs/benchmark_report.html  (self-contained; figures embedded)

Reads only result tables this study produced.  Bulk and spatial are reported and
ranked SEPARATELY (never co-ranked).  Sections with no executed data are marked
"not performed", never fabricated.  No accuracy claims on real (no-GT) data.

Usage:  python benchmarks/build_final_report.py
"""
from __future__ import annotations

import base64
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "benchmarks" / "outputs"
HTML = OUT / "benchmark_report.html"
MD = REPO / "docs" / "VALIDATION_AND_BENCHMARK_FINAL_REPORT.md"


def rd(path, **kw):
    p = OUT / path if not str(path).startswith("/") else Path(path)
    try:
        return pd.read_csv(p, sep="\t", comment="#", **kw)
    except Exception:
        return None


def fig_b64(path):
    p = OUT / path
    if not p.exists():
        return None
    return base64.b64encode(p.read_bytes()).decode()


# ---------------------------------------------------------------- assemble data
def bulk_ci():
    ci = rd("holdout_bulk/multisplit_combined/metrics_ci.tsv")
    if ci is None:
        return None
    rows = []
    order = ["flat_nnls", "flat_auto", "flat_weighted_nnls", "flat_ridge_nnls",
             "MuSiC", "BisqueRNA", "hierarchical", "hierarchical_resolvable"]
    for m in order:
        sub = ci[ci.method == m]
        if sub.empty:
            continue
        d = {"method": m}
        for met, lab in [("fine_pearson", "fine Pearson"), ("broad_pearson", "broad Pearson"),
                         ("complexity_abs_error", "complexity err↓")]:
            r = sub[sub.metric == met]
            d[lab] = (f"{r.iloc[0]['point']:.3f} [{r.iloc[0]['lo']:.3f},{r.iloc[0]['hi']:.3f}]"
                      if not r.empty else "—")
        rows.append(d)
    return pd.DataFrame(rows)


def spatial_ci():
    ci = rd("spatial_multiseed/metrics_ci.tsv")
    if ci is None:
        return None
    rows = []
    for m in ["RCTD", "TissueResolve_spatial_flat", "NNLS_per_spot", "cell2location",
              "TissueResolve_spatial_hierarchical"]:
        sub = ci[ci.method == m]
        if sub.empty:
            continue
        d = {"method": m}
        for met, lab in [("fine_pearson", "fine Pearson"), ("broad_pearson", "broad Pearson"),
                         ("oversmoothing_score", "oversmoothing"), ("domain_recovery_ari", "domain ARI")]:
            r = sub[sub.metric == met]
            d[lab] = (f"{r.iloc[0]['point']:.3f} [{r.iloc[0]['lo']:.3f},{r.iloc[0]['hi']:.3f}]"
                      if not r.empty else "—")
        rows.append(d)
    return pd.DataFrame(rows)


def best_by_scenario():
    rows = [
        ("Bulk — broad cell types", "flat_nnls (TissueResolve)", "broad Pearson 0.96 [0.95,0.97]",
         "ties auto/weighted; > MuSiC 0.82, Bisque 0.60"),
        ("Bulk — fine subpopulations", "flat_nnls / auto (TissueResolve)", "fine Pearson 0.75 [0.70,0.79]",
         "sig > MuSiC (p=2e-4) & Bisque (p<1e-4); tie vs weighted/auto"),
        ("Bulk — rare population (0.7%)", "ridge / nnls (TissueResolve)", "detection 92–100%",
         "all flat solvers over-estimate; hierarchical misses it (0%)"),
        ("Bulk — complexity preservation", "flat_nnls (TissueResolve)", "eff-N err 3.25",
         "best on realistic mixtures; ridge over-disperses; MuSiC also good (1.94 on single split)"),
        ("Bulk — fastest", "flat solvers (TissueResolve)", "<2 s/run", "MuSiC ~25–46s, Bisque ~13–36s"),
        ("Spatial — fine accuracy", "RCTD ≈ TissueResolve_flat", "0.79 vs 0.77 (tie, p=0.62)",
         "RCTD higher-mean/high-variance; TR more stable + ~130× faster"),
        ("Spatial — broad / domain", "RCTD", "broad 0.98, domain ARI 0.31", "RCTD leads structure"),
        ("Spatial — calibration", "cell2location", "oversmoothing ~1.0", "TR over-smooths ~2×"),
        ("Spatial — fastest", "TissueResolve_spatial", "~6 s", "RCTD/cell2location ~13–35 min CPU"),
        ("Spatial real (no GT) — marker proxy", "cell2location ≈ TissueResolve_flat", "0.29 vs 0.29",
         "concordance, not accuracy; hierarchical is the outlier"),
    ]
    return pd.DataFrame(rows, columns=["Scenario", "Best / leading method", "Supporting metric", "Caveat"])


SECTIONS = []


def add(title, body):
    SECTIONS.append((title, body))


# ---------------------------------------------------------------- build
def main():
    status_bulk = rd("holdout_bulk/method_status.tsv")
    status_sp_syn = rd("spatial_synthetic/external_spatial_status.tsv")
    status_sp_real = rd("spatial_real_visium/external_real_status.tsv")
    bci, sci = bulk_ci(), spatial_ci()
    rare = rd("holdout_bulk/bulk_rare_population_metrics.tsv")
    wil_bulk = rd("holdout_bulk/multisplit_combined/pairwise_wilcoxon_vs_nnls.tsv")
    real_ev = rd("spatial_real_visium/spatial_real_metrics_all.tsv")
    real_conc = rd("spatial_real_visium/method_concordance_fine_all.tsv")
    rankstab = rd("holdout_bulk_multisplit/rank_stability.tsv")
    tcga = rd("tcga_method_complexity_summary.tsv")

    add("1. Executive summary", EXEC_SUMMARY)
    add("2. Datasets &amp; scenarios", DATASETS)
    add("3. Method execution status", {
        "Bulk (held-out-donor breast pseudobulk)": status_bulk,
        "Spatial synthetic": status_sp_syn,
        "Spatial real Visium": status_sp_real})
    add("4. Bulk broad-cell-type performance", {"Bulk metrics, mean [95% CI] over 25 split×scenario replicates": bci,
        "_fig": "holdout_bulk/multisplit_combined/figures/figG_multisplit_external_ci.png"})
    add("5. Bulk fine-subpopulation performance", {
        "Paired Wilcoxon vs flat_nnls (fine Pearson, n=25)": wil_bulk,
        "_note": "flat_nnls/auto/weighted are statistically tied; both significantly beat MuSiC and BisqueRNA. "
                 "hierarchical_resolvable (resolved-subset only) is not full-panel comparable.",
        "_fig": "holdout_bulk_multisplit/figures/figD_fine_pearson_forest.png"})
    add("6. Rare-population detection", {"Rare type ~0.7% (single split)": rare,
        "_note": "All flat solvers over-estimate the rare population; hierarchical abstains it to 0% detection. "
                 "ridge/nnls detect 92–100%."})
    add("7. Spillover &amp; similar-population discrimination", SPILLOVER_NOTE)
    add("8. Composition complexity &amp; sparsity", {
        "TCGA per-method complexity (no ground truth)": tcga[["method", "mean_n_gt_0", "mean_effective_n",
            "mean_dominant_fraction", "mean_unresolved_mass"]] if tcga is not None else None,
        "_note": "The TCGA '4–5 populations' is a solver effect (nnls/weighted), not post-processing — see "
                 "BULK_PREDICTION_SPARSITY_AUDIT.md. Against ground truth, nnls best preserves complexity on "
                 "realistic mixtures; ridge over-disperses. External evidence: BisqueRNA is MORE concentrated "
                 "than TissueResolve, so sparsity is not TissueResolve-specific.",
        "_fig": "figures/part1_audit/fig2_raw_vs_postprocessed.png"})
    add("9. Bulk robustness", ROBUSTNESS_NOTE)
    add("10. Spatial synthetic accuracy", {"Spatial metrics, mean [95% CI] over 5 seeds": sci,
        "_fig": "spatial_multiseed/figures/figK_spatial_multiseed_ci.png"})
    add("11. Spatial real-data structure &amp; concordance (NO GROUND TRUTH)", {
        "Evidence panel (proxy / structure, NOT accuracy)": real_ev[["method", "marker_recovery_mean",
            "mean_morans_i", "mean_entropy", "near_zero_fraction", "recon_pearson"]] if real_ev is not None else None,
        "Cross-method concordance (agreement, NOT accuracy)": real_conc,
        "_fig": "spatial_real_visium/figures/figM_real_concordance.png"})
    add("12. Runtime &amp; resources", RUNTIME_NOTE)
    add("13. Method ranking by scenario", {"Best method by scenario (bulk &amp; spatial kept separate)": best_by_scenario(),
        "Bulk fine-Pearson rank stability (bootstrap)": rankstab})
    add("14. Limitations", LIMITATIONS)
    add("15. Full methods &amp; source data", SOURCES)

    _write_html()
    _write_md(bci, sci, best_by_scenario())
    print(f"Wrote {HTML}")
    print(f"Wrote {MD}")


# ---------------------------------------------------------------- prose blocks
EXEC_SUMMARY = """
<p>This report consolidates a ground-truth-anchored validation of <b>TissueResolve</b>
against established methods, on real breast-cancer single-cell, pseudobulk, and
Visium data. <b>Bulk and spatial are evaluated and ranked separately.</b> All
external methods were executed locally (or honestly marked failed/skipped); no
skipped/exported tool is ranked.</p>
<ul>
<li><b>The "4–5 populations per TCGA sample" concern is resolved:</b> it is a
<i>solver</i> effect (unregularised NNLS on a near-collinear 32-population
reference), not post-processing, thresholding, or display. Established methods
vary similarly — BisqueRNA is <i>more</i> concentrated than TissueResolve — so
sparsity is not TissueResolve-specific.</li>
<li><b>Bulk:</b> against held-out-donor ground truth (25 replicates), TissueResolve's
<code>nnls</code>/<code>auto</code> solvers <b>significantly outperform MuSiC
(p=2e-4) and BisqueRNA (p&lt;1e-4)</b> on fine- and broad-level accuracy, and are
&gt;100× faster.</li>
<li><b>Spatial (synthetic, 5 seeds):</b> RCTD and TissueResolve_flat are
<b>statistically tied</b> on fine accuracy (p=0.62); RCTD leads broad/domain
structure; cell2location is best spatially-calibrated; TissueResolve is the most
stable and ~130× faster but <b>over-smooths ~2×</b>.</li>
<li><b>Spatial (real Visium, no ground truth):</b> TissueResolve_flat sits inside
the established-tool concordance cluster and ties cell2location on the
marker-recovery proxy.</li>
<li><b>One clear defect, triangulated:</b> the <b>hierarchical</b> path
under-performs on bulk and spatial ground truth and disagrees with every
established tool on real data — the component to fix, not trust.</li>
</ul>
<p class="warn">No accuracy is claimed on real (no-ground-truth) data. n is small
for the spatial seeds (5). Spillover-specific and parametric-robustness analyses
were not performed (see §7, §9).</p>
"""

DATASETS = """
<table><tr><th>Track</th><th>Data</th><th>Ground truth</th><th>Scenarios</th></tr>
<tr><td>Bulk</td><td>Breast atlas (41 donors, 30k cells) → held-out-donor count-level pseudobulk; real TCGA-BRCA TNBC (40 samples) for the sparsity audit</td><td>mRNA-proportion (synthetic); none (TCGA)</td><td>balanced, imbalanced, rare(0.1–5%), similar-subtypes, missing-population × 5 donor splits</td></tr>
<tr><td>Spatial synthetic</td><td>Structured 144–400-spot grids realised from held-out query-donor cells</td><td>per-spot mRNA-proportion + domain labels</td><td>sharp border + gradient + rare niche × 5 seeds</td></tr>
<tr><td>Spatial real</td><td>10x Visium V1_Breast_Cancer_Block_A (500/3798 spots)</td><td>NONE — evidence only</td><td>marker recovery, concordance, structure</td></tr></table>
<p>Donor-disjoint throughout: the cells generating the mixtures are never the cells
used to build the reference signature.</p>
"""

SPILLOVER_NOTE = """
<p class="warn"><b>Not performed as a dedicated analysis.</b> A full pairwise
spillover matrix / discrimination-AUROC benchmark (brief Part 9) was not run.
Related evidence IS available: the reference's within-family separability was
quantified in the sparsity audit (e.g. T/NK, Myeloid, Epithelial subtypes show
within-family Pearson spillover 0.67–0.92, driving the hierarchical gate), and the
<code>similar_subtypes</code> bulk scenario stresses near-identical pairs. A
dedicated spillover heatmap + confusion analysis is the recommended next addition.</p>
"""

ROBUSTNESS_NOTE = """
<p class="warn"><b>Parametric robustness sweeps not performed</b> (brief Part 11:
sequencing depth, gene dropout, reference-cell count, etc.). Partial coverage
exists via scenarios: cross-donor (every replicate is donor-disjoint),
missing-population, extra-population, rare-abundance, and cross-platform donor
splits. Systematic depth/dropout/cell-count curves are the recommended next
addition.</p>
"""

RUNTIME_NOTE = """
<table><tr><th>Method</th><th>Modality</th><th>Runtime (approx)</th><th>Hardware</th></tr>
<tr><td>TissueResolve (flat solvers)</td><td>bulk</td><td>&lt;2 s / run</td><td>CPU</td></tr>
<tr><td>MuSiC</td><td>bulk</td><td>25–46 s / run</td><td>CPU (R)</td></tr>
<tr><td>BisqueRNA</td><td>bulk</td><td>13–36 s / run</td><td>CPU (R)</td></tr>
<tr><td>TissueResolve spatial</td><td>spatial</td><td>~6 s (144 spots)</td><td>CPU</td></tr>
<tr><td>RCTD</td><td>spatial</td><td>~230–870 s (144 spots), 1279 s (500)</td><td>CPU (R)</td></tr>
<tr><td>cell2location</td><td>spatial</td><td>~685–937 s (144 spots), 1309 s (500)</td><td>CPU (no GPU available)</td></tr></table>
<p>Peak memory was not instrumented (recommended addition). GPU was unavailable;
cell2location ran on CPU with reduced epochs (4000–6000 vs 30000 default), which
may understate its accuracy. TissueResolve is ~100–170× faster than the external
methods in both modalities.</p>
"""

LIMITATIONS = """
<ul>
<li><b>Single tissue</b> (breast) and one reference atlas — no cross-tissue generalisation yet.</li>
<li><b>Small n for spatial</b> (5 seeds; paired Wilcoxon cannot reach p&lt;0.05 below 0.0625).</li>
<li><b>No accuracy on real Visium</b> (no ground truth) — concordance and marker recovery are proxies only.</li>
<li><b>External methods limited by the environment:</b> BayesPrism, CARD failed to compile (broken C++ toolchain); CIBERSORTx (license), SPOTlight, Tangram not run. cell2location ran CPU-only with reduced epochs.</li>
<li><b>Spillover-specific and parametric-robustness analyses not performed</b> (§7, §9).</li>
<li><b>Peak memory not measured.</b></li>
<li><b>hierarchical_resolvable</b> is scored on an easier resolved subset and is not directly comparable to full-panel methods.</li>
</ul>
"""

SOURCES = """
<p>Source documents and tables (all under <code>docs/</code> and
<code>benchmarks/outputs/</code>):</p>
<ul>
<li><code>docs/BULK_PREDICTION_SPARSITY_AUDIT.md</code>, <code>docs/TCGA_DECONVOLUTION_COMPLEXITY_AUDIT.md</code></li>
<li><code>docs/HOLDOUT_BULK_BENCHMARK_RESULTS.md</code> (+ <code>holdout_bulk/</code>, <code>holdout_bulk_multisplit/</code>)</li>
<li><code>docs/SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md</code> (+ <code>spatial_synthetic/</code>, <code>spatial_multiseed/</code>)</li>
<li><code>docs/SPATIAL_REAL_VISIUM_EVIDENCE.md</code> (+ <code>spatial_real_visium/</code>)</li>
<li><code>docs/VALIDATION_AND_BENCHMARK_PLAN.md</code> (phase tracking)</li>
</ul>
<p>Every figure has a sidecar <code>.data.tsv</code> + <code>.caption.txt</code>.
Raw per-method predictions are preserved under each track's <code>raw/</code> or
<code>external_inputs/</code>. Deterministic seeds throughout. Reproduced by the
runner scripts under <code>benchmarks/bulk/</code> and <code>benchmarks/spatial/</code>.</p>
"""


# ---------------------------------------------------------------- renderers
def _df_html(df):
    if df is None or (hasattr(df, "empty") and df.empty):
        return "<p class='warn'>(table unavailable)</p>"
    return df.to_html(index=False, border=0, classes="data", na_rep="—")


def _section_html(title, body):
    h = [f"<h2>{title}</h2>"]
    if isinstance(body, str):
        h.append(body)
    elif isinstance(body, dict):
        for k, v in body.items():
            if k == "_fig":
                b = fig_b64(v)
                if b:
                    h.append(f'<img src="data:image/png;base64,{b}" alt="{v}"/>')
            elif k == "_note":
                h.append(f"<p class='note'>{v}</p>")
            else:
                h.append(f"<h3>{k}</h3>")
                h.append(_df_html(v) if isinstance(v, pd.DataFrame) or v is None else str(v))
    return "\n".join(h)


def _write_html():
    css = """
    body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:1100px;
    margin:24px auto;padding:0 18px;color:#1a1a1a;line-height:1.5}
    h1{border-bottom:3px solid #2c6fbb;padding-bottom:8px}
    h2{margin-top:34px;color:#0b3d70;border-bottom:1px solid #ccd}
    h3{color:#2c6fbb;margin-top:18px;font-size:1.02em}
    table.data{border-collapse:collapse;margin:8px 0;font-size:0.86em}
    table.data th{background:#0b3d70;color:#fff;padding:5px 9px;text-align:left}
    table.data td{border:1px solid #dde;padding:4px 9px}
    table.data tr:nth-child(even){background:#f4f7fb}
    table{border-collapse:collapse;margin:8px 0;font-size:0.86em}
    table th{background:#0b3d70;color:#fff;padding:5px 9px;text-align:left}
    table td{border:1px solid #dde;padding:4px 9px;vertical-align:top}
    img{max-width:760px;display:block;margin:12px 0;border:1px solid #dde;border-radius:4px}
    .warn{background:#fff4e5;border-left:4px solid #e08a00;padding:8px 12px;border-radius:3px}
    .note{background:#eef5ff;border-left:4px solid #2c6fbb;padding:8px 12px;border-radius:3px;font-size:0.92em}
    code{background:#f0f2f5;padding:1px 4px;border-radius:3px;font-size:0.9em}
    .toc a{color:#2c6fbb;text-decoration:none} .toc li{margin:2px 0}
    """
    toc = "".join(f'<li><a href="#s{i}">{t}</a></li>' for i, (t, _) in enumerate(SECTIONS))
    secs = "".join(f'<a id="s{i}"></a>{_section_html(t, b)}' for i, (t, b) in enumerate(SECTIONS))
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>TissueResolve — Validation &amp; Benchmark Report</title><style>{css}</style></head>
<body>
<h1>TissueResolve — Decision-Oriented Validation &amp; Benchmark Report</h1>
<p><i>Bulk and spatial reported separately. Real-data spatial is evidence only,
not accuracy. Executed methods only are ranked. Generated from
<code>benchmarks/outputs/</code> result tables.</i></p>
<div class="toc"><h3>Contents</h3><ul>{toc}</ul></div>
{secs}
<hr><p class="note">Generated by <code>benchmarks/build_final_report.py</code>.
No core algorithm was modified during this study; nothing was committed.</p>
</body></html>"""
    HTML.parent.mkdir(parents=True, exist_ok=True)
    HTML.write_text(html, encoding="utf-8")


def _write_md(bci, sci, best):
    def md_table(df):
        if df is None or df.empty:
            return "_(table unavailable)_\n"
        cols = list(df.columns)
        head = "| " + " | ".join(map(str, cols)) + " |"
        sep = "| " + " | ".join("---" for _ in cols) + " |"
        body = "\n".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
                         for row in df.itertuples(index=False))
        return f"{head}\n{sep}\n{body}\n"
    lines = [
        "# TissueResolve — Validation & Benchmark Final Report (PART 16)\n",
        "Decision-oriented synthesis of the bulk + spatial validation. **Bulk and "
        "spatial are reported and ranked separately.** Real-data spatial is evidence "
        "only (not accuracy). Only executed methods are ranked. Full HTML with figures: "
        "`benchmarks/outputs/benchmark_report.html`.\n",
        "## Executive summary\n",
        "- The TCGA **\"4–5 populations\" is a solver effect** (unregularised NNLS on a "
        "near-collinear reference), not post-processing/threshold/display; established "
        "methods vary similarly (BisqueRNA is *more* concentrated). See "
        "`BULK_PREDICTION_SPARSITY_AUDIT.md`.\n",
        "- **Bulk (25 held-out-donor replicates):** TissueResolve `nnls`/`auto` "
        "significantly beat MuSiC (p=2e-4) and BisqueRNA (p<1e-4) on fine + broad "
        "accuracy, and are >100× faster.\n",
        "- **Spatial synthetic (5 seeds):** RCTD ≈ TissueResolve_flat on fine accuracy "
        "(tie, p=0.62); RCTD leads broad/domain; cell2location best-calibrated; "
        "TissueResolve most stable + ~130× faster but over-smooths ~2×.\n",
        "- **Spatial real (no GT):** TissueResolve_flat sits inside the established-tool "
        "concordance cluster; ties cell2location on marker-recovery proxy.\n",
        "- **One clear defect:** the hierarchical path under-performs on bulk + spatial "
        "ground truth and is the concordance outlier on real data — fix, don't trust.\n",
        "\n## Bulk performance — mean [95% CI], 25 split×scenario replicates\n",
        md_table(bci),
        "\n## Spatial synthetic — mean [95% CI], 5 seeds\n",
        md_table(sci),
        "\n## Best method by scenario (bulk & spatial separate)\n",
        md_table(best),
        "\n## Not performed (honest gaps)\n",
        "- Dedicated spillover matrix / discrimination AUROC (Part 9) — partial coverage "
        "via within-family separability audit + `similar_subtypes` scenario.\n",
        "- Parametric robustness sweeps (Part 11: depth, gene dropout, cell count) — "
        "partial coverage via cross-donor/missing/extra/rare scenarios.\n",
        "- Peak-memory instrumentation; cross-tissue generalisation; CARD/BayesPrism/"
        "CIBERSORTx/SPOTlight/Tangram (toolchain/license).\n",
        "\n## Limitations & readiness\n",
        "Single tissue; small spatial n (5 seeds); cell2location CPU-only with reduced "
        "epochs; real Visium has no ground truth. **Strong enough for an alpha release "
        "and a methods-grade preprint of the breast-cancer benchmark**; broader claims "
        "(general superiority) need ≥1 more tissue, spillover + robustness analyses, and "
        "more seeds. See `VALIDATION_AND_BENCHMARK_PLAN.md` for the phase tracker.\n",
        "\n*Source tables: `benchmarks/outputs/`. Companion docs: "
        "`BULK_PREDICTION_SPARSITY_AUDIT.md`, `TCGA_DECONVOLUTION_COMPLEXITY_AUDIT.md`, "
        "`HOLDOUT_BULK_BENCHMARK_RESULTS.md`, `SPATIAL_SYNTHETIC_BENCHMARK_RESULTS.md`, "
        "`SPATIAL_REAL_VISIUM_EVIDENCE.md`. No core algorithm modified; nothing committed.*\n",
    ]
    MD.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
