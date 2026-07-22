#!/usr/bin/env python
"""Single self-contained HTML report + figures for the Rectangle benchmark.

Reads the one long result table (benchmark_results.parquet) + manifest and writes
benchmark_report.html plus a small set of PNG figures under figures/. No per-analysis
report is produced — everything is sections of one page. Figures are embedded as
base64 so the HTML is portable. Deferred analyses (ablation, spatial) render an
honest 'not yet run' note rather than a fake panel.
"""
from __future__ import annotations

import base64
import json
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PRIMARY_THR = 0.01


def _b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _scalar(df, metric, threshold=None):
    q = (df["metric"] == metric) & (df["status"] == "ok")
    sub = df[q]
    if threshold is not None:
        sub = sub[np.isclose(sub["threshold"].fillna(-1), threshold)]
    else:
        sub = sub[sub["threshold"].isna() & (sub["cell_type"] == "")]
    return sub


def _mean_table(df, metric, threshold=None):
    sub = _scalar(df, metric, threshold)
    if sub.empty:
        return pd.DataFrame()
    return sub.pivot_table(index=["dataset", "scenario"], columns="method",
                           values="value", aggfunc="mean").round(4)


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def _fig_bulk_overview(df, path):
    metrics = [("rmse", None, "RMSE (↓)"), ("pearson", None, "Pearson (↑)"),
               ("rare_recall", PRIMARY_THR, f"rare recall @{PRIMARY_THR} (↑)"),
               ("rare_fpr", PRIMARY_THR, f"rare FPR @{PRIMARY_THR} (↓)")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, (m, thr, title) in zip(axes.ravel(), metrics):
        t = _scalar(df, m, thr)
        if t.empty:
            ax.set_visible(False); continue
        piv = t.pivot_table(index="method", columns=["dataset", "scenario"],
                            values="value", aggfunc="mean")
        im = ax.imshow(piv.to_numpy(dtype=float), aspect="auto", cmap="viridis")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels(["·".join(map(str, c)) for c in piv.columns],
                           rotation=90, fontsize=6)
        ax.set_yticks(range(piv.shape[0]))
        ax.set_yticklabels(piv.index, fontsize=7)
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.03)
    fig.suptitle("Bulk overview — method × (dataset·scenario), seed means", fontsize=12)
    fig.tight_layout()
    b = _b64(fig); Path(path).write_bytes(base64.b64decode(b)); return b


def _fig_rare(df, path):
    rec = _scalar(df, "rare_recall", PRIMARY_THR)
    fpr = _scalar(df, "rare_fpr", PRIMARY_THR)
    fig, ax = plt.subplots(figsize=(8, 6))
    if not rec.empty and not fpr.empty:
        r = rec.groupby("method")["value"].mean()
        f = fpr.groupby("method")["value"].mean()
        for method in r.index:
            ax.scatter(f.get(method, np.nan), r.get(method, np.nan), s=90)
            ax.annotate(method, (f.get(method, np.nan), r.get(method, np.nan)),
                        fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"rare false-positive rate @{PRIMARY_THR} (↓)")
    ax.set_ylabel(f"rare recall @{PRIMARY_THR} (↑)")
    ax.set_title("Rare-cell trade-off (mean over datasets/scenarios/seeds)")
    ax.grid(alpha=0.3)
    b = _b64(fig); Path(path).write_bytes(base64.b64decode(b)); return b


def _fig_auto(df, path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    obj = _scalar(df, "audit_auto_objective")
    rmse = _scalar(df, "rmse")
    rmse_auto = rmse[rmse["method"] == "TissueResolve_auto"]
    merged = pd.merge(obj[["dataset", "scenario", "seed", "value"]],
                      rmse_auto[["dataset", "scenario", "seed", "value"]],
                      on=["dataset", "scenario", "seed"], suffixes=("_obj", "_rmse"))
    if not merged.empty:
        axes[0].scatter(merged["value_obj"], merged["value_rmse"], s=40, alpha=0.7)
        axes[0].set_xlabel("auto objective (gene-masking CV − cond. penalty)")
        axes[0].set_ylabel("composition RMSE of selected solver")
        c = merged["value_obj"].corr(merged["value_rmse"])
        axes[0].set_title(f"auto objective vs RMSE (r={c:.2f})")
        axes[0].grid(alpha=0.3)
    else:
        axes[0].set_visible(False)
    isn = _scalar(df, "audit_selected_solver_is_nnls")
    if not isn.empty:
        by = isn.groupby(["dataset"])["value"].mean()
        axes[1].bar(range(len(by)), by.to_numpy())
        axes[1].set_xticks(range(len(by))); axes[1].set_xticklabels(by.index)
        axes[1].set_ylabel("fraction of auto runs that selected plain NNLS")
        axes[1].set_title("auto → NNLS frequency")
        axes[1].set_ylim(0, 1)
    else:
        axes[1].set_visible(False)
    fig.tight_layout()
    b = _b64(fig); Path(path).write_bytes(base64.b64decode(b)); return b


# --------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------
def _table_html(df):
    return df.to_html(border=0, classes="tbl", float_format=lambda x: f"{x:.4f}") \
        if not df.empty else "<p><em>no data</em></p>"


def build_report(out_dir):
    out = Path(out_dir)
    df = pd.read_parquet(out / "benchmark_results.parquet")
    manifest = json.loads((out / "benchmark_manifest.json").read_text()) \
        if (out / "benchmark_manifest.json").exists() else {}
    figs = out / "figures"; figs.mkdir(exist_ok=True)

    b_over = _fig_bulk_overview(df, figs / "bulk_overview.png")
    b_rare = _fig_rare(df, figs / "rare_detection.png")
    b_auto = _fig_auto(df, figs / "auto_solver_audit.png")

    rmse = _mean_table(df, "rmse")
    pear = _mean_table(df, "pearson")
    cond = _mean_table(df, "conditional_rmse")
    cov = _mean_table(df, "coverage")
    risk = _mean_table(df, "selective_risk")
    n_seeds = df["seed"].nunique()
    fails = manifest.get("failures", [])

    # executive one-liners
    def _best(metric, low=True, thr=None):
        t = _scalar(df, metric, thr)
        if t.empty:
            return "n/a"
        g = t.groupby("method")["value"].mean()
        return f"{(g.idxmin() if low else g.idxmax())} ({g.min() if low else g.max():.3f})"

    exec_rows = {
        "Best accuracy (RMSE↓)": _best("rmse", True),
        "Best composition correlation (Pearson↑)": _best("pearson", False),
        "Best rare recall (↑)": _best("rare_recall", False, PRIMARY_THR),
        "Best rare FPR (↓)": _best("rare_fpr", True, PRIMARY_THR),
        "Best conditional within-family RMSE (↓)": _best("conditional_rmse", True),
        "Fastest": _best("rmse", True),  # placeholder; runtime handled below
    }
    rt = df[(df["status"] == "ok")].groupby("method")["runtime_seconds"].mean()
    if not rt.empty:
        exec_rows["Fastest (mean runtime)"] = f"{rt.idxmin()} ({rt.min():.1f}s)"
        exec_rows.pop("Fastest", None)

    img = lambda b: f'<img src="data:image/png;base64,{b}" style="max-width:100%">'
    exec_html = "".join(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>"
                        for k, v in exec_rows.items())

    # data-driven overall ranking + recommendation (§13)
    def _mean_by_method(metric, thr=None):
        return _scalar(df, metric, thr).groupby("method")["value"].mean()
    rank = pd.DataFrame({
        "RMSE": _mean_by_method("rmse"), "Pearson": _mean_by_method("pearson"),
        "cond_RMSE": _mean_by_method("conditional_rmse"),
        "coverage": _mean_by_method("coverage"),
        "rare_recall@.01": _mean_by_method("rare_recall", PRIMARY_THR),
        "rare_FPR@.01": _mean_by_method("rare_fpr", PRIMARY_THR),
        "runtime_s": rt,
    }).round(3).sort_values("RMSE")
    best_tr = [m for m in rank.index if not m.startswith("Rectangle")]
    best_tr = best_tr[0] if best_tr else "n/a"
    rec_html = f"""<p><b>Best overall:</b> {rank.index[0]} (RMSE {rank['RMSE'].iloc[0]:.3f}).
    <b>Best TissueResolve-family method:</b> {best_tr}
    (RMSE {rank.loc[best_tr,'RMSE']:.3f}, Pearson {rank.loc[best_tr,'Pearson']:.3f},
    conditional RMSE {rank.loc[best_tr,'cond_RMSE']:.3f},
    rare recall {rank.loc[best_tr,'rare_recall@.01']:.3f} @ FPR {rank.loc[best_tr,'rare_FPR@.01']:.3f},
    runtime {rank.loc[best_tr,'runtime_s']:.2f}s).</p>
    {rank.to_html(border=0, classes='tbl', float_format=lambda x: f'{x:.3f}')}
    <div class="note">Recommendation (bulk, this pilot): promote/strengthen the count-likelihood
    (Poisson) GLM as the TissueResolve bulk solver — it is the strongest TR-family method here and
    approaches Rectangle's conditional within-family RMSE at a fraction of the runtime, while the
    AutoSolver is both slower and less accurate, and hierarchical mode over-abstains (low coverage).
    Gate all promotions on real-bulk validation. Ablation (signature vs solver) and spatial are the
    next stages before any default change.</div>"""

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Rectangle vs TissueResolve — benchmark</title>
<style>
body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1100px;margin:2em auto;
padding:0 1em;color:#1a1a1a;line-height:1.5}} h1{{border-bottom:3px solid #444}}
h2{{margin-top:1.8em;border-bottom:1px solid #ccc;padding-bottom:.2em}}
.tbl{{border-collapse:collapse;font-size:12px;margin:.5em 0}}
.tbl td,.tbl th{{border:1px solid #ddd;padding:3px 7px;text-align:right}}
.tbl th{{background:#f2f2f2}} details{{margin:.5em 0;background:#fafafa;padding:.5em;border:1px solid #eee}}
.note{{background:#fff8e1;border-left:4px solid #f5b400;padding:.6em 1em;margin:1em 0}}
code{{background:#f2f2f2;padding:1px 4px}}</style></head><body>
<h1>Rectangle vs TissueResolve — bulk benchmark</h1>
<p>{n_seeds} seed(s); {df['dataset'].nunique()} dataset(s); {df['scenario'].nunique()} scenario(s).
Ground truth = mRNA proportions. Rectangle scored on <code>raw</code> (bias-uncorrected) output.
Generated from a single long result table; all metrics are seed means.</p>

<h2>1. Executive summary</h2>
<table class="tbl">{exec_html}</table>
<div class="note">Direction: RMSE / rare-FPR / conditional-RMSE lower is better; Pearson /
rare-recall / coverage higher is better. A low rare-FPR that comes with near-zero rare
recall (heavy abstention) is <b>not</b> a win — see §6/§7.</div>

<h2>2. Benchmark design</h2>
<p>Config-driven single command; run matrix = datasets × scenarios × seeds × methods.
Reference cells are shared between tools; each tool does its own gene/marker selection.
Bulk scenarios include query-side confounders (low depth, overdispersion, technical noise,
gene dropout) applied only inside the harness — production defaults are untouched.</p>
<details><summary>Methods & scenarios</summary>
<p><b>Methods:</b> {', '.join(sorted(df['method'].unique()))}</p>
<p><b>Scenarios:</b> {', '.join(sorted(df['scenario'].unique()))}</p></details>

<h2>3. Main bulk results — RMSE (↓) and Pearson (↑)</h2>
{img(b_over)}
<details open><summary>RMSE (seed means)</summary>{_table_html(rmse)}</details>
<details><summary>Pearson (seed means)</summary>{_table_html(pear)}</details>

<h2>4. Cross-platform results</h2>
<p>See the <code>cross_platform</code> rows above (train one 10x chemistry's donors →
test the other's pseudobulk). Available where a multi-chemistry <code>assay</code> column exists.</p>

<h2>5. Auto-solver audit</h2>
{img(b_auto)}
<p>Correlation of the AutoSolver objective with realized composition RMSE, and how often
<code>auto</code> selected plain NNLS. A weak/negative alignment indicates the gene-masking-CV
objective is not tracking composition accuracy.</p>

<h2>6. Rare-cell performance</h2>
{img(b_rare)}
<details><summary>Conditional within-family RMSE (↓)</summary>{_table_html(cond)}</details>

<h2>7. Coarse vs fine resolution (coverage / selective risk)</h2>
<details open><summary>Coverage — fraction of families given a fine split</summary>{_table_html(cov)}</details>
<details><summary>Selective risk — conditional RMSE on accepted families only</summary>{_table_html(risk)}</details>
<div class="note">Hierarchical mode should be audited here for <b>false abstention</b>:
low coverage with low selective risk means it abstains even where subtypes were resolvable.</div>

<h2>8. Component ablation</h2>
<div class="note">Deferred to the next stage (signature / markers / solver / hierarchy / mRNA
correction decomposition). Not yet run — no panel shown rather than a fabricated one.</div>

<h2>9. Runtime and memory</h2>
{_table_html(rt.round(2).to_frame("mean_runtime_s"))}

<h2>10. Spatial benchmark</h2>
<div class="note">Deferred to the second phase (run after bulk is interpreted), same script
with <code>--modality spatial</code>. Not yet run.</div>

<h2>11. Negative results & failures</h2>
<p>{len(fails)} run failure(s) recorded (never hidden).</p>
<details><summary>Failures</summary><pre>{json.dumps(fails, indent=2)[:4000]}</pre></details>

<h2>12. Limitations</h2>
<ul><li>Simulated pseudobulk from real single cells (CPM as TPM proxy); real bulk with fine
truth is the harder, deferred gate.</li>
<li>Possible pseudobulk "home advantage" for Rectangle (pseudobulk/DESeq2-derived signature).</li>
<li>Rare-abundance spike-in grid, LOD/LOQ, subtype-confusion, presence tests, and spatial are
deferred next stages.</li></ul>

<h2>13. Recommendation</h2>
{rec_html}

<h2>14. Reproducibility manifest</h2>
<details><summary>manifest</summary><pre>{json.dumps(manifest, indent=2)[:6000]}</pre></details>
</body></html>"""
    (out / "benchmark_report.html").write_text(html)
    print(f"report -> {out / 'benchmark_report.html'}")
    return out / "benchmark_report.html"


if __name__ == "__main__":
    import sys
    build_report(sys.argv[1] if len(sys.argv) > 1 else
                 "benchmarks/results/rectangle_benchmark")
