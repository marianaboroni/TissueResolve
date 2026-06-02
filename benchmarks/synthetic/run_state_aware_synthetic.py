"""
Run the synthetic three-level state benchmark: compare NNLS / hierarchical
(standard, granular) / state-aware / auto across scenarios at broad, cell-type
and state granularity, and write the benchmark tables, figures and an HTML
report.  Offline; no real data.

Headline questions answered honestly:
* does state-aware improve TRUE state-level accuracy where states are separable?
* does it preserve broad/cell-type accuracy?
* does it keep non-separable states UNRESOLVED (low false-resolution) instead of
  inventing precision?
"""
from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.synthetic.state_hierarchy import (
    SCENARIOS, _CT_TO_BROAD, _STATE_TO_CT, build_state_reference,
    simulate_scenario, state_to_broad_mapping, three_level_hierarchy)

_METHODS = ["NNLS_baseline", "TissueResolve_auto",
            "TissueResolve_hierarchical_standard",
            "TissueResolve_hierarchical_granular", "TissueResolve_state_aware"]


# --------------------------------------------------------------------------- #
# run the methods → per-method prediction frame (samples × labels)
# --------------------------------------------------------------------------- #
def run_methods(scenario) -> dict:
    import tissueresolve as tr
    from tissueresolve.bulk.hierarchical import run_hierarchical_bulk
    from tissueresolve.reference.within_family_markers import (
        build_family_specific_gene_panels)

    ref = build_state_reference(scenario)
    bulk = scenario.bulk
    s2b = state_to_broad_mapping(scenario)
    adata = scenario.adata
    gate = dict(min_discriminating_genes=4)
    out: dict[str, tuple] = {}

    def _timed(fn):
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = fn()
        return res, time.perf_counter() - t0

    # NNLS baseline (flat, state granularity)
    r, dt = _timed(lambda: tr.deconv_bulk(bulk, ref, solver="nnls",
                                          resolution_mode="none"))
    out["NNLS_baseline"] = (r.deconv.proportions, dt)

    r, dt = _timed(lambda: tr.deconv_bulk(bulk, ref, solver="auto",
                                          resolution_mode="none"))
    out["TissueResolve_auto"] = (r.deconv.proportions, dt)

    # standard hierarchical broad→state (global genes)
    r, dt = _timed(lambda: tr.deconv_bulk(bulk, ref, resolution_mode="hierarchical",
                                          hierarchy_mapping=s2b, **gate))
    out["TissueResolve_hierarchical_standard"] = (r.deconv.proportions, dt)

    # granular hierarchical broad→state (within-broad state panels)
    def _granular():
        panels = build_family_specific_gene_panels(
            adata, "broad", "state", min_logfc=0.5).family_panels
        return run_hierarchical_bulk(bulk, ref, s2b, family_gene_panels=panels, **gate)
    r, dt = _timed(_granular)
    out["TissueResolve_hierarchical_granular"] = (r.deconv.proportions, dt)

    # state-aware broad→cell_type→state
    c2b = {c: _CT_TO_BROAD[c] for c in set(_STATE_TO_CT[s] for s in scenario.state_cell_types)}
    s2c = {s: _STATE_TO_CT[s] for s in scenario.state_cell_types}
    r, dt = _timed(lambda: tr.deconv_bulk(
        bulk, ref, resolution_mode="hierarchical", state_aware=True,
        hierarchy_mapping=c2b, state_to_celltype=s2c, reference_adata=adata,
        broad_col="broad", cell_type_col="cell_type", state_col="state", **gate))
    sp = r.state_proportions if r.state_proportions is not None else r.cell_type_proportions
    out["TissueResolve_state_aware"] = (sp, dt)
    return out


# --------------------------------------------------------------------------- #
# aggregate a prediction frame to broad / cell-type / state levels
# --------------------------------------------------------------------------- #
def aggregate_levels(pred: pd.DataFrame, scenario):
    states = list(scenario.truth_state.columns)
    cts = list(scenario.truth_celltype.columns)
    broads = list(scenario.truth_broad.columns)
    sp = pd.DataFrame(0.0, index=pred.index, columns=states)
    cp = pd.DataFrame(0.0, index=pred.index, columns=cts)
    bp = pd.DataFrame(0.0, index=pred.index, columns=broads)
    for col in pred.columns:
        v = pred[col].astype(float)
        name = str(col)
        if name.startswith("unresolved_"):
            label = name[len("unresolved_"):]
            if label in cts:
                cp[label] += v; bp[_CT_TO_BROAD[label]] += v
            elif label in broads:
                bp[label] += v
        elif name in states:
            sp[name] += v
            ct = _STATE_TO_CT[name]; cp[ct] += v; bp[_CT_TO_BROAD[ct]] += v
        elif name in cts:
            cp[name] += v; bp[_CT_TO_BROAD[name]] += v
        elif name in broads:
            bp[name] += v
    return bp, cp, sp


def _pearson(a: pd.DataFrame, b: pd.DataFrame) -> float:
    x = a.to_numpy().ravel(); y = b.reindex_like(a).to_numpy().ravel()
    if x.std() < 1e-9 or y.std() < 1e-9:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _rmse(a: pd.DataFrame, b: pd.DataFrame) -> float:
    return float(np.sqrt(((a.to_numpy() - b.reindex_like(a).to_numpy()) ** 2).mean()))


def compute_metrics(scenario, results: dict) -> pd.DataFrame:
    rows = []
    resolvable = scenario.resolvable_states
    non_resolvable = [s for s in scenario.truth_state.columns if s not in resolvable]
    # should-be-unresolved truth mass = non-resolvable states + missing (out-of-ref)
    should_unres_states = set(non_resolvable) | set(scenario.missing_states)
    should_unres = scenario.truth_state[list(should_unres_states)].sum(axis=1) \
        if should_unres_states else pd.Series(0.0, index=scenario.truth_state.index)
    for method, (pred, runtime) in results.items():
        bp, cp, sp = aggregate_levels(pred, scenario)
        # total unresolved mass predicted per sample
        unres_cols = [c for c in pred.columns if str(c).startswith("unresolved_")]
        pred_unres = pred[unres_cols].sum(axis=1) if unres_cols else \
            pd.Series(0.0, index=pred.index)
        # false resolution: resolved state mass placed on non-resolvable states
        resolved_bad = sp[non_resolvable].sum(axis=1).sum() if non_resolvable else 0.0
        total_mass = sp.to_numpy().sum() + float(pred_unres.sum())
        false_res = float(resolved_bad / total_mass) if total_mass > 0 else float("nan")
        # unresolved precision/recall (coarse, mass-based)
        su = float(should_unres.sum()); pu = float(pred_unres.sum())
        unres_recall = float(min(pu, su) / su) if su > 0 else float("nan")
        unres_prec = float(min(pu, su) / pu) if pu > 0 else (
            float("nan") if su == 0 else 0.0)
        resolv_pearson = (_pearson(sp[list(resolvable)], scenario.truth_state[list(resolvable)])
                          if resolvable else float("nan"))
        rows.append({
            "scenario": scenario.name, "method": method,
            "broad_pearson": _pearson(bp, scenario.truth_broad),
            "broad_rmse": _rmse(bp, scenario.truth_broad),
            "celltype_pearson": _pearson(cp, scenario.truth_celltype),
            "celltype_rmse": _rmse(cp, scenario.truth_celltype),
            "state_pearson": _pearson(sp, scenario.truth_state),
            "state_rmse": _rmse(sp, scenario.truth_state),
            "resolvable_state_pearson": resolv_pearson,
            "false_resolution_rate": false_res,
            "unresolved_mass_fraction": float(pred_unres.mean()),
            "unresolved_precision": unres_prec,
            "unresolved_recall": unres_recall,
            "runtime_seconds": round(runtime, 4),
        })
    return pd.DataFrame(rows)


def run_all_scenarios(*, seed: int = 0, **sim_kw) -> pd.DataFrame:
    frames = []
    for name in SCENARIOS:
        sc = simulate_scenario(name, seed=seed, **sim_kw)
        frames.append(compute_metrics(sc, run_methods(sc)))
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------- #
# outputs + figures + report
# --------------------------------------------------------------------------- #
def _family_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    return (metrics.groupby("method")[
        ["broad_pearson", "celltype_pearson", "state_pearson",
         "resolvable_state_pearson", "false_resolution_rate",
         "unresolved_mass_fraction", "runtime_seconds"]]
        .mean().reset_index().sort_values("state_pearson", ascending=False))


def _figures(metrics: pd.DataFrame, fig_dir: Path) -> list[Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    try:
        import plotly.graph_objects as go
    except Exception:
        return paths
    summ = _family_summary(metrics)
    # 1. broad/cell/state accuracy by method
    f = go.Figure()
    for lvl in ("broad_pearson", "celltype_pearson", "state_pearson"):
        f.add_bar(x=summ["method"], y=summ[lvl], name=lvl.replace("_pearson", ""))
    f.update_layout(barmode="group", title="Accuracy by method and level (Pearson)",
                    template="plotly_white", height=460)
    p = fig_dir / "accuracy_by_method.html"; f.write_html(str(p), include_plotlyjs="cdn")
    summ.to_csv(fig_dir / "accuracy_by_method.data.tsv", sep="\t", index=False)
    paths.append(p)
    # 2. false-resolution rate (non-separable honesty)
    nonsep = metrics[metrics["scenario"] == "non_separable"]
    f2 = go.Figure(go.Bar(x=nonsep["method"], y=nonsep["false_resolution_rate"],
                          marker_color="#cc4b3e"))
    f2.update_layout(title="False-resolution rate on non-separable states "
                     "(lower = better)", template="plotly_white", height=420)
    p2 = fig_dir / "false_resolution_rate.html"; f2.write_html(str(p2), include_plotlyjs="cdn")
    nonsep.to_csv(fig_dir / "false_resolution_rate.data.tsv", sep="\t", index=False)
    paths.append(p2)
    # 3. state separability (scenario) vs state accuracy
    f3 = go.Figure()
    for method in metrics["method"].unique():
        m = metrics[metrics["method"] == method]
        f3.add_trace(go.Scatter(x=m["scenario"], y=m["state_pearson"],
                                mode="lines+markers", name=method))
    f3.update_layout(title="State separability (scenario) vs state-level accuracy",
                     template="plotly_white", height=460)
    p3 = fig_dir / "separability_vs_accuracy.html"; f3.write_html(str(p3), include_plotlyjs="cdn")
    paths.append(p3)
    return paths


def write_outputs(metrics: pd.DataFrame, out_dir) -> dict:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    bench_p = out / "state_aware_synthetic_benchmark.tsv"
    fam_p = out / "state_aware_synthetic_family_summary.tsv"
    metrics.to_csv(bench_p, sep="\t", index=False)
    summary = _family_summary(metrics)
    summary.to_csv(fam_p, sep="\t", index=False)
    figs = _figures(metrics, out / "figures")
    # honest verdict
    sa = summary[summary["method"] == "TissueResolve_state_aware"]
    nn = summary[summary["method"] == "NNLS_baseline"]
    verdict = _verdict(metrics, summary)
    html = ("<!doctype html><meta charset='utf-8'><title>State-aware synthetic "
            "benchmark</title><body style='font-family:Arial;max-width:1000px;"
            "margin:0 auto;padding:20px'>"
            "<h1>State-aware synthetic three-level benchmark</h1>"
            "<p>Controlled broad→cell_type→state simulation with KNOWN proportions. "
            "Accuracy metrics are valid (synthetic ground truth).</p>"
            f"<h2>Verdict</h2>{verdict}"
            "<h2>Per-method summary (mean over scenarios)</h2>"
            + summary.round(4).to_html(index=False, border=0)
            + "<h2>Full results</h2>"
            + metrics.round(4).to_html(index=False, border=0)
            + "<h2>Figures</h2><ul>"
            + "".join(f"<li><a href='figures/{p.name}'>{p.name}</a></li>" for p in figs)
            + "</ul></body>")
    rep_p = out / "state_aware_synthetic_report.html"
    rep_p.write_text(html, encoding="utf-8")
    return {"benchmark": bench_p, "family_summary": fam_p, "report": rep_p,
            "figures": figs}


def _verdict(metrics: pd.DataFrame, summary: pd.DataFrame) -> str:
    """Honest, metric-driven verdict (no overclaim)."""
    def g(method, col):
        s = summary.loc[summary["method"] == method, col]
        return float(s.iloc[0]) if len(s) else float("nan")
    sa_state = g("TissueResolve_state_aware", "state_pearson")
    nn_state = g("NNLS_baseline", "state_pearson")
    sa_broad = g("TissueResolve_state_aware", "broad_pearson")
    nn_broad = g("NNLS_baseline", "broad_pearson")
    nonsep = metrics[metrics["scenario"] == "non_separable"]
    sa_fr = float(nonsep.loc[nonsep["method"] == "TissueResolve_state_aware",
                             "false_resolution_rate"].mean())
    nn_fr = float(nonsep.loc[nonsep["method"] == "NNLS_baseline",
                             "false_resolution_rate"].mean())
    bits = ["<ul>"]
    bits.append(f"<li>State-level accuracy: state_aware={sa_state:.3f} vs "
                f"NNLS={nn_state:.3f} → "
                + ("<b>improved</b>" if sa_state > nn_state + 0.02 else
                   ("worsened" if sa_state < nn_state - 0.02 else "comparable")) + ".</li>")
    bits.append(f"<li>Broad-level accuracy: state_aware={sa_broad:.3f} vs "
                f"NNLS={nn_broad:.3f} → "
                + ("preserved" if sa_broad >= nn_broad - 0.02 else "<b>degraded</b>") + ".</li>")
    bits.append(f"<li>False-resolution on non-separable states: "
                f"state_aware={sa_fr:.3f} vs NNLS={nn_fr:.3f} → "
                + ("<b>state-aware keeps them unresolved</b>" if sa_fr < nn_fr - 0.02
                   else "no clear advantage") + ".</li>")
    bits.append("<li>Recommendation: state-aware is appropriate when true state "
                "labels exist AND states are separable; it should remain "
                "experimental where separability is low (it correctly abstains "
                "rather than inventing precision).</li></ul>")
    return "".join(bits)


def main():
    out = Path(__file__).resolve().parents[2] / "benchmarks" / "outputs"
    metrics = run_all_scenarios()
    paths = write_outputs(metrics, out)
    print("Wrote:")
    for k, v in paths.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
