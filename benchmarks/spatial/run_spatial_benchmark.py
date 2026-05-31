#!/usr/bin/env python
"""
Spatial deconvolution benchmark runner.

Compares TissueResolve (flat + hierarchical) and an internal spot-level NNLS
baseline against optional external methods (RCTD, cell2location, stereoscope,
SPOTlight, Tangram) on toy synthetic data (with ground truth) or the existing
real breast-cancer Visium section (no ground truth → concordance/structure).

Usage
-----
    python benchmarks/spatial/run_spatial_benchmark.py --dry-run
    python benchmarks/spatial/run_spatial_benchmark.py --toy
    python benchmarks/spatial/run_spatial_benchmark.py --use-existing-real-data --max-spots 800
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared.method_registry import spatial_methods  # noqa: E402
from benchmarks.shared.method_selection import build_compatibility_table  # noqa: E402
from benchmarks.shared.io import OUTPUTS_DIR, write_tsv, write_json  # noqa: E402
from benchmarks.shared import report as R  # noqa: E402
from benchmarks.shared import environment as ENV  # noqa: E402
from benchmarks.shared import plotting as PL  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RBC = REPO / "examples" / "real_breast_cancer"
OUT = OUTPUTS_DIR / "spatial"


def _toy_scenario():
    from benchmarks.shared.synthetic import toy_reference, toy_spatial
    ref, mapping = toy_reference()
    sc, truth = toy_spatial(ref)
    sc.update({"reference": ref, "hierarchy_mapping": mapping, "truth": truth,
               "normalization_status": "counts", "has_ground_truth": True,
               "reference_library_type": "scrna_10x_3p", "data_source": "toy_synthetic"})
    return sc


def _real_scenario(max_spots: int, seed: int):
    import anndata as ad
    import scipy.sparse as sp
    import warnings
    from tissueresolve.results import ReferenceSignature
    from tissueresolve.reference.hierarchy import (
        load_hierarchy_mapping, build_cell_type_hierarchy,
    )
    ref_dir = RBC / "outputs" / "reference" / "breast_cancer_reference"
    visium = RBC / "data" / "spatial" / "human_breast_cancer_1.h5ad"
    mapping_p = RBC / "config" / "breast_cancer_cell_type_hierarchy.tsv"
    for p in (ref_dir, visium, mapping_p):
        if not p.exists():
            raise FileNotFoundError(f"required real-data file missing: {p}")
    ref = ReferenceSignature.load(ref_dir)
    adata = ad.read_h5ad(visium)
    n_total = adata.n_obs
    if max_spots and n_total > max_spots:
        rng = np.random.default_rng(seed)
        keep = np.sort(rng.choice(n_total, size=max_spots, replace=False))
        adata = adata[keep].copy()
        print(f"  subsampled {max_spots}/{n_total} spots (seed={seed}; logged).")
    Y = adata.X
    gene_names = [str(g) for g in adata.var_names]
    lib = (np.asarray(Y.sum(axis=1)).ravel() if sp.issparse(Y)
           else Y.sum(axis=1)).astype("float32")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mapping = build_cell_type_hierarchy(list(ref.cell_types),
                                            load_hierarchy_mapping(mapping_p))
    return {"reference": ref, "Y": Y, "gene_names": gene_names,
            "array_row": adata.obs["array_row"].to_numpy(),
            "array_col": adata.obs["array_col"].to_numpy(),
            "lib_sizes": lib, "spot_ids": list(adata.obs_names),
            "hierarchy_mapping": mapping, "truth": None, "has_ground_truth": False,
            "normalization_status": "counts", "n_spots_total": int(n_total),
            "reference_library_type": "scrna_10x_3p", "data_source": "real_breast_cancer"}


def main(argv=None) -> int:
    from benchmarks.shared import analysis as AN

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--toy", action="store_true")
    ap.add_argument("--use-existing-real-data", action="store_true")
    ap.add_argument("--no-external", action="store_true")
    ap.add_argument("--include-imported", action="store_true",
                    help="Also compare external results imported via import_external_results.py")
    ap.add_argument("--max-spots", type=int, default=600)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--spatial-auto", dest="spatial_auto", action="store_true",
                    default=True, help="Run lambda_spatial auto-selection (old vs improved).")
    ap.add_argument("--no-spatial-auto", dest="spatial_auto", action="store_false")
    args = ap.parse_args(argv)

    methods = spatial_methods(include_external=not args.no_external)
    if args.include_imported:
        from benchmarks.shared.imported import discover_imported
        methods += discover_imported("spatial")

    if args.dry_run:
        print("Spatial benchmark plan (dry run):")
        for m in methods:
            avail = "available" if m.is_available() else f"SKIP ({m.install_hint()})"
            print(f"  - {m.name:28s} external={m.external!s:5} {avail}")
        print(f"Data source: {'real' if args.use_existing_real_data else 'toy'}")
        print(f"Outputs would be written under {OUT}/")
        return 0

    scenario = (_real_scenario(args.max_spots, args.seed)
                if args.use_existing_real_data else _toy_scenario())
    has_gt = scenario.get("has_ground_truth", False)
    mapping = scenario.get("hierarchy_mapping", {})
    print(f"Spatial benchmark on {scenario['data_source']} (ground truth: {has_gt}).")

    compat = build_compatibility_table(methods, scenario)
    write_tsv(compat, OUT / "method_compatibility.tsv")

    results, per_method_metrics, fair_by_method, fine_by_method, preds = [], [], {}, {}, {}
    family_rows, unresolved_rows = [], []
    for m in methods:
        res = m.run(scenario)
        results.append(res)
        print(f"  {m.name:28s} → {res.status} ({res.runtime_s:.2f}s)")
        if res.status == "success" and res.predictions is not None:
            preds[m.name] = res.predictions
            if has_gt:
                acc = M.accuracy_metrics(scenario["truth"], res.predictions)
                acc["dominant_acc"] = M.dominant_accuracy(scenario["truth"], res.predictions)
                fine_by_method[m.name] = acc
                per_method_metrics.append({"method": m.name, **acc})
                fair = M.hierarchical_fair_metrics(scenario["truth"], res.predictions,
                                                   mapping, set())
                fair_by_method[m.name] = fair
                family_rows.append({"method": m.name, "family_pearson": fair["family_pearson"],
                                    "family_rmse": fair["family_rmse"]})
            ucols = M.unresolved_columns(res.predictions)
            if ucols:
                tot = res.predictions.to_numpy(float).sum()
                unresolved_rows.append({
                    "method": m.name,
                    "unresolved_mass_fraction": float(
                        res.predictions[ucols].to_numpy(float).sum() / tot) if tot else 0.0,
                    "n_unresolved_families": len(ucols)})

    cats = AN.categorize(results)
    cap = AN.capability_matrix(methods)
    exec_table = AN.build_executive_table(results, methods, scenario,
                                          fair_by_method, fine_by_method)
    status_table = AN.build_status_table(results, methods)
    best = AN.best_methods(exec_table)
    # concordance + structure (needed for best-method summary)
    pmc, struct = {}, pd.DataFrame()
    if len(preds) > 1:
        pmc = M.pairwise_method_correlation(preds)
        struct = pd.DataFrame([
            {"method": k, "mean_entropy": float(M.proportion_entropy(v).mean()),
             "near_zero_fraction": M.near_zero_fraction(v)}
            for k, v in preds.items()]).set_index("method")
    best_summary = AN.best_method_summary(exec_table, fair_by_method=fair_by_method,
                                          concordance=pmc, structure=struct,
                                          modality="spatial")
    interp = AN.interpretation_paragraph("spatial", best, has_gt)
    n_nontr = sum(1 for r in cats["executed"] if not r.method.startswith("TissueResolve"))
    conclusion = AN.conclusion_text("spatial", exec_table, best, has_gt, n_nontr)

    write_tsv(exec_table, OUT / "executive_summary.tsv")
    write_tsv(status_table, OUT / "method_status.tsv")
    write_tsv(best_summary, OUT / "best_method_summary.tsv")
    write_tsv(cap, OUT / "capability_matrix.tsv")
    if family_rows:
        write_tsv(pd.DataFrame(family_rows).set_index("method"), OUT / "family_level_metrics.tsv")
        write_tsv(pd.DataFrame([{"method": k, **v} for k, v in fair_by_method.items()]
                               ).set_index("method"), OUT / "hierarchical_fair_metrics.tsv")
    if unresolved_rows:
        write_tsv(pd.DataFrame(unresolved_rows).set_index("method"), OUT / "unresolved_metrics.tsv")

    sections = {"1. Best-method summary": best_summary,
                "2. Interpretation": f"<p>{interp}</p>",
                "3. Method status (executed / imported / exported-only / skipped)": status_table,
                "4. Which methods were compared?": _status_summary_html(cats),
                "5. Executive summary (per-method)": exec_table,
                "8. Method capability matrix": cap,
                "Method compatibility": compat}
    figdir = OUT / "figures"
    if has_gt and per_method_metrics:
        metrics_df = pd.DataFrame(per_method_metrics).set_index("method")
        write_tsv(metrics_df, OUT / "accuracy_metrics.tsv")
        sections["6. Accuracy (synthetic ground truth)"] = metrics_df
    if len(preds) > 1:
        pmc_df = pd.DataFrame({"pair": list(pmc), "pearson": list(pmc.values())})
        write_tsv(pmc_df.set_index("pair"), OUT / "method_concordance.tsv")
        sections["7. Cross-method concordance (no ground truth)"] = pmc_df.set_index("pair")
        write_tsv(struct, OUT / "spatial_structure.tsv")
        sections["7b. Spatial structure / stability"] = struct
        # concordance heatmap
        names = list(preds)
        mat = pd.DataFrame(1.0, index=names, columns=names)
        for k, v in pmc.items():
            a, b = k.split("__vs__")
            mat.loc[a, b] = mat.loc[b, a] = v
        PL.heatmap_figure(mat, title="Spatial cross-method concordance (Pearson)",
                          path=figdir / "spatial_concordance_heatmap.html")
        PL.bar_figure(struct.reset_index()[["method", "near_zero_fraction"]],
                      x="method", y="near_zero_fraction",
                      title="Spatial stability: near-zero fraction by method",
                      path=figdir / "spatial_near_zero.html")
    # spatial auto-parameter selection: old (default lambda) vs improved (auto)
    if args.spatial_auto and "Y" in scenario:
        try:
            from tissueresolve.spatial.auto_params import (
                select_lambda_spatial, save_spatial_auto_outputs)
            import warnings as _w
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                best_lam, comp_lam, selected = select_lambda_spatial(
                    scenario["Y"], scenario["reference"], scenario["array_row"],
                    scenario["array_col"], scenario["lib_sizes"],
                    scenario["gene_names"], scenario.get("spot_ids"),
                    candidates=(0.0, 0.1, 0.5), max_iter=12)
            save_spatial_auto_outputs(comp_lam, selected, OUT)
            default_lam = 0.1
            selected["default_lambda_spatial"] = default_lam
            sections["10. Spatial auto-parameter selection (old vs improved)"] = (
                f"<p>Default λ_spatial={default_lam} → auto-selected "
                f"<b>λ_spatial={best_lam}</b>. {selected['reason']}</p>"
                + comp_lam.to_html(border=0)
                + f"<p class='note'>Graph diagnostics: {selected['graph_diagnostics']}. "
                "Spatial accuracy is not claimed without ground truth — selection "
                "balances reconstruction against over-smoothing.</p>")
            print(f"Spatial auto: default λ=0.1 → selected λ={best_lam}")
        except Exception as exc:  # noqa: BLE001
            print(f"  (spatial auto-parameter selection skipped: {exc})")

    sections["9. Benchmark conclusion"] = f"<p>{conclusion}</p>"

    PL.heatmap_figure(cap.select_dtypes("bool").astype(int) if not cap.empty else cap,
                      title="Method capability matrix",
                      path=figdir / "capability_matrix.html")
    write_json({"environment": ENV.package_versions(),
                "data_source": scenario["data_source"], "has_ground_truth": has_gt,
                "n_spots_total": scenario.get("n_spots_total"),
                "n_nontissueresolve_executed": n_nontr},
               OUT / "benchmark_metadata.json")

    report_path = R.build_modality_report("spatial", sections,
                                          OUT / "spatial_benchmark_report.html")
    print(f"Wrote {report_path}")
    return 0


def _status_summary_html(cats: dict) -> str:
    def names(key):
        return ", ".join(r.method for r in cats[key]) or "—"
    return (f"<ul><li><b>Executed:</b> {names('executed')}</li>"
            f"<li><b>Imported (executed externally):</b> {names('imported')}</li>"
            f"<li><b>Exported only (NOT benchmarked):</b> {names('exported')}</li>"
            f"<li><b>Skipped (not installed):</b> {names('skipped')}</li>"
            f"<li><b>Failed:</b> {names('failed')}</li></ul>")


if __name__ == "__main__":
    raise SystemExit(main())
