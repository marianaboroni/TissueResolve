#!/usr/bin/env python
"""
Bulk deconvolution benchmark runner.

Compares TissueResolve (flat + hierarchical) and an internal NNLS baseline
against optional external methods (MuSiC, Bisque, DWLS, CIBERSORTx-export) on
toy synthetic data (with ground truth) or the existing real breast-cancer
pseudobulk.  External tools are skipped gracefully when not installed.

Usage
-----
    python benchmarks/bulk/run_bulk_benchmark.py --dry-run
    python benchmarks/bulk/run_bulk_benchmark.py --toy
    python benchmarks/bulk/run_bulk_benchmark.py --use-existing-real-data
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# allow running as a script
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd  # noqa: E402

from benchmarks.shared import metrics as M  # noqa: E402
from benchmarks.shared import normalization as norm  # noqa: E402
from benchmarks.shared.method_registry import bulk_methods  # noqa: E402
from benchmarks.shared.method_selection import build_compatibility_table  # noqa: E402
from benchmarks.shared.io import OUTPUTS_DIR, write_tsv, write_json  # noqa: E402
from benchmarks.shared import report as R  # noqa: E402
from benchmarks.shared import environment as ENV  # noqa: E402
from benchmarks.shared import plotting as PL  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RBC = REPO / "examples" / "real_breast_cancer"
OUT = OUTPUTS_DIR / "bulk"


def _toy_scenario():
    from benchmarks.shared.synthetic import toy_reference, toy_bulk
    ref, mapping = toy_reference()
    bulk, truth = toy_bulk(ref)
    return {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
            "truth": truth, "normalization_status": "counts",
            "reference_library_type": "scrna_10x_3p", "data_source": "toy_synthetic"}


def _real_scenario():
    from tissueresolve.results import ReferenceSignature
    from tissueresolve.reference.hierarchy import (
        load_hierarchy_mapping, build_cell_type_hierarchy,
    )
    import warnings
    ref_dir = RBC / "outputs" / "reference" / "breast_cancer_reference"
    counts = RBC / "data" / "derived" / "pseudobulk_counts.tsv"
    truth_p = RBC / "data" / "derived" / "pseudobulk_true_proportions.tsv"
    mapping_p = RBC / "config" / "breast_cancer_cell_type_hierarchy.tsv"
    for p in (ref_dir, counts, truth_p, mapping_p):
        if not p.exists():
            raise FileNotFoundError(
                f"required real-data file missing: {p}. Run the breast-cancer "
                "harness (scripts 01–02) first.")
    ref = ReferenceSignature.load(ref_dir)
    bulk = pd.read_csv(counts, sep="\t", index_col=0, comment="#")
    bulk.index = bulk.index.map(str)
    ref_genes = set(map(str, ref.gene_names))
    bulk = bulk.loc[[g for g in bulk.index if g in ref_genes]]
    truth = pd.read_csv(truth_p, sep="\t", index_col=0, comment="#")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mapping = build_cell_type_hierarchy(list(ref.cell_types),
                                            load_hierarchy_mapping(mapping_p))
    return {"reference": ref, "bulk": bulk, "hierarchy_mapping": mapping,
            "truth": truth, "normalization_status": norm.detect_normalization_status(
                bulk.to_numpy()),
            "reference_library_type": "scrna_10x_3p", "data_source": "real_breast_cancer"}


def _high_risk_families(ref, mapping):
    """Families known (from the reference) to be non-separable → abstention is fair."""
    import warnings
    from tissueresolve.reference.hierarchy import evaluate_within_family_resolvability
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = evaluate_within_family_resolvability(ref, mapping,
                                                       min_discriminating_genes=10)
        return {str(f) for f in res.index[~res["resolvable"].astype(bool)]}
    except Exception:
        return set()


def main(argv=None) -> int:
    from benchmarks.shared import analysis as AN

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--toy", action="store_true")
    ap.add_argument("--use-existing-real-data", action="store_true")
    ap.add_argument("--no-external", action="store_true")
    ap.add_argument("--include-imported", action="store_true",
                    help="Also compare external results imported via import_external_results.py")
    args = ap.parse_args(argv)

    methods = bulk_methods(include_external=not args.no_external)
    if args.include_imported:
        from benchmarks.shared.imported import discover_imported
        methods += discover_imported("bulk")

    if args.dry_run:
        print("Bulk benchmark plan (dry run):")
        for m in methods:
            avail = "available" if m.is_available() else f"SKIP ({m.install_hint()})"
            print(f"  - {m.name:28s} external={m.external!s:5} {avail}")
        print(f"Data source: {'real' if args.use_existing_real_data else 'toy'}")
        print(f"Outputs would be written under {OUT}/")
        return 0

    scenario = (_real_scenario() if args.use_existing_real_data else _toy_scenario())
    truth = scenario["truth"]
    mapping = scenario["hierarchy_mapping"]
    high_risk = _high_risk_families(scenario["reference"], mapping)
    print(f"Bulk benchmark on {scenario['data_source']} "
          f"({truth.shape[0]} samples × {truth.shape[1]} cell types). "
          f"High-risk (non-separable) families: {len(high_risk)}.")

    compat = build_compatibility_table(methods, scenario)
    write_tsv(compat, OUT / "method_compatibility.tsv")

    results = []
    fine_by_method, fair_by_method, family_rows, unresolved_rows = {}, {}, [], []
    preds = {}
    for m in methods:
        res = m.run(scenario)
        results.append(res)
        print(f"  {m.name:28s} → {res.status} ({res.runtime_s:.2f}s)")
        if res.status == "success" and res.predictions is not None:
            preds[m.name] = res.predictions
            fine = M.accuracy_metrics(truth, res.predictions)
            fine["dominant_acc"] = M.dominant_accuracy(truth, res.predictions)
            fine_by_method[m.name] = fine
            fair = M.hierarchical_fair_metrics(truth, res.predictions, mapping, high_risk)
            fair_by_method[m.name] = fair
            family_rows.append({"method": m.name, "family_pearson": fair["family_pearson"],
                                "family_rmse": fair["family_rmse"]})
            if fair.get("abstains"):
                unresolved_rows.append({"method": m.name,
                                        "unresolved_mass_fraction": fair["unresolved_mass_fraction"],
                                        "unresolved_precision": fair["unresolved_precision"],
                                        "unresolved_recall": fair["unresolved_recall"],
                                        "n_unresolved_families": fair["n_unresolved_families"]})

    cats = AN.categorize(results)
    exec_table = AN.build_executive_table(results, methods, scenario,
                                          fair_by_method, fine_by_method)
    cap = AN.capability_matrix(methods)
    best = AN.best_methods(exec_table)
    rec = AN.recommendation_by_use_case(exec_table, best)
    n_nontr = sum(1 for r in cats["executed"] if not r.method.startswith("TissueResolve"))
    conclusion = AN.conclusion_text("bulk", exec_table, best, True, n_nontr)

    # tables
    write_tsv(exec_table, OUT / "executive_summary.tsv")
    write_tsv(cap, OUT / "capability_matrix.tsv")
    if family_rows:
        write_tsv(pd.DataFrame(family_rows).set_index("method"),
                  OUT / "family_level_metrics.tsv")
        write_tsv(pd.DataFrame([{"method": k, **v} for k, v in fair_by_method.items()]
                               ).set_index("method"),
                  OUT / "hierarchical_fair_metrics.tsv")
    if unresolved_rows:
        write_tsv(pd.DataFrame(unresolved_rows).set_index("method"),
                  OUT / "unresolved_metrics.tsv")
    write_json({"environment": ENV.package_versions(),
                "normalization_status": scenario["normalization_status"],
                "data_source": scenario["data_source"],
                "best_methods": best, "recommendation_by_use_case": rec,
                "high_risk_families": sorted(high_risk),
                "n_nontissueresolve_executed": n_nontr},
               OUT / "benchmark_metadata.json")

    # figures (each saves .data.tsv)
    figdir = OUT / "figures"
    ex_exec = exec_table[exec_table["executed_or_exported"].isin(
        ["executed", "executed_imported"])].reset_index()
    if not ex_exec.empty:
        PL.bar_figure(ex_exec[["method", "bulk_fine_pearson"]].dropna(),
                      x="method", y="bulk_fine_pearson",
                      title="Bulk fine-level accuracy (executed methods)",
                      path=figdir / "bulk_fine_leaderboard.html")
        PL.bar_figure(ex_exec[["method", "bulk_family_pearson"]].dropna(),
                      x="method", y="bulk_family_pearson",
                      title="Bulk family-level accuracy (executed methods)",
                      path=figdir / "bulk_family_leaderboard.html")
        PL.bar_figure(ex_exec[["method", "runtime_seconds"]],
                      x="method", y="runtime_seconds",
                      title="Runtime by method (executed)",
                      path=figdir / "bulk_runtime.html")
    PL.heatmap_figure(cap.astype(float, errors="ignore").select_dtypes("bool").astype(int)
                      if not cap.empty else cap,
                      title="Method capability matrix",
                      path=figdir / "capability_matrix.html")

    sections = _bulk_sections(cats, exec_table, cap, compat, best, rec, conclusion,
                              fair_by_method, scenario)
    report_path = R.build_modality_report("bulk", sections, OUT / "bulk_benchmark_report.html")
    print(f"Best fine-level: {best.get('best_bulk_fine_pearson')}; "
          f"best family-level: {best.get('best_bulk_family_pearson')}; "
          f"fastest: {best.get('best_runtime')}.")
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


def _bulk_sections(cats, exec_table, cap, compat, best, rec, conclusion,
                   fair_by_method, scenario):
    best_html = "<ul>" + "".join(
        f"<li>{k.replace('_', ' ')}: <b>{v}</b></li>" for k, v in best.items()) + "</ul>"
    rec_html = "<ul>" + "".join(
        f"<li>{k.replace('_', ' ')}: <b>{v}</b></li>" for k, v in rec.items()) + "</ul>"
    fair_df = (pd.DataFrame([{"method": k, **v} for k, v in fair_by_method.items()]
                            ).set_index("method") if fair_by_method else pd.DataFrame())
    return {
        "1. Executive summary": exec_table,
        "2. Which methods were compared?": _status_summary_html(cats),
        "3. Executed vs skipped/exported": _status_summary_html(cats),
        "4. Best method by criterion": best_html,
        "5. Recommended method by use case": rec_html,
        "6. Fair hierarchical / family-level / unresolved-aware metrics": fair_df
        if not fair_df.empty else "<p>No hierarchical method executed.</p>",
        "7. Method capability matrix": cap,
        "8. Method compatibility": compat,
        "9. Benchmark conclusion": f"<p>{conclusion}</p>",
    }


if __name__ == "__main__":
    raise SystemExit(main())
