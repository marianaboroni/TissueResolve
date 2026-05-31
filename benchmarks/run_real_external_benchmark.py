#!/usr/bin/env python
"""
Real external-method benchmark orchestrator.

Compares TissueResolve against external published tools using harmonised inputs.
External tools run via their own scripts (R/Python) and are recorded with an
honest status; tools that are not installed are **skipped** (never counted as
benchmarked). One failing tool never stops the run.

Modes:
    --dry-run                 list intended tools + environment availability
    --prepare-inputs          write harmonised inputs
    --install-tools           attempt installs, write tool_installation_status.tsv
    --run-bulk / --run-spatial / --run-all
    --fast / --full           fast = reduced iterations for local testing
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from benchmarks.shared.io import OUTPUTS_DIR, write_tsv, write_json  # noqa: E402
from benchmarks.shared import environment as ENV  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
BULK_EXTERNAL = {"MuSiC": ("R", "benchmarks/bulk/methods/run_music.R"),
                 "BayesPrism": ("R", "benchmarks/bulk/methods/run_bayesprism.R"),
                 "BisqueRNA": ("R", "benchmarks/bulk/methods/run_bisque.R")}
SPATIAL_EXTERNAL = {"cell2location": ("py", "benchmarks/spatial/methods/run_cell2location.py"),
                    "CARD": ("R", "benchmarks/spatial/methods/run_card.R"),
                    "SPOTlight": ("R", "benchmarks/spatial/methods/run_spotlight.R")}


def _r_available() -> bool:
    import shutil
    return shutil.which("Rscript") is not None


def _r_pkg_installed(pkg: str) -> bool:
    """Check an R package in the project-local Rlib via Rscript (no rpy2)."""
    if not _r_available():
        return False
    try:
        expr = ('.libPaths(c("benchmarks/envs/Rlib", .libPaths())); '
                f'cat(requireNamespace("{pkg}", quietly=TRUE))')
        p = subprocess.run(["Rscript", "-e", expr], cwd=str(REPO),
                           capture_output=True, text=True, timeout=120)
        return "TRUE" in p.stdout
    except Exception:
        return False


def _invoke(kind: str, script: str, fast: bool) -> dict:
    """Run an external runner script; capture status (never raises)."""
    cmd = (["Rscript", script] if kind == "R"
           else [sys.executable, script] + (["--fast"] if fast else []))
    if kind == "R" and not _r_available():
        return {"status": "skipped", "error": "Rscript not available"}
    try:
        p = subprocess.run(cmd, cwd=str(REPO), capture_output=True, text=True,
                           timeout=1800)
        return {"status": "ran" if p.returncode == 0 else "failed",
                "returncode": p.returncode, "stdout": p.stdout[-500:],
                "stderr": p.stderr[-500:]}
    except Exception as exc:  # noqa: BLE001
        return {"status": "failed", "error": str(exc)}


def _read_metadata(modality: str) -> pd.DataFrame:
    d = OUTPUTS_DIR / modality / "method_metadata"
    rows = []
    if d.exists():
        for p in sorted(d.glob("*.json")):
            try:
                rows.append(json.loads(p.read_text()))
            except Exception:
                pass
    return pd.DataFrame(rows)


def _internal_bulk(prep_summary) -> pd.DataFrame:
    """Compute internal TissueResolve + NNLS bulk predictions/metrics from prepared inputs."""
    import warnings
    from tissueresolve.results import ReferenceSignature
    from tissueresolve.solver import AutoSolver, NNLSSolver
    from benchmarks.shared import metrics as M
    prep = OUTPUTS_DIR / "prepared_inputs"
    sig = pd.read_csv(prep / "reference" / "reference_signature_fine.tsv", sep="\t", index_col=0)
    bulk = pd.read_csv(prep / "bulk" / "bulk_counts_genes_by_samples.tsv", sep="\t", index_col=0)
    bulk.index = bulk.index.map(str)
    truth_p = prep / "bulk" / "bulk_truth_proportions.tsv"
    truth = pd.read_csv(truth_p, sep="\t", index_col=0) if truth_p.exists() else None
    # rebuild a minimal ReferenceSignature from the signature matrix (CPM)
    R = sig.to_numpy(float).T.astype(np.float32)
    ref = ReferenceSignature(gene_names=list(sig.index), cell_types=list(sig.columns),
                             R_cpm=R, R_log=np.log1p(R).astype(np.float32))
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, solver in (("TissueResolve_auto", AutoSolver(n_splits=2)),
                             ("NNLS_baseline", NNLSSolver())):
            t0 = time.perf_counter()
            res = solver.solve(bulk, ref)
            rt = time.perf_counter() - t0
            acc = (M.accuracy_metrics(truth, res.proportions)["pearson"]
                   if truth is not None else np.nan)
            rows.append({"method": name, "modality": "bulk", "status": "executed",
                         "executed": True, "imported": False, "exported_only": False,
                         "failed": False, "runtime_seconds": round(rt, 3),
                         "accuracy": acc})
    return pd.DataFrame(rows).set_index("method")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--prepare-inputs", action="store_true")
    ap.add_argument("--install-tools", action="store_true")
    ap.add_argument("--run-bulk", action="store_true")
    ap.add_argument("--run-spatial", action="store_true")
    ap.add_argument("--run-all", action="store_true")
    ap.add_argument("--use-existing-real-data", action="store_true")
    ap.add_argument("--install-tool", default=None,
                    help="Install one tool (e.g. BisqueRNA) into benchmarks/envs/Rlib")
    ap.add_argument("--methods", default=None,
                    help="Comma-separated subset of external methods to run")
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args(argv)

    if args.install_tool:
        return _install_one(args.install_tool)

    if args.dry_run:
        print("Real external benchmark — plan (dry run):")
        print(f"  R available: {_r_available()}")
        for mod, tools in (("bulk", BULK_EXTERNAL), ("spatial", SPATIAL_EXTERNAL)):
            for name, (kind, script) in tools.items():
                avail = _r_pkg_installed(name) if kind == "R" else ENV.python_module_available(name.lower())
                print(f"  [{mod}] {name:14s} ({kind}) → "
                      f"{'available' if avail else 'NOT installed (will skip)'}")
        print("Internal executed: TissueResolve_auto, TissueResolve_hierarchical, NNLS_baseline")
        print(f"Outputs under {OUTPUTS_DIR}/ ; report: real_external_benchmark_report.html")
        return 0

    if args.install_tools:
        print("Attempting external-tool installation (logs in benchmarks/logs/) …")
        script = REPO / "benchmarks" / "envs" / "install_external_tools.sh"
        try:
            subprocess.run(["bash", str(script)], cwd=str(REPO), timeout=3600)
        except Exception as exc:  # noqa: BLE001
            print(f"  install script error: {exc}")
        _write_install_status()
        return 0

    prep = {}
    if args.prepare_inputs or args.run_bulk or args.run_spatial or args.run_all:
        from benchmarks.shared.prepare_external_inputs import prepare_inputs
        prep = prepare_inputs(use_existing_real_data=args.use_existing_real_data,
                              toy=not args.use_existing_real_data)
        print(f"Prepared inputs: {prep['summary']}")
        if args.prepare_inputs and not (args.run_bulk or args.run_spatial or args.run_all):
            return 0

    methods_filter = set(m.strip() for m in args.methods.split(",")) if args.methods else None
    status_rows = []
    # --- bulk ---
    if args.run_bulk or args.run_all:
        # 1) run external runner scripts (they execute if installed, else skip)
        for name, (kind, script) in BULK_EXTERNAL.items():
            if methods_filter and name not in methods_filter:
                continue
            inv = _invoke(kind, script, args.fast)
            print(f"  [bulk] {name}: {inv['status']}")
        meta = _read_metadata("bulk")
        meta_by = {r["method"]: r for _, r in meta.iterrows()} if not meta.empty else {}
        if args.use_existing_real_data:
            # 2) comprehensive comparison on real data: TissueResolve modes +
            #    baselines + locally-executed/imported externals, full metrics
            from benchmarks.bulk import run_bulk_benchmark as BB
            BB.main(["--use-existing-real-data", "--include-imported"])
            ex = pd.read_csv(OUTPUTS_DIR / "bulk" / "executive_summary.tsv", sep="\t")
            for _, r in ex.iterrows():
                m = r["method"]; st = str(r.get("status", ""))
                is_exec = st == "success"
                md = meta_by.get(m, {})
                status_rows.append({"method": m, "modality": "bulk",
                                    "status": "executed" if is_exec else st,
                                    "executed": is_exec,
                                    "imported": bool(md.get("imported", False)),
                                    "exported_only": st == "exported_not_run",
                                    "failed": st == "failed",
                                    "runtime_seconds": float(r.get("runtime_seconds", 0) or 0),
                                    "accuracy": float(r.get("bulk_fine_pearson"))
                                    if pd.notna(r.get("bulk_fine_pearson")) else np.nan})
            seen = set(ex["method"])
        else:
            # toy / offline: internal compute on the prepared toy inputs
            internal = _internal_bulk(prep)
            for m, r in internal.iterrows():
                status_rows.append({"method": m, **r.to_dict()})
            seen = set(internal.index)
        # include skipped/failed external tools recorded only in metadata
        for name, r in meta_by.items():
            if name not in seen:
                st = str(r.get("status", "skipped"))
                status_rows.append({"method": name, "modality": "bulk", "status": st,
                                    "executed": bool(r.get("executed", False)),
                                    "imported": False, "exported_only": False,
                                    "failed": st == "failed",
                                    "runtime_seconds": float(r.get("runtime_seconds", 0) or 0),
                                    "accuracy": np.nan})
    # --- spatial ---
    if args.run_spatial or args.run_all:
        for name, (kind, script) in SPATIAL_EXTERNAL.items():
            inv = _invoke(kind, script, args.fast)
            print(f"  [spatial] {name}: {inv['status']}")
        meta = _read_metadata("spatial")
        if not meta.empty:
            for _, r in meta.iterrows():
                status_rows.append({"method": r["method"], "modality": "spatial",
                                    "status": r.get("status", "skipped"),
                                    "executed": bool(r.get("executed", False)),
                                    "imported": bool(r.get("imported", False)),
                                    "exported_only": bool(r.get("exported_only", False)),
                                    "failed": str(r.get("status")) == "failed",
                                    "runtime_seconds": float(r.get("runtime_seconds", 0) or 0),
                                    "accuracy": np.nan})

    status = pd.DataFrame(status_rows).drop_duplicates("method").set_index("method") \
        if status_rows else pd.DataFrame()
    if not status.empty:
        write_tsv(status, OUTPUTS_DIR / "real_external_method_status.tsv")
        from benchmarks.shared.composite_score import compute_composite_scores
        comp = compute_composite_scores(status, has_ground_truth=True)
        write_tsv(comp, OUTPUTS_DIR / "composite_scores.tsv")
        _write_report(status, comp)
        print(f"\nReport: {OUTPUTS_DIR / 'real_external_benchmark_report.html'}")
        execed = status[status["executed"]].index.tolist()
        skipped = status[status["status"] == "skipped"].index.tolist()
        print(f"Executed: {execed}")
        print(f"Skipped (not installed): {skipped}")
    return 0


def _install_one(tool: str) -> int:
    """Install a single tool into the project-local R library (logged)."""
    rlib = REPO / "benchmarks" / "envs" / "Rlib"
    rlib.mkdir(parents=True, exist_ok=True)
    logf = REPO / "benchmarks" / "logs" / f"install_{tool.lower()}.log"
    logf.parent.mkdir(parents=True, exist_ok=True)
    print(f"Installing {tool} into {rlib} (log: {logf}) …")
    if tool.lower() in ("bisquerna", "bisque"):
        rexpr = (
            '.libPaths(c("benchmarks/envs/Rlib", .libPaths())); '
            'options(repos="https://cloud.r-project.org"); '
            'for (p in c("quadprog","lpSolve","limSolve")) '
            'if(!requireNamespace(p,quietly=TRUE)) '
            'try(install.packages(p, lib="benchmarks/envs/Rlib", type="source")); '
            'if(!requireNamespace("BisqueRNA",quietly=TRUE)) '
            'try(install.packages('
            '"https://cran.r-project.org/src/contrib/Archive/BisqueRNA/BisqueRNA_1.0.5.tar.gz",'
            ' lib="benchmarks/envs/Rlib", repos=NULL, type="source")); '
            'cat("installed=", requireNamespace("BisqueRNA",quietly=TRUE), "\\n")')
        try:
            with open(logf, "w") as fh:
                subprocess.run(["Rscript", "-e", rexpr], cwd=str(REPO),
                               stdout=fh, stderr=subprocess.STDOUT, timeout=1800)
        except Exception as exc:  # noqa: BLE001
            print(f"  install error: {exc}")
    elif tool.lower() == "music":
        rexpr = (
            '.libPaths(c("benchmarks/envs/Rlib", .libPaths())); '
            'options(repos="https://cloud.r-project.org"); '
            'if(!requireNamespace("remotes",quietly=TRUE)) '
            'install.packages("remotes", lib="benchmarks/envs/Rlib"); '
            'if(!requireNamespace("MuSiC",quietly=TRUE)) '
            'try(remotes::install_github("xuranw/MuSiC", lib="benchmarks/envs/Rlib", '
            'upgrade="never", dependencies=TRUE)); '
            'cat("installed=", requireNamespace("MuSiC",quietly=TRUE), "\\n")')
        try:
            with open(logf, "w") as fh:
                subprocess.run(["Rscript", "-e", rexpr], cwd=str(REPO),
                               stdout=fh, stderr=subprocess.STDOUT, timeout=2400)
        except Exception as exc:  # noqa: BLE001
            print(f"  install error: {exc}")
    else:
        print(f"  {tool}: per-tool install not scripted; see "
              "benchmarks/envs/install_external_tools.sh")
    _write_install_status()
    return 0


def _write_install_status() -> None:
    logdir = REPO / "benchmarks" / "logs"
    rows = []
    tools = [("MuSiC", "bulk", "R"), ("BayesPrism", "bulk", "R"),
             ("BisqueRNA", "bulk", "R"), ("CARD", "spatial", "R"),
             ("SPOTlight", "spatial", "R"), ("cell2location", "spatial", "py")]
    for name, mod, env in tools:
        log = logdir / f"install_{name.lower()}.log"
        installed = (_r_pkg_installed(name) if env == "R"
                     else ENV.python_module_available(name.lower()))
        err = ""
        if log.exists() and not installed:
            err = log.read_text()[-300:].replace("\n", " ")
        rows.append({"tool": name, "modality": mod, "environment": env,
                     "attempted": log.exists(), "installed": installed,
                     "version": "", "install_command": "see install_external_tools.sh",
                     "error_message": err,
                     "install_hint": ENV.INSTALL_HINTS.get(name.lower(), ""),
                     "status": "installed" if installed else "not_installed"})
    write_tsv(pd.DataFrame(rows).set_index("tool"),
              OUTPUTS_DIR / "tool_installation_status.tsv")
    print(f"Wrote {OUTPUTS_DIR / 'tool_installation_status.tsv'}")


def _write_report(status: pd.DataFrame, comp: pd.DataFrame) -> None:
    import html as _h
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    execed = status[status["executed"]].index.tolist()
    imported = status[status.get("imported", False)].index.tolist()
    skipped = status[status["status"] == "skipped"].index.tolist()
    failed = status[status.get("failed", False)].index.tolist()
    ranked = comp[comp["final_score"].notna()].sort_values("final_score", ascending=False)
    best = ranked.index[0] if not ranked.empty else "—"
    css = "body{font-family:Arial,sans-serif;max-width:1000px;margin:0 auto;padding:20px}" \
          "table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ccc;padding:4px 8px}" \
          "th{background:#f0f3f8}.warn{background:#fff4e5;border:1px solid #dd8452;padding:10px}"
    body = [f"<h1>TissueResolve — real external benchmark</h1>",
            "<h2>1. Executive summary</h2>",
            f"<ul><li>Executed: {execed or '—'}</li>"
            f"<li>Imported: {imported or '—'}</li>"
            f"<li>Skipped (not installed): {skipped or '—'}</li>"
            f"<li>Failed: {failed or '—'}</li>"
            f"<li><b>Best overall composite:</b> {best}</li></ul>",
            "<h2>2. Method status</h2>", status.to_html(border=0),
            "<h2>6. Composite score</h2>", comp.to_html(border=0),
            "<h2>7. Best tool by scenario</h2>",
            f"<p>Best overall composite: <b>{best}</b>. Accuracy/runtime per the "
            "composite table; only executed/imported tools are scored.</p>",
            "<h2>8. Metric explanations</h2>",
            "<ul><li><b>Pearson/Spearman</b>: correlation of predicted vs true "
            "proportions.</li><li><b>RMSE/MAE</b>: error magnitude.</li>"
            "<li><b>concordance</b>: agreement between methods (no ground truth).</li>"
            "<li><b>entropy / near-zero fraction</b>: spatial composition spread / sparsity.</li>"
            "<li><b>composite score</b>: weighted accuracy/robustness/usability/"
            "interpretability/resolution-awareness/runtime.</li></ul>",
            "<h2>9. Limitations</h2>",
            "<div class='warn'>External tools were executed only if installed in this "
            "environment; otherwise they are skipped (not benchmarked) with install "
            "hints, or can be imported via import_external_results.py. Real Visium has "
            "no ground truth (concordance/structure only). Fast mode uses reduced "
            "iterations. Imported results may come from a different environment.</div>"]
    (OUTPUTS_DIR / "real_external_benchmark_report.html").write_text(
        f"<!doctype html><html><head><meta charset='utf-8'><style>{css}</style></head>"
        f"<body>{''.join(body)}</body></html>", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
