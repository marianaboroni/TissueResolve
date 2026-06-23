"""Gold-truth external bulk benchmark for TissueResolve.

This script runs external bulk deconvolution tools on the same donor-held-out
pseudobulk mixtures used by the internal gold-truth benchmark. It writes
harmonized inputs, records tool status, scores executed external predictions
with the same metrics as the internal benchmark, and combines internal +
external results.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import nnls

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import benchmarks.diagnostics.gold_truth_performance_benchmark as gtp
from benchmarks.diagnostics.gold_truth_performance_benchmark import (
    DATASETS,
    SCENARIOS,
    DatasetContext,
    MixtureBundle,
    aggregate_truth_to_broad,
    concat_full_overlap,
    conditional_pred_from_estimate,
    conditional_truth_from_fine,
    generate_all_mixtures,
    prepare_dataset,
    realize_mixtures,
    score_bulk_predictions,
    _fine_for_scoring,
    _broad_for_scoring,
    _merge_rows,
    _normalise_rows,
    write_report,
)

# Internal gold-truth benchmark outputs (truth + internal TissueResolve metrics).
GOLD_TRUTH_PERF_DIR = REPO / "benchmarks" / "outputs" / "gold_truth_performance"

# Keys returned by score_bulk_predictions().
SCORE_KEYS = ["broad", "fine", "conditional", "rare", "spillover",
              "resolution", "runtime", "family", "subtype"]

# External tool runner scripts in the existing harness.
EXTERNAL_TOOL_SCRIPTS = {
    "MuSiC": {
        "modality": "bulk",
        "language": "R",
        "script": REPO / "benchmarks" / "bulk" / "methods" / "run_music.R",
        "package_or_command": "MuSiC",
        "outputs_cell_fraction_or_rna_fraction": "RNA_fraction",
        "requires_license_or_web": "no",
    },
    "BisqueRNA": {
        "modality": "bulk",
        "language": "R",
        "script": REPO / "benchmarks" / "bulk" / "methods" / "run_bisque.R",
        "package_or_command": "BisqueRNA",
        "outputs_cell_fraction_or_rna_fraction": "RNA_fraction",
        "requires_license_or_web": "no",
    },
    "BayesPrism": {
        "modality": "bulk",
        "language": "R",
        "script": REPO / "benchmarks" / "bulk" / "methods" / "run_bayesprism.R",
        "package_or_command": "BayesPrism",
        "outputs_cell_fraction_or_rna_fraction": "RNA_fraction",
        "requires_license_or_web": "no",
    },
}

OUTPUT_ROOT = REPO / "benchmarks" / "outputs" / "gold_truth_external"
DEFAULT_OUTPUT_ROOT = REPO / "benchmarks" / "outputs" / "gold_truth_external"
INPUT_ROOT = REPO / "benchmarks" / "external_tools" / "inputs" / "gold_truth_bulk"
RUN_MANIFEST_PATH = REPO / "benchmarks" / "external_tools" / "run_manifest.tsv"
# Curated catalog (human-maintained); intentionally NOT written by this script.
CURATED_TOOL_REGISTRY = REPO / "benchmarks" / "external_tools" / "tool_registry.tsv"
PREPARED_INPUTS_ROOT = REPO / "benchmarks" / "outputs" / "prepared_inputs"
BULK_PREDICTIONS_ROOT = REPO / "benchmarks" / "outputs" / "bulk" / "predictions"
BULK_METADATA_ROOT = REPO / "benchmarks" / "outputs" / "bulk" / "method_metadata"

# -----------------------------------------------------------------------------
# Data export / harmonization helpers
# -----------------------------------------------------------------------------


@dataclass
class ExternalDatasetInputs:
    dataset: str
    reference_counts: Path
    reference_cell_metadata: Path
    bulk_counts: Path
    bulk_sample_metadata: Path
    fine_truth: Path
    broad_truth: Path
    conditional_truth: Path
    cell_fraction_truth: Path
    selected_genes: Path
    label_mapping: Path


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_tsv(df: pd.DataFrame, path: Path, index: bool = True) -> Path:
    path = Path(path)
    _ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=index)
    return path


def _write_tsv_no_index(df: pd.DataFrame, path: Path) -> Path:
    path = Path(path)
    _ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False)
    return path


def _prepare_dataset_inputs(ctx: DatasetContext, bundle: MixtureBundle) -> ExternalDatasetInputs:
    out_dir = INPUT_ROOT / ctx.spec.name
    ref_dir = out_dir / "reference"
    bulk_dir = out_dir / "bulk"
    truth_dir = out_dir / "truth"
    for d in (ref_dir, bulk_dir, truth_dir):
        _ensure_dir(d)

    gene_idx = np.array([ctx.adata.var_names.get_loc(g) for g in ctx.selected_genes], dtype=int)
    ref_rows = ctx.train_cells
    X = ctx.adata.X[ref_rows[:, None], gene_idx]
    if hasattr(X, "toarray"):
        X = np.asarray(X.toarray(), dtype=float)
    else:
        X = np.asarray(X, dtype=float)
    gene_names = [ctx.selected_genes[i] for i in range(len(gene_idx))]
    cell_ids = [f"cell_{ctx.spec.name}_{i:05d}" for i in range(len(ref_rows))]
    ref_counts = pd.DataFrame(X.T, index=gene_names, columns=cell_ids).astype(int)
    _write_tsv(ref_counts, ref_dir / "reference_counts_genes_by_cells.tsv")

    obs = ctx.adata.obs.iloc[ref_rows]
    cell_meta = pd.DataFrame({
        "cell_id": cell_ids,
        "cellType": obs["__fine__"].astype(str).to_numpy(),
        "SubjectName": obs["__donor__"].astype(str).to_numpy(),
    })
    _write_tsv_no_index(cell_meta, ref_dir / "reference_cell_metadata.tsv")

    _write_tsv(bundle.counts, bulk_dir / "bulk_counts_genes_by_samples.tsv")

    sample_meta = bundle.metadata.copy()
    sample_meta = sample_meta.rename_axis("sample").reset_index()
    sample_meta = sample_meta.assign(dataset=ctx.spec.name)
    _write_tsv_no_index(sample_meta, bulk_dir / "bulk_sample_metadata.tsv")

    _write_tsv(bundle.fine_truth, truth_dir / "fine_truth.tsv")
    _write_tsv(bundle.broad_truth, truth_dir / "broad_truth.tsv")
    _write_tsv(bundle.conditional_truth, truth_dir / "conditional_truth.tsv")
    _write_tsv(bundle.cell_fraction_truth, truth_dir / "cell_fraction_truth.tsv")

    selected_genes_df = pd.DataFrame({"gene_id": ctx.selected_genes})
    _write_tsv_no_index(selected_genes_df, out_dir / "selected_genes.tsv")

    label_map = pd.DataFrame({
        "dataset": ctx.spec.name,
        "fine_cell_type": list(ctx.mapping.keys()),
        "broad_family": [ctx.mapping[f] for f in ctx.mapping.keys()],
    })
    _write_tsv_no_index(label_map, out_dir / "label_mapping.tsv")

    return ExternalDatasetInputs(
        dataset=ctx.spec.name,
        reference_counts=ref_dir / "reference_counts_genes_by_cells.tsv",
        reference_cell_metadata=ref_dir / "reference_cell_metadata.tsv",
        bulk_counts=bulk_dir / "bulk_counts_genes_by_samples.tsv",
        bulk_sample_metadata=bulk_dir / "bulk_sample_metadata.tsv",
        fine_truth=truth_dir / "fine_truth.tsv",
        broad_truth=truth_dir / "broad_truth.tsv",
        conditional_truth=truth_dir / "conditional_truth.tsv",
        cell_fraction_truth=truth_dir / "cell_fraction_truth.tsv",
        selected_genes=out_dir / "selected_genes.tsv",
        label_mapping=out_dir / "label_mapping.tsv",
    )


# -----------------------------------------------------------------------------
# Tool runtime and adapter helpers
# -----------------------------------------------------------------------------


@dataclass
class ToolRunResult:
    tool_name: str
    dataset: str
    status: str
    executed_this_session: bool
    runtime_seconds: float | None
    run_status: str
    version: str | None
    reason: str
    output_path: Path | None
    metadata_path: Path | None
    error_log_path: Path | None
    estimate_type: str


def _prepare_runner_inputs(dataset_inputs: ExternalDatasetInputs) -> None:
    prepared_ref = PREPARED_INPUTS_ROOT / "reference"
    prepared_bulk = PREPARED_INPUTS_ROOT / "bulk"
    if prepared_ref.exists():
        shutil.rmtree(prepared_ref)
    if prepared_bulk.exists():
        shutil.rmtree(prepared_bulk)
    shutil.copytree(dataset_inputs.reference_counts.parent, prepared_ref)
    shutil.copytree(dataset_inputs.bulk_counts.parent, prepared_bulk)


def _move_prediction_output(tool_name: str, dataset: str) -> tuple[Path | None, Path | None]:
    out_pred = BULK_PREDICTIONS_ROOT / f"{tool_name}.tsv"
    out_meta = BULK_METADATA_ROOT / f"{tool_name}.json"
    if not out_pred.exists() and not out_meta.exists():
        return None, None
    dataset_pred = OUTPUT_ROOT / "raw_predictions" / f"{tool_name}__{dataset}.tsv"
    dataset_meta = OUTPUT_ROOT / "raw_predictions" / f"{tool_name}__{dataset}.json"
    _ensure_dir(dataset_pred.parent)
    if out_pred.exists():
        shutil.move(str(out_pred), dataset_pred)
    if out_meta.exists():
        shutil.move(str(out_meta), dataset_meta)
    return (dataset_pred if dataset_pred.exists() else None,
            dataset_meta if dataset_meta.exists() else None)


def _tool_is_available(tool_name: str) -> bool:
    if EXTERNAL_TOOL_SCRIPTS[tool_name]["language"] == "R":
        return shutil.which("Rscript") is not None
    return False


def _run_external_tool(tool_name: str, dataset: str) -> ToolRunResult:
    tool = EXTERNAL_TOOL_SCRIPTS[tool_name]
    status = "unavailable"
    executed = False
    runtime = None
    output_path = None
    metadata_path = None
    error_log_path = None
    reason = ""
    version = None
    est_type = tool["outputs_cell_fraction_or_rna_fraction"]
    if not _tool_is_available(tool_name):
        return ToolRunResult(tool_name, dataset, "skipped", False, None, "skipped",
                             None, "Rscript unavailable", None, None, None, est_type)
    script = tool["script"]
    if not script.exists():
        return ToolRunResult(tool_name, dataset, "skipped", False, None, "skipped",
                             None, f"missing runner script {script}", None, None, None,
                             est_type)
    cmd = ["Rscript", str(script)]
    t0 = time.perf_counter()
    log_path = OUTPUT_ROOT / "tool_logs" / f"{tool_name}__{dataset}.log"
    _ensure_dir(log_path.parent)
    try:
        proc = subprocess.run(cmd, cwd=str(REPO), capture_output=True, text=True, timeout=5400)
        runtime = float(time.perf_counter() - t0)
        if proc.returncode != 0:
            status = "failed"
            reason = proc.stderr.strip()[:1000] or proc.stdout.strip()[:1000] or "non-zero return code"
        else:
            status = "executed"
            executed = True
            reason = proc.stdout.strip()[:1000]
    except subprocess.TimeoutExpired as exc:
        runtime = float(time.perf_counter() - t0)
        status = "timeout"
        reason = str(exc)
    except Exception as exc:
        runtime = float(time.perf_counter() - t0)
        status = "failed"
        reason = str(exc)
    with log_path.open("w", encoding="utf-8") as fh:
        fh.write("stdout:\n")
        fh.write(proc.stdout if 'proc' in locals() else "")
        fh.write("\nstderr:\n")
        fh.write(proc.stderr if 'proc' in locals() else "")
    error_log_path = log_path
    pred_path, meta_path = _move_prediction_output(tool_name, dataset)
    # R runners fail gracefully with exit 0; trust their metadata JSON for the
    # true status (executed / skipped / failed / exported_only) and version.
    estimate_type = tool["outputs_cell_fraction_or_rna_fraction"]
    if meta_path and meta_path.exists():
        try:
            rmeta = json.loads(meta_path.read_text(encoding="utf-8"))
            r_status = str(rmeta.get("status", status))
            if r_status:
                status = r_status
                executed = (r_status == "executed")
            version = rmeta.get("version") or version
            if rmeta.get("error_message"):
                reason = str(rmeta["error_message"])[:1000]
        except Exception:  # noqa: BLE001
            pass
    if status == "executed" and pred_path is None:
        status = "failed"
        executed = False
        reason = reason or "executed but no prediction file produced"
    return ToolRunResult(
        tool_name=tool_name,
        dataset=dataset,
        status=status,
        executed_this_session=executed,
        runtime_seconds=runtime,
        run_status=status,
        version=version,
        reason=reason,
        output_path=pred_path if status == "executed" else pred_path,
        metadata_path=meta_path,
        error_log_path=error_log_path,
        estimate_type=estimate_type,
    )


# -----------------------------------------------------------------------------
# Prediction normalization and scoring
# -----------------------------------------------------------------------------


def _normalize_predictions(pred: pd.DataFrame) -> pd.DataFrame:
    pred = pred.copy().astype(float).clip(lower=0.0)
    if pred.shape[1] == 0:
        return pred
    pred = _normalise_rows(pred)
    return pred


def _align_external_predictions(pred: pd.DataFrame, truth: MixtureBundle) -> pd.DataFrame:
    if pred.index.name is None:
        pred.index = pred.index.map(str)
    pred = pred.reindex(index=truth.fine_truth.index).fillna(0.0)
    if truth.fine_truth.columns.difference(pred.columns).any():
        missing = truth.fine_truth.columns.difference(pred.columns).tolist()
        pred = pred.reindex(columns=list(pred.columns) + missing, fill_value=0.0)
    return pred


def _read_prediction_file(path: Path) -> pd.DataFrame:
    if not path or not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.map(str)
    df.columns = df.columns.map(str)
    return df


def score_external_predictions(ctx: DatasetContext, bundle: MixtureBundle,
                               tool_name: str, pred: pd.DataFrame,
                               meta: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    pred = _normalize_predictions(pred)
    pred = _align_external_predictions(pred, bundle)
    return score_bulk_predictions(ctx, pred, meta, bundle, tool_name)


# -----------------------------------------------------------------------------
# NNLS external control (independent of the TissueResolve solver)
# -----------------------------------------------------------------------------


def _reference_signature(ctx: DatasetContext) -> pd.DataFrame:
    """Genes x fine cell-type CPM signature from the same reference used internally.

    Built directly from ReferenceSignature.as_R_cpm() so the control shares the
    reference profiles but NOT the TissueResolve weighted-NNLS / hierarchical
    solver. Columns are restricted to ctx.selected_fine in order.
    """
    ref = ctx.reference
    cpm = np.asarray(ref.as_R_cpm(), dtype=float)  # K x G
    sig = pd.DataFrame(cpm.T, index=list(ref.gene_names), columns=list(ref.cell_types))
    cols = [c for c in ctx.selected_fine if c in sig.columns]
    return sig[cols]


def run_nnls_control(ctx: DatasetContext, bundle: MixtureBundle) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Plain non-negative least squares against the shared reference signature.

    Returns RNA-derived proportions (samples x fine cell types), row-normalised.
    """
    t0 = time.perf_counter()
    sig = _reference_signature(ctx)
    counts = bundle.counts  # genes x samples
    shared = [g for g in sig.index if g in counts.index]
    S = sig.loc[shared].to_numpy(dtype=float)
    Y = counts.loc[shared].to_numpy(dtype=float)  # genes x samples
    sample_ids = list(counts.columns)
    out = np.zeros((len(sample_ids), S.shape[1]), dtype=float)
    for j in range(Y.shape[1]):
        w, _ = nnls(S, Y[:, j])
        out[j] = w
    pred = pd.DataFrame(out, index=sample_ids, columns=list(sig.columns))
    runtime = float(time.perf_counter() - t0)
    meta = {
        "runtime_seconds": runtime,
        "status": "executed",
        "n_shared_genes": len(shared),
        "estimate_type": "RNA_fraction",
    }
    return pred, meta


# -----------------------------------------------------------------------------
# Output assembly
# -----------------------------------------------------------------------------


def _write_status_manifest(rows: list[dict[str, Any]]) -> None:
    df = pd.DataFrame(rows)
    # Public manifest schema uses a `tool` column (see tests/test_benchmark_outputs.py);
    # keep the richer `tool_name`/status fields alongside it.
    if "tool_name" in df.columns and "tool" not in df.columns:
        df.insert(0, "tool", df["tool_name"])
    # Always write a per-run copy under the run's output dir.
    _write_tsv(df, OUTPUT_ROOT / "run_manifest.tsv", index=False)
    # Only a full run (default output dir) refreshes the canonical versioned
    # manifest; subset/scratch runs must not clobber it.
    if OUTPUT_ROOT == DEFAULT_OUTPUT_ROOT:
        _write_tsv(df, RUN_MANIFEST_PATH, index=False)


def _write_tool_registry(rows: list[dict[str, Any]]) -> None:
    # Write a per-run registry under the outputs dir. The shared curated catalog
    # benchmarks/external_tools/tool_registry.tsv is human-maintained and is NOT
    # overwritten here (it documents the full tool catalog incl. spatial tools).
    df = pd.DataFrame(rows)
    _write_tsv(df, OUTPUT_ROOT / "external_tool_registry.tsv", index=False)


def _combine_internal_and_external(metrics_file: Path, external_df: pd.DataFrame,
                                   mode_col: str = "mode") -> pd.DataFrame:
    if not metrics_file.exists():
        return external_df.copy()
    internal = pd.read_csv(metrics_file, sep="\t")
    combined = pd.concat([internal, external_df], ignore_index=True, sort=False)
    if mode_col not in combined.columns and "method" in combined.columns:
        combined = combined.rename(columns={"method": mode_col})
    return combined


def _write_summary_tables(rows: dict[str, list[dict[str, Any]]]) -> None:
    # map: output filename -> rows[] key returned by score_bulk_predictions
    tables = {
        "broad_metrics.tsv": "broad",
        "fine_metrics.tsv": "fine",
        "conditional_family_metrics.tsv": "conditional",
        "rare_subtype_metrics.tsv": "rare",
        "spillover_metrics.tsv": "spillover",
        "resolution_metrics.tsv": "resolution",
        "runtime_metrics.tsv": "runtime",
        "family_summary.tsv": "family",
        "subtype_summary.tsv": "subtype",
    }
    for name, key in tables.items():
        df = pd.DataFrame(rows.get(key, []))
        _write_tsv(df, OUTPUT_ROOT / name, index=False)


def _safe_read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


# Internal TissueResolve mode -> explicit combined method name (STEP 10).
INTERNAL_MODE_RENAME = {
    "flat": "TissueResolve_flat",
    "hierarchical_soft": "TissueResolve_hierarchical_soft",
    "auto": "TissueResolve_auto",
    "broad_only": "TissueResolve_broad_only",
    "hierarchical_hard_legacy": "TissueResolve_hard_legacy",
    "hierarchical_ungated_diagnostic": "TissueResolve_ungated_diagnostic",
}


def _rename_internal_modes(df: pd.DataFrame) -> pd.DataFrame:
    if "mode" in df.columns:
        df = df.copy()
        df["mode"] = df["mode"].map(lambda m: INTERNAL_MODE_RENAME.get(m, m))
    return df


def _write_combined_tables(external_files: dict[str, Path]) -> None:
    combined = {
        "combined_broad_metrics.tsv": (GOLD_TRUTH_PERF_DIR / "broad_metrics.tsv", external_files["broad_metrics"]),
        "combined_fine_metrics.tsv": (GOLD_TRUTH_PERF_DIR / "fine_metrics.tsv", external_files["fine_metrics"]),
        "combined_conditional_family_metrics.tsv": (GOLD_TRUTH_PERF_DIR / "conditional_family_metrics.tsv", external_files["conditional_family_metrics"]),
        "combined_runtime_metrics.tsv": (GOLD_TRUTH_PERF_DIR / "runtime_metrics.tsv", external_files["runtime_metrics"]),
    }
    for out_name, (internal_path, external_path) in combined.items():
        internal = _rename_internal_modes(_safe_read_tsv(internal_path))
        external = _safe_read_tsv(external_path)
        df = pd.concat([internal, external], ignore_index=True, sort=False)
        _write_tsv(df, OUTPUT_ROOT / out_name, index=False)


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------


# Tool selection keys for the --tools flag (lower-case).
TOOL_KEYS = {
    "nnls": "NNLS_external_control",
    "music": "MuSiC",
    "bisque": "BisqueRNA",
    "bayesprism": "BayesPrism",
}


def _load_internal_run_params() -> dict[str, Any]:
    """Read the internal gold-truth run metadata so we regenerate IDENTICAL
    mixtures (same seed / samples / cells / reference downsampling). Falls back
    to the parent module defaults if the metadata file is absent."""
    meta_path = GOLD_TRUTH_PERF_DIR / "run_metadata.json"
    defaults = {"seed": 2026, "samples_per_scenario": 1,
                "cells_per_mixture": 180, "max_ref_cells_per_type": 120}
    if not meta_path.exists():
        return defaults
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return defaults
    return {k: meta.get(k, defaults[k]) for k in defaults}


def _limited_scenarios(max_scenarios: int | None):
    """Context manager that temporarily limits the scenario list (smoke tests).

    Slicing the first N scenarios preserves per-scenario seeds (seed offsets are
    keyed on enumerate index), so the generated mixtures match the internal
    benchmark for those scenarios."""
    class _Ctx:
        def __enter__(self):
            self._orig = list(gtp.SCENARIOS)
            if max_scenarios is not None:
                gtp.SCENARIOS[:] = self._orig[:max_scenarios]
            return self

        def __exit__(self, *exc):
            gtp.SCENARIOS[:] = self._orig
            return False

    return _Ctx()


def _write_input_manifests(per_dataset: list[dict[str, Any]]) -> None:
    """STEP 6 manifests describing the harmonized external inputs."""
    rows_input, rows_gene, rows_truth = [], [], []
    for d in per_dataset:
        rows_input.append({
            "dataset": d["dataset"],
            "reference_counts": str(d["reference_counts"]),
            "reference_cell_metadata": str(d["reference_cell_metadata"]),
            "bulk_counts": str(d["bulk_counts"]),
            "bulk_sample_metadata": str(d["bulk_sample_metadata"]),
            "selected_genes": str(d["selected_genes"]),
            "label_mapping": str(d["label_mapping"]),
            "n_reference_cells": d["n_reference_cells"],
            "n_bulk_samples": d["n_bulk_samples"],
            "n_selected_genes": d["n_selected_genes"],
        })
        rows_gene.append({
            "dataset": d["dataset"],
            "n_selected_genes": d["n_selected_genes"],
            "n_genes_in_bulk": d["n_genes_in_bulk"],
            "n_shared_genes": d["n_shared_genes"],
        })
        rows_truth.append({
            "dataset": d["dataset"],
            "fine_truth": str(d["fine_truth"]),
            "broad_truth": str(d["broad_truth"]),
            "conditional_truth": str(d["conditional_truth"]),
            "truth_type": "RNA_derived_fraction",
            "n_fine_labels": d["n_fine_labels"],
        })
    _write_tsv_no_index(pd.DataFrame(rows_input), INPUT_ROOT / "input_manifest.tsv")
    _write_tsv_no_index(pd.DataFrame(rows_gene), INPUT_ROOT / "gene_overlap_summary.tsv")
    _write_tsv_no_index(pd.DataFrame(rows_truth), INPUT_ROOT / "truth_manifest.tsv")


def _accumulate_scored(rows: dict[str, list], scored: dict[str, list],
                       ds_name: str, method: str) -> None:
    for key in SCORE_KEYS:
        for row in scored.get(key, []):
            row["dataset"] = ds_name
            row["mode"] = method  # explicit method name; matches internal "mode" column
        rows[key].extend(scored.get(key, []))


def run_benchmark(args: argparse.Namespace) -> int:
    global OUTPUT_ROOT
    if getattr(args, "out_dir", None):
        OUTPUT_ROOT = Path(args.out_dir)
    _ensure_dir(OUTPUT_ROOT)
    _ensure_dir(INPUT_ROOT)

    selected_tools = [TOOL_KEYS[t] for t in args.tools]
    rows = {k: [] for k in SCORE_KEYS}
    status_rows: list[dict[str, Any]] = []
    tool_registry_rows: list[dict[str, Any]] = []
    input_records: list[dict[str, Any]] = []

    with _limited_scenarios(args.max_scenarios):
        for ds_name in args.datasets:
            spec = DATASETS[ds_name]
            if not spec.path.exists():
                print(f"[skip] {ds_name}: missing {spec.path}")
                continue
            print(f"[dataset] {ds_name} (spec.name={spec.name})")
            ctx = prepare_dataset(spec, seed=args.seed,
                                  max_ref_cells_per_type=args.max_ref_cells_per_type)
            bundles = generate_all_mixtures(
                ctx, n_samples_per_scenario=args.samples_per_scenario,
                cells_per_sample=args.cells_per_mixture, seed=args.seed)
            full_bundle = concat_full_overlap(bundles)
            dataset_inputs = _prepare_dataset_inputs(ctx, full_bundle)

            sig = _reference_signature(ctx)
            n_shared = len([g for g in sig.index if g in full_bundle.counts.index])
            input_records.append({
                "dataset": spec.name,
                "reference_counts": dataset_inputs.reference_counts,
                "reference_cell_metadata": dataset_inputs.reference_cell_metadata,
                "bulk_counts": dataset_inputs.bulk_counts,
                "bulk_sample_metadata": dataset_inputs.bulk_sample_metadata,
                "selected_genes": dataset_inputs.selected_genes,
                "label_mapping": dataset_inputs.label_mapping,
                "fine_truth": dataset_inputs.fine_truth,
                "broad_truth": dataset_inputs.broad_truth,
                "conditional_truth": dataset_inputs.conditional_truth,
                "n_reference_cells": int(ctx.train_cells.size),
                "n_bulk_samples": int(full_bundle.counts.shape[1]),
                "n_selected_genes": len(ctx.selected_genes),
                "n_genes_in_bulk": int(full_bundle.counts.shape[0]),
                "n_shared_genes": n_shared,
                "n_fine_labels": len(ctx.selected_fine),
            })

            if args.dry_run:
                print(f"  [dry-run] prepared inputs + truth for {spec.name}; "
                      f"would run tools: {selected_tools}")
                continue

            # ---- NNLS external control (smoke test) ----
            if "NNLS_external_control" in selected_tools:
                print("  [tool] NNLS_external_control")
                try:
                    pred, meta = run_nnls_control(ctx, full_bundle)
                    _ensure_dir(OUTPUT_ROOT / "raw_predictions")
                    _write_tsv(pred, OUTPUT_ROOT / "raw_predictions"
                               / f"NNLS_external_control__{spec.name}.tsv")
                    scored = score_external_predictions(
                        ctx, full_bundle, "NNLS_external_control", pred, meta)
                    _accumulate_scored(rows, scored, spec.name, "NNLS_external_control")
                    status_rows.append({
                        "tool_name": "NNLS_external_control", "dataset": spec.name,
                        "executed_this_session": True, "run_status": "executed",
                        "status": "executed", "runtime_seconds": meta["runtime_seconds"],
                        "reason": f"n_shared_genes={meta['n_shared_genes']}",
                        "output_path": str(OUTPUT_ROOT / "raw_predictions"
                                           / f"NNLS_external_control__{spec.name}.tsv"),
                        "metadata_path": "", "error_log_path": "",
                        "estimate_type": meta["estimate_type"],
                    })
                except Exception as exc:  # noqa: BLE001
                    print(f"    [NNLS] FAILED: {exc}")
                    status_rows.append({
                        "tool_name": "NNLS_external_control", "dataset": spec.name,
                        "executed_this_session": False, "run_status": "failed",
                        "status": "failed", "runtime_seconds": None,
                        "reason": str(exc)[:500], "output_path": "",
                        "metadata_path": "", "error_log_path": "",
                        "estimate_type": "RNA_fraction",
                    })
                    return 2  # STEP 7: stop if the control smoke test fails.

            # ---- external R tools ----
            _prepare_runner_inputs(dataset_inputs)
            for tool_name in EXTERNAL_TOOL_SCRIPTS:
                if tool_name not in selected_tools:
                    continue
                if args.run_available and not _tool_is_available(tool_name):
                    print(f"  [tool] {tool_name}: skipped (--run-available; unavailable)")
                    status_rows.append({
                        "tool_name": tool_name, "dataset": spec.name,
                        "executed_this_session": False, "run_status": "skipped",
                        "status": "skipped", "runtime_seconds": None,
                        "reason": "dependency unavailable", "output_path": "",
                        "metadata_path": "", "error_log_path": "",
                        "estimate_type": EXTERNAL_TOOL_SCRIPTS[tool_name][
                            "outputs_cell_fraction_or_rna_fraction"],
                    })
                    continue
                print(f"  [tool] {tool_name}")
                result = _run_external_tool(tool_name, spec.name)
                status_rows.append({
                    "tool_name": tool_name, "dataset": spec.name,
                    "executed_this_session": result.executed_this_session,
                    "run_status": result.run_status, "status": result.status,
                    "runtime_seconds": result.runtime_seconds, "reason": result.reason,
                    "output_path": str(result.output_path) if result.output_path else "",
                    "metadata_path": str(result.metadata_path) if result.metadata_path else "",
                    "error_log_path": str(result.error_log_path) if result.error_log_path else "",
                    "estimate_type": result.estimate_type,
                })
                registry = EXTERNAL_TOOL_SCRIPTS[tool_name]
                tool_registry_rows.append({
                    "tool_name": tool_name, "modality": registry["modality"],
                    "language": registry["language"],
                    "package_or_command": registry["package_or_command"],
                    "install_status": _tool_is_available(tool_name),
                    "run_status": result.status, "version": result.version or "unknown",
                    "environment": "Renv" if registry["language"] == "R" else "python",
                    "supports_bulk": registry["modality"] == "bulk",
                    "supports_spatial": registry["modality"] == "spatial",
                    "supports_reference_sc": True,
                    "outputs_cell_fraction_or_rna_fraction":
                        registry["outputs_cell_fraction_or_rna_fraction"],
                    "requires_license_or_web": registry["requires_license_or_web"],
                    "notes": result.reason[:200],
                })
                if result.status == "executed" and result.output_path:
                    try:
                        pred = _read_prediction_file(result.output_path)
                        meta = {"runtime_seconds": result.runtime_seconds,
                                "status": result.status}
                        scored = score_external_predictions(
                            ctx, full_bundle, tool_name, pred, meta)
                        _accumulate_scored(rows, scored, spec.name, tool_name)
                    except Exception as exc:  # noqa: BLE001
                        print(f"    [score] failed for {tool_name} on {spec.name}: {exc}")

    _write_input_manifests(input_records)

    if args.dry_run:
        print(f"[dry-run] wrote harmonized inputs + manifests under {INPUT_ROOT}")
        return 0

    _write_summary_tables(rows)
    _write_status_manifest(status_rows)
    if tool_registry_rows:
        _write_tool_registry(tool_registry_rows)
    _write_tsv(pd.DataFrame(status_rows), OUTPUT_ROOT / "external_method_status.tsv", index=False)

    external_metric_paths = {
        "broad_metrics": OUTPUT_ROOT / "broad_metrics.tsv",
        "fine_metrics": OUTPUT_ROOT / "fine_metrics.tsv",
        "conditional_family_metrics": OUTPUT_ROOT / "conditional_family_metrics.tsv",
        "runtime_metrics": OUTPUT_ROOT / "runtime_metrics.tsv",
    }
    _write_combined_tables(external_metric_paths)
    # combined method status (internal modes are always 'executed' in the saved run)
    _write_tsv(pd.DataFrame(status_rows), OUTPUT_ROOT / "combined_method_status.tsv", index=False)

    print(f"[done] wrote external benchmark summaries under {OUTPUT_ROOT}")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    internal = _load_internal_run_params()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets", nargs="+", choices=sorted(DATASETS), default=sorted(DATASETS),
                    help="Datasets to benchmark")
    ap.add_argument("--tools", nargs="+", choices=sorted(TOOL_KEYS), default=["nnls"],
                    help="External tools to attempt (default: NNLS control only)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Prepare/export harmonized inputs + manifests but run no tools")
    ap.add_argument("--run-available", action="store_true",
                    help="Skip tools whose dependency is unavailable instead of attempting")
    ap.add_argument("--max-scenarios", type=int, default=None,
                    help="Limit number of scenarios (smoke tests)")
    ap.add_argument("--out-dir", default=str(OUTPUT_ROOT),
                    help="Output directory for external benchmark summaries")
    ap.add_argument("--seed", type=int, default=internal["seed"])
    ap.add_argument("--samples-per-scenario", type=int, default=internal["samples_per_scenario"])
    ap.add_argument("--cells-per-mixture", type=int, default=internal["cells_per_mixture"])
    ap.add_argument("--max-ref-cells-per-type", type=int, default=internal["max_ref_cells_per_type"])
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    return run_benchmark(args)


if __name__ == "__main__":
    raise SystemExit(main())
