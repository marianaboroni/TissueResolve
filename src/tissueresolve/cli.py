"""
TissueResolve command-line interface.

Commands
--------
``tissueresolve bulk``      — bulk RNA-seq deconvolution (CHIMERA algorithm).
``tissueresolve spatial``   — 10x Visium spatial deconvolution (SpatCAR algorithm).
``tissueresolve info``      — print version and dependency information.

Each sub-command group is implemented in its own pipeline module and
registered here.  Sub-commands not yet implemented raise a clear error
rather than silently failing.
"""
from __future__ import annotations

import sys
from typing import Any

import click
import json
from pathlib import Path

from tissueresolve import __version__
from tissueresolve.io import autodetect
from tissueresolve.presets import get_preset
from tissueresolve.protocol import detect as protocol_detect


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="tissueresolve")
def cli() -> None:
    """TissueResolve — unified cell-type deconvolution."""


@cli.command(name="run")
@click.option("--reference", required=True, help="Single-cell reference (.h5ad) or saved reference dir")
@click.option("--query", required=True, help="Bulk counts table or Visium .h5ad or folder")
@click.option("--out", required=True, help="Output directory for analysis bundle")
@click.option("--mode", type=click.Choice(["auto", "bulk", "spatial"]), default="auto")
@click.option("--resolution-mode",
              type=click.Choice(["auto", "hierarchical", "flat", "none", "suggest"]),
              default="auto", show_default=True,
              help="auto = use hierarchical broad→fine when broad/fine labels or a "
                   "mapping are available, else flat (recommended default). "
                   "hierarchical = force broad→fine (requires labels/mapping). "
                   "flat/none = fine-only. suggest = flat + family recommendations.")
@click.option("--broad-cell-type-col", default="auto", show_default=True,
              help="obs column with broad/compartment labels (hierarchical mode). "
                   "'auto' detects a known candidate column.")
@click.option("--fine-cell-type-col", default="auto", show_default=True,
              help="obs column with fine/subpopulation labels (hierarchical mode). "
                   "'auto' detects a known candidate column.")
@click.option("--cell-type-hierarchy", "hierarchy_path", default=None,
              type=click.Path(exists=True),
              help="Optional fine→broad mapping TSV (columns: fine_cell_type, broad_cell_type). "
                   "Use when the reference has only fine labels.")
@click.option("--allow-unresolved/--no-allow-unresolved", default=True,
              help="Keep non-separable families at the broad level as unresolved mass.")
@click.option("--state-aware", is_flag=True, default=False,
              help="EXPERIMENTAL: state-aware broad→cell-type→state hierarchical "
                   "deconvolution (requires --resolution-mode hierarchical). It is "
                   "NOT part of the default v0.1 workflow and has not been "
                   "validated across real datasets. Falls back to a two-level "
                   "broad→cell-type run when no state labels exist.")
@click.option("--solver",
              type=click.Choice(["auto", "nnls", "weighted_nnls", "marker_nnls",
                                 "ridge_nnls", "ensemble_nnls", "pipeline"]),
              default="auto", show_default=True,
              help="Bulk solver backbone. 'auto' picks the best by gene-masking CV; "
                   "'pipeline' uses the protocol-aware weighted pipeline.")
@click.option("--preset", type=click.Choice(["quick", "standard", "publication", "diagnostic"]), default="standard")
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--force", is_flag=True, default=False,
              help="Allow writing into an output directory that already holds a "
                   "run of a DIFFERENT modality (overwrites it).  By default such "
                   "a cross-modality overwrite is refused.")
def run_cli(reference: str, query: str, out: str, mode: str, resolution_mode: str,
            broad_cell_type_col: str, fine_cell_type_col: str,
            hierarchy_path: str | None, allow_unresolved: bool, state_aware: bool,
            solver: str, preset: str, dry_run: bool, force: bool) -> None:
    """User-friendly top-level run: auto-detect inputs, write analysis plan, optionally run pipelines."""
    rc = _run_top_level(
        reference, query, out, mode, preset, resolution_mode,
        broad_cell_type_col=broad_cell_type_col,
        fine_cell_type_col=fine_cell_type_col,
        hierarchy_path=hierarchy_path,
        allow_unresolved=allow_unresolved,
        state_aware=state_aware,
        solver=solver,
        dry_run=dry_run,
        force=force,
    )
    if rc != 0:
        raise click.ClickException(f"tissueresolve run failed with code {rc}")


def _prior_run_modality(outp: Path) -> str | None:
    """Return the modality (``"bulk"``/``"spatial"``) of a previous run in *outp*,
    read from ``analysis_plan.json`` or ``run_metadata.json`` — or ``None`` if the
    directory holds no recognisable prior run."""
    for name in ("analysis_plan.json", "run_metadata.json"):
        p = outp / name
        if not p.exists():
            continue
        try:
            data = json.loads(p.read_text())
        except Exception:
            continue
        mode = data.get("mode")
        if mode is None and isinstance(data.get("analysis_plan"), dict):
            mode = data["analysis_plan"].get("mode")
        if mode in ("bulk", "spatial"):
            return mode
    return None


def _guard_output_modality(outp: Path, resolved_mode: str, force: bool) -> None:
    """Refuse to silently overwrite a different-modality run in *outp*.

    One ``tissueresolve run`` processes one modality.  Writing a spatial run on
    top of a bulk run (or vice-versa) would clobber ``report.html``, ``deconv/``,
    ``qc/`` and the JSON/methods/warnings files.  Same-modality re-runs are
    allowed (they intentionally refresh the directory); cross-modality writes
    require an explicit ``--force``.
    """
    prior = _prior_run_modality(outp)
    if prior is None or prior == resolved_mode or force:
        return
    raise click.ClickException(
        f"Output directory {str(outp)!r} already contains a {prior!r} run, but "
        f"this is a {resolved_mode!r} run.  One `tissueresolve run` processes a "
        "single modality, and writing here would overwrite the existing "
        f"{prior!r} results (report.html, deconv/, qc/, *.json, methods.txt, "
        "warnings.json).\n"
        "Use separate output directories, e.g.:\n"
        "  results/bulk    (--mode bulk)\n"
        "  results/spatial (--mode spatial)\n"
        "then combine them with `tissueresolve combine-report`.\n"
        "Pass --force to intentionally overwrite this directory.")


def _run_top_level(
    reference: str,
    query: str,
    out: str,
    mode: str,
    preset: str,
    resolution_mode: str,
    *,
    broad_cell_type_col: str = "auto",
    fine_cell_type_col: str = "auto",
    hierarchy_path: str | None = None,
    allow_unresolved: bool = True,
    state_aware: bool = False,
    solver: str = "auto",
    dry_run: bool = False,
    force: bool = False,
) -> int:
    outp = Path(out)
    outp.mkdir(parents=True, exist_ok=True)

    detected_ref = autodetect.detect_input_type(reference)
    detected_query = autodetect.detect_input_type(query)
    preset_params = get_preset(preset)

    proto: dict[str, Any] = {"reference": {}, "query": {}}
    try:
        import anndata as ad

        if Path(reference).suffix == ".h5ad":
            adata = ad.read_h5ad(reference)
            proto["reference"] = protocol_detect.detect_reference_protocol(adata)
    except Exception:
        pass

    resolved_mode = mode
    if mode == "auto":
        if detected_query == "bulk_counts_table":
            resolved_mode = "bulk"
        elif detected_query in ("spatial_visium_h5ad", "spatial_visium_folder"):
            resolved_mode = "spatial"
        else:
            raise ValueError(
                "Could not auto-detect query mode; use --mode to specify 'bulk' or 'spatial'."
            )

    # Refuse to silently overwrite a different-modality run in this directory.
    _guard_output_modality(outp, resolved_mode, force)

    # Resolve the requested resolution mode to a concrete one (auto → hierarchical
    # when broad/fine labels/mapping are available; else flat-with-caution).
    requested_resolution_mode = resolution_mode
    resolution_mode, resolution_reason = _select_resolution_mode(
        requested_resolution_mode, reference=reference,
        broad_col=broad_cell_type_col,
        fine_col=fine_cell_type_col,
        hierarchy_path=hierarchy_path, preset=preset)
    print(f"Resolution mode: requested={requested_resolution_mode!r} → "
          f"using {resolution_mode!r}.")
    print(f"  reason: {resolution_reason}")

    plan = {
        "detected_reference": detected_ref,
        "detected_query": detected_query,
        "mode": resolved_mode,
        "requested_resolution_mode": requested_resolution_mode,
        "resolution_mode": resolution_mode,
        "resolution_mode_reason": resolution_reason,
        "preset": preset,
        "preset_params": preset_params,
        "protocol": proto,
        "solver": solver,
    }
    # State-aware is experimental and only meaningful in hierarchical mode.
    state_aware_effective = bool(state_aware and resolution_mode == "hierarchical")
    plan["hierarchy_mode"] = ("state_aware" if state_aware_effective
                              else ("hierarchical" if resolution_mode == "hierarchical"
                                    else "standard"))
    plan["state_aware_enabled"] = state_aware_effective
    plan["state_aware_feature_status"] = "experimental"
    if state_aware and not state_aware_effective:
        plan["state_aware_fallback_reason"] = (
            "--state-aware ignored: requires --resolution-mode hierarchical")
    if resolution_mode == "hierarchical":
        plan["hierarchical"] = {
            "broad_cell_type_col": broad_cell_type_col,
            "fine_cell_type_col": fine_cell_type_col,
            "cell_type_hierarchy": hierarchy_path,
            "allow_unresolved": allow_unresolved,
            "state_aware": state_aware_effective,
        }
    (outp / "analysis_plan.json").write_text(json.dumps(plan, indent=2))

    if dry_run:
        _print_summary_table(
            [
                ("mode", resolved_mode),
                ("resolution_mode", resolution_mode),
                ("resolution_reason", resolution_reason),
                ("preset", preset),
                ("reference type", detected_ref),
                ("query type", detected_query),
                ("analysis_plan", str(outp / "analysis_plan.json")),
            ]
        )
        print("Dry run: no algorithms executed.")
        return 0

    cfg = _configure_from_preset(preset_params)
    if resolution_mode == "hierarchical":
        cfg.hierarchical.broad_cell_type_col = broad_cell_type_col
        cfg.hierarchical.fine_cell_type_col = fine_cell_type_col
        cfg.hierarchical.allow_unresolved = allow_unresolved
    if resolved_mode == "bulk":
        result = _execute_bulk(
            reference, query, outp, cfg, resolution_mode=resolution_mode,
            hierarchy_path=hierarchy_path, solver=solver,
            state_aware=state_aware_effective)
    else:
        result = _execute_spatial(
            reference, query, outp, cfg, resolution_mode=resolution_mode,
            hierarchy_path=hierarchy_path)

    run_metadata = {
        "tissueresolve_version": __version__,
        "requested_resolution_mode": requested_resolution_mode,
        "resolution_mode": resolution_mode,
        "resolution_mode_reason": resolution_reason,
        "analysis_plan": plan,
        "pipeline_run": result.run_metadata,
    }
    # surface the selected mode + reason on the result so the report can state it
    try:
        result.run_metadata.setdefault("resolution_mode", resolution_mode)
        result.run_metadata["resolution_mode_reason"] = resolution_reason
    except Exception:
        pass
    (outp / "run_metadata.json").write_text(json.dumps(run_metadata, indent=2, default=str))

    # Report bundle: methods.txt, warnings.json, report.html (state-aware uses a
    # different result shape and writes its own outputs, so skip it here).
    if not state_aware_effective and hasattr(result, "deconv"):
        _write_run_report_bundle(result, outp, resolved_mode)

    _print_summary_table(
        [
            ("mode", resolved_mode),
            ("resolution_mode", resolution_mode),
            ("preset", preset),
            ("reference type", detected_ref),
            ("query type", detected_query),
            ("results_dir", str(outp)),
            ("n_cell_types", getattr(result.deconv, "n_cell_types", "—")),
            ("n_samples/spots", getattr(result.deconv, "n_samples", getattr(result.deconv, "n_spots", "—"))),
            ("analysis_plan", str(outp / "analysis_plan.json")),
        ]
    )
    print(f"Results written to {outp}/")
    return 0


def _configure_from_preset(preset_params: dict[str, Any]):
    from tissueresolve.config import TissueResolveConfig

    cfg = TissueResolveConfig()
    if preset_params.get("bootstrap"):
        cfg.bootstrap.n_bootstrap = preset_params.get("n_bootstrap", cfg.bootstrap.n_bootstrap)
    hp = preset_params.get("hierarchical") or {}
    for k, v in hp.items():
        if hasattr(cfg.hierarchical, k):
            setattr(cfg.hierarchical, k, v)
    return cfg


def _detect_hierarchy_availability(reference: str, broad_col: str,
                                   fine_col: str, hierarchy_path):
    """Cheaply check whether hierarchical annotations are available.

    Returns ``(available: bool, source: str | None)`` without building the full
    mapping.  A ``--cell-type-hierarchy`` file always counts; otherwise the
    reference ``.h5ad`` ``obs`` is inspected (backed) for broad + fine columns.
    Saved ReferenceSignature directories carry no per-cell annotations, so only
    the mapping-file route applies to them.
    """
    if hierarchy_path:
        return True, f"mapping file: {hierarchy_path}"
    path = Path(reference)
    if path.suffix in {".h5ad", ".h5"}:
        try:
            import anndata as ad
            from tissueresolve.io import validation as v

            obs = ad.read_h5ad(reference, backed="r").obs
            b = (v.detect_broad_cell_type_col(obs)
                 if broad_col in (None, "auto") else
                 (broad_col if broad_col in obs.columns else None))
            f = (v.detect_fine_cell_type_col(obs, exclude=b)
                 if fine_col in (None, "auto") else
                 (fine_col if fine_col in obs.columns else None))
            if b and f and b != f:
                return True, f"obs columns broad={b!r}, fine={f!r}"
        except Exception:
            return False, None
    return False, None


def _select_resolution_mode(requested: str, *, reference: str, broad_col: str,
                            fine_col: str, hierarchy_path, preset: str):
    """Resolve a user-requested resolution mode to a concrete one + a reason.

    ``auto`` (the default) selects **hierarchical** broad→fine when broad/fine
    annotations or a mapping are available, because that reduces spillover and
    yields more reliable interpretation.  When they are not available, ``auto``
    falls back to flat (fine-only) with a caution — or, for the
    ``publication``/``diagnostic`` presets, stops and asks for broad/fine
    labels (those presets imply a publication-grade claim).

    Returns ``(resolved_mode, reason)`` where *resolved_mode* is one of
    ``"hierarchical" | "none" | "suggest"`` (``"flat"`` maps to ``"none"``).
    """
    if requested == "flat":
        return "none", "user explicitly requested flat (fine-only) deconvolution"
    if requested == "none":
        return "none", "user requested none (flat, fine-only)"
    if requested == "suggest":
        return "suggest", "user requested suggest (flat + family recommendations)"
    if requested == "hierarchical":
        return "hierarchical", "user requested hierarchical broad→fine deconvolution"

    # requested == "auto"
    available, source = _detect_hierarchy_availability(
        reference, broad_col, fine_col, hierarchy_path)
    if available:
        return "hierarchical", (
            f"auto: hierarchical broad→fine selected because hierarchical "
            f"annotations are available ({source})")
    if preset in ("publication", "diagnostic"):
        raise click.ClickException(
            "auto resolution-mode with the "
            f"'{preset}' preset requires broad/fine cell-type annotations for "
            "publication-grade hierarchical deconvolution, but none were found. "
            "Add broad/fine columns to the reference (and set "
            "--broad-cell-type-col / --fine-cell-type-col), provide "
            "--cell-type-hierarchy mapping.tsv, or pass --resolution-mode flat "
            "to run fine-only deconvolution explicitly.")
    return "none", (
        "auto: no broad/fine annotations or mapping found; falling back to flat "
        "(fine-only).  Provide broad/fine labels or --cell-type-hierarchy to "
        "enable the recommended hierarchical broad→fine workflow.")


def _resolve_hierarchy_mapping(reference_path: str, ref, cfg, hierarchy_path):
    """Resolve a fine→broad mapping for hierarchical mode (or fail clearly).

    Resolution order:

    1. ``--cell-type-hierarchy`` mapping file (if provided);
    2. broad/fine annotation columns in the reference ``.h5ad`` ``obs``;
    3. otherwise raise an actionable error.

    Returns ``(mapping, source_str)``.
    """
    from tissueresolve.reference.hierarchy import (
        build_cell_type_hierarchy, load_hierarchy_mapping,
    )

    cell_types = list(ref.cell_types)

    if hierarchy_path:
        raw = load_hierarchy_mapping(hierarchy_path)
        mapping = build_cell_type_hierarchy(cell_types, raw)
        return mapping, f"mapping file: {hierarchy_path}"

    # try broad/fine columns from the source h5ad
    path = Path(reference_path)
    if path.suffix in {".h5ad", ".h5"}:
        try:
            import anndata as ad
            from tissueresolve.io import validation as v

            adata = ad.read_h5ad(reference_path, backed="r")
            obs = adata.obs
            hcfg = cfg.hierarchical
            broad_col = (v.detect_broad_cell_type_col(obs)
                         if hcfg.broad_cell_type_col in (None, "auto")
                         else hcfg.broad_cell_type_col)
            fine_col = (v.detect_fine_cell_type_col(obs, exclude=broad_col)
                        if hcfg.fine_cell_type_col in (None, "auto")
                        else hcfg.fine_cell_type_col)
            if broad_col and fine_col:
                info = v.validate_hierarchical_annotations(
                    obs.copy(), broad_col, fine_col)
                click.echo(
                    f"Detected hierarchical annotations: broad={broad_col!r}, "
                    f"fine={fine_col!r} ({info['n_broad']} families, "
                    f"{info['n_fine']} fine types).")
                mapping = build_cell_type_hierarchy(cell_types, info["mapping"])
                return mapping, f"obs columns: broad={broad_col}, fine={fine_col}"
        except (KeyError, ValueError):
            raise
        except Exception:
            pass

    raise click.ClickException(
        "Hierarchical deconvolution requires broad and fine cell-type "
        "annotations.  Add two columns to adata.obs (and set "
        "--broad-cell-type-col / --fine-cell-type-col) or provide "
        "--cell-type-hierarchy mapping.tsv (columns: fine_cell_type, "
        "broad_cell_type)."
    )


def _read_counts_table(path: Path):
    import pandas as pd

    if not path.exists():
        raise FileNotFoundError(f"Bulk counts file not found: {path}")
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    table = pd.read_csv(path, sep=sep, index_col=0)
    if table.empty:
        raise ValueError(f"Bulk counts table is empty: {path}")
    return table


def _load_reference_signature(path: Path, cfg, estimate_overdispersion: bool = False):
    from tissueresolve.api import build_reference
    from tissueresolve.results import ReferenceSignature

    if path.is_dir() and (path / "metadata.json").exists():
        return ReferenceSignature.load(path)
    if path.suffix in {".h5ad", ".h5"}:
        return build_reference(str(path), config=cfg, cell_type_col="cell_type", estimate_overdispersion=estimate_overdispersion)
    raise ValueError(
        f"Cannot interpret reference {path!r}: expected a saved ReferenceSignature directory or an .h5ad file."
    )


def _execute_bulk(reference: str, query: str, outp: Path, cfg, resolution_mode: str,
                  *, hierarchy_path: str | None = None, solver: str = "auto",
                  state_aware: bool = False):
    from tissueresolve.api import deconv_bulk

    ref_path = Path(reference)
    if ref_path.is_dir() and (ref_path / "metadata.json").exists():
        from tissueresolve.results import ReferenceSignature
        ref = ReferenceSignature.load(ref_path)
    else:
        ref = _load_reference_signature(ref_path, cfg, estimate_overdispersion=False)

    hierarchy_mapping = None
    if resolution_mode == "hierarchical":
        hierarchy_mapping, source = _resolve_hierarchy_mapping(
            reference, ref, cfg, hierarchy_path)
        click.echo(f"Hierarchy source: {source}")

    bulk = _read_counts_table(Path(query))
    # solver backbone applies to flat (non-hierarchical) runs; 'pipeline' or
    # hierarchical mode use the protocol-aware weighted pipeline.
    solver_arg = None if (solver in (None, "pipeline") or
                          resolution_mode == "hierarchical") else solver
    state_aware_eff = bool(state_aware and resolution_mode == "hierarchical")
    result = deconv_bulk(
        bulk,
        ref,
        config=cfg,
        resolution_mode=resolution_mode,
        hierarchy_mapping=hierarchy_mapping,
        solver=solver_arg,
        state_aware=state_aware_eff,
        n_bootstrap=cfg.bootstrap.n_bootstrap,
    )

    if state_aware_eff:
        # state-aware returns a StateAwareBulkResult (different shape): write its
        # own outputs and expose run metadata for the report.
        from tissueresolve.bulk.state_aware_hierarchical import (
            write_state_aware_outputs)
        write_state_aware_outputs(result, outp / "deconvolution")
        result.run_metadata = getattr(result, "metadata", {})  # for downstream callers
        return result

    result.deconv.save(outp / "deconv")
    result.qc.save(outp / "qc")
    if resolution_mode == "hierarchical":
        from tissueresolve.bulk.hierarchical import save_hierarchical_bulk_outputs
        save_hierarchical_bulk_outputs(result, outp / "hierarchical")
    return result


def _execute_spatial(reference: str, query: str, outp: Path, cfg, resolution_mode: str,
                     *, hierarchy_path: str | None = None):
    from tissueresolve.api import deconv_spatial
    from tissueresolve.io.spatial import load_visium
    from tissueresolve.results import ReferenceSignature

    ref_path = Path(reference)
    if ref_path.is_dir() and (ref_path / "metadata.json").exists():
        ref = ReferenceSignature.load(ref_path)
    else:
        ref = _load_reference_signature(ref_path, cfg, estimate_overdispersion=True)

    hierarchy_mapping = None
    if resolution_mode == "hierarchical":
        hierarchy_mapping, source = _resolve_hierarchy_mapping(
            reference, ref, cfg, hierarchy_path)
        click.echo(f"Hierarchy source: {source}")

    visium = load_visium(query, min_counts=0, min_genes=0)
    Y = visium.X
    array_row = visium.obs["array_row"].to_numpy()
    array_col = visium.obs["array_col"].to_numpy()
    lib_sizes = visium.obs["total_counts"].to_numpy().astype("float32")
    visium_gene_names = list(visium.var_names)
    spot_ids = list(visium.obs_names)

    result = deconv_spatial(
        Y,
        ref,
        array_row,
        array_col,
        lib_sizes,
        visium_gene_names,
        spot_ids=spot_ids,
        config=cfg,
        resolution_mode=resolution_mode,
        hierarchy_mapping=hierarchy_mapping,
        run_neighbourhood=False,
    )

    result.deconv.save(outp / "deconv")
    result.qc.save(outp / "qc")
    if resolution_mode == "hierarchical":
        from tissueresolve.spatial.hierarchical import save_hierarchical_spatial_outputs
        save_hierarchical_spatial_outputs(result, outp / "hierarchical")
    return result


def _collect_run_warnings(result, modality: str) -> list[dict[str, str]]:
    """Collect human-readable warnings from a finished run for ``warnings.json``.

    Warnings are surfaced, never hidden (CLAUDE.md rule 2): the estimate-type
    caveat is always recorded, plus QC recommendations, non-convergence
    (spatial), and protocol risk when present.
    """
    warns: list[dict[str, str]] = []
    if modality == "bulk":
        warns.append({"severity": "info", "category": "estimate_type",
                      "message": "Estimates are RNA-derived mRNA proportions, "
                                 "not absolute cell fractions."})
    else:
        warns.append({"severity": "info", "category": "estimate_type",
                      "message": "Estimates are spot-level RNA-derived "
                                 "composition, not single-cell counts."})
    qc = getattr(result, "qc", None)
    for rec in (getattr(qc, "recommendations", None) or []):
        warns.append({"severity": "warning", "category": "qc", "message": str(rec)})
    deconv = getattr(result, "deconv", None)
    if modality == "spatial" and deconv is not None and \
            getattr(deconv, "converged", True) is False:
        warns.append({"severity": "error", "category": "convergence",
                      "message": f"Spatial solver did not converge within "
                                 f"{getattr(deconv, 'n_iter', '?')} iterations."})
    risk = getattr(result, "protocol_risk", None)
    if risk is not None and getattr(risk, "risk_level", None) not in (None, "low"):
        warns.append({"severity": "warning", "category": "protocol_risk",
                      "message": f"Protocol risk level: "
                                 f"{getattr(risk, 'risk_level', 'unknown')}."})
    return warns


def _write_run_report_bundle(result, outp: Path, modality: str) -> None:
    """Write methods.txt, warnings.json and report.html for a finished run.

    Uses the in-memory result (no new figures are rendered): predictions, QC,
    methods and warnings are populated from the pipeline result.  Failures here
    never abort a successful deconvolution — they are reported, not hidden.
    """
    from tissueresolve.report import methods_text as _mt

    # methods.txt
    try:
        methods = (_mt.compose_bulk_methods(result) if modality == "bulk"
                   else _mt.compose_spatial_methods(result))
        (outp / "methods.txt").write_text(methods, encoding="utf-8")
    except Exception as exc:  # noqa: BLE001
        click.echo(f"  warning: could not write methods.txt: {exc}", err=True)

    # warnings.json
    try:
        warns = _collect_run_warnings(result, modality)
        (outp / "warnings.json").write_text(
            json.dumps(warns, indent=2), encoding="utf-8")
    except Exception as exc:  # noqa: BLE001
        click.echo(f"  warning: could not write warnings.json: {exc}", err=True)

    # report.html (from the in-memory result → populated predictions + QC)
    try:
        from tissueresolve.report import generate_report
        generate_report(modality, result, out=outp / "report.html")
    except Exception as exc:  # noqa: BLE001
        click.echo(f"  warning: could not generate report.html: {exc}", err=True)


def _print_summary_table(rows: list[tuple[str, Any]]) -> None:
    width = max(len(k) for k, _ in rows)
    for k, v in rows:
        print(f"{k:<{width}}  {v}")


def run(argv: list | None = None) -> int:
    """Programmatic entrypoint used by tests: parse args and write analysis plan.

    Returns 0 on success or raises SystemExit for --help.
    """
    import argparse

    parser = argparse.ArgumentParser(prog="tissueresolve")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--query", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--mode", choices=("auto", "bulk", "spatial"), default="auto")
    parser.add_argument("--resolution-mode",
                        choices=("auto", "hierarchical", "flat", "none", "suggest"),
                        default="auto")
    parser.add_argument("--broad-cell-type-col", default="auto")
    parser.add_argument("--fine-cell-type-col", default="auto")
    parser.add_argument("--cell-type-hierarchy", dest="hierarchy_path", default=None)
    parser.add_argument("--allow-unresolved", dest="allow_unresolved",
                        action="store_true", default=True)
    parser.add_argument("--no-allow-unresolved", dest="allow_unresolved",
                        action="store_false")
    parser.add_argument("--solver",
                        choices=("auto", "nnls", "weighted_nnls", "marker_nnls",
                                 "ridge_nnls", "ensemble_nnls", "pipeline"),
                        default="auto")
    parser.add_argument("--preset", choices=("quick", "standard", "publication", "diagnostic"), default="standard")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true", default=False)
    args = parser.parse_args(argv)

    return _run_top_level(
        args.reference,
        args.query,
        args.out,
        args.mode,
        args.preset,
        args.resolution_mode,
        broad_cell_type_col=args.broad_cell_type_col,
        fine_cell_type_col=args.fine_cell_type_col,
        hierarchy_path=args.hierarchy_path,
        allow_unresolved=args.allow_unresolved,
        solver=args.solver,
        dry_run=args.dry_run,
        force=args.force,
    )


def wizard(argv: list | None = None) -> int:
    """Interactive wizard to write a simple tissueresolve_config.yaml in cwd.

    Returns 0 on success.
    """
    ref = input("Reference path (saved reference dir or .h5ad): ")
    query = input("Query path (bulk counts table or Visium .h5ad/folder): ")
    mode = input("Mode (auto/bulk/spatial): ")
    preset = input("Preset (quick/standard/publication/diagnostic): ")
    contents = [
        "# TissueResolve config written by wizard",
        f"reference: {ref}",
        f"query: {query}",
        f"mode: {mode}",
        f"preset: {preset}",
    ]
    Path("tissueresolve_config.yaml").write_text("\n".join(contents))
    return 0


# ---------------------------------------------------------------------------
# bulk sub-group (populated in Stage 3)
# ---------------------------------------------------------------------------

@cli.group()
def bulk() -> None:
    """Bulk RNA-seq deconvolution (CHIMERA algorithm)."""


@bulk.command(name="run")
def bulk_run() -> None:
    """Run bulk RNA-seq deconvolution.  (Not yet implemented.)"""
    click.echo("tissueresolve bulk run — not yet implemented.", err=True)
    sys.exit(2)


@bulk.command(name="check-compatibility")
def bulk_check() -> None:
    """Check bulk/reference protocol compatibility.  (Not yet implemented.)"""
    click.echo("tissueresolve bulk check-compatibility — not yet implemented.", err=True)
    sys.exit(2)


@bulk.command(name="benchmark")
def bulk_benchmark() -> None:
    """Run bulk deconvolution benchmark.  (Not yet implemented.)"""
    click.echo("tissueresolve bulk benchmark — not yet implemented.", err=True)
    sys.exit(2)


@bulk.command(name="report")
@click.option("--results-dir", required=True, type=click.Path(exists=True),
              help="Directory with tables/ and figures/ from a bulk run.")
@click.option("--out", "out_path", default=None, type=click.Path(),
              help="Output HTML path (default: <results-dir>/report.html).")
def bulk_report(results_dir: str, out_path: str | None) -> None:
    """Generate a bulk HTML report from a results directory."""
    from tissueresolve.report import generate_report

    out = generate_report("bulk", results_dir, out_path)
    click.echo(f"Wrote bulk report -> {out}")


# ---------------------------------------------------------------------------
# spatial sub-group (Stage 4)
# ---------------------------------------------------------------------------

@cli.group()
def spatial() -> None:
    """10x Visium spatial deconvolution (SpatCAR algorithm)."""


@spatial.command(name="run")
@click.option(
    "--visium", "visium_path", required=True,
    type=click.Path(exists=True),
    help="SpaceRanger output directory or Visium .h5ad file.",
)
@click.option(
    "--reference", "reference_path", required=True,
    type=click.Path(exists=True),
    help="Saved ReferenceSignature directory, or a single-cell reference .h5ad.",
)
@click.option(
    "--output", "output_dir", required=True,
    type=click.Path(file_okay=False),
    help="Output directory for results.",
)
@click.option(
    "--config", "config_path", type=click.Path(exists=True), default=None,
    help="Optional TissueResolveConfig YAML.  Defaults are used otherwise.",
)
@click.option(
    "--cell-type-col", default="cell_type", show_default=True,
    help="obs column with cell-type labels (when building a reference from .h5ad).",
)
@click.option(
    "--marker-genes", "marker_genes_path", type=click.Path(exists=True), default=None,
    help="Optional text file of marker gene names (one per line).  "
         "When omitted, markers are selected automatically.",
)
@click.option("--lambda-spatial", type=float, default=None,
              help="Override CAR spatial regularisation strength λ.")
@click.option("--max-iter", type=int, default=None,
              help="Override maximum solver iterations.")
@click.option("--random-state", type=int, default=None,
              help="Override RNG seed.")
@click.option("--min-counts", type=int, default=100, show_default=True,
              help="Minimum UMI per spot (Visium loading filter).")
@click.option("--min-genes", type=int, default=200, show_default=True,
              help="Minimum detected genes per spot.")
@click.option("--neighbourhood/--no-neighbourhood", default=False,
              help="Compute neighbourhood co-occurrence statistics.")
@click.option("--genome", type=click.Choice(["hg38", "mm10"]), default=None,
              help="Reference genome assembly (affects gene-filter lists).")
def spatial_run(
    visium_path: str,
    reference_path: str,
    output_dir: str,
    config_path: str | None,
    cell_type_col: str,
    marker_genes_path: str | None,
    lambda_spatial: float | None,
    max_iter: int | None,
    random_state: int | None,
    min_counts: int,
    min_genes: int,
    neighbourhood: bool,
    genome: str | None,
) -> None:
    """Run 10x Visium spatial deconvolution (SpatCAR NB-CAR model).

    Outputs (written under ``--output``):

    \b
    deconv/   spot × cell-type composition estimates (spot_rna_composition)
    qc/       per-spot QC, Moran's I, model QC
    run_metadata.json   resolved parameters and provenance
    """
    from pathlib import Path

    import numpy as np

    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.results import ReferenceSignature

    # --- config ---
    cfg = (
        TissueResolveConfig.from_yaml(config_path)
        if config_path else TissueResolveConfig()
    )
    if lambda_spatial is not None:
        cfg.spatial_solver.lambda_spatial = lambda_spatial
    if max_iter is not None:
        cfg.spatial_solver.max_iter = max_iter
    if random_state is not None:
        cfg.spatial_solver.random_state = random_state
    if genome is not None:
        cfg.reference.genome = genome

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # --- load Visium ---
    from tissueresolve.io.spatial import load_visium
    import scipy.sparse as sp

    click.echo(f"Loading Visium from {visium_path} …")
    adata = load_visium(visium_path, min_counts=min_counts, min_genes=min_genes)

    Y = adata.X if sp.issparse(adata.X) else np.asarray(adata.X)
    array_row = adata.obs["array_row"].to_numpy()
    array_col = adata.obs["array_col"].to_numpy()
    lib_sizes = adata.obs["total_counts"].to_numpy().astype("float32")
    visium_gene_names = list(adata.var_names)
    spot_ids = list(adata.obs_names)

    # --- reference: saved signature directory or build from .h5ad ---
    ref = _load_or_build_reference(reference_path, cfg, cell_type_col)
    click.echo(f"Reference: {ref!r}")

    # --- optional marker gene list ---
    marker_genes = None
    if marker_genes_path is not None:
        marker_genes = [
            g.strip() for g in Path(marker_genes_path).read_text().splitlines()
            if g.strip()
        ]
        click.echo(f"Using {len(marker_genes)} marker genes from file.")

    # --- run pipeline ---
    from tissueresolve.spatial.pipeline import SpatialPipeline

    click.echo("Running spatial deconvolution …")
    result = SpatialPipeline(cfg).run(
        Y, ref, array_row, array_col, lib_sizes, visium_gene_names, spot_ids,
        marker_genes=marker_genes,
        run_neighbourhood=neighbourhood,
    )

    # --- save ---
    result.deconv.save(out / "deconv")
    result.qc.save(out / "qc")

    import json
    with (out / "run_metadata.json").open("w", encoding="utf-8") as fh:
        json.dump(result.run_metadata, fh, indent=2, default=str)

    m = result.deconv
    click.echo(
        f"Done.  {m.proportions.shape[0]} spots × {m.n_cell_types} cell types, "
        f"converged={m.converged} ({m.n_iter} iters), "
        f"λ={m.lambda_spatial} (smoothing recorded)."
    )
    click.echo(f"Estimate type: {m.ESTIMATE_TYPE} (not single-cell counts).")
    click.echo(f"Results written to {out}/")


@spatial.command(name="report")
@click.option("--results-dir", required=True, type=click.Path(exists=True),
              help="Directory with tables/ and figures/ from a spatial run.")
@click.option("--out", "out_path", default=None, type=click.Path(),
              help="Output HTML path (default: <results-dir>/report.html).")
def spatial_report(results_dir: str, out_path: str | None) -> None:
    """Generate a spatial HTML report from a results directory."""
    from tissueresolve.report import generate_report

    out = generate_report("spatial", results_dir, out_path)
    click.echo(f"Wrote spatial report -> {out}")


@spatial.command(name="benchmark")
@click.option("--output", "output_dir", required=True,
              type=click.Path(file_okay=False),
              help="Output directory for the benchmark results table.")
@click.option("--n-spots", type=int, default=500, show_default=True)
@click.option("--n-types", type=int, default=5, show_default=True)
@click.option("--n-genes", type=int, default=200, show_default=True)
@click.option("--scenarios", default="all", show_default=True,
              help="'all' or a comma-separated subset: basic,mismatch,ablation,shuffle.")
@click.option("--seed", type=int, default=0, show_default=True)
@click.option("--max-iter", type=int, default=100, show_default=True,
              help="SpatCAR solver iterations.")
def spatial_benchmark(
    output_dir: str,
    n_spots: int,
    n_types: int,
    n_genes: int,
    scenarios: str,
    seed: int,
    max_iter: int,
) -> None:
    """Run the synthetic Visium benchmark (SpatCAR vs DWLS/Spatial-NNLS/NNLS).

    Writes ``benchmark_results.csv`` (the underlying data) under ``--output``.
    """
    from tissueresolve.spatial.benchmark import (
        benchmark_results_to_frame,
        run_benchmark_suite,
    )

    click.echo(
        f"Benchmark: {n_spots} spots × {n_types} types × {n_genes} genes, "
        f"scenarios={scenarios}, seed={seed} …"
    )
    results = run_benchmark_suite(
        n_spots=n_spots, n_types=n_types, n_genes=n_genes,
        scenarios=scenarios, output_dir=output_dir,
        seed=seed, spatcar_max_iter=max_iter,
    )
    frame = benchmark_results_to_frame(results)
    click.echo(frame.to_string(index=False))
    click.echo(f"\nSaved benchmark_results.csv to {output_dir}/")


@spatial.command(name="info")
@click.option("--config", "config_path", type=click.Path(exists=True), default=None,
              help="Optional TissueResolveConfig YAML to inspect.")
def spatial_info(config_path: str | None) -> None:
    """Print the resolved spatial solver and QC configuration."""
    from tissueresolve.config import TissueResolveConfig

    cfg = (
        TissueResolveConfig.from_yaml(config_path)
        if config_path else TissueResolveConfig()
    )
    s = cfg.spatial_solver
    click.echo("Spatial solver configuration:")
    click.echo(f"  lambda_spatial         {s.lambda_spatial}")
    click.echo(f"  max_iter               {s.max_iter}")
    click.echo(f"  tol                    {s.tol}")
    click.echo(f"  n_marker_genes         {s.n_marker_genes}")
    click.echo(f"  min_lfc                {s.min_lfc}")
    click.echo(f"  update_mismatch_every  {s.update_mismatch_every}")
    click.echo(f"  n_jobs                 {s.n_jobs}")
    click.echo(f"  random_state           {s.random_state}")
    click.echo(f"  reference genome       {cfg.reference.genome}")


def _load_or_build_reference(reference_path: str, cfg, cell_type_col: str):
    """Load a saved ReferenceSignature directory, or build one from an .h5ad.

    A directory containing ``metadata.json`` is treated as a saved
    :class:`~tissueresolve.results.ReferenceSignature`.  An ``.h5ad`` file is
    aggregated into a new reference (with NB overdispersion, required by the
    spatial model).
    """
    from pathlib import Path

    from tissueresolve.results import ReferenceSignature

    path = Path(reference_path)
    if path.is_dir() and (path / "metadata.json").exists():
        click.echo(f"Loading saved reference signature from {path} …")
        return ReferenceSignature.load(path)

    if path.suffix in {".h5ad", ".h5"}:
        click.echo(f"Building reference from {path} (cell_type_col={cell_type_col!r}) …")
        from tissueresolve.reference.build import ReferenceBuilder

        ref_cfg = cfg.reference
        ref_cfg.celltype_col = cell_type_col
        builder = ReferenceBuilder(ref_cfg)
        # phi_g (NB overdispersion) is required for the spatial NB model.
        return builder.build_from_h5ad(str(path), estimate_overdispersion=True)

    raise click.ClickException(
        f"Cannot interpret reference {reference_path!r}: expected a saved "
        "ReferenceSignature directory or an .h5ad file."
    )


# ---------------------------------------------------------------------------
# top-level report
# ---------------------------------------------------------------------------

@cli.command(name="report")
@click.option("--modality", type=click.Choice(["bulk", "spatial"]), required=True,
              help="Which report to generate.")
@click.option("--results-dir", required=True, type=click.Path(exists=True),
              help="Directory with tables/ and figures/.")
@click.option("--out", "out_path", default=None, type=click.Path(),
              help="Output HTML path (default: <results-dir>/report.html).")
def report(modality: str, results_dir: str, out_path: str | None) -> None:
    """Generate a bulk or spatial HTML report from a results directory."""
    from tissueresolve.report import generate_report

    out = generate_report(modality, results_dir, out_path)
    click.echo(f"Wrote {modality} report -> {out}")


@cli.command(name="combine-report")
@click.option("--bulk-dir", "bulk_dir", required=True, type=click.Path(exists=True),
              help="An existing `tissueresolve run --mode bulk` output directory.")
@click.option("--spatial-dir", "spatial_dir", required=True, type=click.Path(exists=True),
              help="An existing `tissueresolve run --mode spatial` output directory.")
@click.option("--out", "out_dir", required=True, type=click.Path(file_okay=False),
              help="Output directory for the combined report bundle.")
def combine_report(bulk_dir: str, spatial_dir: str, out_dir: str) -> None:
    """Combine an existing bulk run and an existing spatial run into one report.

    Reads the two run directories (it does NOT re-run deconvolution) and writes
    ``report.html`` + ``methods.txt`` + ``warnings.json`` + ``run_metadata.json``
    under ``--out``.  Bulk and spatial sections (and benchmark summaries) are kept
    separate; the report states it summarises two separate runs sharing a
    reference, not a single joint model.
    """
    from tissueresolve.report.combined import generate_combined_report

    try:
        out = generate_combined_report(bulk_dir, spatial_dir, out_dir)
    except ValueError as exc:
        raise click.ClickException(str(exc))
    click.echo(f"Wrote combined report -> {out}")


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

@cli.command()
def info() -> None:
    """Print version and installed dependency information."""
    import importlib.metadata as _meta

    click.echo(f"tissueresolve  {__version__}")

    deps = [
        "numpy", "scipy", "pandas", "anndata", "click",
        "pyyaml", "joblib", "tqdm",
    ]
    optional = ["scanpy", "scikit-learn", "jinja2", "matplotlib", "seaborn"]

    click.echo("\nCore dependencies:")
    for dep in deps:
        try:
            v = _meta.version(dep)
            click.echo(f"  {dep:<20} {v}")
        except _meta.PackageNotFoundError:
            click.echo(f"  {dep:<20} NOT INSTALLED")

    click.echo("\nOptional dependencies:")
    for dep in optional:
        try:
            v = _meta.version(dep)
            click.echo(f"  {dep:<20} {v}")
        except _meta.PackageNotFoundError:
            click.echo(f"  {dep:<20} not installed")
