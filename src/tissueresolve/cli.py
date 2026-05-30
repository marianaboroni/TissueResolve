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

import click

from tissueresolve import __version__


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="tissueresolve")
def cli() -> None:
    """TissueResolve — unified cell-type deconvolution."""


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
def spatial_report() -> None:
    """Generate HTML report.  (Not yet implemented — Stage 5.)"""
    click.echo(
        "tissueresolve spatial report — not yet implemented (report module "
        "arrives in Stage 5).",
        err=True,
    )
    sys.exit(2)


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
