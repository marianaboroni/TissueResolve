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
# spatial sub-group (populated in Stage 4)
# ---------------------------------------------------------------------------

@cli.group()
def spatial() -> None:
    """10x Visium spatial deconvolution (SpatCAR algorithm)."""


@spatial.command(name="run")
def spatial_run() -> None:
    """Run spatial deconvolution.  (Not yet implemented.)"""
    click.echo("tissueresolve spatial run — not yet implemented.", err=True)
    sys.exit(2)


@spatial.command(name="report")
def spatial_report() -> None:
    """Generate HTML report.  (Not yet implemented.)"""
    click.echo("tissueresolve spatial report — not yet implemented.", err=True)
    sys.exit(2)


@spatial.command(name="benchmark")
def spatial_benchmark() -> None:
    """Run spatial deconvolution benchmark.  (Not yet implemented.)"""
    click.echo("tissueresolve spatial benchmark — not yet implemented.", err=True)
    sys.exit(2)


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
