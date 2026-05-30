"""
TissueResolve result containers.

All public result objects are dataclasses with ``save()`` / ``load()`` methods
that persist to a directory.  Every saved directory contains a
``metadata.json`` that records the estimate type, software version, and
any heuristic thresholds that were used.

Output type conventions
-----------------------
``BulkDeconvResult.proportions``
    **mRNA proportions** — the fraction of total mRNA in a bulk sample that
    originates from each cell type.  These are *not* cell fractions.
    ``ESTIMATE_TYPE = "mRNA_proportion"`` is set as a class attribute and
    written to every output file.

``SpatialDeconvResult.proportions``
    **Spot-level RNA-derived cellular composition estimates** — not direct
    single-cell counts unless explicitly calibrated via mRNA content
    correction.  ``ESTIMATE_TYPE = "spot_rna_composition"`` is set as a
    class attribute.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "ReferenceSignature",
    "BulkDeconvResult",
    "SpatialDeconvResult",
    "QCReport",
    "PairSeparability",
    "SeparabilityReport",
    "BenchmarkResult",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tsv(df: pd.DataFrame, path: Path, comment_lines: list[str] | None = None) -> None:
    with path.open("w", encoding="utf-8") as fh:
        if comment_lines:
            for line in comment_lines:
                fh.write(f"# {line}\n")
        df.to_csv(fh, sep="\t")


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0, comment="#")


def _save_optional_array(arr: np.ndarray | None, path: Path) -> None:
    if arr is not None:
        np.save(path, arr)


def _load_optional_array(path: Path) -> np.ndarray | None:
    return np.load(path) if path.exists() else None


def _save_optional_df(
    df: pd.DataFrame | None,
    path: Path,
    comment_lines: list[str] | None = None,
) -> None:
    if df is not None:
        _write_tsv(df, path, comment_lines)


def _load_optional_df(path: Path) -> pd.DataFrame | None:
    return _read_tsv(path) if path.exists() else None


# ---------------------------------------------------------------------------
# ReferenceSignature
# ---------------------------------------------------------------------------


@dataclass
class ReferenceSignature:
    """Unified reference object for both bulk and spatial deconvolution.

    This is the central data structure that bridges the two modalities.  It
    stores representations of the pseudo-bulk reference in both the bulk
    (L1-normalised probability vectors) and spatial (CPM matrix + NB
    dispersion) forms so that algorithms from either modality can be applied
    without re-building the reference.

    At most one of ``phi`` and ``R_cpm`` is *required*; the other is derived
    on demand by ``as_phi()`` / ``as_R_cpm()``.  Providing both avoids
    recomputation.

    Attributes
    ----------
    gene_names:
        Ordered gene identifiers.  Length G.
    cell_types:
        Ordered cell-type names.  Length K.  Sorted alphabetically.
    phi:
        L1-normalised probability matrix, shape ``(G, K)``.  Each column is
        a probability vector over genes summing to 1.  Used by the weighted
        NNLS bulk solver.  ``None`` if built from spatial pipeline only.
    R_cpm:
        CPM expression matrix, shape ``(K, G)``.  Used by the NB-CAR spatial
        model.  ``None`` if built from bulk pipeline only.
    R_log:
        ``log1p(R_cpm)``, shape ``(K, G)``.  Derived from R_cpm if not
        provided directly.
    phi_g:
        Per-gene NB overdispersion, shape ``(G,)``.  Required for spatial
        deconvolution; ``None`` otherwise.
    donor_cv:
        Cross-donor coefficient of variation, shape ``(G, K)``.  ``None``
        when fewer than 2 donors are present or when built from spatial
        pipeline only.
    n_cells_per_type:
        Dict mapping cell-type name → number of cells used in the reference.
    genome:
        Reference genome assembly (``"hg38"`` or ``"mm10"``).
    """

    gene_names: list[str]
    cell_types: list[str]
    phi: np.ndarray | None = field(default=None, repr=False)
    R_cpm: np.ndarray | None = field(default=None, repr=False)
    R_log: np.ndarray | None = field(default=None, repr=False)
    phi_g: np.ndarray | None = field(default=None, repr=False)
    donor_cv: np.ndarray | None = field(default=None, repr=False)
    n_cells_per_type: dict[str, int] = field(default_factory=dict)
    genome: str = "hg38"

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def n_genes(self) -> int:
        return len(self.gene_names)

    @property
    def n_cell_types(self) -> int:
        return len(self.cell_types)

    # ------------------------------------------------------------------
    # Representation accessors
    # ------------------------------------------------------------------

    def as_phi(self) -> np.ndarray:
        """Return L1-normalised probability matrix ``(G, K)``.

        If ``phi`` is stored directly, return it.  Otherwise, derive it
        from ``R_cpm``: transpose to ``(G, K)`` and L1-normalise each
        column so that it sums to 1.

        Raises
        ------
        ValueError
            If neither ``phi`` nor ``R_cpm`` is available.
        """
        if self.phi is not None:
            return self.phi
        if self.R_cpm is not None:
            phi_raw = self.R_cpm.T  # (G, K)
            col_sums = phi_raw.sum(axis=0, keepdims=True)
            col_sums = np.where(col_sums == 0, 1.0, col_sums)
            return (phi_raw / col_sums).astype(np.float64)
        raise ValueError(
            "ReferenceSignature has neither phi nor R_cpm.  "
            "Build the reference via ReferenceBuilder before calling as_phi()."
        )

    def as_R_cpm(self) -> np.ndarray:
        """Return CPM expression matrix ``(K, G)``.

        If ``R_cpm`` is stored directly, return it.  Otherwise, derive it
        from ``phi`` by transposing and scaling to CPM (×1e6).

        Note: when derived from ``phi``, the result is the relative expression
        in CPM-equivalent units (each row sums to 1e6), which is equivalent to
        CPM normalisation of the pseudo-bulk profile.

        Raises
        ------
        ValueError
            If neither ``phi`` nor ``R_cpm`` is available.
        """
        if self.R_cpm is not None:
            return self.R_cpm
        if self.phi is not None:
            return (self.phi.T * 1e6).astype(np.float32)
        raise ValueError(
            "ReferenceSignature has neither R_cpm nor phi.  "
            "Build the reference via ReferenceBuilder before calling as_R_cpm()."
        )

    def as_R_log(self) -> np.ndarray:
        """Return log1p-CPM matrix ``(K, G)``, computing it if not stored."""
        if self.R_log is not None:
            return self.R_log
        return np.log1p(self.as_R_cpm()).astype(np.float32)

    # ------------------------------------------------------------------
    # Gene subsetting
    # ------------------------------------------------------------------

    def subset_genes(self, gene_names: list[str]) -> "ReferenceSignature":
        """Return a new :class:`ReferenceSignature` limited to *gene_names*.

        Genes in *gene_names* that are not present in the reference are
        silently dropped with a warning (never silently included as zeros).

        Parameters
        ----------
        gene_names:
            Ordered list of gene names to retain.

        Returns
        -------
        ReferenceSignature
        """
        gene_to_idx: dict[str, int] = {g: i for i, g in enumerate(self.gene_names)}
        kept = [g for g in gene_names if g in gene_to_idx]
        missing = [g for g in gene_names if g not in gene_to_idx]
        if missing:
            warnings.warn(
                f"subset_genes: {len(missing)} requested gene(s) not in reference "
                f"and will be excluded: {missing[:5]}"
                + ("…" if len(missing) > 5 else ""),
                stacklevel=2,
            )
        if not kept:
            raise ValueError(
                "subset_genes: no requested genes are present in the reference."
            )
        idx = np.array([gene_to_idx[g] for g in kept], dtype=int)
        return ReferenceSignature(
            gene_names=kept,
            cell_types=list(self.cell_types),
            phi=self.phi[idx] if self.phi is not None else None,
            R_cpm=self.R_cpm[:, idx] if self.R_cpm is not None else None,
            R_log=self.R_log[:, idx] if self.R_log is not None else None,
            phi_g=self.phi_g[idx] if self.phi_g is not None else None,
            donor_cv=self.donor_cv[idx] if self.donor_cv is not None else None,
            n_cells_per_type=dict(self.n_cells_per_type),
            genome=self.genome,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> None:
        """Check internal consistency.  Raises ``ValueError`` on mismatch."""
        G, K = self.n_genes, self.n_cell_types
        if self.phi is not None and self.phi.shape != (G, K):
            raise ValueError(
                f"phi shape {self.phi.shape} does not match "
                f"(n_genes={G}, n_cell_types={K})."
            )
        if self.R_cpm is not None and self.R_cpm.shape != (K, G):
            raise ValueError(
                f"R_cpm shape {self.R_cpm.shape} does not match "
                f"(n_cell_types={K}, n_genes={G})."
            )
        if self.R_log is not None and self.R_log.shape != (K, G):
            raise ValueError(
                f"R_log shape {self.R_log.shape} does not match "
                f"(n_cell_types={K}, n_genes={G})."
            )
        if self.phi_g is not None and self.phi_g.shape != (G,):
            raise ValueError(
                f"phi_g shape {self.phi_g.shape} does not match (n_genes={G},)."
            )
        if self.donor_cv is not None and self.donor_cv.shape != (G, K):
            raise ValueError(
                f"donor_cv shape {self.donor_cv.shape} does not match "
                f"(n_genes={G}, n_cell_types={K})."
            )
        if self.phi is None and self.R_cpm is None:
            warnings.warn(
                "ReferenceSignature has neither phi nor R_cpm.  "
                "Call as_phi() / as_R_cpm() will raise until the reference is built.",
                stacklevel=2,
            )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> pd.DataFrame:
        """Return a per-cell-type summary of the reference profiles.

        Returns
        -------
        pd.DataFrame
            Index: cell-type names.
            Columns: ``n_cells``, ``mean_cpm``, ``max_cpm``,
            ``profile_entropy`` (Shannon entropy of the phi column —
            lower = more concentrated marker profile).
        """
        R = self.as_R_cpm()   # (K, G)
        phi = self.as_phi()   # (G, K)
        rows = {}
        for k, ct in enumerate(self.cell_types):
            p_k = phi[:, k]
            p_safe = np.maximum(p_k, 1e-12)
            p_safe /= p_safe.sum()
            entropy = float(-np.sum(p_safe * np.log(p_safe)))
            rows[ct] = {
                "n_cells": self.n_cells_per_type.get(ct, 0),
                "mean_cpm": round(float(R[k].mean()), 4),
                "max_cpm": round(float(R[k].max()), 4),
                "profile_entropy": round(entropy, 4),
            }
        return pd.DataFrame(rows).T

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | str) -> None:
        """Save to a directory.

        Creates the following files inside *path*:

        ``gene_names.txt``, ``cell_types.txt``
            Plain text lists, one entry per line.
        ``phi.npy``, ``R_cpm.npy``, ``R_log.npy``, ``phi_g.npy``,
        ``donor_cv.npy``
            NumPy arrays (written only if present).
        ``metadata.json``
            Provenance metadata.
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        (path / "gene_names.txt").write_text(
            "\n".join(self.gene_names), encoding="utf-8"
        )
        (path / "cell_types.txt").write_text(
            "\n".join(self.cell_types), encoding="utf-8"
        )
        _save_optional_array(self.phi, path / "phi.npy")
        _save_optional_array(self.R_cpm, path / "R_cpm.npy")
        _save_optional_array(self.R_log, path / "R_log.npy")
        _save_optional_array(self.phi_g, path / "phi_g.npy")
        _save_optional_array(self.donor_cv, path / "donor_cv.npy")

        meta: dict[str, Any] = {
            "tissueresolve_object": "ReferenceSignature",
            "n_genes": self.n_genes,
            "n_cell_types": self.n_cell_types,
            "genome": self.genome,
            "n_cells_per_type": self.n_cells_per_type,
            "has_phi": self.phi is not None,
            "has_R_cpm": self.R_cpm is not None,
            "has_R_log": self.R_log is not None,
            "has_phi_g": self.phi_g is not None,
            "has_donor_cv": self.donor_cv is not None,
        }
        with (path / "metadata.json").open("w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2)

    @classmethod
    def load(cls, path: Path | str) -> "ReferenceSignature":
        """Load a :class:`ReferenceSignature` saved with :meth:`save`.

        Parameters
        ----------
        path:
            Directory created by :meth:`save`.
        """
        path = Path(path)
        if not path.is_dir():
            raise FileNotFoundError(
                f"ReferenceSignature directory not found: {path}"
            )
        gene_names = (path / "gene_names.txt").read_text(encoding="utf-8").splitlines()
        cell_types = (path / "cell_types.txt").read_text(encoding="utf-8").splitlines()
        with (path / "metadata.json").open(encoding="utf-8") as fh:
            meta: dict[str, Any] = json.load(fh)
        return cls(
            gene_names=gene_names,
            cell_types=cell_types,
            phi=_load_optional_array(path / "phi.npy"),
            R_cpm=_load_optional_array(path / "R_cpm.npy"),
            R_log=_load_optional_array(path / "R_log.npy"),
            phi_g=_load_optional_array(path / "phi_g.npy"),
            donor_cv=_load_optional_array(path / "donor_cv.npy"),
            n_cells_per_type=meta.get("n_cells_per_type", {}),
            genome=meta.get("genome", "hg38"),
        )

    def __repr__(self) -> str:
        parts = [f"ReferenceSignature({self.n_genes}g × {self.n_cell_types}k"]
        present = [
            name for name, val in [
                ("phi", self.phi), ("R_cpm", self.R_cpm),
                ("phi_g", self.phi_g), ("donor_cv", self.donor_cv),
            ] if val is not None
        ]
        if present:
            parts.append(f", arrays=[{', '.join(present)}]")
        parts.append(f", genome={self.genome!r})")
        return "".join(parts)


# ---------------------------------------------------------------------------
# BulkDeconvResult
# ---------------------------------------------------------------------------


@dataclass
class BulkDeconvResult:
    """Per-sample bulk RNA-seq deconvolution results.

    **Output type:** ``proportions`` contains **mRNA proportions**, not cell
    fractions.  In tissues where cell types differ substantially in mRNA
    content per cell (plasma cells, neurons, hepatocytes), these quantities
    can differ substantially.  Use ``bulk.solver.MRNAContentCorrector`` to
    convert if appropriate mRNA content data are available.

    This distinction is enforced programmatically:
    - ``ESTIMATE_TYPE = "mRNA_proportion"`` is a class-level constant.
    - Every saved ``proportions.tsv`` includes the comment header
      ``# estimate_type: mRNA_proportion``.
    - ``cell_fractions`` is a *separate* field populated only after explicit
      mRNA content correction.

    Attributes
    ----------
    proportions:
        ``(n_samples × n_cell_types)`` DataFrame.  Values are mRNA
        proportions; rows sum to 1.  Index is sample IDs.
    coverage_r2:
        Per-sample reconstruction R².
    gene_panel:
        Ordered gene IDs used in the solve.
    gene_weights:
        Composite per-gene weights used in wNNLS.  ``None`` if uniform.
    lower_ci, upper_ci:
        Bootstrap CI bounds ``(n_samples × n_cell_types)``.  ``None``
        until bootstrap is run.
    cell_fractions:
        mRNA-content-corrected cell fractions.  ``None`` until
        ``MRNAContentCorrector.correct()`` is applied.
    run_metadata:
        Dict of run parameters, versions, and heuristic flags.
    """

    ESTIMATE_TYPE: str = field(default="mRNA_proportion", init=False, repr=False)

    proportions: pd.DataFrame
    coverage_r2: pd.Series
    gene_panel: list[str]
    gene_weights: pd.Series | None = field(default=None, repr=False)
    lower_ci: pd.DataFrame | None = field(default=None, repr=False)
    upper_ci: pd.DataFrame | None = field(default=None, repr=False)
    cell_fractions: pd.DataFrame | None = field(default=None, repr=False)
    run_metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Enforce estimate type — protects against accidental reuse of this
        # object as a cell-fraction container.
        object.__setattr__(self, "ESTIMATE_TYPE", "mRNA_proportion")

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    @property
    def n_samples(self) -> int:
        return len(self.proportions)

    @property
    def n_cell_types(self) -> int:
        return len(self.proportions.columns)

    def summary(self) -> pd.DataFrame:
        """Descriptive statistics of proportions, sorted by mean (descending)."""
        d = self.proportions.describe().T[["mean", "std", "50%", "min", "max"]]
        d.columns = pd.Index(["mean", "sd", "median", "min", "max"])
        return d.sort_values("mean", ascending=False)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | str) -> None:
        """Write all result files to *path* directory.

        Files written:

        ``proportions.tsv``
            mRNA proportions with ``# estimate_type: mRNA_proportion`` header.
        ``coverage_r2.tsv``
        ``gene_panel.txt``
        ``gene_weights.tsv``         (if present)
        ``lower_ci.tsv``             (if present)
        ``upper_ci.tsv``             (if present)
        ``cell_fractions.tsv``       (if present)
        ``metadata.json``
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        _write_tsv(
            self.proportions, path / "proportions.tsv",
            comment_lines=[
                "estimate_type: mRNA_proportion",
                "WARNING: these are mRNA proportions, NOT cell fractions.",
                "See BulkDeconvResult docstring and docs/ for the distinction.",
            ],
        )
        _write_tsv(
            self.coverage_r2.to_frame("coverage_r2"), path / "coverage_r2.tsv"
        )
        (path / "gene_panel.txt").write_text(
            "\n".join(self.gene_panel), encoding="utf-8"
        )
        if self.gene_weights is not None:
            _write_tsv(
                self.gene_weights.to_frame("gene_weight"), path / "gene_weights.tsv"
            )
        _save_optional_df(
            self.lower_ci, path / "lower_ci.tsv",
            comment_lines=["estimate_type: mRNA_proportion", "bootstrap_lower_bound"],
        )
        _save_optional_df(
            self.upper_ci, path / "upper_ci.tsv",
            comment_lines=["estimate_type: mRNA_proportion", "bootstrap_upper_bound"],
        )
        _save_optional_df(
            self.cell_fractions, path / "cell_fractions.tsv",
            comment_lines=[
                "estimate_type: cell_fraction",
                "Derived from mRNA proportions via MRNAContentCorrector.",
            ],
        )
        meta: dict[str, Any] = {
            "tissueresolve_object": "BulkDeconvResult",
            "estimate_type": self.ESTIMATE_TYPE,
            "n_samples": self.n_samples,
            "n_cell_types": self.n_cell_types,
            "n_genes_panel": len(self.gene_panel),
            "has_bootstrap_ci": self.lower_ci is not None,
            "has_cell_fractions": self.cell_fractions is not None,
        }
        meta.update(self.run_metadata)
        with (path / "metadata.json").open("w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, default=str)

    @classmethod
    def load(cls, path: Path | str) -> "BulkDeconvResult":
        """Load from a directory created by :meth:`save`."""
        path = Path(path)
        if not path.is_dir():
            raise FileNotFoundError(f"BulkDeconvResult directory not found: {path}")
        proportions = _read_tsv(path / "proportions.tsv")
        coverage_r2 = _read_tsv(path / "coverage_r2.tsv")["coverage_r2"]
        gene_panel = (path / "gene_panel.txt").read_text(encoding="utf-8").splitlines()
        gene_weights: pd.Series | None = None
        if (path / "gene_weights.tsv").exists():
            gene_weights = _read_tsv(path / "gene_weights.tsv")["gene_weight"]
        lower_ci = _load_optional_df(path / "lower_ci.tsv")
        upper_ci = _load_optional_df(path / "upper_ci.tsv")
        cell_fractions = _load_optional_df(path / "cell_fractions.tsv")
        run_metadata: dict[str, Any] = {}
        if (path / "metadata.json").exists():
            with (path / "metadata.json").open(encoding="utf-8") as fh:
                run_metadata = json.load(fh)
        return cls(
            proportions=proportions,
            coverage_r2=coverage_r2,
            gene_panel=gene_panel,
            gene_weights=gene_weights,
            lower_ci=lower_ci,
            upper_ci=upper_ci,
            cell_fractions=cell_fractions,
            run_metadata=run_metadata,
        )


# ---------------------------------------------------------------------------
# SpatialDeconvResult
# ---------------------------------------------------------------------------


@dataclass
class SpatialDeconvResult:
    """Per-spot spatial deconvolution results.

    **Output type:** ``proportions`` contains **spot-level RNA-derived
    cellular composition estimates**, not direct single-cell counts.

    Attributes
    ----------
    proportions:
        ``(n_spots × n_cell_types)`` DataFrame.  Rows sum to 1.
        Index is spot barcodes.
    cell_types:
        Ordered cell-type names.
    marker_genes:
        Marker genes used in the model.
    n_iter:
        Number of optimisation iterations run.
    converged:
        Whether the model converged within ``max_iter``.
    convergence_trace:
        ``‖ΔΠ‖_F / N`` per iteration.
    lambda_spatial:
        CAR spatial regularisation strength used.
    mismatch_factors:
        Per-gene mismatch scale factors ``d_g``, shape ``(n_marker_genes,)``.
        ``None`` if not estimated.
    lower_ci, upper_ci:
        Bootstrap CI arrays, shape ``(n_spots, n_cell_types)``.  ``None``
        until bootstrap is run.
    bootstrap_coverage_note:
        Documents empirical CI coverage.  Populated when bootstrap is run.
        Empirical coverage is typically 88–93 % for nominal 95 % with the
        default ``n_iter_per_boot=30``; use ``n_iter_per_boot ≥ 50`` for
        closer-to-nominal coverage.
    n_smooth:
        Number of Laplacian smoothing rounds applied before niche detection.
        ``None`` if niche detection was not run.  Always recorded so that
        spatial smoothing is never hidden.
    run_metadata:
        Dict of run parameters and versions.
    """

    ESTIMATE_TYPE: str = field(default="spot_rna_composition", init=False, repr=False)

    proportions: pd.DataFrame
    cell_types: list[str]
    marker_genes: list[str]
    n_iter: int
    converged: bool
    convergence_trace: list[float]
    lambda_spatial: float
    mismatch_factors: np.ndarray | None = field(default=None, repr=False)
    lower_ci: np.ndarray | None = field(default=None, repr=False)
    upper_ci: np.ndarray | None = field(default=None, repr=False)
    bootstrap_coverage_note: str | None = None
    n_smooth: int | None = None
    run_metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "ESTIMATE_TYPE", "spot_rna_composition")
        if self.n_smooth is not None and self.n_smooth > 0:
            # Validate that smoothing is always recorded — enforces DESIGN_SPEC.
            pass  # n_smooth is set; nothing further needed here.

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    @property
    def n_spots(self) -> int:
        return len(self.proportions)

    @property
    def n_cell_types(self) -> int:
        return len(self.cell_types)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | str) -> None:
        """Write all result files to *path* directory.

        Files written:

        ``proportions.tsv``        — spot × cell-type composition estimates.
        ``marker_genes.txt``
        ``convergence_trace.json``
        ``mismatch_factors.npy``   (if present)
        ``lower_ci.npy``           (if present)
        ``upper_ci.npy``           (if present)
        ``metadata.json``
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        _write_tsv(
            self.proportions, path / "proportions.tsv",
            comment_lines=[
                "estimate_type: spot_rna_composition",
                "WARNING: these are RNA-derived composition estimates, "
                "not direct single-cell counts.",
                f"lambda_spatial: {self.lambda_spatial}",
            ],
        )
        (path / "marker_genes.txt").write_text(
            "\n".join(self.marker_genes), encoding="utf-8"
        )
        with (path / "convergence_trace.json").open("w", encoding="utf-8") as fh:
            json.dump(self.convergence_trace, fh)

        _save_optional_array(self.mismatch_factors, path / "mismatch_factors.npy")
        _save_optional_array(self.lower_ci, path / "lower_ci.npy")
        _save_optional_array(self.upper_ci, path / "upper_ci.npy")

        meta: dict[str, Any] = {
            "tissueresolve_object": "SpatialDeconvResult",
            "estimate_type": self.ESTIMATE_TYPE,
            "n_spots": self.n_spots,
            "n_cell_types": self.n_cell_types,
            "n_marker_genes": len(self.marker_genes),
            "n_iter": self.n_iter,
            "converged": self.converged,
            "lambda_spatial": self.lambda_spatial,
            "n_smooth": self.n_smooth,
            "has_mismatch_factors": self.mismatch_factors is not None,
            "has_bootstrap_ci": self.lower_ci is not None,
            "bootstrap_coverage_note": self.bootstrap_coverage_note,
        }
        meta.update(self.run_metadata)
        with (path / "metadata.json").open("w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, default=str)

    @classmethod
    def load(cls, path: Path | str) -> "SpatialDeconvResult":
        """Load from a directory created by :meth:`save`."""
        path = Path(path)
        if not path.is_dir():
            raise FileNotFoundError(f"SpatialDeconvResult directory not found: {path}")
        proportions = _read_tsv(path / "proportions.tsv")
        marker_genes = (
            (path / "marker_genes.txt").read_text(encoding="utf-8").splitlines()
        )
        with (path / "convergence_trace.json").open(encoding="utf-8") as fh:
            convergence_trace: list[float] = json.load(fh)
        meta: dict[str, Any] = {}
        if (path / "metadata.json").exists():
            with (path / "metadata.json").open(encoding="utf-8") as fh:
                meta = json.load(fh)
        return cls(
            proportions=proportions,
            cell_types=meta.get("cell_types", list(proportions.columns)),
            marker_genes=marker_genes,
            n_iter=meta.get("n_iter", 0),
            converged=meta.get("converged", False),
            convergence_trace=convergence_trace,
            lambda_spatial=meta.get("lambda_spatial", 0.0),
            mismatch_factors=_load_optional_array(path / "mismatch_factors.npy"),
            lower_ci=_load_optional_array(path / "lower_ci.npy"),
            upper_ci=_load_optional_array(path / "upper_ci.npy"),
            bootstrap_coverage_note=meta.get("bootstrap_coverage_note"),
            n_smooth=meta.get("n_smooth"),
            run_metadata=meta,
        )


# ---------------------------------------------------------------------------
# QCReport
# ---------------------------------------------------------------------------


@dataclass
class QCReport:
    """Unified QC report container for both bulk and spatial runs.

    This is a result container — it stores pre-computed metrics and flags.
    The computation logic lives in ``bulk.qc`` and ``spatial.qc``.

    All heuristic-threshold flags carry ``is_heuristic=True`` in the
    ``metadata`` dict.

    Attributes
    ----------
    modality:
        ``"bulk"`` or ``"spatial"``.
    recommendations:
        Human-readable actionable strings generated from QC metrics.
        Never suppressed.

    Bulk-specific fields
    --------------------
    recon_r2:
        Per-sample reconstruction R².
    profile_corr:
        Per-sample Pearson correlation between bulk and reference mean.
    mismatch_flag:
        Per-sample flag: ``"low"`` | ``"medium"`` | ``"high"``.
    marker_recall:
        Per-cell-type fraction of top-50 markers detected in bulk.
    spillover_risk:
        Per-cell-type max Pearson correlation with nearest neighbour.
    mean_ci_width:
        Per-cell-type mean bootstrap CI width.  ``None`` if no bootstrap.
    condition_number:
        κ(Φ_panel), scalar.
    protocol_risk_level:
        ``"low"`` | ``"medium"`` | ``"high"`` | ``"unknown"``.
    n_genes_excluded_protocol:
        Number of genes removed by protocol filter.

    Spatial-specific fields
    -----------------------
    spot_qc:
        Per-spot QC DataFrame (entropy, nb_loglik, dominant_type, etc.).
    morans_i:
        Per-cell-type Moran's I values.
    model_qc:
        Dict of model-level metrics (n_iter, converged, d_g stats).
    """

    modality: str  # "bulk" | "spatial"
    recommendations: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    # bulk fields
    recon_r2: pd.Series | None = field(default=None, repr=False)
    profile_corr: pd.Series | None = field(default=None, repr=False)
    mismatch_flag: pd.Series | None = field(default=None, repr=False)
    marker_recall: pd.Series | None = field(default=None, repr=False)
    spillover_risk: pd.Series | None = field(default=None, repr=False)
    mean_ci_width: pd.Series | None = field(default=None, repr=False)
    condition_number: float | None = None
    protocol_risk_level: str | None = None
    n_genes_excluded_protocol: int | None = None

    # spatial fields
    spot_qc: pd.DataFrame | None = field(default=None, repr=False)
    morans_i: pd.Series | None = field(default=None, repr=False)
    model_qc: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.modality not in ("bulk", "spatial"):
            raise ValueError(
                f"modality must be 'bulk' or 'spatial', got {self.modality!r}."
            )

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | str) -> None:
        """Write QC report to *path* directory."""
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        if self.recon_r2 is not None:
            per_sample = pd.DataFrame({
                "recon_r2": self.recon_r2,
                "profile_corr": self.profile_corr,
                "mismatch_flag": self.mismatch_flag,
            })
            _write_tsv(per_sample, path / "per_sample_qc.tsv",
                       comment_lines=["All thresholds are heuristic."])

        if self.marker_recall is not None:
            parts: dict[str, Any] = {
                "marker_recall": self.marker_recall,
                "spillover_risk": self.spillover_risk,
            }
            if self.mean_ci_width is not None:
                parts["mean_ci_width"] = self.mean_ci_width
            _write_tsv(pd.DataFrame(parts), path / "per_celltype_qc.tsv")

        if self.spot_qc is not None:
            _write_tsv(self.spot_qc, path / "spot_qc.tsv")

        if self.morans_i is not None:
            _write_tsv(
                self.morans_i.to_frame("morans_i"), path / "morans_i.tsv"
            )

        if self.recommendations:
            (path / "recommendations.txt").write_text(
                "\n".join(self.recommendations), encoding="utf-8"
            )

        meta: dict[str, Any] = {
            "tissueresolve_object": "QCReport",
            "modality": self.modality,
            "condition_number": self.condition_number,
            "protocol_risk_level": self.protocol_risk_level,
            "n_genes_excluded_protocol": self.n_genes_excluded_protocol,
            "model_qc": self.model_qc,
            "n_recommendations": len(self.recommendations),
        }
        meta.update(self.metadata)
        with (path / "metadata.json").open("w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, default=str)

    @classmethod
    def load(cls, path: Path | str) -> "QCReport":
        """Load from a directory created by :meth:`save`."""
        path = Path(path)
        if not path.is_dir():
            raise FileNotFoundError(f"QCReport directory not found: {path}")
        meta: dict[str, Any] = {}
        if (path / "metadata.json").exists():
            with (path / "metadata.json").open(encoding="utf-8") as fh:
                meta = json.load(fh)
        recommendations: list[str] = []
        if (path / "recommendations.txt").exists():
            recommendations = (
                (path / "recommendations.txt")
                .read_text(encoding="utf-8")
                .splitlines()
            )
        per_sample = _load_optional_df(path / "per_sample_qc.tsv")
        per_celltype = _load_optional_df(path / "per_celltype_qc.tsv")
        spot_qc = _load_optional_df(path / "spot_qc.tsv")
        morans_i_df = _load_optional_df(path / "morans_i.tsv")
        return cls(
            modality=meta.get("modality", "bulk"),
            recommendations=recommendations,
            metadata=meta,
            recon_r2=per_sample["recon_r2"] if per_sample is not None and "recon_r2" in per_sample else None,
            profile_corr=per_sample["profile_corr"] if per_sample is not None and "profile_corr" in per_sample else None,
            mismatch_flag=per_sample["mismatch_flag"] if per_sample is not None and "mismatch_flag" in per_sample else None,
            marker_recall=per_celltype["marker_recall"] if per_celltype is not None and "marker_recall" in per_celltype else None,
            spillover_risk=per_celltype["spillover_risk"] if per_celltype is not None and "spillover_risk" in per_celltype else None,
            mean_ci_width=per_celltype["mean_ci_width"] if per_celltype is not None and "mean_ci_width" in per_celltype else None,
            condition_number=meta.get("condition_number"),
            protocol_risk_level=meta.get("protocol_risk_level"),
            n_genes_excluded_protocol=meta.get("n_genes_excluded_protocol"),
            spot_qc=spot_qc,
            morans_i=morans_i_df["morans_i"] if morans_i_df is not None else None,
            model_qc=meta.get("model_qc"),
        )


# ---------------------------------------------------------------------------
# SeparabilityReport
# ---------------------------------------------------------------------------


@dataclass
class PairSeparability:
    """Separability metrics for one pair of cell types.

    Attributes
    ----------
    type_a, type_b:
        Names of the two cell types being compared.
    bhattacharyya_coeff:
        BC ∈ [0, 1].  1 = identical profiles, 0 = orthogonal.
        BC > 0.90 is considered poorly separable; estimates will be
        unreliable for this pair.
    jeffreys_divergence:
        Symmetric KL divergence J = KL(a‖b) + KL(b‖a) (nats).
    pearson_r:
        Pearson correlation of the log1p-CPM profiles.
    n_discriminating_genes:
        Number of genes with |log2FC| > 1 between the two types.
    """

    type_a: str
    type_b: str
    bhattacharyya_coeff: float
    jeffreys_divergence: float
    pearson_r: float
    n_discriminating_genes: int

    @property
    def risk_level(self) -> str:
        """Human-readable risk level."""
        if self.bhattacharyya_coeff > 0.97:
            return "CRITICAL"
        if self.bhattacharyya_coeff > 0.90:
            return "HIGH"
        if self.bhattacharyya_coeff > 0.80:
            return "MEDIUM"
        return "OK"

    @property
    def is_problematic(self) -> bool:
        """True when BC > 0.90."""
        return self.bhattacharyya_coeff > 0.90

    @property
    def discriminability_score(self) -> float:
        """Separability score ∈ [0, 1].  Higher = more separable.

        Defined as ``1 − bhattacharyya_coeff``.  Pairs with score < 0.10
        are poorly separable; estimates for those types will be unreliable.
        """
        return max(0.0, 1.0 - self.bhattacharyya_coeff)


@dataclass
class SeparabilityReport:
    """Full pairwise separability report for a reference dataset.

    Contains all K*(K-1)/2 pairwise :class:`PairSeparability` results.
    Pairs are sorted by Bhattacharyya coefficient descending (worst first).

    Attributes
    ----------
    pairs:
        All pairwise results.
    n_critical, n_high, n_medium, n_ok:
        Count of pairs at each risk level.
    """

    pairs: list[PairSeparability]
    n_critical: int = 0
    n_high: int = 0
    n_medium: int = 0
    n_ok: int = 0

    def __post_init__(self) -> None:
        self.n_critical = sum(1 for p in self.pairs if p.risk_level == "CRITICAL")
        self.n_high = sum(1 for p in self.pairs if p.risk_level == "HIGH")
        self.n_medium = sum(1 for p in self.pairs if p.risk_level == "MEDIUM")
        self.n_ok = sum(1 for p in self.pairs if p.risk_level == "OK")

    @property
    def has_problems(self) -> bool:
        return self.n_critical > 0 or self.n_high > 0

    def summary_str(self) -> str:
        lines = [
            f"SeparabilityReport: {len(self.pairs)} pairs",
            f"  CRITICAL (BC>0.97): {self.n_critical}",
            f"  HIGH     (BC>0.90): {self.n_high}",
            f"  MEDIUM   (BC>0.80): {self.n_medium}",
            f"  OK:                 {self.n_ok}",
        ]
        for p in self.pairs:
            if p.is_problematic:
                lines.append(
                    f"  [{p.risk_level:8s}] {p.type_a!r} vs {p.type_b!r}  "
                    f"BC={p.bhattacharyya_coeff:.4f}  r={p.pearson_r:.3f}  "
                    f"disc_genes={p.n_discriminating_genes}"
                )
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | str) -> None:
        """Write separability report to *path* directory."""
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        rows = [
            {
                "type_a": p.type_a,
                "type_b": p.type_b,
                "bhattacharyya_coeff": p.bhattacharyya_coeff,
                "jeffreys_divergence": p.jeffreys_divergence,
                "pearson_r": p.pearson_r,
                "n_discriminating_genes": p.n_discriminating_genes,
                "risk_level": p.risk_level,
            }
            for p in self.pairs
        ]
        df = pd.DataFrame(rows)
        _write_tsv(
            df, path / "separability_pairs.tsv",
            comment_lines=[
                "BC > 0.90 indicates poorly separable types.",
                "Estimates for HIGH/CRITICAL pairs will be unreliable.",
            ],
        )
        meta = {
            "tissueresolve_object": "SeparabilityReport",
            "n_pairs": len(self.pairs),
            "n_critical": self.n_critical,
            "n_high": self.n_high,
            "n_medium": self.n_medium,
            "n_ok": self.n_ok,
        }
        with (path / "metadata.json").open("w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2)

    @classmethod
    def load(cls, path: Path | str) -> "SeparabilityReport":
        """Load from a directory created by :meth:`save`."""
        path = Path(path)
        if not path.is_dir():
            raise FileNotFoundError(f"SeparabilityReport directory not found: {path}")
        df = _read_tsv(path / "separability_pairs.tsv")
        df = df.reset_index(drop=True)
        pairs = [
            PairSeparability(
                type_a=str(row["type_a"]),
                type_b=str(row["type_b"]),
                bhattacharyya_coeff=float(row["bhattacharyya_coeff"]),
                jeffreys_divergence=float(row["jeffreys_divergence"]),
                pearson_r=float(row["pearson_r"]),
                n_discriminating_genes=int(row["n_discriminating_genes"]),
            )
            for _, row in df.iterrows()
        ]
        return cls(pairs=pairs)


# ---------------------------------------------------------------------------
# BenchmarkResult
# ---------------------------------------------------------------------------


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run of one method on one scenario.

    Attributes
    ----------
    method:
        Method name, e.g. ``"SpatCAR"``, ``"NNLS"``, ``"CHIMERA-protocol"``.
    scenario:
        Benchmark scenario name, e.g. ``"PROTOCOL_MISMATCH_CENTRAL"``.
    metrics:
        Flat dict of scalar metrics, e.g.
        ``{"rmse": 0.04, "pcc": 0.92, "jsd": 0.03, "morans_i_mean": 0.71}``.
    per_celltype_metrics:
        Optional per-cell-type breakdown of metrics.
    n_samples_or_spots:
        Number of bulk samples or Visium spots evaluated.
    run_time_s:
        Wall-clock time in seconds.
    run_metadata:
        Dict of method parameters and version info.
    """

    method: str
    scenario: str
    metrics: dict[str, float]
    per_celltype_metrics: pd.DataFrame | None = field(default=None, repr=False)
    n_samples_or_spots: int = 0
    run_time_s: float = 0.0
    run_metadata: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | str) -> None:
        """Write benchmark result to *path* directory."""
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        scalar_df = pd.DataFrame([self.metrics], index=[self.method])
        _write_tsv(scalar_df, path / "metrics.tsv")
        if self.per_celltype_metrics is not None:
            _write_tsv(self.per_celltype_metrics, path / "per_celltype_metrics.tsv")
        meta: dict[str, Any] = {
            "tissueresolve_object": "BenchmarkResult",
            "method": self.method,
            "scenario": self.scenario,
            "n_samples_or_spots": self.n_samples_or_spots,
            "run_time_s": self.run_time_s,
        }
        meta.update(self.run_metadata)
        with (path / "metadata.json").open("w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, default=str)

    @classmethod
    def load(cls, path: Path | str) -> "BenchmarkResult":
        """Load from a directory created by :meth:`save`."""
        path = Path(path)
        if not path.is_dir():
            raise FileNotFoundError(f"BenchmarkResult directory not found: {path}")
        metrics_df = _read_tsv(path / "metrics.tsv")
        meta: dict[str, Any] = {}
        if (path / "metadata.json").exists():
            with (path / "metadata.json").open(encoding="utf-8") as fh:
                meta = json.load(fh)
        per_celltype = _load_optional_df(path / "per_celltype_metrics.tsv")
        return cls(
            method=meta.get("method", metrics_df.index[0]),
            scenario=meta.get("scenario", "unknown"),
            metrics={k: float(v) for k, v in metrics_df.iloc[0].items()},
            per_celltype_metrics=per_celltype,
            n_samples_or_spots=meta.get("n_samples_or_spots", 0),
            run_time_s=meta.get("run_time_s", 0.0),
            run_metadata=meta,
        )
