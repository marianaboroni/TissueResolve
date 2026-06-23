"""
TissueResolve configuration.

Single source of truth for every tunable parameter.  All thresholds that
are heuristic are explicitly labelled as such in docstrings and in
``run_metadata.json`` under ``"is_heuristic": true``.

Hierarchy
---------
TissueResolveConfig
  ├── reference: ReferenceConfig
  ├── genes: GeneConfig
  ├── discordance: DiscordanceConfig   (bulk Tier-3 only)
  ├── bulk_solver: BulkSolverConfig
  ├── spatial_solver: SpatialSolverConfig
  ├── bootstrap: BootstrapConfig
  ├── bulk_qc: BulkQCConfig
  └── spatial_qc: SpatialQCConfig
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

__all__ = [
    "TissueResolveConfig",
    "ReferenceConfig",
    "GeneConfig",
    "DiscordanceConfig",
    "BulkSolverConfig",
    "SpatialSolverConfig",
    "BootstrapConfig",
    "BulkQCConfig",
    "SpatialQCConfig",
    "HierarchicalConfig",
]


# ---------------------------------------------------------------------------
# Sub-configs
# ---------------------------------------------------------------------------


@dataclass
class ReferenceConfig:
    """Settings for reference matrix construction.

    Attributes
    ----------
    celltype_col:
        obs column holding cell-type labels.
    donor_col:
        obs column holding donor IDs.  When present, profiles are computed
        as the mean of per-donor averages (donor-aware aggregation).
        ``None`` → simple cell mean.
    min_cells:
        Minimum cells required to retain a cell type.
    genome:
        Reference genome assembly (``"hg38"`` or ``"mm10"``).  Determines
        which gene filter lists are loaded.
    """

    celltype_col: str = "cell_type"
    donor_col: str | None = "donor"
    min_cells: int = 10
    genome: str = "hg38"


@dataclass
class GeneConfig:
    """Gene-panel selection settings.

    Composite weight for gene *g* (bulk)::

        weight_g = w_spec  × specificity_g
                 × w_stab  × stability_g
                 × w_prot  × (1 − protocol_risk_g)
                 × w_conc  × bulk_concordance_g

    Attributes
    ----------
    n_genes:
        Maximum panel size.
    min_log2fc:
        Minimum pairwise log₂FC for a gene to be considered a marker.
    min_mean_cpm:
        Minimum CPM-like expression in the top cell type.
    max_donor_cv:
        Maximum cross-donor CV.  High-CV genes are unstable markers.
        Set to ``1e9`` to disable.
    weight_specificity, weight_stability, weight_protocol, weight_concordance:
        Relative contributions of each component in the composite score.
    n_per_type_spatial:
        Maximum marker genes per cell type for spatial (log2FC selection).
    min_lfc_spatial:
        Minimum log2FC for spatial marker selection.
    """

    n_genes: int = 500
    min_log2fc: float = 1.0
    min_mean_cpm: float = 0.5
    max_donor_cv: float = 1.0
    weight_specificity: float = 1.0
    weight_stability: float = 0.5
    weight_protocol: float = 1.0
    weight_concordance: float = 0.5
    n_per_type_spatial: int = 100
    min_lfc_spatial: float = 1.0


@dataclass
class DiscordanceConfig:
    """Paired bulk↔pseudobulk discordance filter (bulk Tier-3 only).

    Requires ≥ 2 paired donors with both bulk and pseudobulk data.

    Attributes
    ----------
    fdr:
        Benjamini-Hochberg FDR threshold.
    log2fc_threshold:
        Minimum |mean log₂FC| between bulk and pseudobulk.
    require_both:
        If True, a gene is flagged only when both FDR and log₂FC criteria
        are met simultaneously.
    """

    fdr: float = 0.05
    log2fc_threshold: float = 1.0
    require_both: bool = True


@dataclass
class BulkSolverConfig:
    """Weighted NNLS solver settings for bulk deconvolution.

    Attributes
    ----------
    seed:
        Random seed (passed to bootstrap; the solver itself is deterministic).
    """

    seed: int = 42


@dataclass
class SpatialSolverConfig:
    """NB-CAR model settings for spatial deconvolution.

    Attributes
    ----------
    lambda_spatial:
        CAR spatial regularisation strength λ ≥ 0.  ``0`` disables spatial
        smoothing (pure NB-MAP per spot).
    max_iter:
        Maximum outer optimisation iterations.
    tol:
        Convergence tolerance on ``‖ΔΠ‖_F / N``.
    n_marker_genes:
        Maximum marker genes per cell type for auto-selection.
    min_lfc:
        Minimum log2FC for auto marker-gene selection.
    update_mismatch_every:
        Update per-gene mismatch scale factors every this many iterations.
    n_jobs:
        Parallel workers for NNLS warm-start.  ``-1`` uses all CPUs.
    random_state:
        Integer seed for reproducibility.
    """

    lambda_spatial: float = 0.1
    max_iter: int = 200
    tol: float = 5e-5
    n_marker_genes: int = 100
    min_lfc: float = 1.0
    update_mismatch_every: int = 5
    n_jobs: int = 1
    random_state: int = 42


@dataclass
class BootstrapConfig:
    """Bootstrap confidence interval settings (both bulk and spatial).

    Attributes
    ----------
    n_bootstrap:
        Number of bootstrap iterations.  ``0`` → skip CIs.
    bootstrap_frac:
        Fraction of genes resampled per iteration (bulk).
    ci_level:
        Nominal CI level (0.95 → 2.5%–97.5% percentiles).
    seed:
        Random seed.
    n_iter_per_boot_spatial:
        Re-fitting iterations per spatial bootstrap resample.
        Fewer iterations → faster but wider empirical under-coverage.
        See ``SpatialDeconvResult.bootstrap_coverage_note``.
    """

    n_bootstrap: int = 200
    bootstrap_frac: float = 0.8
    ci_level: float = 0.95
    seed: int = 42
    n_iter_per_boot_spatial: int = 30


@dataclass
class BulkQCConfig:
    """QC flag thresholds for bulk deconvolution.

    All thresholds are **heuristic** values.  They are not universal
    calibrated decision boundaries.  Every threshold is written to
    ``run_metadata.json`` with ``"is_heuristic": true``.

    Sources are documented in ``docs/qc_thresholds.md``.

    Attributes
    ----------
    r2_warn, r2_fail:
        Reconstruction R² thresholds.
    profile_corr_warn, profile_corr_fail:
        Pearson correlation between bulk sample and reference column mean.
    marker_recall_warn, marker_recall_fail:
        Fraction of top-50 markers per cell type expressed ≥ 0.5 CPM.
    condition_number_warn, condition_number_fail:
        κ(Φ_panel).  Sensitive to linear methods; regularised methods are
        less affected.
    bootstrap_cv_warn, bootstrap_cv_fail:
        CV of bootstrap proportion estimates per cell type.
    """

    r2_warn: float = 0.50
    r2_fail: float = 0.40
    profile_corr_warn: float = 0.60
    profile_corr_fail: float = 0.45
    marker_recall_warn: float = 0.55
    marker_recall_fail: float = 0.35
    condition_number_warn: float = 500.0
    condition_number_fail: float = 2000.0
    bootstrap_cv_warn: float = 0.35
    bootstrap_cv_fail: float = 0.60


@dataclass
class SpatialQCConfig:
    """QC flag thresholds for spatial deconvolution.

    ``None`` values trigger automatic threshold derivation from the data
    (percentile-based) at run time.

    Attributes
    ----------
    entropy_threshold:
        Spots with Shannon entropy of proportions **above** this are flagged
        as uncertain.  ``None`` → 95th-percentile of the section.
    loglik_threshold:
        Spots with NB log-likelihood **below** this are flagged as poor fit.
        ``None`` → 5th-percentile of the section.
    spatial_resid_threshold:
        Spots with spatial residual ``‖π_s − π̄_neighbours‖₂`` **above** this
        are flagged.  ``None`` → 3× median of the section.
    morans_i_warn:
        Cell types with Moran's I below this value are flagged as lacking
        spatial structure (expected for biologically real types).
    """

    entropy_threshold: float | None = None
    loglik_threshold: float | None = None
    spatial_resid_threshold: float | None = None
    morans_i_warn: float = 0.05


@dataclass
class HierarchicalConfig:
    """Broad-to-fine (hierarchical) deconvolution settings.

    Hierarchical mode first estimates broad cell-type families, then estimates
    fine subpopulations *within* each family.  Fine subtype splits are only
    trusted when the subtypes are demonstrably separable within their family;
    otherwise the family's mass is reported as ``unresolved_<family>`` rather
    than split into subtypes for which there is no evidence.

    All thresholds below are **heuristic**.

    Attributes
    ----------
    broad_cell_type_col:
        obs column with broad/compartment labels.  ``"auto"`` → detect from a
        list of known candidates; ``None`` → not provided.
    fine_cell_type_col:
        obs column with fine/subpopulation labels.  ``"auto"`` → detect.
    allow_unresolved:
        When True (default), families whose subtypes are not separable keep
        their mass at the broad level (``unresolved_<family>``); when False,
        fine splits are always produced (and a warning is emitted).
    unresolved_threshold:
        A family is treated as unresolved when its mean within-family
        separability score (``1 − Bhattacharyya``) is **below** this value.
    min_discriminating_genes:
        A family is treated as unresolved when any within-family pair has fewer
        than this many discriminating genes (|log2FC| > 1).
    within_family_spillover_threshold:
        A family is treated as unresolved when its mean within-family spillover
        (max correlation to a family sibling) is **above** this value.
    within_family_marker_selection:
        ``"auto"`` (default), ``"all"`` (use the shared panel), or ``"pairwise"``
        (augment with pairwise within-family discriminative genes).
    hierarchy_level:
        ``"fine"``, ``"family"``, or ``"both"`` (default) — which estimates to
        emphasise in outputs/reports.  All levels are always saved.
    """

    broad_cell_type_col: str | None = "auto"
    fine_cell_type_col: str | None = "auto"
    allow_unresolved: bool = True
    unresolved_threshold: float = 0.10
    min_discriminating_genes: int = 10
    within_family_spillover_threshold: float = 0.30
    within_family_marker_selection: str = "auto"
    hierarchy_level: str = "both"
    # partial resolution: assign confident subtype mass and keep only the
    # ambiguous remainder as unresolved_<family> (not all-or-nothing).
    allow_partial_resolution: bool = True
    subtype_confidence_threshold: float = 0.10
    # Gating mode for hierarchical inference (validated on breast + lung benchmarks):
    #   "soft"    — DEFAULT. Partial confidence-weighted unresolved mass: each
    #               subtype's mass is multiplied by a calibrated confidence in
    #               [0,1] and the residual family mass goes to unresolved_<family>.
    #   "hard"    — LEGACY. Binary threshold gate (confident subtypes keep full
    #               mass, the rest are zeroed). Over-abstains in collinear families.
    #   "ungated" — DIAGNOSTIC only. No abstention; not calibrated.
    hierarchical_gating: str = "soft"
    gating_version: str = "soft_gating-1.0"


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------


@dataclass
class TissueResolveConfig:
    """Top-level TissueResolve configuration.

    Serialises to / deserialises from YAML.  The output directory and
    verbosity flag are the only top-level scalars; all algorithm parameters
    live in sub-configs.
    """

    reference: ReferenceConfig = field(default_factory=ReferenceConfig)
    genes: GeneConfig = field(default_factory=GeneConfig)
    discordance: DiscordanceConfig = field(default_factory=DiscordanceConfig)
    bulk_solver: BulkSolverConfig = field(default_factory=BulkSolverConfig)
    spatial_solver: SpatialSolverConfig = field(default_factory=SpatialSolverConfig)
    bootstrap: BootstrapConfig = field(default_factory=BootstrapConfig)
    bulk_qc: BulkQCConfig = field(default_factory=BulkQCConfig)
    spatial_qc: SpatialQCConfig = field(default_factory=SpatialQCConfig)
    hierarchical: HierarchicalConfig = field(default_factory=HierarchicalConfig)
    output_dir: Path = field(default_factory=lambda: Path("tissueresolve_out"))
    verbose: bool = True

    # ------------------------------------------------------------------
    # YAML I/O
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: "Path | str") -> "TissueResolveConfig":
        """Load from a YAML file.

        Unknown keys are silently ignored so that configs written by a
        newer version can be loaded by an older one.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")
        with path.open("r", encoding="utf-8") as fh:
            data: dict[str, Any] = yaml.safe_load(fh) or {}
        return cls._from_dict(data)

    def to_yaml(self, path: "Path | str") -> None:
        """Persist to a YAML file.  Parent directories are created."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as fh:
            yaml.safe_dump(self._to_dict(), fh, sort_keys=False)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        # Ensure output_dir is always a Path object even after YAML load.
        self.output_dir = Path(self.output_dir)

    @classmethod
    def _from_dict(cls, d: dict[str, Any]) -> "TissueResolveConfig":
        sub_map: dict[str, type] = {
            "reference": ReferenceConfig,
            "genes": GeneConfig,
            "discordance": DiscordanceConfig,
            "bulk_solver": BulkSolverConfig,
            "spatial_solver": SpatialSolverConfig,
            "bootstrap": BootstrapConfig,
            "bulk_qc": BulkQCConfig,
            "spatial_qc": SpatialQCConfig,
            "hierarchical": HierarchicalConfig,
        }
        kwargs: dict[str, Any] = {}
        for key, sub_cls in sub_map.items():
            raw = d.get(key, {})
            if isinstance(raw, dict):
                valid = {
                    k: v for k, v in raw.items()
                    if k in sub_cls.__dataclass_fields__  # type: ignore[attr-defined]
                }
                kwargs[key] = sub_cls(**valid)
            else:
                kwargs[key] = sub_cls()
        if "output_dir" in d:
            kwargs["output_dir"] = Path(d["output_dir"])
        if "verbose" in d:
            kwargs["verbose"] = bool(d["verbose"])
        return cls(**kwargs)

    def _to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["output_dir"] = str(self.output_dir)
        return d

    def heuristic_thresholds(self) -> dict[str, dict]:
        """Return all heuristic QC thresholds with provenance metadata.

        Suitable for embedding in ``run_metadata.json``.
        """
        qb = self.bulk_qc
        return {
            "r2_warn":               {"value": qb.r2_warn,            "is_heuristic": True,
                                      "source": "Finotello 2019; BAL 2026"},
            "r2_fail":               {"value": qb.r2_fail,            "is_heuristic": True,
                                      "source": "same"},
            "profile_corr_warn":     {"value": qb.profile_corr_warn,  "is_heuristic": True,
                                      "source": "empirical; see docs/qc_thresholds.md"},
            "profile_corr_fail":     {"value": qb.profile_corr_fail,  "is_heuristic": True,
                                      "source": "same"},
            "marker_recall_warn":    {"value": qb.marker_recall_warn, "is_heuristic": True,
                                      "source": "empirical"},
            "marker_recall_fail":    {"value": qb.marker_recall_fail, "is_heuristic": True,
                                      "source": "empirical"},
            "condition_number_warn": {"value": qb.condition_number_warn, "is_heuristic": True,
                                      "source": "numerical analysis; not tissue-validated"},
            "condition_number_fail": {"value": qb.condition_number_fail, "is_heuristic": True,
                                      "source": "same"},
        }
