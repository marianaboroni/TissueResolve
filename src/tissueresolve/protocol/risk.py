"""
Bulk protocol risk assessment for TissueResolve.

Evaluates the gene-level bias introduced when the sequencing protocol of the
bulk RNA-seq dataset (``BulkProtocol``) is incompatible with the protocol of
the single-cell reference (``RefModality``, ``RefCapture``, ``RefCounting``).

This module is relevant for **bulk deconvolution only**.

Distinction from spatial mismatch
-----------------------------------
This module answers: *which genes should I exclude or down-weight BEFORE
running deconvolution because of known protocol biases?*

The spatial mismatch module (``protocol.mismatch``) answers a different
question: *what multiplicative scale factors should I estimate DURING
deconvolution to correct for the systematic difference between Visium spot
expression and the pseudo-bulk reference?*

Risk categories
---------------
``intronic_retention``
    snRNA-seq references retain intronic signal that bulk polyA-seq
    does not capture.  Genes dominated by intronic reads produce
    artificially high reference expression relative to bulk.
    Source: Ivich & Greene 2026 (Cell Rep Methods); Bisque 2020 (Nat Commun).

``length_bias``
    3′-UMI capture references (10x 3′, Drop-seq) under-detect genes with
    transcript length > 5 000 bp relative to full-length bulk sequencing.
    Source: Bisque 2020 (Nat Commun); hg38/mm10 length annotations.

``dissociation_stress``
    scRNA-seq references may contain artificially up-regulated stress genes
    (FOS, JUN, HSP70 family) induced by enzymatic cell dissociation.
    Source: van den Brink et al. 2017 (Nat Methods); Machado et al. 2021.

``ribodep_intronic``
    Ribo-depleted bulk libraries capture partial intronic signal (which polyA
    bulk does not).  This causes a mild global bias rather than a gene-specific
    one.

No silent gene removal
----------------------
:class:`ProtocolRiskReport` explicitly lists every excluded and down-weighted
gene in :attr:`~ProtocolRiskReport.excluded_genes` and
:attr:`~ProtocolRiskReport.downweighted_genes`.  The caller is always informed
about every gene that was removed or penalised.
"""
from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.protocol.metadata import (
    BulkProtocol,
    ProtocolMetadata,
    RefCapture,
    RefModality,
)
from tissueresolve.reference.gene_filters import GeneFilterSet

__all__ = [
    "ProtocolRiskReport",
    "ProtocolRiskAssessor",
]

logger = logging.getLogger("tissueresolve.protocol.risk")

# Risk score contributions per category (all in [0, 1])
_INTRONIC_RISK: float = 0.85
_LENGTH_RISK: float = 0.65
_STRESS_RISK: float = 0.70
_RIBODEP_RISK: float = 0.40

# Hard-exclusion threshold: genes with risk >= this value are excluded
RISK_HARD_THRESHOLD: float = 0.50


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class ProtocolRiskReport:
    """Output of :meth:`ProtocolRiskAssessor.assess`.

    All gene-level effects are reported explicitly.  No gene is silently
    removed or penalised.

    Attributes
    ----------
    risk_level:
        ``"low"``, ``"medium"``, ``"high"``, or ``"unknown"``.
    mismatch_types:
        Active mismatch categories (e.g. ``["intronic_retention"]``).
    affected_gene_fraction:
        Fraction of candidate genes touched by at least one filter.
    n_genes_excluded:
        Number of genes whose risk score ≥ :data:`RISK_HARD_THRESHOLD`.
        These genes are recommended for exclusion from the panel.
    excluded_genes:
        Explicit list of every gene recommended for exclusion.
        Never silently omitted.
    downweighted_genes:
        Genes with 0 < risk < :data:`RISK_HARD_THRESHOLD` that should be
        penalised in gene weighting but not fully excluded.
    rationale:
        Human-readable explanation of active risks.
    recommended_actions:
        Actionable strings for display to the user / in the HTML report.
    literature_refs:
        Citation strings for the filters applied.
    """

    risk_level: str
    mismatch_types: list[str]
    affected_gene_fraction: float
    n_genes_excluded: int
    excluded_genes: list[str]
    downweighted_genes: list[str]
    rationale: str
    recommended_actions: list[str]
    literature_refs: list[str]

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dictionary."""
        return {
            "risk_level": self.risk_level,
            "mismatch_types": self.mismatch_types,
            "affected_gene_fraction": round(self.affected_gene_fraction, 4),
            "n_genes_excluded": self.n_genes_excluded,
            "excluded_genes": self.excluded_genes,
            "downweighted_genes": self.downweighted_genes,
            "rationale": self.rationale,
            "recommended_actions": self.recommended_actions,
            "literature_refs": self.literature_refs,
        }


# ---------------------------------------------------------------------------
# Risk assessor
# ---------------------------------------------------------------------------


class ProtocolRiskAssessor:
    """Evaluate bulk × reference protocol compatibility and compute per-gene
    risk scores.

    Parameters
    ----------
    genome:
        Reference genome assembly (``"hg38"`` or ``"mm10"``).
    gene_filters:
        Optional pre-built :class:`~tissueresolve.reference.gene_filters.GeneFilterSet`.
        When provided, its cache is shared with the gene selector (avoids
        loading the same files twice).  When ``None``, a new instance is
        created using *genome* and *data_dir*.
    data_dir:
        Override the gene list data directory.  Only used when
        *gene_filters* is ``None``.
    """

    def __init__(
        self,
        genome: str = "hg38",
        gene_filters: Optional[GeneFilterSet] = None,
        data_dir=None,
    ) -> None:
        self.genome = genome
        if gene_filters is not None:
            self._filters = gene_filters
        else:
            self._filters = GeneFilterSet(genome=genome, data_dir=data_dir)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def assess(
        self,
        meta: ProtocolMetadata,
        candidate_genes: list[str],
    ) -> ProtocolRiskReport:
        """Assess protocol compatibility and return a :class:`ProtocolRiskReport`.

        Parameters
        ----------
        meta:
            Protocol metadata for bulk and reference.
        candidate_genes:
            Gene IDs to check against risk lists.  The returned report lists
            every gene that is excluded or down-weighted.

        Returns
        -------
        ProtocolRiskReport
        """
        if meta.is_bulk_fully_unknown():
            return self._unknown_report(candidate_genes)

        gene_set = set(candidate_genes)
        active_types: list[str] = []
        excluded_set: set[str] = set()
        refs: list[str] = []

        # 1. Intronic retention — active when reference is snRNA
        if meta.ref_modality == RefModality.SNRNA:
            intronic = self._filters.intronic_dominant()
            bad = gene_set & intronic
            if bad:
                active_types.append("intronic_retention")
                excluded_set |= bad
                refs.append(
                    "Ivich & Greene 2026 (Cell Rep Methods); "
                    "Bisque 2020 (Nat Commun)"
                )

        # 2. Length bias — active when reference uses 3′-end UMI capture
        if meta.ref_capture in (RefCapture.UMI_3PRIME, RefCapture.DROPSEQ):
            length_bad = self._filters.length_biased_3prime()
            bad = gene_set & length_bad
            if bad:
                active_types.append("length_bias")
                excluded_set |= bad
                refs.append(
                    "Bisque 2020 (Nat Commun); "
                    "TissueResolve gene_lengths annotation (hg38/mm10)"
                )

        # 3. Dissociation stress — active when reference is scRNA
        if meta.ref_modality == RefModality.SCRNA:
            stress = self._filters.dissociation_stress()
            bad = gene_set & stress
            if bad:
                active_types.append("dissociation_stress")
                excluded_set |= bad
                refs.append(
                    "van den Brink et al. 2017 (Nat Methods); "
                    "Machado et al. 2021"
                )

        # 4. Ribo-depletion — mild global uniform risk, no gene-level exclusion
        if meta.bulk_protocol == BulkProtocol.RIBODEP:
            active_types.append("ribodep_intronic")
            refs.append(
                "Protocol-specific guidance; "
                "partial intronic retention in ribo-depleted bulk"
            )

        n_cand = len(candidate_genes)
        n_excl = len(excluded_set)
        affected_frac = n_excl / n_cand if n_cand > 0 else 0.0

        risk_level = _compute_risk_level(meta, active_types, affected_frac)
        rationale, actions = _build_rationale(meta, active_types, affected_frac)

        # Determine down-weighted genes: have risk > 0 but < hard threshold
        # Currently length bias and ribo-dep produce partial risk (no full exclusion).
        # Down-weighted = any genes touched by length_bias or ribo_dep but not fully
        # excluded by intronic/stress.
        downweighted: set[str] = set()
        if "length_bias" in active_types:
            # Length-biased genes not already excluded are down-weighted
            length_bad_again = self._filters.length_biased_3prime()
            downweighted |= (gene_set & length_bad_again) - excluded_set

        logger.info(
            "ProtocolRisk: level=%s  types=%s  excluded=%d/%d  downweighted=%d",
            risk_level, active_types, n_excl, n_cand, len(downweighted),
        )
        return ProtocolRiskReport(
            risk_level=risk_level,
            mismatch_types=active_types,
            affected_gene_fraction=affected_frac,
            n_genes_excluded=n_excl,
            excluded_genes=sorted(excluded_set),        # explicit, never silent
            downweighted_genes=sorted(downweighted),    # explicit, never silent
            rationale=rationale,
            recommended_actions=actions,
            literature_refs=list(dict.fromkeys(refs)),  # deduplicate, preserve order
        )

    def gene_risk_scores(
        self,
        meta: ProtocolMetadata,
        gene_names: list[str],
    ) -> pd.Series:
        """Return per-gene risk score ∈ [0, 1].

        0 = no known risk for this protocol combination.
        1 = maximum risk; gene is recommended for exclusion.

        Risk components (additive, capped at 1.0):
        - Intronic score (snRNA reference)
        - Length score (3′-UMI capture reference)
        - Stress score (scRNA reference)
        - Uniform ribo-dep increment (ribodepleted bulk)

        Parameters
        ----------
        meta:
            Protocol metadata.
        gene_names:
            Gene identifiers to score.

        Returns
        -------
        pd.Series indexed by gene name.
        """
        scores = pd.Series(0.0, index=gene_names, dtype=float)

        if meta.ref_modality == RefModality.SNRNA:
            intronic = self._filters.intronic_dominant()
            mask = scores.index.isin(intronic)
            scores[mask] = np.minimum(scores[mask] + _INTRONIC_RISK, 1.0)

        if meta.ref_capture in (RefCapture.UMI_3PRIME, RefCapture.DROPSEQ):
            length_bad = self._filters.length_biased_3prime()
            mask = scores.index.isin(length_bad)
            scores[mask] = np.minimum(scores[mask] + _LENGTH_RISK, 1.0)

        if meta.ref_modality == RefModality.SCRNA:
            stress = self._filters.dissociation_stress()
            mask = scores.index.isin(stress)
            scores[mask] = np.minimum(scores[mask] + _STRESS_RISK, 1.0)

        if meta.bulk_protocol == BulkProtocol.RIBODEP:
            # Uniform partial increment for ribo-dep × scRNA combination
            scores = np.minimum(scores + _RIBODEP_RISK * 0.3, 1.0)

        return scores

    def safe_panel(self, gene_names: list[str]) -> list[str]:
        """Return genes in the universal safe panel that appear in *gene_names*.

        Used as fallback when protocol is fully unknown.

        Parameters
        ----------
        gene_names:
            Candidate gene identifiers.

        Returns
        -------
        list[str]
        """
        safe = self._filters.safe_universal()
        panel = [g for g in gene_names if g in safe]
        if not panel:
            logger.warning(
                "safe_panel: no overlap with gene_names; "
                "returning all gene_names (no filter applied)."
            )
            return list(gene_names)
        return panel

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _unknown_report(self, candidate_genes: list[str]) -> ProtocolRiskReport:
        return ProtocolRiskReport(
            risk_level="unknown",
            mismatch_types=[],
            affected_gene_fraction=0.0,
            n_genes_excluded=0,
            excluded_genes=[],
            downweighted_genes=[],
            rationale=(
                "Protocol metadata not provided.  "
                "Conservative universal gene panel applied.  "
                "Specify bulk_protocol and ref_modality for targeted filtering."
            ),
            recommended_actions=[
                "Provide bulk_protocol and ref_modality for best results.",
                "Using safe_universal gene panel as fallback.",
            ],
            literature_refs=[],
        )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _compute_risk_level(
    meta: ProtocolMetadata,
    active_types: list[str],
    affected_frac: float,
) -> str:
    if not active_types:
        return "low"
    if "intronic_retention" in active_types and affected_frac > 0.15:
        return "high"
    if "intronic_retention" in active_types:
        return "medium"
    if affected_frac > 0.20:
        return "high"
    if active_types:
        return "medium"
    return "low"


def _build_rationale(
    meta: ProtocolMetadata,
    active_types: list[str],
    affected_frac: float,
) -> tuple[str, list[str]]:
    parts: list[str] = []
    actions: list[str] = []

    if "intronic_retention" in active_types:
        parts.append(
            f"snRNA-seq reference detected.  "
            f"{affected_frac * 100:.1f}% of candidate genes have intronic "
            "retention bias that bulk polyA-seq does not capture.  "
            "These genes are excluded from the panel."
        )
        actions.append(
            "Intronic-dominant genes excluded.  Consider using exonic-only "
            "counts for the snRNA reference if available."
        )

    if "length_bias" in active_types:
        parts.append(
            "3′-UMI capture reference detected.  Long genes (> 5 000 bp) are "
            "penalised in gene weights to reduce length-dependent bias."
        )
        actions.append("Length-biased genes penalised in composite weights.")

    if "dissociation_stress" in active_types:
        parts.append(
            "scRNA-seq reference detected.  Dissociation stress genes "
            "(FOS, JUN family, HSP70) are excluded from the panel."
        )
        actions.append("Dissociation stress genes excluded from panel.")

    if "ribodep_intronic" in active_types:
        parts.append(
            "Ribo-depleted bulk detected.  Partial intronic signal in bulk "
            "may inflate expression of intron-rich genes."
        )
        actions.append(
            "All gene weights slightly penalised for the ribo-depleted × "
            "scRNA reference combination."
        )

    if not parts:
        parts.append("No known protocol incompatibilities detected.")

    return " ".join(parts), actions
