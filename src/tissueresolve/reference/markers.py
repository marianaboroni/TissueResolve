"""
Marker gene selection for TissueResolve.

Combines CHIMERA's composite log-additive scoring with SpatCAR's pairwise
discriminability augmentation.

Selection strategy
------------------
Phase 1 — per-type composite scoring (from CHIMERA)
    For each gene, a composite score is computed from four components:

    * **specificity** — pairwise max log₂FC across all type pairs.
    * **stability** — ``1 / (1 + mean_donor_cv)`` (0 when donor_cv is absent).
    * **protocol_weight** — ``1 − protocol_risk_score`` (1.0 by default; filled
      in during Stage 2 when the protocol module is available).
    * **concordance** — similarity between gene expression in bulk and the
      reference (1.0 by default; filled in during Stage 3 pipeline build).

    Log-additive combination::

        score = w_spec  × log(specificity + ε)
              + w_stab  × log(stability    + ε)
              + w_prot  × log(prot_weight  + ε)
              + w_conc  × log(concordance  + ε)

    Hard filters (genes failing these are excluded before scoring):
    * Pairwise log₂FC < ``min_log2fc``
    * Max CPM < ``min_mean_cpm``
    * Mean donor CV > ``max_donor_cv``
    * Blacklist prefix match
    * Not protein-coding (when protein-coding set is available)
    * Hypervariable inflammatory gene

Phase 2 — pairwise discriminability augmentation (from SpatCAR)
    For each pair of cell types, the top ``n_per_pair`` genes ranked by
    |log₂FC(type_a, type_b)| are added to the panel.  This is critical when
    two types are similar overall but differ on a small gene set.

The final panel is the union of Phase 1 and Phase 2 genes, optionally
filtered to ``query_genes`` (e.g. genes present in a Visium dataset).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.config import GeneConfig
from tissueresolve.results import ReferenceSignature

__all__ = [
    "MarkerSelectionResult",
    "GeneSelector",
]

logger = logging.getLogger("tissueresolve.reference.markers")

_LOG_EPS = 1e-12
_PROTOCOL_RISK_HARD_THRESHOLD = 0.50


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class MarkerSelectionResult:
    """Result of a :class:`GeneSelector` run.

    Attributes
    ----------
    selected_genes:
        Ordered list of selected marker gene names.
    rejected_genes:
        Dict mapping each rejected gene → primary rejection reason.
        Only genes that failed a hard filter are included.
        Genes that simply ranked below the top-N are not included.
    score_components:
        DataFrame (selected genes × components) with the individual
        score components for each selected gene.
        Columns: ``specificity``, ``stability``, ``prot_weight``,
        ``concordance``, ``composite_score``.
    n_per_type:
        Dict mapping cell-type name → number of genes for which this type
        was the primary marker (i.e. the type with highest expression).
    gene_panel_overlap:
        Intersection of ``selected_genes`` with ``query_genes`` when
        ``query_genes`` was provided.  ``None`` otherwise.
    """

    selected_genes: list[str]
    rejected_genes: dict[str, str]
    score_components: pd.DataFrame
    n_per_type: dict[str, int]
    gene_panel_overlap: Optional[list[str]] = None

    @property
    def n_selected(self) -> int:
        return len(self.selected_genes)

    @property
    def n_rejected(self) -> int:
        return len(self.rejected_genes)

    def rejection_summary(self) -> pd.Series:
        """Return counts of rejections per reason."""
        if not self.rejected_genes:
            return pd.Series(dtype=int)
        return pd.Series(list(self.rejected_genes.values())).value_counts()


# ---------------------------------------------------------------------------
# GeneSelector
# ---------------------------------------------------------------------------


class GeneSelector:
    """Composite marker gene selector.

    Parameters
    ----------
    config:
        :class:`~tissueresolve.config.GeneConfig` with scoring weights and
        thresholds.
    """

    def __init__(self, config: Optional[GeneConfig] = None) -> None:
        self.cfg = config or GeneConfig()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def select(
        self,
        ref: ReferenceSignature,
        *,
        query_genes: Optional[list[str]] = None,
        protocol_risk_scores: Optional[pd.Series] = None,
        bulk_mean: Optional[pd.Series] = None,
        discordant_genes: Optional[set] = None,
        gene_filters: Optional["GeneFilterSet"] = None,  # type: ignore[name-defined]
        n_per_pair: int = 15,
    ) -> MarkerSelectionResult:
        """Select marker genes from *ref*.

        Parameters
        ----------
        ref:
            Reference built by
            :class:`~tissueresolve.reference.build.ReferenceBuilder`.
        query_genes:
            When provided, the final panel is restricted to genes present in
            this list (e.g. genes available in the Visium dataset or bulk
            counts).
        protocol_risk_scores:
            Per-gene protocol risk score ∈ [0, 1] (from
            ``protocol.risk.ProtocolRiskAssessor``).  Defaults to 0 for all
            genes (no penalty) when not provided.
        bulk_mean:
            Per-gene mean expression in the bulk dataset.  Used for the
            concordance component.
        discordant_genes:
            Genes to exclude due to paired bulk↔pseudobulk discordance.
        gene_filters:
            :class:`~tissueresolve.reference.gene_filters.GeneFilterSet`
            instance.  When ``None``, hard filters based on gene lists are
            skipped (protein-coding, hypervariable, blacklist).
        n_per_pair:
            Number of pairwise discriminability genes to add per cell-type pair
            (Phase 2).  Set to 0 to disable pairwise augmentation.

        Returns
        -------
        MarkerSelectionResult
        """
        phi = ref.as_phi()          # (G, K)
        genes = np.asarray(ref.gene_names)
        G, K = phi.shape

        # --- Phase 1: hard filters ----------------------------------------
        rejected: dict[str, str] = {}
        candidate_mask = np.ones(G, dtype=bool)

        # Pairwise max log₂FC filter
        pw_fc = _pairwise_max_fc(phi)
        cpm_max = phi.max(axis=1) * 1e6

        fc_mask = pw_fc >= self.cfg.min_log2fc
        expr_mask = cpm_max >= self.cfg.min_mean_cpm

        for i, g in enumerate(genes):
            if not fc_mask[i]:
                rejected[g] = "below_min_log2fc"
                candidate_mask[i] = False
            elif not expr_mask[i]:
                rejected[g] = "below_min_mean_cpm"
                candidate_mask[i] = False

        # Donor CV filter
        if ref.donor_cv is not None:
            mean_cv = ref.donor_cv.mean(axis=1)
            for i, g in enumerate(genes):
                if candidate_mask[i] and mean_cv[i] > self.cfg.max_donor_cv:
                    rejected[g] = "above_max_donor_cv"
                    candidate_mask[i] = False

        # Gene list filters (applied only when gene_filters provided)
        if gene_filters is not None:
            bl_prefixes = gene_filters.blacklist_prefixes()
            hv_genes = gene_filters.hypervariable_inflammatory()
            pc_genes = gene_filters.protein_coding()

            for i, g in enumerate(genes):
                if not candidate_mask[i]:
                    continue
                if any(g.startswith(p) for p in bl_prefixes):
                    rejected[g] = "blacklist_prefix"
                    candidate_mask[i] = False
                elif pc_genes and g not in pc_genes:
                    rejected[g] = "not_protein_coding"
                    candidate_mask[i] = False
                elif g in hv_genes:
                    rejected[g] = "hypervariable"
                    candidate_mask[i] = False

        # Protocol risk hard filter
        if protocol_risk_scores is not None:
            prot_risk = protocol_risk_scores.reindex(list(genes)).fillna(0.0)
            for i, g in enumerate(genes):
                if candidate_mask[i] and prot_risk.iloc[i] >= _PROTOCOL_RISK_HARD_THRESHOLD:
                    rejected[g] = "protocol_risk_high"
                    candidate_mask[i] = False

        # Discordance filter
        if discordant_genes:
            for i, g in enumerate(genes):
                if candidate_mask[i] and g in discordant_genes:
                    rejected[g] = "discordant"
                    candidate_mask[i] = False

        # Relaxation fallback: if too few candidates, loosen FC threshold
        candidate_idx = np.where(candidate_mask)[0]
        min_required = max(10, self.cfg.n_genes // 10)
        if len(candidate_idx) < min_required:
            relaxed = self.cfg.min_log2fc / 2.0
            logger.warning(
                "Only %d candidates pass filters; relaxing min_log2fc %.2f→%.2f.",
                len(candidate_idx), self.cfg.min_log2fc, relaxed,
            )
            fc_mask2 = pw_fc >= relaxed
            for i, g in enumerate(genes):
                if not candidate_mask[i] and rejected.get(g) == "below_min_log2fc" and fc_mask2[i]:
                    del rejected[g]
                    candidate_mask[i] = True
            candidate_idx = np.where(candidate_mask)[0]

        if len(candidate_idx) == 0:
            raise RuntimeError(
                "GeneSelector: gene panel is empty after all filters.  "
                "Lower GeneConfig.min_log2fc or min_mean_cpm, or check that "
                "the reference contains usable marker genes."
            )

        # --- Phase 1: composite scoring -----------------------------------
        prot_w = np.ones(len(candidate_idx))
        if protocol_risk_scores is not None:
            prs = protocol_risk_scores.reindex(genes[candidate_idx]).fillna(0.0)
            prot_w = np.clip(1.0 - prs.to_numpy(), 0.0, 1.0)

        stab_w = np.ones(len(candidate_idx))
        if ref.donor_cv is not None:
            cv_c = ref.donor_cv[candidate_idx].mean(axis=1)
            stab_w = 1.0 / (1.0 + cv_c)

        conc_w = np.ones(len(candidate_idx))
        if bulk_mean is not None:
            conc_w = _bulk_concordance(
                phi, candidate_idx, bulk_mean, genes
            )

        sep_raw = pw_fc[candidate_idx] * np.log1p(cpm_max[candidate_idx])
        sep_raw = np.clip(sep_raw, 0, None)

        cfg = self.cfg
        eps = _LOG_EPS
        log_score = (
            cfg.weight_specificity  * np.log(sep_raw  + eps)
            + cfg.weight_stability  * np.log(stab_w   + eps)
            + cfg.weight_protocol   * np.log(prot_w   + eps)
            + cfg.weight_concordance * np.log(conc_w  + eps)
        )
        log_score -= log_score.min()
        composite = np.exp(log_score)
        mu = composite.mean()
        if mu > 0:
            composite /= mu

        order = np.argsort(composite)[::-1]
        top_idx = candidate_idx[order[:self.cfg.n_genes]]
        phase1_genes = set(genes[top_idx].tolist())

        # Score component DataFrame
        score_df = pd.DataFrame({
            "specificity":   pw_fc[top_idx],
            "stability":     stab_w[order[:self.cfg.n_genes]],
            "prot_weight":   prot_w[order[:self.cfg.n_genes]],
            "concordance":   conc_w[order[:self.cfg.n_genes]],
            "composite_score": composite[order[:self.cfg.n_genes]],
        }, index=genes[top_idx])

        # Per-type dominant gene count
        dom_type = np.argmax(phi, axis=1)   # (G,) index of dominant cell type
        n_per_type = {
            ct: int((dom_type[top_idx] == k).sum())
            for k, ct in enumerate(ref.cell_types)
        }

        # --- Phase 2: pairwise discriminability augmentation --------------
        phase2_genes: set[str] = set()
        if n_per_pair > 0 and K >= 2:
            pairwise = _compute_pairwise_discriminability(
                ref, top_n_per_pair=n_per_pair
            )
            for gene_list in pairwise.values():
                phase2_genes.update(gene_list)
        # Hard-rejected genes must never enter the panel via Phase 2.
        phase2_genes -= set(rejected.keys())

        # n_genes caps the TOTAL panel (phase1 + phase2 combined).
        # Phase 1 genes take priority; phase 2 fills remaining budget.
        all_selected: set[str]
        if len(phase1_genes) >= self.cfg.n_genes:
            all_selected = phase1_genes
        else:
            extra_budget = self.cfg.n_genes - len(phase1_genes)
            extra_p2 = sorted(phase2_genes - phase1_genes)[:extra_budget]
            all_selected = phase1_genes | set(extra_p2)

        # --- Query filter -------------------------------------------------
        overlap: Optional[list[str]] = None
        if query_genes is not None:
            query_set = set(query_genes)
            not_in_query = [g for g in all_selected if g not in query_set]
            for g in not_in_query:
                if g not in rejected:
                    rejected[g] = "not_in_query"
            all_selected = all_selected & query_set
            overlap = sorted(all_selected)

        final = sorted(all_selected)

        logger.info(
            "Marker selection: %d genes selected "
            "(phase1=%d, phase2_added=%d, rejected=%d).",
            len(final), len(phase1_genes),
            len(phase2_genes - phase1_genes),
            len(rejected),
        )

        return MarkerSelectionResult(
            selected_genes=final,
            rejected_genes=rejected,
            score_components=score_df,
            n_per_type=n_per_type,
            gene_panel_overlap=overlap,
        )

    def compute_gene_weights(
        self,
        ref: ReferenceSignature,
        panel: list[str],
        *,
        protocol_risk_scores: Optional[pd.Series] = None,
        bulk_mean: Optional[pd.Series] = None,
        bulk_full: Optional[pd.DataFrame] = None,
    ) -> pd.Series:
        """Compute per-gene composite weights for the NNLS solver.

        Parameters
        ----------
        ref:
            Reference (full, not subsetted).
        panel:
            Gene panel (subset of ref.gene_names).
        protocol_risk_scores:
            Optional per-gene risk ∈ [0, 1].
        bulk_mean:
            Per-gene mean bulk expression.
        bulk_full:
            Full (genes × samples) bulk matrix; used for noise weighting.

        Returns
        -------
        pd.Series
            Log-additive composite weights indexed by gene name.
        """
        phi_df = pd.DataFrame(
            ref.as_phi(), index=ref.gene_names, columns=ref.cell_types
        ).loc[panel]
        phi = phi_df.to_numpy()
        genes_p = np.asarray(panel)

        pw_fc = _pairwise_max_fc(phi)
        cpm_max = phi.max(axis=1) * 1e6
        sep_raw = np.clip(pw_fc * np.log1p(cpm_max), 0, None)

        stab = np.ones(len(panel))
        if ref.donor_cv is not None:
            cv_df = pd.DataFrame(
                ref.donor_cv, index=ref.gene_names, columns=ref.cell_types
            ).reindex(panel)
            stab = 1.0 / (1.0 + cv_df.to_numpy().mean(axis=1))

        prot_w = np.ones(len(panel))
        if protocol_risk_scores is not None:
            prot_w = np.clip(
                1.0 - protocol_risk_scores.reindex(panel).fillna(0.0).to_numpy(),
                0.0, 1.0,
            )

        conc = np.ones(len(panel))
        if bulk_mean is not None:
            conc = _bulk_concordance(phi, np.arange(len(panel)), bulk_mean, genes_p)

        noise_w = np.ones(len(panel))
        if bulk_full is not None:
            sub = bulk_full.reindex(panel).fillna(0.0)
            mu_bulk = sub.mean(axis=1).to_numpy() + 1.0
            noise_w = np.exp(-0.5 * sub.std(axis=1).to_numpy() / mu_bulk)

        cfg = self.cfg
        eps = _LOG_EPS
        logw = (
            cfg.weight_specificity   * np.log(sep_raw   + eps)
            + cfg.weight_stability   * np.log(stab      + eps)
            + cfg.weight_protocol    * np.log(prot_w    + eps)
            + cfg.weight_concordance * np.log(conc      + eps)
            + 0.5                    * np.log(noise_w   + eps)
        )
        logw -= logw.min()
        w = np.exp(logw)
        mu = w.mean()
        return pd.Series(w / mu if mu > 0 else w, index=panel)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _pairwise_max_fc(phi: np.ndarray) -> np.ndarray:
    """Max absolute log₂FC across all cell-type pairs for each gene.

    Parameters
    ----------
    phi:
        (G, K) L1-normalised reference matrix.

    Returns
    -------
    np.ndarray shape (G,)
    """
    eps = 1e-12
    G, K = phi.shape
    if K < 2:
        return np.zeros(G)
    log_phi = np.log2(phi + eps)
    best = np.zeros(G)
    for j in range(K):
        for k in range(j + 1, K):
            diff = np.abs(log_phi[:, j] - log_phi[:, k])
            np.maximum(best, diff, out=best)
    return best


def _bulk_concordance(
    phi: np.ndarray,
    candidate_idx: np.ndarray,
    bulk_mean: pd.Series,
    gene_names: np.ndarray,
) -> np.ndarray:
    """Concordance between reference CPM and bulk mean expression."""
    panel_genes = gene_names[candidate_idx]
    bv = bulk_mean.reindex(panel_genes).fillna(0.0).to_numpy()
    rm = phi[candidate_idx].mean(axis=1) * 1e6
    with np.errstate(divide="ignore", invalid="ignore"):
        lb = np.log2(bv + 1.0)
        lr = np.log2(rm + 1.0)
    return np.exp(-np.abs(lb - lr) / 2.0)


def _compute_pairwise_discriminability(
    ref: ReferenceSignature,
    *,
    top_n_per_pair: int = 20,
) -> dict[tuple, list[str]]:
    """Top discriminating genes for each pair of cell types.

    Returns
    -------
    dict
        ``(type_a, type_b) → list[gene_name]``.
    """
    K = ref.n_cell_types
    ct = ref.cell_types
    R_log = ref.as_R_log().astype(np.float64)  # (K, G)
    gnames = np.asarray(ref.gene_names)
    result: dict[tuple, list[str]] = {}
    for i in range(K):
        for j in range(i + 1, K):
            lfc = np.abs(R_log[i] - R_log[j])
            top_idx = np.argsort(-lfc)[:top_n_per_pair]
            result[(ct[i], ct[j])] = gnames[top_idx].tolist()
    return result
