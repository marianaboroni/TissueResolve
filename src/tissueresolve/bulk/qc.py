"""
Per-sample and per-cell-type quality control for bulk deconvolution.

All thresholds are heuristic.  Every threshold is written to output JSON with
``"is_heuristic": true``.  Provenance is documented in ``docs/qc_thresholds.md``.

Metrics
-------
Per sample:
    recon_r2      — Reconstruction R².  From ``BulkDeconvResult.coverage_r2``.
    recon_rmse    — Per-sample RMSE between L1-normalised bulk and reconstruction.
    profile_corr  — Pearson(L1-norm(bulk_n)[panel], mean_column(phi[panel])).
                    Scale-independent; comparable across protocols.
    mismatch_flag — "low" | "medium" | "high".  Heuristic combination of r2
                    and profile_corr.

Per cell type:
    marker_recall  — Fraction of top-50 specificity markers that are expressed
                     ≥ 0.5 CPM-equivalent in the bulk dataset.
    spillover_risk — Max Pearson(φ_k, φ_j) across j ≠ k on the panel.
                     High values indicate collinear types; estimates unreliable.
    mean_ci_width  — Mean bootstrap CI width.  None if bootstrap was not run.

Global:
    condition_number          — κ(Φ_panel).  Heuristic collinearity indicator.
    n_genes_excluded_protocol — From ProtocolRiskReport.
    protocol_risk_level       — From ProtocolRiskReport.
    recommendations           — Non-suppressible human-readable strings.
"""
from __future__ import annotations

import logging
from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from tissueresolve.config import BulkQCConfig
from tissueresolve.results import BulkDeconvResult, QCReport, ReferenceSignature

__all__ = ["BulkQC"]

logger = logging.getLogger("tissueresolve.bulk.qc")

_TOP_MARKERS = 50
_MIN_BULK_CPM = 0.5


class BulkQC:
    """Compute QC metrics for a completed bulk deconvolution run.

    Parameters
    ----------
    config:
        :class:`~tissueresolve.config.BulkQCConfig` with heuristic thresholds.
    """

    def __init__(self, config: Optional[BulkQCConfig] = None) -> None:
        self.cfg = config or BulkQCConfig()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def compute(
        self,
        result: BulkDeconvResult,
        ref: ReferenceSignature,
        bulk: pd.DataFrame,
        protocol_risk_report: Optional[Any] = None,
    ) -> QCReport:
        """Compute all QC metrics and return a :class:`~tissueresolve.results.QCReport`.

        Parameters
        ----------
        result:
            Completed deconvolution result.
        ref:
            Reference used in the deconvolution.
        bulk:
            Original genes × samples bulk count matrix.
        protocol_risk_report:
            Optional :class:`~tissueresolve.protocol.risk.ProtocolRiskReport`.
            If provided, ``n_genes_excluded`` and ``risk_level`` are used in QC.

        Returns
        -------
        QCReport
        """
        panel = result.gene_panel
        phi_mat = ref.as_phi()
        phi_df = pd.DataFrame(phi_mat, index=ref.gene_names, columns=ref.cell_types)
        phi_panel = phi_df.reindex(panel).dropna(how="all")
        panel = list(phi_panel.index)

        recon_r2 = result.coverage_r2.copy()
        recon_rmse = self._recon_rmse(bulk, phi_panel, result.proportions, panel)
        profile_corr = self._profile_correlation(bulk, phi_panel, panel)

        mismatch_flag = pd.Series(
            [
                self._mismatch_flag(r2, pc)
                for r2, pc in zip(
                    recon_r2.reindex(profile_corr.index),
                    profile_corr,
                )
            ],
            index=profile_corr.index,
            name="mismatch_flag",
        )

        marker_recall = self._marker_recall(ref, bulk, panel)
        spillover_risk = self._spillover_risk(phi_panel)

        mean_ci_width: Optional[pd.Series] = None
        if result.lower_ci is not None and result.upper_ci is not None:
            mean_ci_width = (result.upper_ci - result.lower_ci).mean(axis=0)
            mean_ci_width.name = "mean_ci_width"

        phi_arr = phi_panel.to_numpy()
        try:
            cond_num = float(np.linalg.cond(phi_arr))
        except np.linalg.LinAlgError:
            cond_num = float("nan")

        n_excl = 0
        risk_lv = "unknown"
        if protocol_risk_report is not None:
            n_excl = protocol_risk_report.n_genes_excluded
            risk_lv = protocol_risk_report.risk_level

        recommendations = self._generate_recommendations(
            recon_r2, profile_corr, mismatch_flag,
            marker_recall, spillover_risk, cond_num,
            risk_lv, n_excl,
        )

        logger.info(
            "BulkQC: mean_R²=%.3f  mean_profile_corr=%.3f  κ(Φ)=%.0f  "
            "high_mismatch=%d  mean_RMSE=%.4f",
            float(recon_r2.mean()),
            float(profile_corr.dropna().mean()) if not profile_corr.isna().all() else float("nan"),
            cond_num,
            int((mismatch_flag == "high").sum()),
            float(recon_rmse.mean()),
        )

        meta: dict[str, Any] = {
            "is_heuristic": True,
            "recon_rmse_mean": round(float(recon_rmse.mean()), 6),
            "condition_number": cond_num,
        }
        return QCReport(
            modality="bulk",
            recommendations=recommendations,
            metadata=meta,
            recon_r2=recon_r2,
            profile_corr=profile_corr,
            mismatch_flag=mismatch_flag,
            marker_recall=marker_recall,
            spillover_risk=spillover_risk,
            mean_ci_width=mean_ci_width,
            condition_number=cond_num,
            protocol_risk_level=risk_lv,
            n_genes_excluded_protocol=n_excl,
        )

    # ------------------------------------------------------------------
    # Per-sample helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _recon_rmse(
        bulk: pd.DataFrame,
        phi_panel: pd.DataFrame,
        proportions: pd.DataFrame,
        panel: list[str],
    ) -> pd.Series:
        """Per-sample RMSE between L1-normalised bulk and reconstruction."""
        shared = [g for g in panel if g in bulk.index]
        if not shared:
            return pd.Series(np.nan, index=bulk.columns, name="recon_rmse")

        B = bulk.loc[shared].to_numpy(dtype=np.float64)
        col_sums = B.sum(axis=0, keepdims=True)
        col_sums = np.where(col_sums == 0, 1.0, col_sums)
        B_norm = B / col_sums

        R = phi_panel.loc[shared].to_numpy(dtype=np.float64)
        col_sums_r = R.sum(axis=0, keepdims=True)
        col_sums_r = np.where(col_sums_r == 0, 1.0, col_sums_r)
        R_norm = R / col_sums_r

        W = proportions.reindex(bulk.columns).to_numpy(dtype=np.float64)  # (N, K)
        recon = R_norm @ W.T   # (G, N)

        rmse_vals = np.sqrt(np.mean((B_norm - recon) ** 2, axis=0))
        return pd.Series(rmse_vals, index=list(bulk.columns), name="recon_rmse")

    @staticmethod
    def _profile_correlation(
        bulk: pd.DataFrame,
        phi_panel: pd.DataFrame,
        panel: list[str],
    ) -> pd.Series:
        """Pearson(L1-norm(bulk_n)[panel], mean_reference_profile[panel])."""
        shared = [g for g in panel if g in bulk.index]
        if not shared:
            return pd.Series(np.nan, index=bulk.columns, name="profile_corr")

        phi_mean = phi_panel.loc[shared].mean(axis=1).values
        results: dict[str, float] = {}
        for sid in bulk.columns:
            b_raw = bulk.loc[shared, sid].to_numpy(dtype=float)
            b_sum = b_raw.sum()
            b = b_raw / b_sum if b_sum > 0 else b_raw
            if np.std(b) < 1e-9 or np.std(phi_mean) < 1e-9:
                results[sid] = float("nan")
            else:
                r, _ = pearsonr(b, phi_mean)
                results[sid] = float(r)
        return pd.Series(results, name="profile_corr")

    def _mismatch_flag(self, r2: float, profile_corr: float) -> str:
        """Heuristic flag: 'high' | 'medium' | 'low'.

        Applied in order of severity.  Thresholds from ``BulkQCConfig``.
        All thresholds are heuristic and are written to output JSON.
        """
        r2_v = float(r2) if np.isfinite(r2) else 0.0
        pc_v = float(profile_corr) if np.isfinite(profile_corr) else 0.0
        if r2_v < self.cfg.r2_fail or pc_v < self.cfg.profile_corr_fail:
            return "high"
        if r2_v < self.cfg.r2_warn or pc_v < self.cfg.profile_corr_warn:
            return "medium"
        return "low"

    # ------------------------------------------------------------------
    # Per-cell-type helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _marker_recall(
        ref: ReferenceSignature,
        bulk: pd.DataFrame,
        panel: list[str],
    ) -> pd.Series:
        """Fraction of top-50 specificity markers detected in bulk (≥ 0.5 CPM)."""
        phi = ref.as_phi()  # (G, K)
        gene_arr = np.asarray(ref.gene_names)
        eps = 1e-12
        log_phi = np.log2(phi + eps)

        bulk_in_panel = bulk.reindex(list(gene_arr)).fillna(0.0)
        total_bulk = bulk_in_panel.sum().sum()
        bulk_mean_cpm = (
            bulk_in_panel.sum(axis=1) / max(total_bulk, 1) * 1e6
        )

        results: dict[str, float] = {}
        for k, ct in enumerate(ref.cell_types):
            fc_vs_mean = log_phi[:, k] - np.log2(phi.mean(axis=1) + eps)
            top_idx = np.argsort(fc_vs_mean)[::-1][:_TOP_MARKERS]
            top_genes = gene_arr[top_idx]
            detected = sum(
                1 for g in top_genes
                if g in bulk_mean_cpm.index and bulk_mean_cpm[g] >= _MIN_BULK_CPM
            )
            results[ct] = detected / max(_TOP_MARKERS, 1)
        return pd.Series(results, name="marker_recall")

    @staticmethod
    def _spillover_risk(phi_panel: pd.DataFrame) -> pd.Series:
        """Max Pearson correlation between each cell type and its nearest neighbour."""
        phi = phi_panel.to_numpy()   # (G, K)
        labels = list(phi_panel.columns)
        K = phi.shape[1]
        results: dict[str, float] = {}
        for k in range(K):
            max_corr = -1.0
            for j in range(K):
                if j == k:
                    continue
                v1, v2 = phi[:, k], phi[:, j]
                if np.std(v1) < 1e-9 or np.std(v2) < 1e-9:
                    continue
                r, _ = pearsonr(v1, v2)
                max_corr = max(max_corr, float(r))
            results[labels[k]] = max_corr
        return pd.Series(results, name="spillover_risk")

    # ------------------------------------------------------------------
    # Recommendations
    # ------------------------------------------------------------------

    def _generate_recommendations(
        self,
        recon_r2: pd.Series,
        profile_corr: pd.Series,
        mismatch_flag: pd.Series,
        marker_recall: pd.Series,
        spillover_risk: pd.Series,
        cond_num: float,
        protocol_risk: str,
        n_excl: int,
    ) -> list[str]:
        recs: list[str] = []

        n_high = int((mismatch_flag == "high").sum())
        n_med = int((mismatch_flag == "medium").sum())
        if n_high > 0:
            recs.append(
                f"{n_high} sample(s) have mismatch_flag=high "
                f"(recon_R² < {self.cfg.r2_fail} or profile_corr "
                f"< {self.cfg.profile_corr_fail}).  "
                "Consider excluding these or using a more appropriate reference."
            )
        if n_med > 0:
            recs.append(
                f"{n_med} sample(s) have mismatch_flag=medium.  "
                "Use caution in downstream analysis."
            )

        low_recall = marker_recall[marker_recall < self.cfg.marker_recall_fail]
        for ct, val in low_recall.items():
            recs.append(
                f"Low marker recall for '{ct}' ({val:.2f}).  "
                "This cell type may be absent from bulk samples or its markers "
                "may be incompatible with the bulk protocol."
            )

        high_spill = spillover_risk[spillover_risk > 0.85]
        for ct, val in high_spill.items():
            recs.append(
                f"High spillover risk for '{ct}' (r={val:.2f}).  "
                "Collinear cell types: subtype-level estimates may be unreliable."
            )

        if np.isfinite(cond_num):
            if cond_num > self.cfg.condition_number_fail:
                recs.append(
                    f"κ(Φ_panel) = {cond_num:.0f} "
                    f"(threshold: {self.cfg.condition_number_fail:.0f}).  "
                    "Estimates may be unstable for linear methods.  "
                    "Consider coarser cell-type annotations."
                )
            elif cond_num > self.cfg.condition_number_warn:
                recs.append(
                    f"κ(Φ_panel) = {cond_num:.0f}.  Moderate collinearity detected."
                )

        if protocol_risk in ("medium", "high"):
            recs.append(
                f"Protocol risk level: {protocol_risk}.  "
                f"{n_excl} genes excluded by protocol filter.  "
                "Results may be affected by residual protocol mismatch."
            )

        return recs
