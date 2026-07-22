"""Separability-aware adaptive per-cell-type gene budget (Stage 1.6, experimental).

Estimates the gene count PER fine cell type from how hard that population is to separate —
easy types get compact panels, collinear types get more genes *only while they help*, and
structurally non-resolvable types stop expanding and are flagged. It does NOT fix a global or
per-type constant; those are only benchmark comparators.

Reuses existing pieces (no new distance metric, no default change):
  * `gene_selection.donor_aware_de` — donor-aware one-vs-rest candidate ranking.
  * `separability.compute_separability` — pairwise Bhattacharyya → closest confounder.
  * a leakage-safe **donor-held-out** per-type incremental-gain curve: a 2-component
    (target vs aggregate-rest) recovery on mixtures simulated from held-out TRAIN donors,
    evaluated over growing panels; stop at diminishing returns.

Bulk-only for now (validated_modality="bulk_experimental"; spatial unvalidated).
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from tissueresolve.reference import gene_selection as GS

__all__ = ["SeparabilityAwareBudgetAllocator", "AdaptiveSignature", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "separability_aware_allocator-0.1.0"


@dataclass
class AdaptiveSignature:
    genes: list                       # global deduped signature
    per_type: dict                    # type -> {genes, gene_count, status, best_rmse, ...}
    curves: pd.DataFrame              # long per-(type,size) curve
    manifest: dict = field(default_factory=dict)

    def save(self, out_dir):
        out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"gene": self.genes}).to_csv(out / "adaptive_signature.tsv", sep="\t", index=False)
        rows = [{"cell_type": k, **{kk: vv for kk, vv in v.items() if kk != "genes"},
                 "n_genes": v["gene_count"]} for k, v in self.per_type.items()]
        pd.DataFrame(rows).to_csv(out / "adaptive_gene_budget_by_celltype.tsv", sep="\t", index=False)
        self.curves.to_csv(out / "celltype_budget_curve.tsv", sep="\t", index=False)
        scores = [{"cell_type": k, "gene": g} for k, v in self.per_type.items() for g in v["genes"]]
        pd.DataFrame(scores).to_csv(out / "adaptive_signature_scores.tsv", sep="\t", index=False)
        (out / "adaptive_signature_manifest.json").write_text(json.dumps(self.manifest, indent=2))


class SeparabilityAwareBudgetAllocator:
    def __init__(self, celltype_col: str, donor_col: Optional[str], mapping: dict, *,
                 sizes: Sequence[int] = (5, 10, 15, 20, 30, 40, 60, 80, 100),
                 perf_tol: float = 0.03, gain_thr: float = 0.01,
                 min_genes: int = 5, max_genes: int = 100,
                 unresolvable_rmse: float = 0.30, easy_rmse: float = 0.10,
                 min_cells: int = 10, seed: int = 0):
        self.ctc = celltype_col; self.dc = donor_col; self.mapping = mapping
        self.sizes = sizes; self.perf_tol = perf_tol; self.gain_thr = gain_thr
        self.min_genes = min_genes; self.max_genes = max_genes
        self.unresolvable_rmse = unresolvable_rmse; self.easy_rmse = easy_rmse
        self.min_cells = min_cells; self.seed = seed

    # --- per-(donor,type) linear CPM profiles -----------------------------
    def _cpm_by_donor_type(self, adata):
        import scipy.sparse as sp
        obs = adata.obs
        ct = obs[self.ctc].astype(str).to_numpy(); dn = obs[self.dc].astype(str).to_numpy()
        M = adata.X
        genes = [str(g) for g in adata.var_names]
        prof = {}   # (donor, type) -> CPM vector
        for t in pd.unique(ct):
            for d in pd.unique(dn[ct == t]):
                idx = np.where((ct == t) & (dn == d))[0]
                if len(idx) < self.min_cells:
                    continue
                sub = M[idx]
                v = np.asarray(sub.sum(0)).ravel() if sp.issparse(sub) else np.asarray(sub).sum(0).ravel()
                tot = v.sum()
                prof[(d, t)] = (v / tot * 1e6) if tot > 0 else v
        return prof, genes

    # --- donor-held-out target-vs-CLOSEST-CONFOUNDER recovery curve ------
    def _type_curve(self, prof, gidx, k, confounder, target_ranked, conf_ranked):
        """Difficulty = donor-held-out recovery of the target's fraction against its CLOSEST
        confounder. The fit basis is the UNION of BOTH types' markers, so the two columns are
        non-proportional exactly when discriminating genes exist — a near-identical confounder
        (no discriminating genes) correctly reads as unresolvable. Panel grows with `size`."""
        from scipy.optimize import nnls
        rng = np.random.default_rng(self.seed)
        donors_k = sorted({d for (d, t) in prof if t == k})
        if confounder is None or not any(t == confounder for (_, t) in prof) or len(donors_k) < 2:
            return None
        donors_c = sorted({d for (d, t) in prof if t == confounder})
        wgrid = np.array([0.05, 0.1, 0.25, 0.5, 0.75])
        curve = []
        for size in self.sizes:
            panel = list(dict.fromkeys(list(target_ranked[:size]) + list(conf_ranked[:size])))
            sel = np.array([gidx[g] for g in panel if g in gidx])
            if sel.size < 2:
                continue
            fold_rmse = []
            for h in donors_k:                       # leave-one-donor-out (within train)
                tr_k = [prof[(d, k)][sel] for d in donors_k if d != h]
                tr_c = [prof[(d, confounder)][sel] for d in donors_c if d != h] \
                    or [prof[(d, confounder)][sel] for d in donors_c]
                if not tr_k or not tr_c or (h, k) not in prof:
                    continue
                pk_tr = np.mean(tr_k, axis=0); pc_tr = np.mean(tr_c, axis=0)
                pk_h = prof[(h, k)][sel]
                pc_h = prof[(h, confounder)][sel] if (h, confounder) in prof else pc_tr
                R = np.vstack([pk_tr, pc_tr]).T
                if not np.isfinite(R).all() or R.sum() == 0:
                    continue
                errs = []
                for w in wgrid:
                    mix = rng.poisson(np.clip(w * pk_h + (1 - w) * pc_h, 0, None)).astype(float)
                    est = nnls(R, mix)[0]
                    we = est[0] / est.sum() if est.sum() > 0 else 0.0
                    errs.append((we - w) ** 2)
                if errs:
                    fold_rmse.append(np.sqrt(np.mean(errs)))
            if fold_rmse:
                curve.append((size, float(np.mean(fold_rmse))))
        return curve

    def _choose(self, curve):
        """Smallest panel within perf_tol of best, respecting diminishing returns/status."""
        if not curve:
            return self.min_genes, float("nan"), "INSUFFICIENT_REFERENCE"
        sizes = [c[0] for c in curve]; rmse = [c[1] for c in curve]
        best = min(rmse)
        if best > self.unresolvable_rmse:
            return sizes[int(np.argmin(rmse))], best, "UNRESOLVABLE"
        # smallest size that is good enough: within perf_tol of best OR already below the
        # absolute easy floor (avoids chasing a near-zero best for an easily-separable type)
        thr = max(best * (1 + self.perf_tol), self.easy_rmse)
        chosen = next((s for s, r in curve if r <= thr), sizes[-1])
        chosen = int(np.clip(chosen, self.min_genes, self.max_genes))
        # status reflects achievable DIFFICULTY (best rmse) and the panel it needed
        if best <= self.easy_rmse and chosen <= self.sizes[len(self.sizes) // 2]:
            status = "EASY_COMPACT"
        elif best <= self.easy_rmse:
            status = "MODERATE"          # easy to resolve but needed a wider panel
        else:
            status = "DIFFICULT_EXPANDED"
        return chosen, best, status

    def allocate(self, adata) -> AdaptiveSignature:
        from tissueresolve.reference.separability import compute_separability
        obs = adata.obs
        types = [t for t in pd.unique(obs[self.ctc].astype(str))
                 if (obs[self.ctc].astype(str) == t).sum() >= self.min_cells]
        if self.dc is None or self.dc not in obs.columns or obs[self.dc].astype(str).nunique() < 2:
            warnings.warn("adaptive allocator: no usable donor column — cannot estimate "
                          "donor-held-out difficulty; returning INSUFFICIENT_REFERENCE.", stacklevel=2)
            return AdaptiveSignature([], {t: {"genes": [], "gene_count": 0,
                                              "status": "INSUFFICIENT_REFERENCE"} for t in types},
                                     pd.DataFrame(), {"algorithm_version": ALGORITHM_VERSION,
                                                      "feature_status": FEATURE_STATUS})
        if len(types) < 2:
            return AdaptiveSignature([], {types[0]: {"genes": [], "gene_count": 0,
                                          "status": "SINGLE_SUBTYPE"}} if types else {},
                                     pd.DataFrame(), {"algorithm_version": ALGORITHM_VERSION})
        prof, genes = self._cpm_by_donor_type(adata)
        gidx = {g: i for i, g in enumerate(genes)}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sep = compute_separability(self._ref_for_sep(adata, genes), warn_threshold=2.0)
        confounder = self._closest_confounders(sep)
        # per-type one-vs-rest ranked candidates, computed ONCE
        ranked_by_type = {}
        for k in types:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                de = GS.donor_aware_de(adata, self.ctc, self.dc, k, mode="one_vs_rest",
                                       top_n=self.max_genes, min_cells=self.min_cells)
            ranked_by_type[k] = [g for g in de.genes if g in gidx]
        per_type, curve_rows, all_genes = {}, [], []
        for k in types:
            ranked = ranked_by_type[k]
            c = confounder.get(k, {}).get("other")
            conf_ranked = ranked_by_type.get(c, []) if c else []
            curve = self._type_curve(prof, gidx, k, c, ranked, conf_ranked)
            n, best, status = self._choose(curve or [])
            n = min(n, len(ranked))
            sel_genes = ranked[:n]
            per_type[k] = {"genes": sel_genes, "gene_count": len(sel_genes),
                           "status": status, "best_rmse": round(best, 4) if best == best else None,
                           "separability_bc": round(confounder.get(k, {}).get("bc", float("nan")), 4),
                           "closest_confounder": confounder.get(k, {}).get("other", None),
                           "n_donors": int(len({d for (d, t) in prof if t == k})),
                           "n_candidates": len(ranked)}
            all_genes.extend(sel_genes)
            for s, r in (curve or []):
                curve_rows.append({"cell_type": k, "size": s, "held_out_rmse": r})
        signature = list(dict.fromkeys(all_genes))
        manifest = {"algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
                    "validated_modality": "bulk_experimental", "spatial_status": "unvalidated",
                    "n_types": len(types), "n_signature_genes": len(signature),
                    "perf_tol": self.perf_tol, "gain_thr": self.gain_thr,
                    "min_genes": self.min_genes, "max_genes": self.max_genes,
                    "unresolvable_rmse": self.unresolvable_rmse, "seed": self.seed,
                    "status_counts": pd.Series([v["status"] for v in per_type.values()]).value_counts().to_dict()}
        return AdaptiveSignature(signature, per_type, pd.DataFrame(curve_rows), manifest)

    def _ref_for_sep(self, adata, genes):
        """Build a minimal mean-profile ReferenceSignature for compute_separability."""
        import scipy.sparse as sp
        from tissueresolve.results import ReferenceSignature
        obs = adata.obs; ct = obs[self.ctc].astype(str).to_numpy(); M = adata.X
        cts = [t for t in pd.unique(ct) if (ct == t).sum() >= self.min_cells]
        R = []
        for t in cts:
            idx = np.where(ct == t)[0]; sub = M[idx]
            v = np.asarray(sub.sum(0)).ravel() if sp.issparse(sub) else np.asarray(sub).sum(0).ravel()
            R.append(v / v.sum() * 1e6 if v.sum() > 0 else v)
        return ReferenceSignature(gene_names=genes, cell_types=list(cts),
                                  R_cpm=np.vstack(R).astype("float32"))

    @staticmethod
    def _closest_confounders(sep):
        out = {}
        for p in sep.pairs:
            for a, b in ((p.type_a, p.type_b), (p.type_b, p.type_a)):
                if a not in out or p.bhattacharyya_coeff > out[a]["bc"]:
                    out[a] = {"other": b, "bc": float(p.bhattacharyya_coeff)}
        return out
