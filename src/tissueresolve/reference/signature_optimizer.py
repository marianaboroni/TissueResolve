"""Multi-panel donor-aware signature optimizer (Stage 1, experimental, opt-in).

Different signatures have different statistical jobs; one panel should not be assumed to
optimise all of them. This optimizer exposes **four explicit panels**, each with an
intended use and statistical contrast (see `docs/dev/STAGE1_MULTIPANEL_IMPLEMENTATION_NOTES.md`):

  * **broad**            — family-vs-family; estimates broad-family mass / reference suitability.
  * **fine_global**      — fine-type vs ALL other fine types; the **primary fine deconvolution**
                           panel. Wraps the validated `donor_de` (one-vs-rest). NOT sibling-only.
  * **sibling**          — fine-type vs its within-family siblings; supports **resolution
                           decisions / soft gating**, not necessarily the deconvolution matrix.
  * **rare_confirmation** — specificity-first per subtype; **presence/absence evidence and
                           false-positive control** (Stage 2), not abundance estimation.

It **composes** the validated donor-aware primitives in `reference.gene_selection` (no new DE).
It changes **no default**: it emits panels + an optional `ReferenceSignature.selected_genes`;
the production pipeline (GeneSelector + wNNLS) is untouched unless the user opts in. Panels
carry an `intended_use`; `panel_for(task)` returns the validated panel and misuse emits a warning.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from tissueresolve.reference import gene_selection as GS

__all__ = ["ReferenceSignatureOptimizer", "SignatureModel", "PanelSpec", "SignaturePanel",
           "validate_panel_use", "deconvolve_with_signature_panel", "warn_if_panel_misused",
           "PanelRoleError", "FEATURE_STATUS"]

FEATURE_STATUS = "experimental"
ALGORITHM_VERSION = "signature_optimizer-0.2.0-multipanel"

_MIN_STABLE_SIBLING_MARKERS = 3
_MIN_DONORS_PASS = 5
_MIN_CELLS_PER_TYPE = 20

#: validated statistical contrast + intended use per panel (misuse -> warning)
PANEL_SPECS = {
    "broad": dict(contrast="broad_family_vs_other_families",
                  intended_use="broad_family_deconvolution"),
    "fine_global": dict(contrast="fine_type_vs_all_other_fine_types",
                        intended_use="primary_fine_deconvolution"),
    "sibling": dict(contrast="fine_type_vs_within_family_siblings",
                    intended_use="resolution_decisions_and_gating"),
    "rare_confirmation": dict(contrast="subtype_specificity_first",
                              intended_use="rare_presence_absence_evidence"),
}


@dataclass
class PanelSpec:
    name: str
    contrast: str
    intended_use: str
    n_genes: int
    gene_budget: Optional[int]
    donor_aware: bool
    status: str
    validated_modality: str = "bulk_experimental"   # NOT spatial-validated yet (§3.7)
    spatial_status: str = "spatial_unvalidated"
    limitations: str = ""


@dataclass
class SignatureModel:
    broad_genes: list
    broad_minimal_genes: list
    fine_global_genes: list        # primary fine deconvolution panel
    sibling_genes: dict            # family -> list (resolution/gating only)
    sibling_genes_by_subtype: dict  # subtype -> list (for per-subtype evidence scores)
    rare_confirmation: dict        # subtype -> list (evidence only)
    status: str
    family_status: dict
    panels: dict                   # name -> PanelSpec
    gene_scores: pd.DataFrame
    selection_trace: pd.DataFrame
    metadata: dict = field(default_factory=dict)

    def as_panel(self, name: str) -> "SignaturePanel":
        """Typed SignaturePanel (genes + validated role/contrast/modality) for safeguards."""
        genes = {"broad": self.broad_genes, "fine_global": self.fine_global_genes,
                 "sibling": self.panel_for("resolution"),
                 "rare_confirmation": self.panel_for("rare_evidence")}[name]
        spec = self.panels[name]
        return SignaturePanel(genes=list(genes), role=spec.intended_use, contrast=spec.contrast,
                              modality="bulk", experimental_status=spec.status,
                              spatial_status=spec.spatial_status)

    def panel_for(self, task: str) -> list:
        """Return the validated panel for a task. Guards against panel misuse.

        task: 'broad_deconvolution' -> broad; 'fine_deconvolution' -> fine_global;
        'resolution' -> sibling (flattened); 'rare_evidence' -> rare_confirmation (flattened).
        """
        if task == "broad_deconvolution":
            return list(self.broad_genes)
        if task == "fine_deconvolution":
            return list(self.fine_global_genes)
        if task == "resolution":
            return list(dict.fromkeys(g for gl in self.sibling_genes.values() for g in gl))
        if task == "rare_evidence":
            return list(dict.fromkeys(g for gl in self.rare_confirmation.values() for g in gl))
        raise ValueError(f"unknown task {task!r}")

    def save(self, out_dir):
        out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"gene": self.broad_genes}).to_csv(out / "broad_signature.tsv", sep="\t", index=False)
        pd.DataFrame({"gene": self.fine_global_genes}).to_csv(
            out / "fine_global_signature.tsv", sep="\t", index=False)
        sib = [{"family": f, "gene": g} for f, gl in self.sibling_genes.items() for g in gl]
        pd.DataFrame(sib or [{"family": None, "gene": None}]).to_csv(
            out / "sibling_signature.tsv", sep="\t", index=False)
        rare = [{"subtype": s, "gene": g} for s, gl in self.rare_confirmation.items() for g in gl]
        pd.DataFrame(rare or [{"subtype": None, "gene": None}]).to_csv(
            out / "rare_confirmation_signature.tsv", sep="\t", index=False)
        self.gene_scores.to_csv(out / "signature_gene_scores.tsv", sep="\t", index=False)
        self.selection_trace.to_csv(out / "signature_selection_trace.tsv", sep="\t", index=False)
        (out / "signature_manifest.json").write_text(json.dumps({
            "status": self.status, "family_status": self.family_status,
            "panels": {k: vars(v) for k, v in self.panels.items()},
            "algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
            **self.metadata}, indent=2))


class PanelRoleError(ValueError):
    """Raised (strict mode) when a signature panel is used outside its validated role."""


@dataclass
class SignaturePanel:
    """A typed signature panel carrying its genes AND its validated role/contrast/modality,
    so downstream inference can refuse or flag misuse (§3.5 safeguards, integrated)."""
    genes: list
    role: str                      # e.g. "primary_fine_deconvolution"
    contrast: str
    modality: str = "bulk"
    experimental_status: str = "bulk_experimental"
    spatial_status: str = "spatial_unvalidated"


def validate_panel_use(panel: "SignaturePanel", expected_role: str, *, mode: str = "permissive",
                       modality: str = "bulk") -> list:
    """Validate a panel against the role/modality the caller expects.

    strict → raise PanelRoleError on incompatibility; permissive → return warnings (and emit
    them) so the caller can record them in metadata. Never silently accepts misuse.
    """
    issues = []
    if panel.role != expected_role:
        issues.append(f"panel role '{panel.role}' (contrast {panel.contrast}) used where "
                      f"'{expected_role}' is expected — not validated for this use")
    if modality == "spatial" and panel.spatial_status == "spatial_unvalidated":
        issues.append("panel is spatial_unvalidated but is being used for spatial inference")
    if issues:
        msg = "; ".join(issues)
        if mode == "strict":
            raise PanelRoleError(msg)
        warnings.warn(f"signature panel misuse: {msg}", stacklevel=2)
    return issues


def deconvolve_with_signature_panel(query, ref, panel: "SignaturePanel", *,
                                    expected_role: str = "primary_fine_deconvolution",
                                    mode: str = "permissive", modality: str = "bulk",
                                    solver: str = "poisson"):
    """Opt-in inference entry that ENFORCES panel-role safeguards before solving.

    Returns the solver's ``SolverResult`` with ``diagnostics['panel_use_warnings']`` recording
    any (permissive) role/modality issues. Backward-compatible: existing default paths are
    untouched; this is a new, explicit, safeguarded entry point.
    """
    from tissueresolve.solver import get_solver
    issues = validate_panel_use(panel, expected_role, mode=mode, modality=modality)
    res = get_solver(solver, genes=list(panel.genes)).solve(query, ref)
    res.diagnostics["panel_use_warnings"] = issues
    res.diagnostics["panel_role"] = panel.role
    return res


def warn_if_panel_misused(panel_name: str, task: str) -> None:
    """Warn if a panel is used outside its validated intended use."""
    ok = {"broad_deconvolution": "broad", "fine_deconvolution": "fine_global",
          "resolution": "sibling", "rare_evidence": "rare_confirmation"}
    expected = ok.get(task)
    if expected and panel_name != expected:
        warnings.warn(
            f"signature panel '{panel_name}' used for task '{task}' but the validated panel "
            f"is '{expected}'. Its contrast ({PANEL_SPECS.get(panel_name, {}).get('contrast')}) "
            f"is not optimised for this task; results may be unreliable.", stacklevel=2)


class ReferenceSignatureOptimizer:
    def __init__(self, celltype_col: str, donor_col: Optional[str],
                 mapping: Optional[dict] = None, broad_col: Optional[str] = None, *,
                 batch_col: Optional[str] = None, tolerance: float = 0.05, seed: int = 0,
                 broad_top_n: int = 20, fine_global_top_n: int = 15, sibling_top_n: int = 15,
                 rare_size: int = 8, gene_budget: Optional[int] = None,
                 min_cells: int = 10, minimal_sizes=(20, 40, 60, 100, 150, 200)):
        if mapping is None and broad_col is None:
            raise ValueError("provide either mapping (fine->broad) or broad_col")
        self.ctc = celltype_col
        self.dc = donor_col
        self.mapping = mapping
        self.broad_col = broad_col
        self.batch_col = batch_col
        self.tolerance = tolerance
        self.seed = seed
        self.broad_top_n = broad_top_n
        self.fine_global_top_n = fine_global_top_n
        self.sibling_top_n = sibling_top_n
        self.rare_size = rare_size
        self.gene_budget = gene_budget        # optional per-panel cap (fairness)
        self.min_cells = min_cells
        self.minimal_sizes = minimal_sizes

    def _broad_adata(self, adata):
        obs = adata.obs.copy()
        if self.broad_col and self.broad_col in obs.columns:
            broad = obs[self.broad_col].astype(str)
        else:
            broad = obs[self.ctc].astype(str).map(lambda c: str(self.mapping.get(c, c)))
        a = adata.copy(); a.obs["__broad__"] = pd.Categorical(broad.values)
        return a

    def _families(self, adata):
        if self.broad_col and self.broad_col in adata.obs.columns:
            df = adata.obs[[self.ctc, self.broad_col]].astype(str).drop_duplicates()
            fam = {}
            for _, r in df.iterrows():
                fam.setdefault(str(r[self.broad_col]), []).append(str(r[self.ctc]))
        else:
            fam = {}
            for c in pd.unique(adata.obs[self.ctc].astype(str)):
                fam.setdefault(str(self.mapping.get(c, c)), []).append(str(c))
        return fam

    def _cap(self, genes):
        """Flat cap — used ONLY for panels that are not per-type grouped (broad)."""
        return genes[: self.gene_budget] if self.gene_budget else genes

    def _round_robin(self, per_type: dict, budget: int, min_per_type: int):
        """Stratified allocation: guarantee `min_per_type` per eligible type, then
        round-robin the remaining budget. Avoids starving later-listed cell types
        (the flat-truncation bug). Returns (genes, shortfall_records)."""
        types = list(per_type)
        selected, seen = [], set()
        # phase 1: minimum coverage per type
        for t in types:
            for g in per_type[t][:min_per_type]:
                if g not in seen:
                    seen.add(g); selected.append(g)
        # phase 2: round-robin fill up to budget
        pos = min_per_type
        max_len = max((len(v) for v in per_type.values()), default=0)
        while len(selected) < budget and pos < max_len:
            for t in types:
                if len(selected) >= budget:
                    break
                if pos < len(per_type[t]):
                    g = per_type[t][pos]
                    if g not in seen:
                        seen.add(g); selected.append(g)
            pos += 1
        selected = selected[:budget]
        sel_set = set(selected)
        shortfall = [{"cell_type": t, "requested_genes": min_per_type,
                      "available_genes": len(per_type[t]),
                      "selected_genes": sum(g in sel_set for g in per_type[t]),
                      "budget_shortfall": max(0, min_per_type - len(per_type[t]))}
                     for t in types]
        return selected, shortfall

    def _minimal_broad(self, badata, ranked_genes):
        if self.dc is None:
            return ranked_genes[: self.minimal_sizes[0]], pd.DataFrame()
        X, labels, donors, genes = GS.donor_pseudobulk(
            badata, "__broad__", self.dc, min_cells=self.min_cells)
        gidx = {g: i for i, g in enumerate(genes)}
        ranked_idx = [gidx[g] for g in ranked_genes if g in gidx]
        uniq_don = sorted(set(donors))
        if len(uniq_don) < 3 or len(ranked_idx) < min(self.minimal_sizes):
            return ranked_genes[: self.minimal_sizes[0]], pd.DataFrame()
        rng = np.random.default_rng(self.seed)
        classes = sorted(set(labels))
        rows = []
        from scipy.optimize import nnls
        for k in [s for s in self.minimal_sizes if s <= len(ranked_idx)]:
            sel = ranked_idx[:k]
            fold_rmse = []
            for d in uniq_don:
                tr = donors != d
                cc = [c for c in classes
                      if (tr & (labels == c)).sum() > 0 and ((donors == d) & (labels == c)).sum() > 0]
                if len(cc) < 2:
                    continue
                R = np.vstack([X[np.ix_(tr & (labels == c), sel)].mean(0) for c in cc])
                held = np.vstack([X[np.ix_((donors == d) & (labels == c), sel)].mean(0) for c in cc])
                if not (np.isfinite(R).all() and np.isfinite(held).all()):
                    continue
                w = rng.dirichlet(np.ones(len(cc)), size=8)
                B = w @ held
                errs = []
                for s in range(B.shape[0]):
                    est = nnls(R.T, B[s])[0]
                    est = est / est.sum() if est.sum() > 0 else est
                    errs.append(np.sqrt(np.mean((est - w[s]) ** 2)))
                fold_rmse.append(float(np.mean(errs)))
            if fold_rmse:
                rows.append({"size": k, "donor_held_out_broad_rmse": float(np.mean(fold_rmse))})
        trace = pd.DataFrame(rows)
        if trace.empty:
            return ranked_genes[: self.minimal_sizes[0]], trace
        best = trace["donor_held_out_broad_rmse"].min()
        ok = trace[trace["donor_held_out_broad_rmse"] <= best * (1 + self.tolerance)]
        return ranked_genes[: int(ok["size"].min())], trace

    def _status(self, adata, families, fine_marker_counts, prov):
        n_types = adata.obs[self.ctc].astype(str).nunique()
        n_grp = adata.obs[self.dc].astype(str).nunique() if self.dc else 0
        fam_status = {}
        for fam, members in families.items():
            present = [m for m in members
                       if (adata.obs[self.ctc].astype(str) == m).sum() >= _MIN_CELLS_PER_TYPE]
            if len(present) < 2:
                fam_status[fam] = "BROAD_ONLY"; continue
            fam_status[fam] = ("PASS" if fine_marker_counts.get(fam, 0) >= _MIN_STABLE_SIBLING_MARKERS
                               else "BROAD_ONLY")
        # explicit fallback/degenerate states
        if n_types < 2:
            return "SINGLE_SUBTYPE", fam_status
        if prov.get("actual_method") == "none" or self.dc is None:
            return "REFERENCE_INADEQUATE", fam_status
        big = [t for t in pd.unique(adata.obs[self.ctc].astype(str))
               if (adata.obs[self.ctc].astype(str) == t).sum() >= _MIN_CELLS_PER_TYPE]
        if len(big) < 2:
            return "REFERENCE_INADEQUATE", fam_status
        # non-donor groupings (batch / pooled) never claim donor-validated PASS
        if not prov.get("donor_aware"):
            return "EXPERIMENTAL_FINE", fam_status
        passed = [f for f, s in fam_status.items() if s == "PASS"]
        if not passed:
            return "BROAD_ONLY", fam_status
        if len(passed) == len(fam_status) and n_grp >= _MIN_DONORS_PASS:
            return "PASS", fam_status
        return "PASS_WITH_RESTRICTIONS", fam_status

    def _resolve_grouping(self, adata):
        """Explicit fallback chain (§3): donor → batch → pooled-cell → none.
        Returns (grouping_col, provenance, adata) without mutating self permanently."""
        obs = adata.obs
        prov = {"requested_method": "donor_aware_pseudobulk_de", "donor_aware": False,
                "batch_aware": False, "confidence_level": "none", "actual_method": "none",
                "fallback_reason": "", "limitations": ""}
        if self.dc and self.dc in obs.columns and obs[self.dc].astype(str).nunique() >= 2:
            prov.update(actual_method="donor_aware_pseudobulk_de", donor_aware=True,
                        confidence_level="standard")
            return self.dc, prov, adata
        if self.batch_col and self.batch_col in obs.columns \
                and obs[self.batch_col].astype(str).nunique() >= 2:
            warnings.warn("no usable donor column — falling back to BATCH-aware pseudobulk DE; "
                          "donor stability is NOT validated.", stacklevel=3)
            prov.update(actual_method="batch_aware_pseudobulk_de", batch_aware=True,
                        confidence_level="reduced", fallback_reason="no donor column; batch grouping",
                        limitations="donor stability NOT validated (batch used as grouping)")
            return self.batch_col, prov, adata
        if obs.shape[0] >= self.min_cells * 2:
            warnings.warn("CRITICAL: no donor/batch column — POOLED-CELL DE (random cell bins as "
                          "pseudo-replicates); results are NOT donor-validated.", stacklevel=3)
            rng = np.random.default_rng(self.seed)
            a = adata.copy()
            a.obs["__pooled_bin__"] = pd.Categorical(
                [f"bin{int(x)}" for x in rng.integers(0, 3, size=a.n_obs)])
            prov.update(actual_method="pooled_cell_de", confidence_level="low",
                        fallback_reason="no donor/batch column; pooled-cell DE",
                        limitations="CRITICAL: cells treated as pseudo-replicates (random bins); "
                                    "NOT donor/batch validated")
            return "__pooled_bin__", prov, a
        prov.update(fallback_reason="insufficient cells/structure for any DE",
                    limitations="reference inadequate for donor-aware signatures")
        return None, prov, adata

    def optimize(self, adata) -> SignatureModel:
        _orig_dc = self.dc                      # never permanently mutate config (§3.1)
        try:
            grp, prov, adata = self._resolve_grouping(adata)
            self.dc = grp
            return self._optimize_inner(adata, prov)
        finally:
            self.dc = _orig_dc                  # restore so the object stays donor-aware

    def _optimize_inner(self, adata, prov) -> SignatureModel:
        badata = self._broad_adata(adata)
        families = self._families(adata)
        score_rows = []

        # A. broad panel (family-vs-family donor-aware one-vs-rest)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            broad_genes = []
            if self.dc:
                for fam in pd.unique(badata.obs["__broad__"].astype(str)):
                    de = GS.donor_aware_de(badata, "__broad__", self.dc, fam,
                                           top_n=self.broad_top_n, min_cells=self.min_cells)
                    for g in de.genes:
                        broad_genes.append(g)
                        r = de.stats.loc[g] if (not de.stats.empty and g in de.stats.index) else None
                        score_rows.append({"gene": g, "level": "broad", "target": fam,
                                           "log2fc": float(r["log2fc"]) if r is not None else np.nan,
                                           "support_frac": float(r["support_frac"]) if r is not None else np.nan})
            broad_genes = list(dict.fromkeys(broad_genes))

        gs = pd.DataFrame(score_rows)
        if not gs.empty:
            gs["rank"] = gs["support_frac"].fillna(0) * np.abs(gs["log2fc"].fillna(0))
            ranked = gs.sort_values("rank", ascending=False).drop_duplicates("gene")["gene"].tolist()
        else:
            ranked = broad_genes
        broad_minimal, trace = self._minimal_broad(badata, ranked)
        broad_genes = self._cap(broad_genes)

        # B. fine-global panel (fine-vs-ALL-other-fine; the validated donor_de) — primary.
        # Stratified per-type allocation when a budget is set (never flat-truncate a
        # type-grouped panel — §3.3), else the natural top_n_per_type union.
        fine_shortfall = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if not self.dc:
                fine_global = []
            elif self.gene_budget:
                fine_types = [t for t in pd.unique(adata.obs[self.ctc].astype(str))
                              if (adata.obs[self.ctc].astype(str) == t).sum() >= self.min_cells]
                per_type = {t: GS.donor_aware_de(adata, self.ctc, self.dc, t,
                                                 mode="one_vs_rest",
                                                 top_n=max(self.fine_global_top_n, 60),
                                                 min_cells=self.min_cells).genes
                            for t in fine_types}
                n_types = max(len(fine_types), 1)
                min_pt = max(2, min(self.fine_global_top_n, self.gene_budget // n_types))
                fine_global, fine_shortfall = self._round_robin(
                    per_type, self.gene_budget, min_pt)
            else:
                fine_global = GS.select_donor_aware_genes(
                    adata, self.ctc, self.dc, top_n_per_type=self.fine_global_top_n,
                    min_cells=self.min_cells)

        # C. sibling panel + D. rare-confirmation (per family)
        sibling_genes, sib_by_sub, rare_conf, fine_marker_counts = {}, {}, {}, {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for fam, members in families.items():
                members = [m for m in members
                           if (adata.obs[self.ctc].astype(str) == m).sum() >= self.min_cells]
                if len(members) < 2 or self.dc is None:
                    continue
                fam_genes, stable_ct = [], []
                for m in members:
                    de = GS.donor_aware_de(adata, self.ctc, self.dc, m, mode="sibling",
                                           siblings=members, top_n=self.sibling_top_n,
                                           min_cells=self.min_cells)
                    fam_genes.extend(de.genes)
                    sib_by_sub[m] = list(de.genes)          # per-subtype (for evidence)
                    if not de.stats.empty:
                        supported = de.stats[(de.stats["log2fc"] > 0) &
                                             (~de.stats["single_donor_penalised"])]
                        stable_ct.append(int(supported.head(self.sibling_top_n).shape[0]))
                        # rare-confirmation = sibling-supported AND globally specific
                        # (member of the fine-global marker set) — combined evidence, not
                        # sibling-only. Falls back to sibling-supported if no overlap.
                        sib_conf = (de.stats[(de.stats["supporting_donors"] >= 2) &
                                             (de.stats["log2fc"] > 0)]
                                    .sort_values(["support_frac", "log2fc"], ascending=False)
                                    .index.tolist())
                        fine_set = set(fine_global)
                        combined = [g for g in sib_conf if g in fine_set][: self.rare_size]
                        conf = combined if combined else sib_conf[: self.rare_size]
                        if conf:
                            rare_conf[m] = conf
                sibling_genes[fam] = list(dict.fromkeys(fam_genes))
                fine_marker_counts[fam] = int(np.median(stable_ct)) if stable_ct else 0

        status, fam_status = self._status(adata, families, fine_marker_counts, prov)

        def _spec(name, genes, donor_aware):
            return PanelSpec(name=name, contrast=PANEL_SPECS[name]["contrast"],
                             intended_use=PANEL_SPECS[name]["intended_use"],
                             n_genes=int(len(genes)), gene_budget=self.gene_budget,
                             donor_aware=bool(donor_aware), status="bulk_experimental")
        panels = {
            "broad": _spec("broad", broad_genes, self.dc is not None),
            "fine_global": _spec("fine_global", fine_global, self.dc is not None),
            "sibling": _spec("sibling", [g for gl in sibling_genes.values() for g in gl], self.dc is not None),
            "rare_confirmation": _spec("rare_confirmation",
                                       [g for gl in rare_conf.values() for g in gl], self.dc is not None),
        }
        meta = {"n_donors": int(adata.obs[self.dc].astype(str).nunique()) if self.dc else 0,
                "n_cell_types": int(adata.obs[self.ctc].astype(str).nunique()),
                "n_broad_families": len(families), "tolerance": self.tolerance,
                "seed": self.seed, "gene_budget": self.gene_budget,
                "budget_strategy": ("stratified_round_robin" if self.gene_budget else "natural_per_type"),
                "fine_global_shortfall": fine_shortfall,
                "selection_provenance": prov,
                "eligible_panels": [k for k, v in panels.items() if v.n_genes > 0],
                "has_donor_col": self.dc is not None}
        return SignatureModel(
            broad_genes=broad_genes, broad_minimal_genes=broad_minimal,
            fine_global_genes=fine_global, sibling_genes=sibling_genes,
            sibling_genes_by_subtype=sib_by_sub,
            rare_confirmation=rare_conf, status=status, family_status=fam_status,
            panels=panels, gene_scores=pd.DataFrame(score_rows),
            selection_trace=trace, metadata=meta)
