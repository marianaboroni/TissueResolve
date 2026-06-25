"""
Within-family marker / HVG selection for fine-subtype resolution.

The global signature (``reference/markers.py``) is good at separating *broad*
families but not *fine subtypes within* a family.  This module recomputes
highly-variable and discriminative genes **inside each broad family**, using the
raw per-cell AnnData (only available at reference-build time), and produces
family-specific gene panels + per-gene scores that the hierarchical solver can
use for the fine-level step.

It is **additive**: the existing global selection is untouched.  Nothing here
overclaims fine-subtype accuracy — a family is only useful for sub-resolution if
its within-family panel actually discriminates its subtypes (reported, not
assumed).

Design note
-----------
HVG and differential expression need per-cell data, so these functions take an
AnnData (``adata.X`` = raw counts, ``adata.var_names`` = genes,
``adata.obs[broad_col]`` / ``[fine_col]`` = labels).  The resulting panels are
plain gene-name lists that can be stored on / threaded through the aggregated
``ReferenceSignature`` used at solve time.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Mapping, Optional

import numpy as np
import pandas as pd

__all__ = [
    "select_within_family_hvgs",
    "select_within_family_de_genes",
    "select_pairwise_discriminative_genes_within_family",
    "build_family_specific_gene_panels",
    "within_family_resolution_summary",
    "FamilyPanels",
]

_EPS = 1e-9


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _dense_counts(adata) -> np.ndarray:
    X = adata.X
    try:
        import scipy.sparse as sp
        if sp.issparse(X):
            X = X.toarray()
    except Exception:
        pass
    return np.asarray(X, dtype=np.float64)


def _normalize_log(counts: np.ndarray, target: float = 1e4) -> np.ndarray:
    """Library-size normalise to *target* then log1p (cells × genes)."""
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(counts / lib * target)


def _is_mito(g: str) -> bool:
    s = str(g).upper()
    return s.startswith("MT-") or s.startswith("MT.") or s.startswith("MT_")


def _is_ribo(g: str) -> bool:
    s = str(g).upper()
    return (s.startswith("RPS") or s.startswith("RPL")
            or s.startswith("MRPS") or s.startswith("MRPL"))


def _stress_genes(genome: str = "hg38") -> frozenset:
    try:
        from tissueresolve.reference.gene_filters import GeneFilterSet
        gf = GeneFilterSet(genome=genome)
        return frozenset({str(g).upper() for g in
                          (gf.dissociation_stress() | gf.hypervariable_inflammatory())})
    except Exception:
        # heuristic fallback (heat-shock / immediate-early / IFN response)
        return frozenset({"FOS", "JUN", "JUNB", "EGR1", "HSPA1A", "HSPA1B",
                          "HSP90AA1", "DNAJB1", "DUSP1", "IER2", "FOSB"})


def _family_members(adata, broad_col: str, fine_col: str,
                    family: str, min_cells: int) -> list[str]:
    obs = adata.obs
    mask = obs[broad_col].astype(str) == str(family)
    sub = obs.loc[mask, fine_col].astype(str)
    vc = sub.value_counts()
    return sorted(str(s) for s, n in vc.items() if n >= min_cells)


def _subtype_stats(Xln: np.ndarray, counts: np.ndarray, labels: np.ndarray,
                   subtypes: list[str]) -> tuple[dict, dict]:
    """Per-subtype mean (log-normalised) and detection fraction."""
    means, detect = {}, {}
    for s in subtypes:
        rows = labels == s
        if rows.sum() == 0:
            continue
        means[s] = Xln[rows].mean(axis=0)
        detect[s] = (counts[rows] > 0).mean(axis=0)
    return means, detect


# --------------------------------------------------------------------------- #
# 1. within-family HVGs
# --------------------------------------------------------------------------- #
def select_within_family_hvgs(
    adata_ref,
    broad_col: str,
    fine_col: str,
    family: str,
    *,
    sample_col: Optional[str] = None,
    n_hvgs: int = 500,
    min_cells_per_subtype: int = 20,
    min_detection: float = 0.05,
    max_genes: Optional[int] = None,
) -> pd.DataFrame:
    """Highly-variable genes computed from cells of one broad *family* only.

    Returns a DataFrame indexed by gene with columns ``mean``, ``dispersion``,
    ``detection``, ``hvg_score`` (0–1), sorted by ``hvg_score`` descending and
    truncated to ``n_hvgs`` (or ``max_genes``).  Raises ``ValueError`` if the
    family has fewer than two qualifying subtypes.
    """
    genes = np.asarray([str(g) for g in adata_ref.var_names])
    members = _family_members(adata_ref, broad_col, fine_col, family,
                              min_cells_per_subtype)
    if len(members) < 2:
        raise ValueError(f"family '{family}' has <2 subtypes with "
                         f">={min_cells_per_subtype} cells; cannot sub-resolve.")
    mask = (adata_ref.obs[broad_col].astype(str) == str(family)).to_numpy()
    counts = _dense_counts(adata_ref)[mask]
    Xln = _normalize_log(counts)
    detection = (counts > 0).mean(axis=0)
    keep = detection >= min_detection
    mean = Xln.mean(axis=0)
    var = Xln.var(axis=0)
    dispersion = var / (mean + _EPS)
    # Highly variable = high within-family (log-normalised) variance: genes that
    # vary across cells of the family.  Genes flat within the family (e.g. broad
    # between-family markers) score low, which is exactly what we want for the
    # fine-level panel.
    score_metric = np.where(keep, var, -np.inf)
    order = np.argsort(score_metric)[::-1]
    n_keep = int(keep.sum())
    score = np.zeros(len(genes))
    if n_keep > 0:
        ranks = np.empty(len(genes))
        ranks[order] = np.arange(len(genes))
        score = np.clip(1.0 - ranks / max(n_keep, 1), 0.0, 1.0)
        score[~keep] = 0.0
    df = pd.DataFrame({"gene": genes, "mean": mean, "dispersion": dispersion,
                       "detection": detection, "hvg_score": score})
    df = df[keep].sort_values("hvg_score", ascending=False)
    top = max_genes if max_genes else n_hvgs
    return df.head(int(top)).set_index("gene")


# --------------------------------------------------------------------------- #
# 2. within-family one-vs-rest DE genes
# --------------------------------------------------------------------------- #
def select_within_family_de_genes(
    adata_ref,
    broad_col: str,
    fine_col: str,
    family: str,
    *,
    sample_col: Optional[str] = None,
    top_n_per_pair: int = 50,
    min_logfc: float = 0.25,
    min_detection: float = 0.05,
    min_cells_per_subtype: int = 20,
) -> pd.DataFrame:
    """One-vs-rest marker genes for each fine subtype within a *family*.

    Returns rows ``gene, subtype, log2fc, detection, de_score`` for genes
    up-regulated in a subtype vs the rest of the family
    (``log2fc >= min_logfc``, detection in the subtype ``>= min_detection``),
    keeping the top ``top_n_per_pair`` per subtype.
    """
    genes = np.asarray([str(g) for g in adata_ref.var_names])
    members = _family_members(adata_ref, broad_col, fine_col, family,
                              min_cells_per_subtype)
    if len(members) < 2:
        raise ValueError(f"family '{family}' has <2 qualifying subtypes.")
    mask = (adata_ref.obs[broad_col].astype(str) == str(family)).to_numpy()
    counts = _dense_counts(adata_ref)[mask]
    Xln = _normalize_log(counts) / np.log(2)   # log2 scale
    labels = adata_ref.obs.loc[mask, fine_col].astype(str).to_numpy()

    rows = []
    for s in members:
        in_s = labels == s
        out_s = np.isin(labels, [m for m in members if m != s])
        if in_s.sum() == 0 or out_s.sum() == 0:
            continue
        mean_in = Xln[in_s].mean(axis=0)
        mean_out = Xln[out_s].mean(axis=0)
        det_in = (counts[in_s] > 0).mean(axis=0)
        lfc = mean_in - mean_out
        ok = (lfc >= min_logfc) & (det_in >= min_detection)
        idx = np.where(ok)[0]
        idx = idx[np.argsort(lfc[idx])[::-1][:top_n_per_pair]]
        for g in idx:
            rows.append({"gene": genes[g], "subtype": s,
                         "log2fc": float(lfc[g]), "detection": float(det_in[g]),
                         "de_score": float(min(lfc[g] / 2.0, 1.0))})
    return pd.DataFrame(rows, columns=["gene", "subtype", "log2fc",
                                       "detection", "de_score"])


# --------------------------------------------------------------------------- #
# 3. within-family pairwise discriminative genes
# --------------------------------------------------------------------------- #
def select_pairwise_discriminative_genes_within_family(
    adata_ref,
    broad_col: str,
    fine_col: str,
    family: str,
    *,
    sample_col: Optional[str] = None,
    top_n_per_pair: int = 50,
    min_logfc: float = 0.25,
    min_detection: float = 0.05,
    min_cells_per_subtype: int = 20,
    stability_across_samples: bool = True,
) -> pd.DataFrame:
    """Genes that discriminate each *pair* of fine subtypes within a family.

    Returns ``gene, type_a, type_b, log2fc, stability, pair_score``.  When
    ``stability_across_samples`` and ``sample_col`` are given, ``stability`` is
    the fraction of samples in which the gene keeps the same up/down direction
    (a sample-robustness check); otherwise it is 1.0.
    """
    genes = np.asarray([str(g) for g in adata_ref.var_names])
    members = _family_members(adata_ref, broad_col, fine_col, family,
                              min_cells_per_subtype)
    if len(members) < 2:
        raise ValueError(f"family '{family}' has <2 qualifying subtypes.")
    mask = (adata_ref.obs[broad_col].astype(str) == str(family)).to_numpy()
    counts = _dense_counts(adata_ref)[mask]
    Xln = _normalize_log(counts) / np.log(2)
    labels = adata_ref.obs.loc[mask, fine_col].astype(str).to_numpy()
    samples = (adata_ref.obs.loc[mask, sample_col].astype(str).to_numpy()
               if (sample_col and sample_col in adata_ref.obs.columns) else None)

    def _mean(rows):
        return Xln[rows].mean(axis=0)

    rows = []
    for a, b in combinations(members, 2):
        ra, rb = labels == a, labels == b
        if ra.sum() == 0 or rb.sum() == 0:
            continue
        lfc = _mean(ra) - _mean(rb)
        det = np.maximum((counts[ra] > 0).mean(axis=0),
                         (counts[rb] > 0).mean(axis=0))
        ok = (np.abs(lfc) >= min_logfc) & (det >= min_detection)
        idx = np.where(ok)[0]
        idx = idx[np.argsort(np.abs(lfc[idx]))[::-1][:top_n_per_pair]]
        # sample stability: same sign of lfc across samples
        stab = np.ones(len(genes))
        if stability_across_samples and samples is not None:
            uniq = [s for s in np.unique(samples)]
            if len(uniq) >= 2:
                signs = []
                for sm in uniq:
                    rasm = ra & (samples == sm)
                    rbsm = rb & (samples == sm)
                    if rasm.sum() and rbsm.sum():
                        signs.append(np.sign(_mean(rasm) - _mean(rbsm)))
                if signs:
                    S = np.vstack(signs)
                    ref_sign = np.sign(lfc)
                    stab = (S == ref_sign).mean(axis=0)
        for g in idx:
            rows.append({"gene": genes[g], "type_a": a, "type_b": b,
                         "log2fc": float(lfc[g]), "stability": float(stab[g]),
                         "pair_score": float(min(abs(lfc[g]) / 2.0, 1.0) * stab[g])})
    return pd.DataFrame(rows, columns=["gene", "type_a", "type_b", "log2fc",
                                       "stability", "pair_score"])


# --------------------------------------------------------------------------- #
# 4. build family-specific panels + scoring
# --------------------------------------------------------------------------- #
@dataclass
class FamilyPanels:
    """Result of :func:`build_family_specific_gene_panels`."""
    family_panels: dict[str, list[str]] = field(default_factory=dict)
    gene_scores: pd.DataFrame = field(default_factory=pd.DataFrame)
    selected_genes: pd.DataFrame = field(default_factory=pd.DataFrame)
    pairwise_markers: pd.DataFrame = field(default_factory=pd.DataFrame)
    gene_weights: pd.DataFrame = field(default_factory=pd.DataFrame)
    excluded_genes: pd.DataFrame = field(default_factory=pd.DataFrame)
    marker_support_by_family: pd.DataFrame = field(default_factory=pd.DataFrame)
    marker_support_by_subtype: pd.DataFrame = field(default_factory=pd.DataFrame)
    marker_support_by_pair: pd.DataFrame = field(default_factory=pd.DataFrame)
    families_skipped: dict[str, str] = field(default_factory=dict)

    def save(self, out_dir) -> dict[str, Path]:
        """Write the Part-3 output tables; returns ``{name: path}``."""
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        tables = {
            "within_family_marker_scores": self.gene_scores,
            "within_family_selected_genes": self.selected_genes,
            "within_family_pairwise_markers": self.pairwise_markers,
            "within_family_gene_weights": self.gene_weights,
            "within_family_excluded_genes": self.excluded_genes,
            "marker_support_by_family": self.marker_support_by_family,
            "marker_support_by_pair": self.marker_support_by_pair,
        }
        paths = {}
        for name, df in tables.items():
            p = out / f"{name}.tsv"
            df.to_csv(p, sep="\t", index=False)
            paths[name] = p
        return paths


def within_family_resolution_summary(
    fine_ref,
    mapping: Mapping[str, str],
    family_gene_panels: Mapping[str, list],
    *,
    unresolved_threshold: float = 0.10,
    min_discriminating_genes: int = 10,
    within_family_spillover_threshold: float = 0.30,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare global vs within-family separability/resolvability per family.

    Computes :func:`evaluate_within_family_resolvability` twice — once on the
    global gene set and once restricted to each family's panel — and returns
    ``(global_resolvability, within_family_resolvability, summary)``.  The
    *summary* flags families whose mean separability improved and whose
    resolvable verdict changed.  Improvement is reported, never assumed.
    """
    from tissueresolve.reference.hierarchy import evaluate_within_family_resolvability
    kw = dict(unresolved_threshold=unresolved_threshold,
              min_discriminating_genes=min_discriminating_genes,
              within_family_spillover_threshold=within_family_spillover_threshold)
    glob = evaluate_within_family_resolvability(fine_ref, dict(mapping), **kw)
    within = evaluate_within_family_resolvability(
        fine_ref, dict(mapping), family_gene_panels=dict(family_gene_panels), **kw)
    rows = []
    for fam in glob.index:
        g_sep = float(glob.loc[fam, "mean_separability"]) if pd.notna(
            glob.loc[fam, "mean_separability"]) else float("nan")
        w_sep = (float(within.loc[fam, "mean_separability"])
                 if fam in within.index and pd.notna(
                     within.loc[fam, "mean_separability"]) else float("nan"))
        g_res = bool(glob.loc[fam, "resolvable"])
        w_res = bool(within.loc[fam, "resolvable"]) if fam in within.index else g_res
        improved = (not np.isnan(g_sep) and not np.isnan(w_sep) and w_sep > g_sep)
        rows.append({
            "family": fam,
            "global_mean_separability": round(g_sep, 4) if not np.isnan(g_sep) else np.nan,
            "within_family_mean_separability": round(w_sep, 4) if not np.isnan(w_sep) else np.nan,
            "global_resolvable": g_res,
            "within_family_resolvable": w_res,
            "separability_improved": bool(improved),
            "verdict_changed": bool(g_res != w_res),
            "n_panel_genes": int(within.loc[fam, "n_panel_genes"])
            if fam in within.index and pd.notna(within.loc[fam, "n_panel_genes"]) else np.nan,
        })
    summary = pd.DataFrame(rows)
    return glob, within, summary


def _specificity(mean_by_sub: np.ndarray) -> np.ndarray:
    """Tau specificity (0 = ubiquitous, 1 = single-subtype) per gene.

    *mean_by_sub*: (n_subtypes, n_genes) non-negative mean expression.
    """
    mx = mean_by_sub.max(axis=0)
    mx_safe = np.where(mx > 0, mx, 1.0)
    ratios = mean_by_sub / mx_safe
    k = mean_by_sub.shape[0]
    tau = (1.0 - ratios).sum(axis=0) / max(k - 1, 1)
    tau[mx <= 0] = 0.0
    return tau


def build_family_specific_gene_panels(
    adata_ref,
    broad_col: str,
    fine_col: str,
    *,
    sample_col: Optional[str] = None,
    query_adata=None,
    n_hvgs_per_family: int = 500,
    top_n_de_per_pair: int = 50,
    min_cells_per_subtype: int = 20,
    min_query_detection: float = 0.01,
    min_logfc: float = 0.25,
    min_detection: float = 0.05,
    remove_mito: bool = True,
    remove_ribo: bool = True,
    remove_stress: bool = True,
    genome: str = "hg38",
    max_genes_per_family: Optional[int] = None,
) -> FamilyPanels:
    """Build one within-family gene panel per multi-subtype broad family.

    For each family: within-family HVGs ∪ pairwise discriminative genes,
    minus mito/ribo/stress and minus genes undetected in the query (when a
    *query_adata* is given), scored and weighted (0–1 within the family).

    Families with <2 qualifying subtypes are skipped (recorded in
    ``families_skipped``) and **not** sub-resolved.
    """
    genes_all = np.asarray([str(g) for g in adata_ref.var_names])
    families = sorted(adata_ref.obs[broad_col].astype(str).unique())
    stress = _stress_genes(genome) if remove_stress else frozenset()

    # query detection (optional)
    query_det = None
    if query_adata is not None:
        qc = _dense_counts(query_adata)
        qgenes = [str(g) for g in query_adata.var_names]
        qser = pd.Series((qc > 0).mean(axis=0), index=qgenes)
        query_det = qser

    panels: dict[str, list[str]] = {}
    skipped: dict[str, str] = {}
    score_rows, sel_rows, weight_rows, excl_rows = [], [], [], []
    pair_frames, support_fam, support_sub, support_pair = [], [], [], []

    for fam in families:
        members = _family_members(adata_ref, broad_col, fine_col, fam,
                                  min_cells_per_subtype)
        if len(members) < 2:
            skipped[fam] = (f"<2 subtypes with >={min_cells_per_subtype} cells "
                            f"({len(members)} found)")
            support_fam.append({"family": fam, "n_subtypes": len(members),
                                "n_panel_genes": 0, "n_pairwise_markers": 0,
                                "resolvable_candidate": False})
            continue

        # per-cell matrices for this family
        mask = (adata_ref.obs[broad_col].astype(str) == str(fam)).to_numpy()
        counts = _dense_counts(adata_ref)[mask]
        Xln = _normalize_log(counts)
        labels = adata_ref.obs.loc[mask, fine_col].astype(str).to_numpy()
        means, _det = _subtype_stats(Xln, counts, labels, members)
        mean_by_sub = np.vstack([means[s] for s in members])
        detection = (counts > 0).mean(axis=0)
        var = Xln.var(axis=0)
        spec = _specificity(mean_by_sub)

        hvg = select_within_family_hvgs(
            adata_ref, broad_col, fine_col, fam, sample_col=sample_col,
            n_hvgs=n_hvgs_per_family, min_cells_per_subtype=min_cells_per_subtype,
            min_detection=min_detection, max_genes=max_genes_per_family)
        pair = select_pairwise_discriminative_genes_within_family(
            adata_ref, broad_col, fine_col, fam, sample_col=sample_col,
            top_n_per_pair=top_n_de_per_pair, min_logfc=min_logfc,
            min_detection=min_detection, min_cells_per_subtype=min_cells_per_subtype)
        if not pair.empty:
            pair = pair.assign(family=fam)
            pair_frames.append(pair)

        # candidate genes = HVG ∪ pairwise-discriminative
        cand = set(hvg.index) | set(pair["gene"]) if not pair.empty else set(hvg.index)

        # discriminatory score: fraction of within-family pairs each gene discriminates
        n_pairs = max(len(list(combinations(members, 2))), 1)
        disc_counts = (pair.groupby("gene").size() if not pair.empty
                       else pd.Series(dtype=float))
        gidx = {g: i for i, g in enumerate(genes_all)}

        for g in sorted(cand):
            i = gidx.get(g)
            if i is None:
                continue
            reasons = []
            if remove_mito and _is_mito(g):
                reasons.append("mito")
            if remove_ribo and _is_ribo(g):
                reasons.append("ribo")
            if remove_stress and g.upper() in stress:
                reasons.append("stress")
            # genes absent from the query are treated as detection 0 (excluded)
            q_det = float(query_det.get(g, 0.0)) if query_det is not None else np.nan
            if query_det is not None and q_det < min_query_detection:
                reasons.append(f"query_detection<{min_query_detection}")
            var_score = float(hvg.loc[g, "hvg_score"]) if g in hvg.index else 0.0
            spec_score = float(spec[i])
            disc_score = float(disc_counts.get(g, 0)) / n_pairs
            det_ref = float(detection[i])
            redundancy_penalty = 0.0   # light: handled by tau specificity
            stress_pen = 1.0 if (g.upper() in stress) else 0.0
            row = {
                "family": fam, "gene": g,
                "within_family_variance_score": var_score,
                "fine_subtype_specificity_score": spec_score,
                "pairwise_discriminatory_score": disc_score,
                "detection_in_reference": det_ref,
                "detection_in_query": q_det,
                "donor_or_sample_stability_score": float(
                    pair.loc[pair["gene"] == g, "stability"].mean()) if not pair.empty
                    and (pair["gene"] == g).any() else np.nan,
                "batch_penalty": 0.0,
                "redundancy_penalty": redundancy_penalty,
                "stress_gene_penalty": stress_pen,
                "excluded": bool(reasons),
                "exclusion_reason": ";".join(reasons),
            }
            score_rows.append(row)
            if reasons:
                excl_rows.append({"family": fam, "gene": g,
                                  "reason": ";".join(reasons)})

        # final weight (0-1 within family) over NON-excluded candidate genes
        fam_scores = pd.DataFrame([r for r in score_rows if r["family"] == fam])
        kept = fam_scores[~fam_scores["excluded"]].copy()
        if kept.empty:
            skipped[fam] = "no within-family genes survived filtering"
            support_fam.append({"family": fam, "n_subtypes": len(members),
                                "n_panel_genes": 0, "n_pairwise_markers":
                                int(0 if pair.empty else pair["gene"].nunique()),
                                "resolvable_candidate": False})
            continue
        raw = (0.35 * kept["pairwise_discriminatory_score"]
               + 0.30 * kept["fine_subtype_specificity_score"]
               + 0.25 * kept["within_family_variance_score"]
               + 0.10 * kept["detection_in_reference"]
               - 0.25 * kept["stress_gene_penalty"])
        raw = raw.clip(lower=0.0)
        rng = raw.max() - raw.min()
        weight = (raw - raw.min()) / rng if rng > _EPS else pd.Series(
            np.ones(len(raw)), index=raw.index)
        kept = kept.assign(final_within_family_weight=weight.values)
        kept = kept.sort_values("final_within_family_weight", ascending=False)
        if max_genes_per_family:
            kept = kept.head(int(max_genes_per_family))
        panel = kept["gene"].tolist()
        panels[fam] = panel
        for _, r in kept.iterrows():
            weight_rows.append({"family": fam, "gene": r["gene"],
                                "final_within_family_weight":
                                float(r["final_within_family_weight"])})
            sel_rows.append({"family": fam, "gene": r["gene"]})

        # marker support tables
        support_fam.append({
            "family": fam, "n_subtypes": len(members),
            "n_panel_genes": len(panel),
            "n_pairwise_markers": int(0 if pair.empty else pair["gene"].nunique()),
            "mean_pairwise_score": float(0.0 if pair.empty else pair["pair_score"].mean()),
            "resolvable_candidate": len(panel) >= 10,
        })
        for s in members:
            n_de = 0
            if not pair.empty:
                n_de = int(((pair["type_a"] == s) | (pair["type_b"] == s)).sum())
            support_sub.append({"family": fam, "subtype": s,
                                "n_pairwise_markers": n_de})
        if not pair.empty:
            for (a, b), grp in pair.groupby(["type_a", "type_b"]):
                support_pair.append({"family": fam, "type_a": a, "type_b": b,
                                     "n_markers": int(grp["gene"].nunique()),
                                     "mean_score": float(grp["pair_score"].mean())})

    scores = pd.DataFrame(score_rows)
    return FamilyPanels(
        family_panels=panels,
        gene_scores=scores,
        selected_genes=pd.DataFrame(sel_rows, columns=["family", "gene"]),
        pairwise_markers=(pd.concat(pair_frames, ignore_index=True)
                          if pair_frames else pd.DataFrame(
                              columns=["gene", "type_a", "type_b", "log2fc",
                                       "stability", "pair_score", "family"])),
        gene_weights=pd.DataFrame(weight_rows,
                                  columns=["family", "gene",
                                           "final_within_family_weight"]),
        excluded_genes=pd.DataFrame(excl_rows, columns=["family", "gene", "reason"]),
        marker_support_by_family=pd.DataFrame(support_fam),
        marker_support_by_subtype=pd.DataFrame(support_sub,
                                               columns=["family", "subtype",
                                                        "n_pairwise_markers"]),
        marker_support_by_pair=pd.DataFrame(support_pair,
                                            columns=["family", "type_a", "type_b",
                                                     "n_markers", "mean_score"]),
        families_skipped=skipped,
    )
