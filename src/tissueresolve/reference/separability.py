"""
Pairwise cell-type separability analysis.

The single most common silent failure mode in deconvolution is feeding a
reference containing two or more cell types whose transcriptomic profiles are
nearly identical.  When that happens the model has no signal to split them and
produces noise-driven, unstable estimates for the pair — while aggregate
metrics may look fine.

This module:

1. **Screens** all K*(K-1)/2 type pairs using the Bhattacharyya coefficient
   (BC) on normalised CPM profiles.  BC = 1 means identical; BC = 0 means
   orthogonal.  Jeffreys divergence and Pearson r are reported alongside.

2. **Warns** via :class:`SeparabilityWarning` when any pair has BC > threshold
   (default 0.90).  The warning is a proper Python warning so that downstream
   code can filter or escalate it.

3. **Merges** non-separable pairs into a single merged type, updating the
   reference.

Ported from SpatCAR v1.1 with the following fix
-------------------------------------------------
``merge_nonseparable_types`` in the original SpatCAR code contained a loop
that effectively set ``phi_new = ref_data.phi_g`` (global dispersion) by
re-applying ``np.maximum(phi_new, ref_data.phi_g)`` on every iteration
regardless of which group was being processed.  In TissueResolve the fix is
explicit: since ``phi_g`` is a global per-gene dispersion (not per-cell-type),
the merged reference inherits it directly via a single ``phi_g.copy()`` call.

References
----------
Bhattacharyya 1943; Kailath 1967; Ma & Zhou 2022 (CARD Supplementary).
"""
from __future__ import annotations

import logging
import warnings

import numpy as np

from tissueresolve.results import PairSeparability, SeparabilityReport, ReferenceSignature

__all__ = [
    "SeparabilityWarning",
    "compute_separability",
    "merge_nonseparable_types",
    "separability_heatmap_data",
]

logger = logging.getLogger("tissueresolve.reference.separability")


class SeparabilityWarning(UserWarning):
    """Raised when a reference contains poorly separable cell-type pairs.

    Always emit; never suppress via the pipeline.  Downstream code may
    catch or re-raise this warning, but it must never be silenced without
    explicit user opt-in.
    """


def compute_separability(
    ref: ReferenceSignature,
    *,
    warn_threshold: float = 0.90,
    raise_on_critical: bool = False,
) -> SeparabilityReport:
    """Compute pairwise cell-type separability for *ref*.

    Parameters
    ----------
    ref:
        Reference built by :class:`~tissueresolve.reference.build.ReferenceBuilder`.
    warn_threshold:
        Emit a :class:`SeparabilityWarning` for pairs with BC > this value.
        Set to ``1.1`` to suppress all warnings (not recommended).
    raise_on_critical:
        If ``True``, raise :exc:`ValueError` when any pair has BC > 0.97.
        Useful in automated pipelines where fitting would waste compute.

    Returns
    -------
    SeparabilityReport
        Pairwise results sorted worst-first (highest BC first).
    """
    K = ref.n_cell_types
    cell_types = ref.cell_types
    R_log = ref.as_R_log().astype(np.float64)   # (K, G)
    R_cpm = ref.as_R_cpm().astype(np.float64)   # (K, G)

    eps = 1e-10
    # Normalised probability vectors for Bhattacharyya coefficient
    P_norm = R_cpm / (R_cpm.sum(axis=1, keepdims=True) + eps)  # (K, G)

    pairs: list[PairSeparability] = []

    for i in range(K):
        for j in range(i + 1, K):
            p = P_norm[i] + eps
            q = P_norm[j] + eps

            # Bhattacharyya coefficient
            bc = float(np.sum(np.sqrt(p * q)))
            bc = float(np.clip(bc, 0.0, 1.0))

            # Jeffreys divergence (symmetric KL, nats)
            jd = float(np.sum(p * np.log(p / q)) + np.sum(q * np.log(q / p)))
            jd = max(0.0, jd)

            # Pearson r on log1p-CPM profiles
            ri, rj = R_log[i], R_log[j]
            ri_c = ri - ri.mean()
            rj_c = rj - rj.mean()
            denom = float(np.sqrt((ri_c ** 2).sum() * (rj_c ** 2).sum())) + eps
            pearson_r = float(np.dot(ri_c, rj_c) / denom)

            # Discriminating genes: |log2FC| > 1 between the two types
            lfc = R_log[i] - R_log[j]
            n_disc = int((np.abs(lfc) > 1.0).sum())

            pairs.append(PairSeparability(
                type_a=cell_types[i],
                type_b=cell_types[j],
                bhattacharyya_coeff=bc,
                jeffreys_divergence=jd,
                pearson_r=pearson_r,
                n_discriminating_genes=n_disc,
            ))

    # Sort worst-first (highest BC)
    pairs.sort(key=lambda p: p.bhattacharyya_coeff, reverse=True)
    report = SeparabilityReport(pairs=pairs)

    logger.info("Separability:\n%s", report.summary_str())

    problematic = [p for p in pairs if p.bhattacharyya_coeff > warn_threshold]
    if problematic:
        n_crit = sum(1 for p in problematic if p.bhattacharyya_coeff > 0.97)
        n_high = len(problematic) - n_crit
        worst = sorted(problematic, key=lambda p: -p.bhattacharyya_coeff)[:5]
        top_lines = "\n".join(
            f"  {p.type_a!r} vs {p.type_b!r}: "
            f"BC={p.bhattacharyya_coeff:.4f}, "
            f"separability_score={p.discriminability_score:.4f}, "
            f"J={p.jeffreys_divergence:.3f}, "
            f"{p.n_discriminating_genes} discriminating genes"
            for p in worst
        )
        more = (f"\n  … and {len(problematic) - len(worst)} more pair(s)."
                if len(problematic) > len(worst) else "")
        msg = (
            f"Reference contains {len(problematic)} poorly separable cell-type "
            f"pair(s) (BC > {warn_threshold}): {n_crit} CRITICAL (BC>0.97), "
            f"{n_high} HIGH.  Subtype-level proportions for these types may be "
            "unreliable.\nTop pairs:\n" + top_lines + more
            + "\n\nRecommended action — use the resolution-aware recommender for a "
            "machine-readable merge plan instead of merging by hand:\n"
            "  from tissueresolve.reference.resolution import (\n"
            "      build_resolution_report, recommend_cell_type_merges,\n"
            "      assign_resolution_families, write_recommended_merges)\n"
            "It writes recommended_merges.tsv / cell_type_families.tsv under the "
            "run's resolution/ outputs, or run the analysis with "
            "--resolution-mode {suggest|auto|hierarchical}.\n"
            "Consider interpreting confusable types at the family level rather "
            "than as confident subtypes."
        )
        warnings.warn(msg, SeparabilityWarning, stacklevel=2)
        logger.warning("SEPARABILITY WARNING:\n%s", msg)

    critical = [p for p in pairs if p.bhattacharyya_coeff > 0.97]
    if critical and raise_on_critical:
        raise ValueError(
            f"Reference has {len(critical)} critically non-separable pair(s): "
            + ", ".join(f"{p.type_a}/{p.type_b}" for p in critical)
        )

    return report


def merge_nonseparable_types(
    ref: ReferenceSignature,
    report: SeparabilityReport,
    *,
    merge_threshold: float = 0.90,
    merged_label_sep: str = "+",
) -> tuple[ReferenceSignature, dict[str, str]]:
    """Merge pairs of non-separable types into single combined types.

    Parameters
    ----------
    ref:
        Original reference.
    report:
        Output of :func:`compute_separability`.
    merge_threshold:
        Pairs with BC > this value are merged (default 0.90).
    merged_label_sep:
        String separator for merged type names (default ``"+"``).

    Returns
    -------
    tuple (merged_ref, merge_map)
        ``merged_ref`` — new :class:`~tissueresolve.results.ReferenceSignature`
        with merged types.
        ``merge_map`` — maps original type names → merged label, e.g.
        ``{"ExcL2": "ExcL2+ExcL3", "ExcL3": "ExcL2+ExcL3"}``.

    Bug fix (P0-3 in AUDIT_AND_MIGRATION_PLAN.md)
    -----------------------------------------------
    The original SpatCAR implementation ran a loop over all original cell
    types and applied ``np.maximum(phi_new, ref_data.phi_g)`` on every
    iteration, effectively setting ``phi_new = ref_data.phi_g`` globally
    regardless of group membership.  In TissueResolve ``phi_g`` is a single
    per-gene dispersion vector (not per-cell-type), so the merged reference
    inherits it via a single ``phi_g.copy()`` call — no loop required.
    """
    # --- union-find for grouping ----------------------------------------
    parent: dict[str, str] = {ct: ct for ct in ref.cell_types}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(x: str, y: str) -> None:
        parent[find(x)] = find(y)

    for p in report.pairs:
        if p.bhattacharyya_coeff > merge_threshold:
            union(p.type_a, p.type_b)

    # Build groups
    groups: dict[str, list[str]] = {}
    for ct in ref.cell_types:
        root = find(ct)
        groups.setdefault(root, []).append(ct)

    merge_map: dict[str, str] = {}
    for root, members in groups.items():
        label = merged_label_sep.join(sorted(members))
        for m in members:
            merge_map[m] = label

    merged_types = sorted(set(merge_map.values()))
    K_new = len(merged_types)
    G = ref.n_genes
    type_to_new_idx = {ct: merged_types.index(ct) for ct in merged_types}

    # --- build merged R_cpm (cell-count-weighted average) ---------------
    R_cpm_src = ref.as_R_cpm()    # (K, G)
    R_cpm_new = np.zeros((K_new, G), dtype=np.float32)
    n_cells_new: dict[str, int] = {ct: 0 for ct in merged_types}

    for orig_type in ref.cell_types:
        merged_label = merge_map[orig_type]
        new_idx = type_to_new_idx[merged_label]
        n_k = ref.n_cells_per_type.get(orig_type, 1)
        k_orig = ref.cell_types.index(orig_type)
        R_cpm_new[new_idx] += R_cpm_src[k_orig] * n_k
        n_cells_new[merged_label] += n_k

    for new_idx, merged_label in enumerate(merged_types):
        n = max(n_cells_new[merged_label], 1)
        R_cpm_new[new_idx] /= n

    R_log_new = np.log1p(R_cpm_new).astype(np.float32)

    # --- re-derive phi (L1-normalised) from merged R_cpm ----------------
    phi_raw = R_cpm_new.T.copy()     # (G, K_new)
    col_sums = phi_raw.sum(axis=0, keepdims=True)
    col_sums = np.where(col_sums == 0, 1.0, col_sums)
    phi_new = (phi_raw / col_sums).astype(np.float64)

    # --- phi_g: global per-gene dispersion — BUG FIX --------------------
    # phi_g is a single (G,) vector shared across all cell types.  The correct
    # action is a direct copy — NOT iterating over all types and applying a
    # global maximum (which was the SpatCAR bug, P0-3).
    phi_g_new = ref.phi_g.copy() if ref.phi_g is not None else None

    # --- donor_cv: cannot be meaningfully merged; set to None -----------
    # Future work: pool donor_cv within each merged group if needed.
    donor_cv_new = None

    merged_ref = ReferenceSignature(
        gene_names=list(ref.gene_names),
        cell_types=merged_types,
        phi=phi_new,
        R_cpm=R_cpm_new,
        R_log=R_log_new,
        phi_g=phi_g_new,
        donor_cv=donor_cv_new,
        n_cells_per_type=n_cells_new,
        genome=ref.genome,
    )

    n_merged = len(ref.cell_types) - K_new
    logger.info(
        "Merged %d type(s) into %d combined type(s): %s",
        n_merged, K_new,
        ", ".join(
            f"{ct!r}" for ct in merged_types if merged_label_sep in ct
        ),
    )
    return merged_ref, merge_map


def separability_heatmap_data(
    report: SeparabilityReport,
    cell_types: list[str],
) -> np.ndarray:
    """Return a K × K Bhattacharyya-coefficient matrix for plotting.

    Diagonal is 1.0 (type with itself).  Off-diagonal is BC for that pair.

    Parameters
    ----------
    report:
        Output of :func:`compute_separability`.
    cell_types:
        Ordered list of type names (defines row/column order).

    Returns
    -------
    np.ndarray
        Shape ``(K, K)`` float32.
    """
    K = len(cell_types)
    idx = {ct: i for i, ct in enumerate(cell_types)}
    M = np.eye(K, dtype=np.float32)
    for p in report.pairs:
        if p.type_a in idx and p.type_b in idx:
            i, j = idx[p.type_a], idx[p.type_b]
            M[i, j] = M[j, i] = float(p.bhattacharyya_coeff)
    return M
