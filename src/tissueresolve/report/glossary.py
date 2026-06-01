"""
Variable / metric dictionary for TissueResolve reports.

Plain-language definitions of every metric shown in figures and tables, so the
report explains what each number means without over-interpreting the biology.
"""
from __future__ import annotations

__all__ = ["GLOSSARY", "define", "subset"]

GLOSSARY: dict[str, str] = {
    "Pearson correlation": "Linear agreement between predicted and true "
        "proportions. Higher is better. Only meaningful when ground truth exists.",
    "Spearman correlation": "Rank agreement between predicted and true "
        "proportions. Higher is better; robust to non-linear scaling.",
    "RMSE": "Root mean squared error between predicted and true proportions. "
        "Lower is better; sensitive to large errors.",
    "MAE": "Mean absolute error between predicted and true proportions. Lower is "
        "better; less sensitive to outliers than RMSE.",
    "signed bias": "Mean (predicted − true). Positive = systematic over-estimation.",
    "dominant fraction": "Proportion assigned to the most abundant predicted "
        "cell type or family in a sample/spot.",
    "dominant-type accuracy": "Fraction of samples/spots whose top predicted "
        "type matches the top true type (needs ground truth).",
    "entropy": "How mixed a sample/spot prediction is. Higher entropy means the "
        "signal is spread across more cell types.",
    "near-zero fraction": "Fraction of predicted values close to zero. Very high "
        "values may indicate sparse or overly conservative predictions.",
    "Moran's I": "Spatial autocorrelation. Higher values indicate stronger "
        "spatial structure for a predicted population (spatial only).",
    "unresolved mass": "RNA-derived signal assigned to a broad family but not "
        "confidently attributable to a fine subpopulation.",
    "spillover": "Potential leakage of signal between similar cell types whose "
        "expression profiles are hard to distinguish.",
    "separability": "How distinguishable two cell types are from the reference "
        "expression profiles (here, 1 − Bhattacharyya coefficient).",
    "gene overlap": "Number of genes shared between the query data and the "
        "reference; more shared genes generally means more reliable estimates.",
    "coverage R²": "Goodness of fit of the reconstructed bulk profile to the "
        "observed bulk profile per sample (diagnostic, not accuracy).",
    "concordance": "Agreement between two methods' predictions when no ground "
        "truth exists; high concordance ≠ correctness.",
    "composite score": "Weighted blend of accuracy, robustness, usability, "
        "interpretability, resolution-awareness and runtime (executed/imported "
        "tools only).",
    "mRNA-derived proportion": "Fraction of RNA attributed to a cell type — NOT "
        "an absolute cell-count fraction.",
}


def define(term: str) -> str:
    """Return the definition for *term* (empty string if unknown)."""
    return GLOSSARY.get(term, "")


def subset(terms) -> dict[str, str]:
    """Return a {term: definition} dict for the given *terms* (known ones only)."""
    return {t: GLOSSARY[t] for t in terms if t in GLOSSARY}
