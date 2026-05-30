"""
Shared, importable logic for the real-data validation harness.

This module holds the reusable (and unit-tested) building blocks so the
numbered ``NN_*.py`` scripts stay thin CLI orchestrators.  Nothing here
touches the network; downloads live in ``00_download_data.py`` and run only
behind an explicit flag.

It also discovers the TissueResolve public API at call time rather than
assuming it — every TissueResolve entry point used here was confirmed against
``src/`` and the tests.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HARNESS_DIR = Path(__file__).resolve().parent.parent  # examples/real_breast_cancer/
DATA_DIR = HARNESS_DIR / "data"
REFERENCE_DATA_DIR = DATA_DIR / "reference"
SPATIAL_DATA_DIR = DATA_DIR / "spatial"
DERIVED_DIR = DATA_DIR / "derived"

OUTPUTS_DIR = HARNESS_DIR / "outputs"
OUT_REFERENCE_DIR = OUTPUTS_DIR / "reference"
OUT_BULK_DIR = OUTPUTS_DIR / "bulk"
OUT_SPATIAL_DIR = OUTPUTS_DIR / "spatial"
OUT_SUMMARY_DIR = OUTPUTS_DIR / "validation_summary"

MANIFEST_PATH = DATA_DIR / "download_manifest.json"

REFERENCE_H5AD = REFERENCE_DATA_DIR / "breast_cancer_sc_reference.h5ad"
SPATIAL_H5AD = SPATIAL_DATA_DIR / "human_breast_cancer_1.h5ad"

SAVED_REFERENCE_DIR = OUT_REFERENCE_DIR / "breast_cancer_reference"

PSEUDOBULK_COUNTS = DERIVED_DIR / "pseudobulk_counts.tsv"
PSEUDOBULK_TRUE_PROPS = DERIVED_DIR / "pseudobulk_true_proportions.tsv"
PSEUDOBULK_META = DERIVED_DIR / "pseudobulk_metadata.tsv"

ALL_DIRS = [
    REFERENCE_DATA_DIR, SPATIAL_DATA_DIR, DERIVED_DIR,
    OUT_REFERENCE_DIR, OUT_BULK_DIR, OUT_SPATIAL_DIR, OUT_SUMMARY_DIR,
]

# Candidate cell-type annotation columns, in priority order.
CELLTYPE_COL_CANDIDATES = [
    "cell_type_major",
    "cell_type",
    "celltype",
    "annotation",
    "author_cell_type",
    "cell_type_fine",
    "cell_type_ontology_term_id",
]

# Downsampling caps for an oversized reference (see README / PLAN).
MAX_REFERENCE_CELLS = 30_000
MAX_REFERENCE_GENES = 5_000
MIN_CELLS_PER_TYPE = 50


# ---------------------------------------------------------------------------
# Environment / flags
# ---------------------------------------------------------------------------

REAL_DATA_ENV = "TISSUERESOLVE_RUN_REAL_DATA"


def real_data_enabled(flag: bool = False) -> bool:
    """True iff the ``--run-real-data`` flag or the env var is set.

    The default (no flag, env unset) keeps every script and test offline.
    """
    if flag:
        return True
    return os.environ.get(REAL_DATA_ENV, "").strip() in {"1", "true", "True", "yes"}


def ensure_dirs() -> None:
    """Create all harness data/output directories (idempotent)."""
    for d in ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Manifest schema
# ---------------------------------------------------------------------------

MANIFEST_SCHEMA_VERSION = 1
_DATASET_REQUIRED_KEYS = {
    "dataset_name", "source", "url_or_id", "access_date",
    "package_versions", "filters_applied", "downsampling", "file", "downloaded",
}


def package_versions() -> dict[str, str]:
    """Record versions of packages relevant to provenance."""
    import importlib.metadata as md

    out: dict[str, str] = {}
    for pkg in ("tissueresolve", "anndata", "scanpy", "numpy", "scipy",
                "pandas", "cellxgene-census", "scikit-learn"):
        try:
            out[pkg] = md.version(pkg)
        except Exception:
            out[pkg] = "not installed"
    return out


def default_manifest(*, access_date: Optional[str] = None) -> dict[str, Any]:
    """A schema-valid manifest skeleton with both dataset entries un-downloaded."""
    if access_date is None:
        from datetime import datetime, timezone

        access_date = datetime.now(timezone.utc).isoformat()
    versions = package_versions()

    def _rel(p: Path) -> str:
        try:
            return str(p.relative_to(HARNESS_DIR))
        except ValueError:
            return str(p)
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "generated_by": "00_download_data.py",
        "datasets": {
            "reference": {
                "dataset_name": "Human breast cancer single-cell atlas (CZ CELLxGENE)",
                "source": "CZ CELLxGENE Census (cellxgene-census)",
                "url_or_id": "https://cellxgene.cziscience.com/  (census collection: breast cancer)",
                "access_date": access_date,
                "package_versions": versions,
                "filters_applied": [
                    "human (Homo sapiens)",
                    "tissue: breast",
                    "primary cells with cell-type annotation",
                ],
                "downsampling": {
                    "max_cells": MAX_REFERENCE_CELLS,
                    "max_genes": MAX_REFERENCE_GENES,
                    "applied": False,
                },
                "file": _rel(REFERENCE_H5AD),
                "downloaded": False,
            },
            "spatial": {
                "dataset_name": "10x Visium Human Breast Cancer 1",
                "source": "OpenProblems / 10x Genomics public Visium",
                "url_or_id": "openproblems: human_breast_cancer_1 "
                             "(≈4263 spots × 14906 genes; counts in layers['counts'])",
                "access_date": access_date,
                "package_versions": versions,
                "filters_applied": [],
                "downsampling": {"applied": False},
                "file": _rel(SPATIAL_H5AD),
                "downloaded": False,
            },
        },
    }


def validate_manifest(manifest: dict[str, Any]) -> None:
    """Raise ``ValueError`` if *manifest* does not match the expected schema."""
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a dict")
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            f"manifest schema_version must be {MANIFEST_SCHEMA_VERSION}, "
            f"got {manifest.get('schema_version')!r}"
        )
    datasets = manifest.get("datasets")
    if not isinstance(datasets, dict) or not {"reference", "spatial"} <= set(datasets):
        raise ValueError("manifest.datasets must contain 'reference' and 'spatial'")
    for name, entry in datasets.items():
        missing = _DATASET_REQUIRED_KEYS - set(entry)
        if missing:
            raise ValueError(f"manifest dataset {name!r} missing keys: {sorted(missing)}")


def write_manifest(manifest: dict[str, Any], path: Path = MANIFEST_PATH) -> Path:
    validate_manifest(manifest)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path


def read_manifest(path: Path = MANIFEST_PATH) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# AnnData helpers
# ---------------------------------------------------------------------------


#: Sentinel override values that mean "auto-detect" rather than a literal column.
AUTO_DETECT_TOKENS = {None, "", "auto"}


def detect_cell_type_col(obs: "pd.DataFrame", override: Optional[str] = None) -> str:
    """Return the cell-type annotation column to use.

    *override* is honoured as a literal column name **unless** it is an
    auto-detect sentinel (``None``, ``""``, or ``"AUTO"`` in any case), in which
    case the documented candidate list is searched in priority order.  Raises a
    clear error listing what was searched and what is available.
    """
    is_auto = override is None or (
        isinstance(override, str) and override.strip().lower() in {"", "auto"}
    )
    if not is_auto:
        if override not in obs.columns:
            raise KeyError(
                f"--cell-type-col {override!r} not in obs columns: {list(obs.columns)}"
            )
        return override
    for cand in CELLTYPE_COL_CANDIDATES:
        if cand in obs.columns:
            return cand
    raise KeyError(
        "Could not auto-detect a cell-type column.  Searched "
        f"{CELLTYPE_COL_CANDIDATES}; available obs columns: {list(obs.columns)}.  "
        "Pass --cell-type-col explicitly."
    )


def resolve_counts_matrix(adata) -> tuple[Any, str]:
    """Return ``(matrix, source)`` of raw counts from an AnnData.

    Prefers ``layers['counts']``, then ``raw.X`` (when gene-compatible), then
    ``X``.  Never silently rescales; the chosen source is returned for logging.
    """
    if "counts" in getattr(adata, "layers", {}):
        return adata.layers["counts"], "layers['counts']"
    raw = getattr(adata, "raw", None)
    if raw is not None and getattr(raw, "X", None) is not None \
            and raw.X.shape[1] == adata.n_vars:
        return raw.X, "raw.X"
    return adata.X, "X"


def adata_with_raw_counts(adata):
    """Return an AnnData whose ``.X`` is the best available raw-count matrix.

    Used before :meth:`ReferenceBuilder.build_from_adata`, which reads
    ``adata.X`` directly and expects integer UMI counts.
    """
    mat, source = resolve_counts_matrix(adata)
    if source == "X":
        return adata, source
    out = adata.copy()
    out.X = mat
    return out, source


# ---------------------------------------------------------------------------
# Reference preparation (script 01 logic)
# ---------------------------------------------------------------------------


@dataclass
class ReferencePrep:
    reference: Any            # tissueresolve.results.ReferenceSignature
    cell_type_col: str
    cell_type_counts: "pd.Series"
    summary: "pd.DataFrame"
    counts_source: str


def prepare_reference(
    adata,
    *,
    cell_type_col: Optional[str] = None,
    min_cells: int = MIN_CELLS_PER_TYPE,
    estimate_overdispersion: bool = True,
) -> ReferencePrep:
    """Build a TissueResolve reference from a single-cell AnnData.

    Uses the discovered public API:
    ``ReferenceBuilder(ReferenceConfig(...)).build_from_adata(adata,
    estimate_overdispersion=...)``.  ``estimate_overdispersion=True`` so the
    same reference serves the spatial NB model.
    """
    from tissueresolve.config import ReferenceConfig
    from tissueresolve.reference.build import ReferenceBuilder

    col = detect_cell_type_col(adata.obs, cell_type_col)
    counts = adata.obs[col].astype(str).value_counts()

    counts_adata, source = adata_with_raw_counts(adata)
    cfg = ReferenceConfig(celltype_col=col, min_cells=min_cells)
    builder = ReferenceBuilder(cfg)
    ref = builder.build_from_adata(
        counts_adata, estimate_overdispersion=estimate_overdispersion
    )

    summary = pd.DataFrame({
        "metric": ["n_cells", "n_genes", "n_cell_types", "cell_type_col",
                   "min_cells", "counts_source"],
        "value": [int(adata.n_obs), int(ref.n_genes), int(ref.n_cell_types),
                  col, int(min_cells), source],
    })
    return ReferencePrep(
        reference=ref, cell_type_col=col,
        cell_type_counts=counts, summary=summary, counts_source=source,
    )


# ---------------------------------------------------------------------------
# Pseudobulk generation (script 02 logic)
# ---------------------------------------------------------------------------

_REGIMES = ("easy", "medium", "hard")


def _regime_weights(regime: str, cell_types: list[str], rng: np.random.Generator) -> np.ndarray:
    """Target cell-mixing weights over cell types for a difficulty regime."""
    K = len(cell_types)
    if regime == "easy":
        # one dominant type (~0.7–0.9), rest small
        w = rng.dirichlet(np.full(K, 0.15))
    elif regime == "medium":
        # 2–3 types share most mass
        w = rng.dirichlet(np.full(K, 0.7))
    else:  # hard — near-uniform / many types
        w = rng.dirichlet(np.full(K, 4.0))
    return w / w.sum()


def generate_pseudobulk(
    adata,
    cell_type_col: str,
    *,
    n_per_regime: int = 4,
    n_cells: int = 200,
    seed: int = 0,
) -> tuple["pd.DataFrame", "pd.DataFrame", "pd.DataFrame"]:
    """Build pseudobulk mixtures with known ground-truth proportions.

    Cells are sampled (with replacement) per mixture according to target
    cell-mixing weights, and their raw counts are summed.  The recorded
    **ground-truth proportions are mRNA proportions** (fraction of total
    counts contributed by each cell type), because that is what the bulk
    deconvolver estimates — bulk outputs are mRNA proportions, not cell
    fractions.  Cell-count proportions are also recorded in the metadata for
    transparency.

    Returns
    -------
    (counts_df, true_props_df, meta_df)
        ``counts_df``: genes × samples integer counts.
        ``true_props_df``: samples × cell types mRNA proportions (rows sum to 1).
        ``meta_df``: per-sample regime, n_cells, total_counts, dominant_type,
        and per-type cell-count proportions.
    """
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    col = cell_type_col
    cell_types = sorted(adata.obs[col].astype(str).unique().tolist())
    gene_names = list(adata.var_names)

    mat, _ = resolve_counts_matrix(adata)
    dense = mat.toarray() if sp.issparse(mat) else np.asarray(mat)
    dense = np.asarray(dense, dtype=np.float64)

    labels = adata.obs[col].astype(str).to_numpy()
    idx_by_type = {ct: np.where(labels == ct)[0] for ct in cell_types}

    counts_cols: dict[str, np.ndarray] = {}
    prop_rows: list[dict[str, float]] = []
    meta_rows: list[dict[str, Any]] = []
    sample_ids: list[str] = []

    for regime in _REGIMES:
        for r in range(n_per_regime):
            sid = f"{regime}_{r:02d}"
            sample_ids.append(sid)
            w = _regime_weights(regime, cell_types, rng)
            # assign each sampled cell to a type, then to a specific cell
            type_choices = rng.choice(len(cell_types), size=n_cells, p=w)
            pseudobulk = np.zeros(dense.shape[1], dtype=np.float64)
            cell_count = np.zeros(len(cell_types), dtype=np.int64)
            type_counts = np.zeros(len(cell_types), dtype=np.float64)
            for k in type_choices:
                pool = idx_by_type[cell_types[k]]
                if len(pool) == 0:
                    continue
                ci = pool[rng.integers(len(pool))]
                vec = dense[ci]
                pseudobulk += vec
                cell_count[k] += 1
                type_counts[k] += vec.sum()

            total = pseudobulk.sum()
            mrna_prop = type_counts / total if total > 0 else np.zeros_like(type_counts)
            cell_prop = cell_count / max(cell_count.sum(), 1)

            counts_cols[sid] = np.rint(pseudobulk).astype(np.int64)
            prop_rows.append(dict(zip(cell_types, mrna_prop)))
            dom = cell_types[int(np.argmax(mrna_prop))] if total > 0 else "NA"
            meta_rows.append({
                "sample": sid, "regime": regime, "n_cells": int(cell_count.sum()),
                "total_counts": float(total), "dominant_type": dom,
                **{f"cellprop_{ct}": float(p) for ct, p in zip(cell_types, cell_prop)},
            })

    counts_df = pd.DataFrame(counts_cols, index=gene_names)
    counts_df.index.name = "gene"
    true_props_df = pd.DataFrame(prop_rows, index=sample_ids)[cell_types]
    true_props_df.index.name = "sample"
    meta_df = pd.DataFrame(meta_rows).set_index("sample")
    return counts_df, true_props_df, meta_df


# ---------------------------------------------------------------------------
# Validation metrics (script 03 logic)
# ---------------------------------------------------------------------------


def align_proportions(
    true_df: "pd.DataFrame", est_df: "pd.DataFrame"
) -> tuple["pd.DataFrame", "pd.DataFrame"]:
    """Restrict both frames to common samples (rows) and cell types (cols)."""
    samples = [s for s in true_df.index if s in set(est_df.index)]
    types = [c for c in true_df.columns if c in set(est_df.columns)]
    if not samples or not types:
        raise ValueError(
            "No overlap between true and estimated proportions "
            f"(samples={len(samples)}, cell_types={len(types)})."
        )
    return true_df.loc[samples, types], est_df.loc[samples, types]


def _safe_corr(a: np.ndarray, b: np.ndarray, kind: str) -> float:
    from scipy.stats import pearsonr, spearmanr

    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    r = pearsonr(a, b)[0] if kind == "pearson" else spearmanr(a, b)[0]
    return float(r)


def compute_bulk_metrics(true_df: "pd.DataFrame", est_df: "pd.DataFrame") -> dict[str, float]:
    """Overall agreement metrics over the flattened sample × cell-type matrix."""
    t, e = align_proportions(true_df, est_df)
    tv = t.to_numpy(dtype=float).ravel()
    ev = e.to_numpy(dtype=float).ravel()
    diff = ev - tv
    return {
        "pearson": _safe_corr(tv, ev, "pearson"),
        "spearman": _safe_corr(tv, ev, "spearman"),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "mae": float(np.mean(np.abs(diff))),
        "signed_bias": float(np.mean(diff)),
        "n_samples": int(t.shape[0]),
        "n_cell_types": int(t.shape[1]),
    }


def per_celltype_metrics(true_df: "pd.DataFrame", est_df: "pd.DataFrame") -> "pd.DataFrame":
    """Per-cell-type RMSE, MAE, signed bias and Pearson r."""
    t, e = align_proportions(true_df, est_df)
    rows = []
    for ct in t.columns:
        tv = t[ct].to_numpy(dtype=float)
        ev = e[ct].to_numpy(dtype=float)
        diff = ev - tv
        rows.append({
            "cell_type": ct,
            "rmse": float(np.sqrt(np.mean(diff ** 2))),
            "mae": float(np.mean(np.abs(diff))),
            "signed_bias": float(np.mean(diff)),
            "pearson": _safe_corr(tv, ev, "pearson"),
        })
    return pd.DataFrame(rows).set_index("cell_type")


def gene_overlap(query_genes, reference_genes) -> dict[str, int]:
    """Report gene-overlap counts between a query and the reference."""
    q, r = set(map(str, query_genes)), set(map(str, reference_genes))
    return {
        "n_query": len(q),
        "n_reference": len(r),
        "n_shared": len(q & r),
        "n_query_only": len(q - r),
        "n_reference_only": len(r - q),
    }


# ---------------------------------------------------------------------------
# Small TSV helpers
# ---------------------------------------------------------------------------


def write_tsv(df: "pd.DataFrame", path: Path, comment: Optional[list[str]] = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for line in comment or []:
            fh.write(f"# {line}\n")
        df.to_csv(fh, sep="\t")
    return path


def write_json(obj: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")
    return path
