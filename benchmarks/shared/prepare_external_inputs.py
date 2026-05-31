"""
Harmonised input preparation for external tools.

Writes standardized reference / bulk / spatial files (both orientations, plus
CPM/logCPM and the hierarchy mapping) so every external tool consumes the *same*
inputs as TissueResolve.  Sources the existing real breast-cancer outputs or a
toy synthetic scenario; never downloads data.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.shared.io import OUTPUTS_DIR, write_tsv, write_json
from benchmarks.shared import normalization as norm

PREP = OUTPUTS_DIR / "prepared_inputs"


def _signature_frames(ref, mapping=None):
    sig_fine = pd.DataFrame(ref.as_R_cpm().T, index=list(ref.gene_names),
                            columns=list(ref.cell_types))
    sig_broad = None
    if mapping:
        from tissueresolve.reference.hierarchy import aggregate_reference_by_family
        fam = aggregate_reference_by_family(ref, mapping)
        sig_broad = pd.DataFrame(fam.as_R_cpm().T, index=list(fam.gene_names),
                                 columns=list(fam.cell_types))
    return sig_fine, sig_broad


def _export_single_cell_reference(h5ad_path, ref_genes, cell_types, out_dir, *,
                                  max_cells_per_type: int = 60, seed: int = 0) -> dict:
    """Export subsampled single-cell counts (genes×cells) + per-cell metadata.

    Restricts to the reference's cell types and genes; subsamples up to
    *max_cells_per_type* cells per type (recorded).  Writes
    ``reference_counts_genes_by_cells.tsv`` and ``reference_cell_metadata.tsv``
    (columns ``cell_id, cellType, SubjectName``)."""
    import anndata as ad
    import scipy.sparse as sp
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        return {"exported": False, "error": f"h5ad not found: {h5ad_path}"}
    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs
    ct_col = "cell_type" if "cell_type" in obs else None
    subj_col = "donor_id" if "donor_id" in obs else (
        "donor" if "donor" in obs else None)
    if ct_col is None:
        return {"exported": False, "error": "no cell_type column in obs"}
    rng = np.random.default_rng(seed)
    keep_types = set(map(str, cell_types))
    idx = []
    obs_ct = obs[ct_col].astype(str)
    for ct in sorted(keep_types):
        cells = np.where(obs_ct.to_numpy() == ct)[0]
        if len(cells) == 0:
            continue
        take = cells if len(cells) <= max_cells_per_type else \
            np.sort(rng.choice(cells, size=max_cells_per_type, replace=False))
        idx.extend(take.tolist())
    idx = np.sort(np.array(idx, dtype=int))
    sub = adata[idx]
    # gene IDs must match the reference/bulk (symbols). Pick the var column with
    # the best overlap to the reference gene set; default to a 'feature_name'
    # symbol column when present, else var_names.
    ref_set = set(map(str, ref_genes))
    genes = list(map(str, sub.var_names))
    best_overlap = len(set(genes) & ref_set)
    for col in sub.var.columns:
        vals = list(map(str, sub.var[col]))
        ov = len(set(vals) & ref_set)
        if ov > best_overlap:
            best_overlap, genes = ov, vals
    X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    cell_ids = [f"cell{i}" for i in range(X.shape[0])]
    counts = pd.DataFrame(X.T, index=genes, columns=cell_ids)   # genes × cells
    counts.to_csv(out_dir / "reference_counts_genes_by_cells.tsv", sep="\t")
    meta = pd.DataFrame({
        "cell_id": cell_ids,
        "cellType": sub.obs[ct_col].astype(str).to_numpy(),
        "SubjectName": (sub.obs[subj_col].astype(str).to_numpy()
                        if subj_col else "donor1")})
    meta.to_csv(out_dir / "reference_cell_metadata.tsv", sep="\t", index=False)
    return {"exported": True, "n_cells": int(X.shape[0]), "n_genes": len(genes),
            "max_cells_per_type": max_cells_per_type,
            "subject_column": subj_col, "n_subjects": int(meta["SubjectName"].nunique())}


def prepare_inputs(*, use_existing_real_data: bool = True, toy: bool = False) -> dict:
    """Write harmonised inputs; return a dict of written paths + summaries."""
    ref_dir = PREP / "reference"; bulk_dir = PREP / "bulk"; sp_dir = PREP / "spatial"
    for d in (ref_dir, bulk_dir, sp_dir):
        d.mkdir(parents=True, exist_ok=True)

    if toy or not use_existing_real_data:
        from benchmarks.shared.synthetic import toy_reference, toy_bulk, toy_spatial
        ref, mapping = toy_reference()
        bulk, truth = toy_bulk(ref)
        sc, sp_truth = toy_spatial(ref)
        spatial_h5ad = ""
        obs = None
    else:
        from tissueresolve.results import ReferenceSignature
        from tissueresolve.reference.hierarchy import (
            load_hierarchy_mapping, build_cell_type_hierarchy)
        import warnings
        RBC = Path(__file__).resolve().parents[2] / "examples" / "real_breast_cancer"
        ref = ReferenceSignature.load(RBC / "outputs" / "reference" / "breast_cancer_reference")
        bulk = pd.read_csv(RBC / "data" / "derived" / "pseudobulk_counts.tsv",
                           sep="\t", index_col=0, comment="#"); bulk.index = bulk.index.map(str)
        rg = set(map(str, ref.gene_names)); bulk = bulk.loc[[g for g in bulk.index if g in rg]]
        truth = pd.read_csv(RBC / "data" / "derived" / "pseudobulk_true_proportions.tsv",
                            sep="\t", index_col=0, comment="#")
        mpath = RBC / "config" / "breast_cancer_cell_type_hierarchy.tsv"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mapping = build_cell_type_hierarchy(list(ref.cell_types),
                                                load_hierarchy_mapping(mpath)) if mpath.exists() else None
        sc, sp_truth = None, None
        spatial_h5ad = str(RBC / "data" / "spatial" / "human_breast_cancer_1.h5ad")

    written = {}
    # reference: single-cell counts + per-cell metadata for methods that need
    # individual cells (Bisque, MuSiC).  Exported from the source h5ad,
    # subsampled per cell type and restricted to the reference genes — the
    # subsampling is recorded (never silent).
    sc_info = {}
    if use_existing_real_data and not toy:
        try:
            h5ad = (Path(__file__).resolve().parents[2] / "examples" /
                    "real_breast_cancer" / "data" / "reference" /
                    "breast_cancer_sc_reference.h5ad")
            sc_info = _export_single_cell_reference(
                h5ad, ref_genes=list(ref.gene_names),
                cell_types=list(ref.cell_types), out_dir=ref_dir,
                max_cells_per_type=60, seed=0)
            written["reference_counts_genes_by_cells"] = ref_dir / "reference_counts_genes_by_cells.tsv"
        except Exception as exc:  # noqa: BLE001
            sc_info = {"exported": False, "error": str(exc)}

    # reference signatures
    sig_fine, sig_broad = _signature_frames(ref, mapping)
    written["reference_signature_fine"] = write_tsv(sig_fine, ref_dir / "reference_signature_fine.tsv")
    if sig_broad is not None:
        written["reference_signature_broad"] = write_tsv(sig_broad, ref_dir / "reference_signature_broad.tsv")
    pd.Series(list(ref.gene_names), name="gene_id").to_csv(ref_dir / "gene_ids.tsv", sep="\t", index=False)
    meta = pd.DataFrame({"cell_type": list(ref.cell_types),
                         "n_cells": [ref.n_cells_per_type.get(c, 0) for c in ref.cell_types]})
    written["reference_metadata"] = write_tsv(meta.set_index("cell_type"), ref_dir / "reference_metadata.tsv")
    if mapping:
        from tissueresolve.reference.hierarchy import hierarchy_to_frame
        hierarchy_to_frame(mapping).to_csv(ref_dir / "hierarchy_mapping.tsv", sep="\t", index=False)
        written["hierarchy_mapping"] = ref_dir / "hierarchy_mapping.tsv"

    # bulk (both orientations + normalization variants)
    written["bulk_genes_by_samples"] = write_tsv(bulk, bulk_dir / "bulk_counts_genes_by_samples.tsv")
    write_tsv(bulk.T, bulk_dir / "bulk_counts_samples_by_genes.tsv")
    write_tsv(norm.convert_counts_to_cpm(bulk), bulk_dir / "bulk_cpm.tsv")
    write_tsv(norm.convert_counts_to_log_cpm(bulk), bulk_dir / "bulk_logcpm.tsv")
    if truth is not None:
        written["bulk_truth"] = write_tsv(truth, bulk_dir / "bulk_truth_proportions.tsv")

    # spatial
    (sp_dir / "spatial_h5ad_path.txt").write_text(spatial_h5ad, encoding="utf-8")
    if sc is not None:
        Y = np.asarray(sc["Y"]); spots = sc["spot_ids"]; genes = sc["gene_names"]
        write_tsv(pd.DataFrame(Y, index=spots, columns=genes),
                  sp_dir / "spatial_counts_spots_by_genes.tsv")
        write_tsv(pd.DataFrame({"array_row": sc["array_row"], "array_col": sc["array_col"]},
                               index=spots), sp_dir / "spatial_coordinates.tsv")
        if sp_truth is not None:
            write_tsv(sp_truth, sp_dir / "spatial_truth_if_synthetic.tsv")

    # summaries
    summaries = {
        "n_reference_genes": ref.n_genes, "n_cell_types": ref.n_cell_types,
        "n_families": (len(set(mapping.values())) if mapping else None),
        "bulk_samples": int(bulk.shape[1]),
        "bulk_normalization": norm.detect_normalization_status(bulk.to_numpy()),
        "data_source": "toy" if (toy or not use_existing_real_data) else "real_breast_cancer",
        "spatial_h5ad": spatial_h5ad or None,
        "single_cell_export": sc_info,
    }
    write_tsv(pd.DataFrame([summaries]).T.rename(columns={0: "value"}),
              PREP / "input_summary.tsv")
    write_json(summaries, PREP / "input_summary.json")
    return {"written": {k: str(v) for k, v in written.items()}, "summary": summaries}
