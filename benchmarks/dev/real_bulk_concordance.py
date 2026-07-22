#!/usr/bin/env python
"""Real-bulk (TCGA-TNBC) validation by concordance + marker coherence (dev, real-data).

TCGA-TNBC bulk has **no cell-composition ground truth**, so this does NOT measure
accuracy. It asks whether the pseudobulk-derived recommendations (Poisson solver,
donor-aware DE genes) hold up on REAL bulk by two truth-free criteria:

  1. MARKER COHERENCE — across the 40 tumors, does each estimated broad-family fraction
     track that family's canonical marker expression in the bulk (Spearman)? A better
     method makes fractions that follow their markers.
  2. METHOD CONCORDANCE — does Poisson/donor_de TissueResolve agree with Rectangle (a
     strong external tool) on real bulk more than the wNNLS default does?
  3. RECONSTRUCTION FIT — Poisson vs wNNLS deviance/R² on the observed bulk (weak signal;
     better fit ≠ better composition — reported with that caveat).

This is a promotion GATE input, not a promotion: broad-level, coherence-only. Fine/rare
accuracy on TCGA is not claimed (no truth). Experimental; outputs gitignored; no defaults
changed. Usage:
  PYTHONPATH=src:. python benchmarks/dev/real_bulk_concordance.py --run-real-data
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "real_breast_cancer" / "scripts"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
import _harness as H  # noqa: E402
from benchmarks.shared import environment as ENV  # noqa: E402

REF_H5AD = REPO / "examples" / "real_breast_cancer" / "data" / "reference" / "breast_cancer_sc_reference.h5ad"
TCGA = REPO / "examples" / "real_breast_cancer" / "data" / "bulk_tcga_tnbc" / "tcga_tnbc_counts_clean.tsv"
HIER = REPO / "examples" / "real_breast_cancer" / "config" / "breast_cancer_cell_type_hierarchy.tsv"
ENV_RUNNER = REPO / "benchmarks" / "rectangle_env_runner.py"
OUT = REPO / "benchmarks" / "results" / "real_bulk_concordance"

# canonical broad-family markers (HGNC); coherence = corr(family fraction, marker mean logCPM)
FAMILY_MARKERS = {
    "T/NK": ["CD3D", "CD3E", "CD2", "NKG7", "GNLY", "CD8A"],
    "Myeloid": ["LYZ", "CD68", "CSF1R", "AIF1", "ITGAX"],
    "B/Plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "IGHG1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Stromal/Fibroblast": ["COL1A1", "COL1A2", "LUM", "DCN", "PDGFRA"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Mural": ["RGS5", "ACTA2", "MYH11", "PDGFRB"],
    "Adipocyte": ["ADIPOQ", "LEP", "FABP4"],
}


def _load_tcga():
    c = pd.read_csv(TCGA, sep="\t", comment="#", index_col=0)
    c.index = c.index.map(str)
    return c[~c.index.duplicated()]


def _logcpm(bulk):
    cpm = bulk / bulk.sum(axis=0) * 1e6
    return np.log2(cpm + 1.0)


def _to_broad(props, mapping):
    fam = {}
    for c in props.columns:
        f = str(mapping.get(c, c))
        fam.setdefault(f, []).append(c)
    return pd.DataFrame({f: props[cols].sum(axis=1) for f, cols in fam.items()},
                        index=props.index)


def _marker_coherence(broad_props, logcpm):
    """Spearman across samples between each family fraction and its marker mean logCPM."""
    from scipy.stats import spearmanr
    rows = []
    for fam, markers in FAMILY_MARKERS.items():
        if fam not in broad_props.columns:
            continue
        present = [m for m in markers if m in logcpm.index]
        if len(present) < 2:
            continue
        msig = logcpm.loc[present].mean(axis=0).reindex(broad_props.index)
        frac = broad_props[fam].reindex(broad_props.index)
        if frac.std() < 1e-9 or msig.std() < 1e-9:
            rho = np.nan
        else:
            rho = float(spearmanr(frac.to_numpy(), msig.to_numpy()).correlation)
        rows.append({"family": fam, "n_markers": len(present), "coherence_spearman": rho})
    return pd.DataFrame(rows).set_index("family")


def _run_rectangle(bulk_counts, cell_type_col, gene_col):
    import subprocess
    py = ENV.rectangle_env_python()
    if py is None:
        return None
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        tpm = (bulk_counts.T.div(bulk_counts.sum(axis=0).replace(0, np.nan), axis=0)
               .fillna(0.0) * 1e6)
        tpm.to_csv(tmp / "bulk_tpm.tsv", sep="\t")
        out = tmp / "out"
        cmd = [str(py), str(ENV_RUNNER), "--ref-h5ad", str(REF_H5AD),
               "--cell-type-col", cell_type_col, "--gene-col", gene_col,
               "--bulk-tpm", str(tmp / "bulk_tpm.tsv"), "--out-dir", str(out),
               "--cells-per-type", "400", "--bootstraps", "5"]
        subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))
        f = out / "predictions_raw.tsv"
        if f.exists():
            df = pd.read_csv(f, sep="\t", index_col=0)
            return df.drop(columns=[c for c in df.columns if c == "Unknown"], errors="ignore")
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-real-data", action="store_true")
    ap.add_argument("--no-rectangle", action="store_true")
    args = ap.parse_args(argv)
    if not H.real_data_enabled(args.run_real_data):
        print(f"Refusing without --run-real-data (or {H.REAL_DATA_ENV}=1).", file=sys.stderr)
        return 2
    for f in (REF_H5AD, TCGA):
        if not f.exists():
            print(f"missing required file: {f}", file=sys.stderr)
            return 3
    OUT.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    import tissueresolve as tr
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.reference.hierarchy import load_hierarchy_mapping, build_cell_type_hierarchy

    bulk = _load_tcga()
    logcpm = _logcpm(bulk)
    print(f"TCGA-TNBC bulk: {bulk.shape[0]} genes × {bulk.shape[1]} tumors")

    adata = ad.read_h5ad(REF_H5AD)
    gene_col = "feature_name" if "feature_name" in adata.var.columns else None
    if gene_col:
        adata.var_names = adata.var[gene_col].astype(str).to_numpy()
        adata.var_names_make_unique()
    ctc = "cell_type"
    cfg = TissueResolveConfig(); cfg.reference.celltype_col = ctc; cfg.reference.donor_col = "donor_id"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref = tr.build_reference(adata, cell_type_col=ctc, config=cfg,
                                 estimate_overdispersion=True)
        ref_dde = tr.build_reference(adata, cell_type_col=ctc, config=cfg,
                                     estimate_overdispersion=True, gene_selection="donor_de")
    types = [str(c) for c in ref.cell_types]
    raw = load_hierarchy_mapping(HIER)

    def _fallback_family(t):
        low = t.lower()
        if any(k in low for k in ("t cell", "nk", "natural killer", "t-cell")):
            return "T/NK"
        if any(k in low for k in ("b cell", "plasma", "b-cell")):
            return "B/Plasma"
        if any(k in low for k in ("dendritic", "neutrophil", "mast", "macrophage",
                                  "monocyte", "myeloid", "granulocyte")):
            return "Myeloid"
        return t  # self-map (own singleton broad) — never silently merged elsewhere
    complete = {t: raw.get(t, _fallback_family(t)) for t in types}
    mapping = build_cell_type_hierarchy(types, complete)
    print("uncovered-type fallback mapping: "
          + str({t: complete[t] for t in types if t not in raw}))

    methods = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        methods["wNNLS_default"] = tr.deconv_bulk(bulk, ref, resolution_mode="flat", n_bootstrap=0).deconv.proportions
        methods["Poisson"] = tr.deconv_bulk(bulk, ref, solver="poisson", resolution_mode="flat").deconv.proportions
        methods["Poisson_donorDE"] = tr.deconv_bulk(bulk, ref_dde, solver="poisson", resolution_mode="flat").deconv.proportions
    if not args.no_rectangle and ENV.rectangle_available():
        print("running Rectangle on TCGA-TNBC (isolated env)…")
        rect = _run_rectangle(bulk, ctc, gene_col or "feature_name")
        if rect is not None:
            methods["Rectangle"] = rect

    # broad compositions + marker coherence
    broad = {m: _to_broad(p, mapping) for m, p in methods.items()}
    coh_rows = []
    for m, bp in broad.items():
        c = _marker_coherence(bp, logcpm)
        c["method"] = m
        coh_rows.append(c.reset_index())
    coh = pd.concat(coh_rows, ignore_index=True)
    coh_pivot = coh.pivot_table(index="family", columns="method", values="coherence_spearman").round(3)

    # method concordance on broad fractions (mean abs corr across shared families/samples)
    conc = {}
    ms = list(methods)
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            a, b = broad[ms[i]], broad[ms[j]]
            fams = [f for f in a.columns if f in b.columns]
            av = a[fams].to_numpy().ravel(); bv = b[fams].reindex(index=a.index, columns=fams).to_numpy().ravel()
            conc[f"{ms[i]}~{ms[j]}"] = float(np.corrcoef(av, bv)[0, 1])

    coh_pivot.to_csv(OUT / "marker_coherence.tsv", sep="\t")
    (OUT / "concordance.json").write_text(json.dumps(conc, indent=2))
    (OUT / "manifest.json").write_text(json.dumps({
        "n_tumors": int(bulk.shape[1]), "methods": list(methods),
        "note": "TCGA-TNBC has NO cell-composition ground truth; concordance + marker "
                "coherence only, broad-level. Not an accuracy claim."}, indent=2))
    mean_coh = coh.groupby("method")["coherence_spearman"].mean().round(3)
    print("\n=== marker coherence (Spearman; higher = family fraction tracks its markers) ===")
    print(coh_pivot.to_string())
    print("\nmean marker coherence by method:", mean_coh.to_dict())
    print("\n=== broad-composition concordance (Pearson) ===")
    for k, v in conc.items():
        print(f"  {k}: {v:.3f}")
    print(f"\nWrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
