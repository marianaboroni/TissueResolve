"""Level-specific spatial smoothing (Task D, experimental).

Broad and fine resolution levels get DIFFERENT CAR strengths (λ_broad, λ_fine) by
running the EXISTING spatial solver twice — once on a family-merged reference at
λ_broad, once on the fine reference at λ_fine — then combining:

    broad_f       = deconv_spatial(family_ref,  λ=λ_broad)          # spots × families
    fine_k        = deconv_spatial(fine_ref,     λ=λ_fine)           # spots × subtypes
    conditional_k = fine_k / Σ_{k∈f} fine_k                          # within-family share
    combined_k    = broad_{f(k)} · conditional_k                     # mass-consistent

No new model, no NB-CAR VI, no PyTorch/JAX (rules). The default pipeline is
untouched; this is opt-in. λ=0 at a level ⇒ pure NB-MAP (no smoothing) there.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

ALGORITHM_VERSION = "level_specific_spatial-0.1.0"
FEATURE_STATUS = "experimental"


@dataclass
class LevelSpecificSpatialResult:
    broad_proportions: pd.DataFrame          # spots × families
    fine_proportions: pd.DataFrame           # spots × subtypes (combined = broad × conditional)
    conditional_proportions: pd.DataFrame    # spots × subtypes (within-family share)
    lambda_broad: float
    lambda_fine: float
    runtime: float
    metadata: dict = field(default_factory=dict)


def _config(lam):
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    cfg.spatial_solver.lambda_spatial = float(lam)
    return cfg


def fit_level_specific_spatial(Y, ref, array_row, array_col, lib_sizes, gene_names,
                               mapping, *, lambda_broad=0.1, lambda_fine=0.02,
                               spot_ids=None) -> LevelSpecificSpatialResult:
    """Run the existing NB-CAR solver at λ_broad (family level) and λ_fine (subtype
    level) and combine into mass-consistent fine proportions."""
    import time
    import warnings
    import tissueresolve as tr
    from tissueresolve.reference.hierarchy import merge_reference_cell_types
    t0 = time.perf_counter()
    family_ref = merge_reference_cell_types(ref, dict(mapping))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        broad = tr.deconv_spatial(Y, family_ref, array_row, array_col, lib_sizes, gene_names,
                                  spot_ids=spot_ids, resolution_mode="none",
                                  config=_config(lambda_broad), run_neighbourhood=False).deconv.proportions
        fine = tr.deconv_spatial(Y, ref, array_row, array_col, lib_sizes, gene_names,
                                 spot_ids=spot_ids, resolution_mode="none",
                                 config=_config(lambda_fine), run_neighbourhood=False).deconv.proportions
    subtypes = [str(c) for c in fine.columns]
    fam_of = {c: str(mapping.get(c, c)) for c in subtypes}
    families = [str(c) for c in broad.columns]
    cond = fine.copy().astype(float)
    combined = fine.copy().astype(float) * 0.0
    for fam in families:
        mem = [c for c in subtypes if fam_of[c] == fam]
        if not mem:
            continue
        s = fine[mem].sum(axis=1)
        share = fine[mem].div(s.replace(0, np.nan), axis=0).fillna(1.0 / len(mem))
        for c in mem:
            cond[c] = share[c]
            combined[c] = broad[fam] * share[c]
    return LevelSpecificSpatialResult(
        broad_proportions=broad, fine_proportions=combined, conditional_proportions=cond,
        lambda_broad=float(lambda_broad), lambda_fine=float(lambda_fine),
        runtime=time.perf_counter() - t0,
        metadata={"algorithm_version": ALGORITHM_VERSION, "feature_status": FEATURE_STATUS,
                  "n_families": len(families), "n_subtypes": len(subtypes),
                  "estimate_type": "spot_rna_composition"})
