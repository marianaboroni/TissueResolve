"""Standard-report integration for the identifiability-aware diagnostics.

Emits, as first-class report tables (best-effort; never breaks a run):
  * ``adaptive_resolution.tsv``   — per-family best supported resolution (P4)
  * ``reference_uncertainty.tsv`` / ``state_reliability.tsv`` /
    ``family_reliability.tsv``    — reference reliability (P2c)

These are **reporting / QC only** — they do not change any estimate. They are
written whenever the necessary inputs are available (a fine→broad mapping for
adaptive resolution; the raw reference AnnData for reference uncertainty); otherwise
they are silently skipped.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Mapping, Optional

import pandas as pd

logger = logging.getLogger("tissueresolve.report.resolution_diagnostics")

__all__ = ["write_resolution_diagnostics", "mean_unresolved_by_family"]


def mean_unresolved_by_family(proportions: pd.DataFrame) -> dict:
    """Mean mass of each ``unresolved_<family>`` column → ``{family: mass}``."""
    out = {}
    for c in proportions.columns:
        cs = str(c)
        if cs.startswith("unresolved_"):
            out[cs[len("unresolved_"):]] = float(proportions[c].mean())
    return out


def write_resolution_diagnostics(
    out_dir,
    ref,
    mapping: Optional[Mapping[str, str]] = None,
    *,
    reference_path: Optional[str] = None,
    cell_type_col: str = "cell_type",
    donor_col: Optional[str] = "donor",
    unresolved_mass: Optional[Mapping[str, float]] = None,
    subtype_confidence: Optional[Mapping[str, float]] = None,
    rare_protection=None,
) -> dict:
    """Write the adaptive-resolution + reference-uncertainty tables under *out_dir*.

    Best-effort and fully defensive: any failure logs a debug note and is skipped so
    a run never breaks. Returns ``{name: path}`` for what was written.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict = {}

    # identity mapping fallback (flat runs) so the table still lists the types
    eff_mapping = dict(mapping) if mapping else {str(c): str(c) for c in ref.cell_types}

    # --- P4: adaptive resolution (needs ref + mapping only) ---
    try:
        from tissueresolve.experimental.adaptive_resolution import (
            build_adaptive_resolution_report)
        rep = build_adaptive_resolution_report(
            ref, eff_mapping, subtype_confidence=subtype_confidence,
            unresolved_mass=unresolved_mass, rare_protection=rare_protection)
        written.update(rep.write(out))
        logger.info("Wrote adaptive_resolution.tsv (%d families).", len(rep.families))
    except Exception as exc:  # noqa: BLE001
        logger.debug("adaptive_resolution skipped: %s", exc)

    # --- P2c: reference uncertainty (needs the raw reference AnnData) ---
    try:
        if reference_path and str(reference_path).endswith(".h5ad"):
            import anndata as ad
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                adata = ad.read_h5ad(reference_path)
            obs_cols = set(adata.obs.columns)
            ctc = cell_type_col if cell_type_col in obs_cols else next(
                (c for c in ("cell_type", "celltype", "cell_type_fine") if c in obs_cols), None)
            dcol = donor_col if (donor_col and donor_col in obs_cols) else next(
                (c for c in ("donor", "donor_id", "Subject", "SubjectName") if c in obs_cols), None)
            if ctc is not None:
                from tissueresolve.experimental.reference_uncertainty import (
                    estimate_reference_uncertainty)
                ru = estimate_reference_uncertainty(
                    adata, ctc, donor_col=dcol, genes=list(ref.gene_names),
                    mapping=eff_mapping)
                written.update(ru.write(out))
                logger.info("Wrote reference_uncertainty/state/family reliability "
                            "(%d states, donor=%s).", ru.metadata["n_states"], dcol)
    except Exception as exc:  # noqa: BLE001
        logger.debug("reference_uncertainty skipped: %s", exc)

    return written
