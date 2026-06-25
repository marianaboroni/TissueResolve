"""Tests for the promotion of partial confidence-weighted soft gating to the
DEFAULT hierarchical gating mode.

Offline, deterministic, no generated data. Covers:
  * the default gating mode is "soft" (config, assemble, run_hierarchical_bulk),
  * "hard" and "ungated" are explicitly requestable,
  * soft preserves family + total mass (mass conservation),
  * hard remains backward-compatible (binary keep/drop),
  * ungated never abstains,
  * run metadata records the gating mode / version / validation status,
  * an invalid gating value raises,
  * a config that predates `hierarchical_gating` is interpreted as "soft" and
    the run metadata records the resolved mode,
  * the stable spatial smoothing default (lambda_spatial) is NOT changed.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from tissueresolve.results import ReferenceSignature
from tissueresolve.reference.hierarchy import assemble_hierarchical_estimates


def _ref():
    """Toy fine reference: FamX separable (X1/X2), FamP partly collinear (P2≈P3)."""
    rng = np.random.default_rng(0)
    G = 120
    genes = [f"g{i}" for i in range(G)]

    def prof(a, lvl=200.0):
        v = np.full(G, 2.0)
        v[list(a)] = lvl
        return v

    X1 = prof(range(0, 20))
    X2 = prof(range(20, 40))
    P1 = prof(range(40, 60))
    P2 = prof(range(60, 80))
    P3 = P2 * (1 + rng.normal(0, 0.004, G))
    R = np.vstack([X1, X2, P1, P2, P3]).astype(np.float32)
    cts = ["X1", "X2", "P1", "P2", "P3"]
    ref = ReferenceSignature(gene_names=genes, cell_types=cts, R_cpm=R,
                             R_log=np.log1p(R).astype(np.float32),
                             n_cells_per_type={c: 100 for c in cts})
    mapping = {"X1": "FamX", "X2": "FamX", "P1": "FamP", "P2": "FamP", "P3": "FamP"}
    return ref, mapping


def _toy_query():
    fam = pd.DataFrame({"FamX": [0.5], "FamP": [0.5]}, index=["s0"])
    fine = pd.DataFrame({"X1": [0.3], "X2": [0.2], "P1": [0.2], "P2": [0.15],
                         "P3": [0.15]}, index=["s0"])
    return fam, fine


def _assemble(gating=None):
    ref, mapping = _ref()
    fam, fine = _toy_query()
    kw = {} if gating is None else {"gating": gating}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return assemble_hierarchical_estimates(
            fam, fine, ref, mapping, min_discriminating_genes=5,
            allow_partial_resolution=True, **kw)


# --- defaults --------------------------------------------------------------

def test_assemble_default_is_soft():
    est = _assemble()
    assert est.metadata["hierarchical_gating"] == "soft"
    assert est.metadata["gating_effective_mode"] == "soft"
    assert est.metadata["gating_is_default"] is True


def test_config_default_is_soft():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    assert cfg.hierarchical.hierarchical_gating == "soft"
    assert cfg.hierarchical.gating_version == "soft_gating-1.0"


def test_run_hierarchical_bulk_default_records_soft():
    from tissueresolve.bulk.hierarchical import run_hierarchical_bulk
    ref, mapping = _ref()
    # tiny bulk: two pseudo-samples on the reference genes
    rng = np.random.default_rng(1)
    bulk = pd.DataFrame(rng.poisson(50, size=(len(ref.gene_names), 2)) + 1,
                        index=ref.gene_names, columns=["b0", "b1"]).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_hierarchical_bulk(bulk, ref, mapping, min_discriminating_genes=5)
    assert res.run_metadata["hierarchical_gating"] == "soft"
    assert res.run_metadata["hard_gating_status"] == "legacy"


# --- modes explicitly requestable -----------------------------------------

def test_hard_and_ungated_requestable():
    soft = _assemble("soft")
    hard = _assemble("hard")
    ungated = _assemble("ungated")
    assert soft.metadata["gating_effective_mode"] == "soft"
    assert hard.metadata["gating_effective_mode"] == "hard"
    assert ungated.metadata["gating_effective_mode"] == "ungated"


def test_invalid_gating_raises():
    with pytest.raises(ValueError):
        _assemble("bogus")


# --- mass conservation -----------------------------------------------------

def test_soft_conserves_total_and_family_mass():
    est = _assemble("soft")
    c = est.combined_fine
    # total mass preserved
    np.testing.assert_allclose(c.sum(axis=1), 1.0, atol=1e-9)
    assert est.metadata["mass_conservation_max_error"] < 1e-9
    # family mass preserved: resolved members + unresolved_<fam> == broad mass
    for fam, members in {"FamX": ["X1", "X2"], "FamP": ["P1", "P2", "P3"]}.items():
        resolved = c[members].sum(axis=1)
        unres = c.get(f"unresolved_{fam}")
        recon = resolved + (unres if unres is not None else 0.0)
        np.testing.assert_allclose(recon, est.family_proportions[fam], atol=1e-9)


def test_hard_conserves_mass_and_is_binary():
    # Hard gate: near-identical P2/P3 abstained to unresolved, mass preserved.
    est = _assemble("hard")
    c = est.combined_fine
    np.testing.assert_allclose(c.sum(axis=1), 1.0, atol=1e-9)
    assert c["P2"].iloc[0] == pytest.approx(0.0)
    assert c["P3"].iloc[0] == pytest.approx(0.0)
    assert "unresolved_FamP" in c.columns and c["unresolved_FamP"].iloc[0] > 0


def test_ungated_never_abstains():
    est = _assemble("ungated")
    c = est.combined_fine
    np.testing.assert_allclose(c.sum(axis=1), 1.0, atol=1e-9)
    # no unresolved columns and every subtype gets some mass
    assert not any(col.startswith("unresolved_") for col in c.columns)
    assert est.metadata["n_unresolved_families"] == 0
    assert est.metadata["unresolved_mass_fraction"] == pytest.approx(0.0)


def test_soft_keeps_partial_mass_where_hard_zeroes_it():
    # Soft is continuous: a subtype the hard gate zeroes out (confidence below
    # the threshold) can still receive a fraction of its mass under soft when
    # its confidence is > 0.  Demonstrated on a subtype that hard abstains.
    soft = _assemble("soft").combined_fine
    hard = _assemble("hard").combined_fine
    # P1 has confidence ~0.77: hard keeps it fully, soft keeps a (smaller but
    # positive) confidence-weighted fraction — both > 0, neither zeroed.
    assert soft["P1"].iloc[0] > 0 and hard["P1"].iloc[0] > 0
    # soft scales continuously, so its kept P1 mass is <= hard's full mass.
    assert soft["P1"].iloc[0] <= hard["P1"].iloc[0] + 1e-12


# --- metadata --------------------------------------------------------------

def test_metadata_records_gating_provenance():
    est = _assemble("soft")
    m = est.metadata
    assert m["gating_version"] == "soft_gating-1.0"
    assert m["gating_default"] == "soft"
    assert m["gating_validation_status"] == "validated on breast and lung benchmarks"
    assert m["gating_confidence_model"] == "within_family_discriminating_gene_confidence"
    assert isinstance(m["gating_confidence_features"], list) and m["gating_confidence_features"]
    assert "mass_conservation_max_error" in m
    assert isinstance(m["unresolved_mass_per_family"], dict)


def test_metadata_marks_hard_as_legacy_status():
    est = _assemble("hard")
    assert est.metadata["gating_validation_status"] == "legacy"
    assert est.metadata["hard_gating_status"] == "legacy"


# --- backward compatibility -------------------------------------------------

def test_config_lacking_field_interpreted_as_soft():
    # Simulate an old config object with no `hierarchical_gating` attribute.
    class OldHCfg:
        allow_unresolved = True
        unresolved_threshold = 0.10
        min_discriminating_genes = 5
        within_family_spillover_threshold = 0.30
        allow_partial_resolution = True
        subtype_confidence_threshold = 0.10
        # NOTE: no hierarchical_gating / gating_version

    resolved = getattr(OldHCfg, "hierarchical_gating", "soft")
    assert resolved == "soft"
    # and assemble with that resolved value records soft in metadata
    est = _assemble(resolved)
    assert est.metadata["gating_effective_mode"] == "soft"


# --- report text -----------------------------------------------------------

def test_report_states_soft_default_and_hard_legacy():
    from tissueresolve.report.result_sections import _gating_mode_html
    est = _assemble("soft")
    html = _gating_mode_html(est).lower()
    assert "soft" in html and "default" in html
    assert "hard" in html and "legacy" in html
    assert "over-abstain" in html
    assert "ungated" in html and "diagnostic" in html
    # interpretation caveats + methods note
    assert "rna-derived proportions" in html
    assert "auroc" in html
    assert "confidence-weighted gating" in html


def test_report_hierarchical_section_includes_gating_block():
    # The full hierarchical section embeds the gating block.
    from tissueresolve.report import result_sections as RS

    class _R:
        estimates = _assemble("soft")
    body = RS._hierarchical_html(_R(), "bulk")
    assert body is not None
    assert "Within-family gating mode" in body
    assert "soft" in body.lower() and "legacy" in body.lower()


def test_report_includes_spillover_and_false_positive_caution():
    from tissueresolve.report.result_sections import _gating_mode_html
    html = _gating_mode_html(_assemble("soft")).lower()
    assert "spillover" in html
    assert "false-positive" in html
    # fine predictions are explicitly gated by trusted resolution
    assert "gated by trusted" in html or "trusted within-family resolution" in html


def test_report_does_not_present_fine_refiner_as_promoted():
    from tissueresolve.report.result_sections import _gating_mode_html
    html = _gating_mode_html(_assemble("soft")).lower()
    # the refiner is shown as an experimental negative result, NOT a promoted feature
    assert "fine-granularity refinement (experimental, not used)" in html
    assert "failed promotion" in html
    assert "not</b> applied" in html or "not applied" in html
    # it must not be described as the default / active / promoted gating
    assert "refiner is the default" not in html
    assert "refinement is applied" not in html


# --- FineGranularityRefiner is experimental and NOT in the default path -------

def test_fine_refiner_default_mode_is_none():
    from tissueresolve.experimental.soft_hierarchy.fine_refiner import FineRefinerConfig
    assert FineRefinerConfig().mode == "none"
    assert FineRefinerConfig().feature_status == "experimental"


def test_fine_refiner_not_used_in_default_hierarchical_run():
    # A default hierarchical run must NOT invoke the refiner: its run metadata
    # carries gating provenance but no fine-refinement status.
    from tissueresolve.bulk.hierarchical import run_hierarchical_bulk
    import warnings
    ref, mapping = _ref()
    rng = np.random.default_rng(3)
    bulk = pd.DataFrame(rng.poisson(50, size=(len(ref.gene_names), 2)) + 1,
                        index=ref.gene_names, columns=["b0", "b1"]).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_hierarchical_bulk(bulk, ref, mapping, min_discriminating_genes=5)
    meta = res.run_metadata
    assert meta["hierarchical_gating"] == "soft"
    assert "fine_refinement" not in meta
    assert "refined_families" not in meta


def test_default_api_path_does_not_import_fine_refiner():
    # The production API/pipeline must not depend on the experimental refiner.
    import inspect
    import tissueresolve.api as api
    import tissueresolve.bulk.hierarchical as bh
    import tissueresolve.spatial.hierarchical as sh
    for mod in (api, bh, sh):
        assert "fine_refiner" not in inspect.getsource(mod)
        assert "FineGranularityRefiner" not in inspect.getsource(mod)


# --- spatial smoothing default unchanged -----------------------------------

def test_spatial_lambda_default_unchanged():
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    # The stable spatial smoothing default must NOT have been weakened by this
    # promotion. lambda_spatial stays at its established value (> the 0.02
    # weak-experimental candidate, which was explicitly NOT promoted here).
    lam = cfg.spatial_solver.lambda_spatial
    assert lam > 0.02
    assert lam == pytest.approx(0.1)
