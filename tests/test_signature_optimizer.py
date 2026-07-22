"""Offline tests for the multi-panel ReferenceSignatureOptimizer (experimental, opt-in)."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


@pytest.fixture
def toy_adata():
    """3 fine types (A,B in fam1; C in fam2), 6 donors, per-type marker blocks."""
    rng = np.random.default_rng(0)
    G, types, donors = 60, ["A", "B", "C"], [f"d{i}" for i in range(6)]
    genes = [f"g{i}" for i in range(G)]
    blocks = {"A": range(0, 10), "B": range(10, 20), "C": range(20, 30)}
    X, ct, dn = [], [], []
    for d in donors:
        for t in types:
            for _ in range(25):
                base = rng.poisson(2.0, size=G).astype(float)
                base[list(blocks[t])] += rng.poisson(40.0, size=10)
                X.append(base); ct.append(t); dn.append(d)
    ad = anndata.AnnData(np.vstack(X).astype("float32"))
    ad.var_names = genes
    ad.obs["cell_type"] = pd.Categorical(ct)
    ad.obs["donor_id"] = pd.Categorical(dn)
    return ad, blocks


MAPPING = {"A": "fam1", "B": "fam1", "C": "fam2"}
KW = dict(mapping=MAPPING, minimal_sizes=(5, 10, 15, 20), seed=0)


def _opt(**extra):
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    return ReferenceSignatureOptimizer("cell_type", "donor_id", **{**KW, **extra})


def test_four_panels_present(toy_adata):
    ad, blocks = toy_adata
    m = _opt().optimize(ad)
    assert m.broad_genes and set(m.broad_genes) <= set(ad.var_names)
    assert m.fine_global_genes                            # primary fine panel
    assert "fam1" in m.sibling_genes                      # 2 fine members -> sibling panel
    assert "fam2" not in m.sibling_genes                  # 1 member -> none
    assert set(m.panels) == {"broad", "fine_global", "sibling", "rare_confirmation"}
    # fine-global recovers markers of both fine types (fine-vs-all-rest)
    allmark = {f"g{i}" for b in blocks.values() for i in b}
    assert sum(g in allmark for g in m.fine_global_genes) >= 15


def test_panel_specs_have_intended_use(toy_adata):
    ad, _ = toy_adata
    m = _opt().optimize(ad)
    assert m.panels["fine_global"].intended_use == "primary_fine_deconvolution"
    assert m.panels["broad"].contrast == "broad_family_vs_other_families"
    assert m.panels["sibling"].intended_use == "resolution_decisions_and_gating"


def test_panel_for_returns_correct_panel(toy_adata):
    ad, _ = toy_adata
    m = _opt().optimize(ad)
    assert m.panel_for("fine_deconvolution") == list(m.fine_global_genes)
    assert m.panel_for("broad_deconvolution") == list(m.broad_genes)
    with pytest.raises(ValueError):
        m.panel_for("nonsense")


def test_misuse_warning():
    from tissueresolve.reference.signature_optimizer import warn_if_panel_misused
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        warn_if_panel_misused("broad", "fine_deconvolution")   # wrong panel for task
    assert any("validated panel is 'fine_global'" in str(x.message) for x in w)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        warn_if_panel_misused("fine_global", "fine_deconvolution")  # correct -> no warn
    assert not w


def test_panel_role_safeguard_strict_and_permissive(toy_adata):
    """§3.5: role safeguard integrated at a real inference entry (not just a helper)."""
    from tissueresolve.reference.signature_optimizer import (
        validate_panel_use, PanelRoleError)
    ad, _ = toy_adata
    m = _opt().optimize(ad)
    broad = m.as_panel("broad")            # role = broad_family_deconvolution
    fine = m.as_panel("fine_global")       # role = primary_fine_deconvolution
    # strict: using broad panel for fine deconvolution must raise
    with pytest.raises(PanelRoleError):
        validate_panel_use(broad, "primary_fine_deconvolution", mode="strict")
    # permissive: warns + returns issues (no raise)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        issues = validate_panel_use(broad, "primary_fine_deconvolution", mode="permissive")
    assert issues and any("misuse" in str(x.message) for x in w)
    # correct role: no issue
    assert validate_panel_use(fine, "primary_fine_deconvolution", mode="strict") == []
    # spatial_unvalidated panel used for spatial: strict raises
    with pytest.raises(PanelRoleError):
        validate_panel_use(fine, "primary_fine_deconvolution", mode="strict", modality="spatial")


def test_deconvolve_with_signature_panel_records_warnings(toy_adata):
    import numpy as np
    from tissueresolve.reference.signature_optimizer import deconvolve_with_signature_panel
    import tissueresolve as tr
    ad, _ = toy_adata
    ref = tr.build_reference(ad, cell_type_col="cell_type")
    m = _opt().optimize(ad)
    rng = np.random.default_rng(4)
    true = rng.dirichlet(np.ones(len(ref.cell_types)), size=3)
    bulk = pd.DataFrame((true @ ref.as_R_cpm()).T, index=list(ref.gene_names),
                        columns=[f"s{i}" for i in range(3)])
    # correct role -> no warnings recorded
    res = deconvolve_with_signature_panel(bulk, ref, m.as_panel("fine_global"),
                                          expected_role="primary_fine_deconvolution", mode="permissive")
    np.testing.assert_allclose(res.proportions.sum(axis=1), 1.0, atol=1e-6)
    assert res.diagnostics["panel_use_warnings"] == []
    # wrong role (broad panel for fine) -> warning recorded in diagnostics, still runs
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res2 = deconvolve_with_signature_panel(bulk, ref, m.as_panel("broad"),
                                               expected_role="primary_fine_deconvolution", mode="permissive")
    assert res2.diagnostics["panel_use_warnings"]


def test_gene_budget_caps_panels(toy_adata):
    ad, _ = toy_adata
    m = _opt(gene_budget=12).optimize(ad)
    assert len(m.broad_genes) <= 12
    assert len(m.fine_global_genes) <= 12
    assert m.panels["fine_global"].gene_budget == 12


def test_optimize_does_not_mutate_donor_config(toy_adata):
    """§3.1: running on a no-donor reference must NOT permanently disable donor-awareness."""
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    ad, _ = toy_adata
    ad_nodonor = ad.copy(); del ad_nodonor.obs["donor_id"]
    opt = ReferenceSignatureOptimizer("cell_type", "donor_id", **KW)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m_nodonor = opt.optimize(ad_nodonor)          # no donor -> capped
        m_donor = opt.optimize(ad)                    # SAME object, donor present again
    assert m_nodonor.status == "EXPERIMENTAL_FINE"
    assert opt.dc == "donor_id"                        # config restored, not mutated
    assert m_donor.status in ("PASS", "PASS_WITH_RESTRICTIONS")
    assert m_donor.fine_global_genes                   # donor-aware again


def test_stratified_budget_covers_all_types(toy_adata):
    """§3.3: budgeted fine-global must not starve later cell types (no flat truncation)."""
    ad, _ = toy_adata
    m = _opt(gene_budget=12).optimize(ad)
    assert len(m.fine_global_genes) <= 12
    short = m.metadata["fine_global_shortfall"]
    assert short and all(r["selected_genes"] > 0 for r in short)   # every eligible type covered
    assert m.metadata["budget_strategy"] == "stratified_round_robin"
    # natural (no budget) records the other strategy
    assert _opt().optimize(ad).metadata["budget_strategy"] == "natural_per_type"


def test_batch_aware_fallback(toy_adata):
    """No donor but batch present → batch-aware DE, reduced confidence, capped status."""
    import warnings
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    ad, _ = toy_adata
    a = ad.copy(); del a.obs["donor_id"]
    a.obs["batch"] = pd.Categorical(["bx" if i % 2 else "by" for i in range(a.n_obs)])
    opt = ReferenceSignatureOptimizer("cell_type", "donor_id", mapping=MAPPING,
                                      batch_col="batch", minimal_sizes=(5, 10))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        m = opt.optimize(a)
    prov = m.metadata["selection_provenance"]
    assert prov["actual_method"] == "batch_aware_pseudobulk_de" and prov["batch_aware"]
    assert prov["confidence_level"] == "reduced"
    assert m.status == "EXPERIMENTAL_FINE"                 # non-donor never claims PASS
    assert m.fine_global_genes                              # not empty (fallback produced a panel)
    assert any("batch" in str(x.message).lower() for x in w)


def test_pooled_cell_fallback_critical(toy_adata):
    """No donor and no batch → pooled-cell DE with CRITICAL warning + low confidence."""
    import warnings
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    ad, _ = toy_adata
    a = ad.copy(); del a.obs["donor_id"]
    opt = ReferenceSignatureOptimizer("cell_type", "donor_id", mapping=MAPPING,
                                      minimal_sizes=(5, 10))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        m = opt.optimize(a)
    prov = m.metadata["selection_provenance"]
    assert prov["actual_method"] == "pooled_cell_de" and prov["confidence_level"] == "low"
    assert "CRITICAL" in prov["limitations"]
    assert any("CRITICAL" in str(x.message) for x in w)
    assert m.status == "EXPERIMENTAL_FINE"


def test_single_subtype_status(toy_adata):
    ad, _ = toy_adata
    a = ad[ad.obs["cell_type"].astype(str) == "A"].copy()   # only one cell type
    m = _opt(mapping={"A": "fam1"}).optimize(a)
    assert m.status == "SINGLE_SUBTYPE"


def test_reference_inadequate_status(toy_adata):
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    ad, _ = toy_adata
    a = ad[:8].copy(); del a.obs["donor_id"]                # too few cells for any DE
    m = ReferenceSignatureOptimizer("cell_type", "donor_id", mapping=MAPPING,
                                    min_cells=10, minimal_sizes=(5,)).optimize(a)
    assert m.status in ("REFERENCE_INADEQUATE", "SINGLE_SUBTYPE")
    assert m.metadata["selection_provenance"]["actual_method"] in ("none", "pooled_cell_de")


def test_spatial_status_labels(toy_adata):
    ad, _ = toy_adata
    m = _opt().optimize(ad)
    assert m.panels["fine_global"].validated_modality == "bulk_experimental"
    assert m.panels["fine_global"].spatial_status == "spatial_unvalidated"


def test_rare_confirmation_combines_global_specificity(toy_adata):
    """§3.4: rare-confirmation genes should be members of the fine-global marker set
    when overlap exists (combined evidence, not sibling-only)."""
    ad, _ = toy_adata
    m = _opt().optimize(ad)
    fine_set = set(m.fine_global_genes)
    overlaps = [any(g in fine_set for g in gl) for gl in m.rare_confirmation.values()]
    assert m.rare_confirmation and all(overlaps)


def test_rare_confirmation_is_specific(toy_adata):
    ad, blocks = toy_adata
    m = _opt().optimize(ad)
    for sub in ("A", "B"):
        if sub in m.rare_confirmation:
            expected = {f"g{i}" for i in blocks[sub]}
            assert sum(g in expected for g in m.rare_confirmation[sub]) >= 2


def test_status_and_no_donor(toy_adata):
    from tissueresolve.reference.signature_optimizer import ReferenceSignatureOptimizer
    ad, _ = toy_adata
    m = _opt().optimize(ad)
    assert m.status in ("PASS", "PASS_WITH_RESTRICTIONS")
    assert m.family_status["fam2"] == "BROAD_ONLY"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        m2 = ReferenceSignatureOptimizer("cell_type", None, mapping=MAPPING,
                                         minimal_sizes=(5, 10)).optimize(ad)
    assert m2.status == "EXPERIMENTAL_FINE"
    assert any("donor" in str(x.message).lower() for x in w)


def test_deterministic_and_save(tmp_path, toy_adata):
    ad, _ = toy_adata
    a = _opt().optimize(ad); b = _opt().optimize(ad)
    assert a.fine_global_genes == b.fine_global_genes
    assert a.broad_genes == b.broad_genes and a.sibling_genes == b.sibling_genes
    a.save(tmp_path / "sig")
    for f in ("broad_signature.tsv", "fine_global_signature.tsv", "sibling_signature.tsv",
              "rare_confirmation_signature.tsv", "signature_manifest.json"):
        assert (tmp_path / "sig" / f).exists(), f"missing {f}"
    import json
    meta = json.loads((tmp_path / "sig" / "signature_manifest.json").read_text())
    assert meta["panels"]["fine_global"]["intended_use"] == "primary_fine_deconvolution"
