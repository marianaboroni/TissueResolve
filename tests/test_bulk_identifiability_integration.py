"""Stage 2C — the calibrated certificate wired into the bulk pipeline as an opt-in reported output.

Guarantees: (a) opt-in attaches a certificate; (b) estimates are byte-identical with/without it;
(c) default is OFF (backward compatible); (d) recommended_merge guidance for confounded clusters;
(e) a warning when the reference lacks donor_cv.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import tissueresolve as tr
from tissueresolve.config import ReferenceConfig
from tissueresolve.reference.build import ReferenceBuilder


def _toy_reference(with_donor=True, seed=0):
    """4 types over 60 genes; A and C share an identical block (confounded); 3 donors."""
    rng = np.random.default_rng(seed)
    import anndata as adann
    blocks = {"A": range(0, 15), "B": range(15, 30), "C": range(0, 15), "D": range(30, 45)}
    G, per = 60, 20
    X, obs_ct, obs_donor = [], [], []
    for donor in ("d1", "d2", "d3"):
        for ct, blk in blocks.items():
            base = np.full(G, 1.0); base[list(blk)] = 30.0
            for _ in range(per):
                X.append(rng.poisson(base)); obs_ct.append(ct); obs_donor.append(donor)
    ad = adann.AnnData(np.array(X, dtype=float),
                       obs=pd.DataFrame({"cell_type": obs_ct, "donor_id": obs_donor}),
                       var=pd.DataFrame(index=[f"g{i}" for i in range(G)]))
    cfg = ReferenceConfig(celltype_col="cell_type", min_cells=5,
                          donor_col=("donor_id" if with_donor else None))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ReferenceBuilder(cfg).build_from_adata(ad)


def _toy_bulk(ref, seed=1):
    rng = np.random.default_rng(seed)
    R = ref.as_R_cpm()                                  # (K, G)
    genes = list(ref.gene_names)
    cols = {}
    for s in range(3):
        theta = rng.dirichlet(np.ones(R.shape[0]))
        mu = 5e5 * (theta @ R) / 1e6
        cols[f"s{s}"] = rng.poisson(np.clip(mu, 0, None))
    return pd.DataFrame(cols, index=genes)


def test_default_off_backward_compatible():
    ref = _toy_reference(); bulk = _toy_bulk(ref)
    res = tr.deconv_bulk(bulk, ref, resolution_mode="none")
    assert getattr(res, "identifiability", None) is None


def test_optin_attaches_certificate_without_changing_estimates():
    ref = _toy_reference(); bulk = _toy_bulk(ref)
    base = tr.deconv_bulk(bulk, ref, resolution_mode="none")
    with_id = tr.deconv_bulk(bulk, ref, resolution_mode="none", identifiability=True)
    assert with_id.identifiability is not None
    assert list(with_id.identifiability.per_type.cell_type) == list(ref.cell_types)
    # estimates must be byte-identical — the certificate is a diagnostic only
    pd.testing.assert_frame_equal(base.deconv.proportions, with_id.deconv.proportions)


def test_recommended_merge_for_confounded_cluster():
    ref = _toy_reference(); bulk = _toy_bulk(ref)
    cert = tr.bulk_identifiability(bulk, ref)
    per = cert.per_type.set_index("cell_type")
    assert "recommended_merge" in cert.per_type.columns
    # A and C are identical → they must land in the same confusable cluster
    assert set(per.loc["A", "cluster"].split(";")) >= {"A", "C"}
    # any confounded type carries a recommended_merge; resolvable types do not
    for c in ref.cell_types:
        if per.loc[c, "recoverability"] in ("GROUP_ONLY", "UNRESOLVABLE"):
            assert per.loc[c, "recommended_merge"]
        elif per.loc[c, "recoverability"] == "RESOLVABLE":
            assert per.loc[c, "recommended_merge"] == ""


def test_default_shift_scale_is_two():
    from tissueresolve.reference.identifiability_calibration import DEFAULT_SHIFT_SCALE
    assert DEFAULT_SHIFT_SCALE == 2.0
    ref = _toy_reference(); bulk = _toy_bulk(ref)
    cert = tr.bulk_identifiability(bulk, ref)
    assert cert.metadata["shift_scale"] == 2.0
    assert cert.metadata["donor_uncertainty"] is True


def test_warns_when_reference_has_no_donor_cv():
    ref = _toy_reference(with_donor=False); bulk = _toy_bulk(ref)
    assert ref.donor_cv is None
    with pytest.warns(UserWarning, match="donor_cv"):
        cert = tr.bulk_identifiability(bulk, ref)
    assert cert.metadata["donor_uncertainty"] is False
