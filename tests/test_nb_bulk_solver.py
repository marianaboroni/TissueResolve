"""Tests for the experimental Poisson / NB GLM bulk solver (P1, opt-in).

Numerical optimiser tests (non-negative, valid compositions, no NaN, loss
decreases, fallback recorded, no silent gene drop) + config/CLI/pipeline
integration (default wNNLS unchanged, schema parity, unknown solver raises).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from tissueresolve.experimental.nb_bulk_solver import (
    fit_bulk_nb_glm, NBGLMBulkSolver, NBGLMResult)


def _toy(G=60, K=4, N=6, seed=0):
    """Reference (G,K) proportion-scale + bulk counts (G,N) from known truth."""
    rng = np.random.default_rng(seed)
    Phi = rng.random((G, K)) * rng.choice([0.2, 1.0, 5.0], size=(G, K))
    Phi = Phi / Phi.sum(0, keepdims=True)            # columns ~ per-type gene dist
    truth = rng.dirichlet(np.ones(K), size=N).T      # (K, N)
    lib = rng.integers(2000, 8000, N).astype(float)
    mu = lib[None, :] * (Phi @ truth)
    Y = rng.poisson(mu).astype(float)                # (G, N) counts
    return Y, Phi, lib, truth


# ---- array-level solver ----

def test_nonneg_valid_compositions_no_nan():
    Y, Phi, lib, _ = _toy()
    r = fit_bulk_nb_glm(Y, Phi, library_sizes=lib, loss="poisson")
    assert isinstance(r, NBGLMResult)
    assert (r.theta >= -1e-12).all()
    assert np.allclose(r.theta.sum(0), 1.0, atol=1e-6)
    assert not np.any(np.isnan(r.theta)) and np.all(np.isfinite(r.theta))


def test_poisson_loss_decreases():
    Y, Phi, lib, _ = _toy()
    r = fit_bulk_nb_glm(Y, Phi, library_sizes=lib, loss="poisson")
    assert r.final_loss <= r.loss_trace[0] + 1e-6
    assert r.loss_monotonic


def test_nb_loss_decreases_with_dispersion():
    Y, Phi, lib, _ = _toy()
    phi_g = np.full(Phi.shape[0], 10.0)
    r = fit_bulk_nb_glm(Y, Phi, library_sizes=lib, gene_dispersion=phi_g, loss="nb")
    assert r.final_loss <= r.loss_trace[0] + 1e-6
    assert r.metadata["dispersion_source"] == "provided"
    assert r.metadata["likelihood"] == "nb"


def test_nb_without_dispersion_falls_back_to_poisson_recorded():
    Y, Phi, lib, _ = _toy()
    with pytest.warns(UserWarning):
        r = fit_bulk_nb_glm(Y, Phi, library_sizes=lib, gene_dispersion=None, loss="nb")
    assert r.metadata["fallback_to_poisson"] is True
    assert r.metadata["dispersion_source"] == "fallback_poisson"
    assert r.metadata["likelihood"] == "poisson"


def test_recovers_truth_reasonably():
    Y, Phi, lib, truth = _toy(G=120, seed=3)
    r = fit_bulk_nb_glm(Y, Phi, library_sizes=lib, loss="poisson", max_iter=800)
    corr = np.corrcoef(r.theta.ravel(), truth.ravel())[0, 1]
    assert corr > 0.9  # well-conditioned toy → strong recovery


def test_unknown_optimizer_and_loss_raise():
    Y, Phi, lib, _ = _toy()
    with pytest.raises(ValueError):
        fit_bulk_nb_glm(Y, Phi, library_sizes=lib, optimizer="lbfgs")
    with pytest.raises(ValueError):
        fit_bulk_nb_glm(Y, Phi, library_sizes=lib, loss="gaussian")


def test_library_size_default_is_column_sums():
    Y, Phi, _, _ = _toy()
    r = fit_bulk_nb_glm(Y, Phi, library_sizes=None, loss="poisson")
    assert r.metadata["library_size_normalization"] is True
    assert np.allclose(r.theta.sum(0), 1.0, atol=1e-6)


# ---- pipeline-compatible wrapper / schema parity ----

def _ref_and_bulk():
    from tissueresolve.results import ReferenceSignature
    rng = np.random.default_rng(1)
    G, K, N = 50, 3, 5
    genes = [f"g{i}" for i in range(G)]
    cts = [f"CT{k}" for k in range(K)]
    R_cpm = rng.random((K, G)) * 100
    R_cpm = R_cpm / R_cpm.sum(1, keepdims=True) * 1e6     # CPM rows sum 1e6
    ref = ReferenceSignature(gene_names=genes, cell_types=cts,
                             phi=(R_cpm.T / R_cpm.T.sum(0, keepdims=True)).astype(float),
                             R_cpm=R_cpm.astype(np.float32),
                             R_log=np.log1p(R_cpm).astype(np.float32),
                             phi_g=np.full(G, 8.0, dtype=np.float32),
                             n_cells_per_type={c: 50 for c in cts}, genome="hg38")
    truth = rng.dirichlet(np.ones(K), size=N)
    lib = rng.integers(3000, 6000, N)
    counts = (truth @ (R_cpm / 1e6)).T * lib[None, :]
    bulk = pd.DataFrame(rng.poisson(counts), index=genes, columns=[f"s{i}" for i in range(N)])
    return ref, bulk, genes


def test_wrapper_returns_bulk_schema():
    from tissueresolve.results import BulkDeconvResult
    ref, bulk, genes = _ref_and_bulk()
    res = NBGLMBulkSolver(loss="nb").solve(bulk, ref, genes)
    assert isinstance(res, BulkDeconvResult)
    assert res.ESTIMATE_TYPE == "mRNA_proportion"
    P = res.proportions.to_numpy()
    assert np.allclose(P.sum(1), 1.0, atol=1e-5) and (P >= -1e-9).all()
    assert list(res.proportions.columns) == list(ref.cell_types)
    assert res.run_metadata["solver"] == "nb_glm_experimental"


def test_no_silent_gene_drop_reports_overlap():
    ref, bulk, genes = _ref_and_bulk()
    # add genes to bulk not in panel → must not be silently used; overlap is explicit
    res = NBGLMBulkSolver(loss="poisson").solve(bulk, ref, genes[:40])
    assert res.run_metadata["n_genes_panel"] == 40
    assert len(res.gene_panel) == 40


# ---- config / pipeline integration ----

def test_default_bulk_solver_is_wNNLS():
    from tissueresolve.config import TissueResolveConfig
    assert TissueResolveConfig().bulk_solver.method == "wNNLS"


def test_pipeline_default_unchanged_uses_wNNLS():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.bulk.pipeline import BulkPipeline
    ref, bulk, genes = _ref_and_bulk()
    res = BulkPipeline(TissueResolveConfig()).run(bulk, ref, gene_panel=genes,
                                                  n_bootstrap=0, run_qc=False)
    assert res.deconv.run_metadata["solver"] == "WNNLSSolver"


def test_pipeline_opt_in_nb_glm():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.bulk.pipeline import BulkPipeline
    ref, bulk, genes = _ref_and_bulk()
    cfg = TissueResolveConfig()
    cfg.bulk_solver.method = "nb_glm_experimental"
    res = BulkPipeline(cfg).run(bulk, ref, gene_panel=genes, n_bootstrap=0, run_qc=False)
    md = res.deconv.run_metadata
    assert md["solver"] == "nb_glm_experimental"
    P = res.deconv.proportions.to_numpy()
    assert np.allclose(P.sum(1), 1.0, atol=1e-5) and (P >= -1e-9).all()


def test_pipeline_unknown_solver_raises():
    from tissueresolve.config import TissueResolveConfig
    from tissueresolve.bulk.pipeline import BulkPipeline
    ref, bulk, genes = _ref_and_bulk()
    cfg = TissueResolveConfig()
    cfg.bulk_solver.method = "totally_unknown"
    with pytest.raises(ValueError):
        BulkPipeline(cfg).run(bulk, ref, gene_panel=genes, n_bootstrap=0, run_qc=False)


def test_cli_accepts_bulk_solver_option():
    from tissueresolve.cli import run_cli
    opt = next((p for p in run_cli.params if p.name == "bulk_solver"), None)
    assert opt is not None
    for c in ("wNNLS", "poisson_glm_experimental", "nb_glm_experimental"):
        assert c in opt.type.choices


def test_config_yaml_roundtrip_bulk_method(tmp_path):
    from tissueresolve.config import TissueResolveConfig
    cfg = TissueResolveConfig()
    cfg.bulk_solver.method = "poisson_glm_experimental"
    p = tmp_path / "cfg.yaml"
    cfg.to_yaml(p)
    assert TissueResolveConfig.from_yaml(p).bulk_solver.method == "poisson_glm_experimental"
