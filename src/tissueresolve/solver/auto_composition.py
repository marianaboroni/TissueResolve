"""Composition-calibrated auto solver (experimental, opt-in) — Etapa 6.

The default :class:`AutoSolver` selects a backbone by **gene-masking CV** (reconstruction
of held-out genes). The Rectangle benchmark showed this objective does not track
**composition** accuracy: on lung it picked plain NNLS (worse than flat) and on breast a
backbone worse than NNLS, and it never considers the count-likelihood (Poisson) solver.

``AutoCompositionSolver`` instead selects the backbone that best recovers **known
proportions** on synthetic mixtures simulated *from the reference itself*, perturbed by
the reference's per-type donor CV (a donor-variation proxy, since the aggregated reference
carries no per-cell data at predict time). It also includes the Poisson GLM as a candidate.
Selection therefore optimises the metric we actually care about (composition RMSE), not
gene reconstruction.

Opt-in / non-default (registered as ``"auto_composition"``); the production ``"auto"``
path is unchanged. Calibration is on **reference-simulated** mixtures — a documented
limitation (not real donor-held-out bulk); it is a within-reference selection heuristic,
benchmark-gated before any default change.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tissueresolve.solver.base import BaseSolver, SolverResult, _normalize


def _composition_candidates():
    from tissueresolve.solver.nnls import NNLSSolver
    from tissueresolve.solver.weighted_nnls import WeightedNNLSSolver
    from tissueresolve.solver.marker_nnls import MarkerNNLSSolver
    from tissueresolve.solver.ridge_nnls import RidgeNNLSSolver
    from tissueresolve.solver.poisson_glm import PoissonGLMSolver
    return [NNLSSolver(), WeightedNNLSSolver(), MarkerNNLSSolver(),
            RidgeNNLSSolver(), PoissonGLMSolver()]


class AutoCompositionSolver(BaseSolver):
    name = "auto_composition"

    def __init__(self, candidates=None, n_sim: int = 40, seed: int = 0,
                 cv_cap: float = 1.0, default_cv: float = 0.3):
        self.candidates = candidates or _composition_candidates()
        self.n_sim = n_sim
        self.seed = seed
        self.cv_cap = cv_cap
        self.default_cv = default_cv

    # -- synthetic, donor-CV-perturbed mixtures from the reference ---------
    def _simulate(self, ref, lib_size: float):
        rng = np.random.default_rng(self.seed)
        R = np.asarray(ref.as_R_cpm(), dtype=float)          # (K, G)
        K, G = R.shape
        # per-(gene,type) donor CV → lognormal sigma (mean-preserving)
        if ref.donor_cv is not None and np.asarray(ref.donor_cv).shape == (G, K):
            sigma = np.clip(np.asarray(ref.donor_cv, float).T, 0, self.cv_cap)  # (K, G)
        else:
            sigma = np.full((K, G), self.default_cv)
        W = rng.dirichlet(np.ones(K), size=self.n_sim)       # (S, K), known truth
        cols = []
        for s in range(self.n_sim):
            # perturb each type's profile by its donor CV, then mix by W
            noise = np.exp(rng.normal(0.0, sigma) - 0.5 * sigma ** 2)   # (K, G)
            Rp = R * noise
            profile = W[s] @ Rp                              # (G,)
            p = profile / profile.sum() if profile.sum() > 0 else np.full(G, 1.0 / G)
            counts = rng.poisson(p * lib_size)               # count noise at query depth
            cols.append(counts.astype(float))
        synth = pd.DataFrame(np.array(cols).T, index=list(ref.gene_names),
                             columns=[f"sim{s}" for s in range(self.n_sim)])
        truth = pd.DataFrame(W, index=synth.columns, columns=list(ref.cell_types))
        return synth, truth

    def select(self, query, ref):
        """Return (best_solver, comparison_df, reason)."""
        col_sums = np.asarray(query.to_numpy(float)).sum(axis=0)
        lib = float(np.median(col_sums[col_sums > 0])) if (col_sums > 0).any() else 1e6
        lib = float(np.clip(lib, 1e4, 5e7))
        synth, truth = self._simulate(ref, lib)
        rows = []
        for sv in self.candidates:
            try:
                est = sv.solve(synth, ref).proportions.reindex(
                    index=truth.index, columns=truth.columns).fillna(0.0)
                rmse = float(np.sqrt(np.mean((est.to_numpy(float) - truth.to_numpy(float)) ** 2)))
            except Exception:  # noqa: BLE001 — a broken candidate must not abort selection
                rmse = float("inf")
            rows.append({"solver": sv.name, "sim_composition_rmse": rmse})
        comp = pd.DataFrame(rows).set_index("solver").sort_values("sim_composition_rmse")
        best_name = comp.index[0]
        best = next(sv for sv in self.candidates if sv.name == best_name)
        reason = (f"selected '{best_name}' by reference-simulated composition RMSE "
                  f"({comp.loc[best_name, 'sim_composition_rmse']:.4f}) over "
                  f"{self.n_sim} donor-CV-perturbed mixtures; candidates include the "
                  f"Poisson GLM (unlike gene-masking 'auto').")
        return best, comp, reason

    def solve(self, query, ref) -> SolverResult:
        best, comp, reason = self.select(query, ref)
        res = best.solve(query, ref)
        res.diagnostics.update({
            "solver": self.name, "selected_solver": best.name,
            "selection_reason": reason, "feature_status": "experimental",
            "sim_composition_comparison": comp.reset_index().to_dict("records"),
        })
        return res
