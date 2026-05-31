"""
Tests for the ``tissueresolve spatial`` CLI commands (Stage 4 wiring).

Fast, deterministic, offline:  synthetic Visium written to a temp .h5ad and a
saved ReferenceSignature directory drive an end-to-end ``spatial run``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from click.testing import CliRunner

from tissueresolve.cli import cli
from tissueresolve.io.spatial import make_synthetic_visium
from tissueresolve.results import ReferenceSignature

N_SPOTS = 40
N_GENES = 30
N_TYPES = 3
GENE_NAMES = [f"GENE_{i:04d}" for i in range(N_GENES)]
CELL_TYPES = [f"CT{k}" for k in range(N_TYPES)]


def _write_reference(tmp_path):
    """Build and save a ReferenceSignature matching the synthetic Visium genes."""
    rng = np.random.default_rng(0)
    block = N_GENES // N_TYPES
    R_cpm = np.full((N_TYPES, N_GENES), 10.0, dtype=np.float32)
    for k in range(N_TYPES):
        R_cpm[k, k * block:(k + 1) * block] = rng.uniform(500, 1500, block)
    phi_g = np.full(N_GENES, 5.0, dtype=np.float32)
    ref = ReferenceSignature(
        gene_names=GENE_NAMES,
        cell_types=CELL_TYPES,
        R_cpm=R_cpm,
        R_log=np.log1p(R_cpm).astype(np.float32),
        phi_g=phi_g,
        n_cells_per_type={ct: 100 for ct in CELL_TYPES},
    )
    ref_dir = tmp_path / "reference"
    ref.save(ref_dir)
    return ref_dir


def _write_visium(tmp_path):
    adata = make_synthetic_visium(
        n_spots=N_SPOTS, n_genes=N_GENES, n_cell_types=N_TYPES, seed=1
    )
    h5ad = tmp_path / "visium.h5ad"
    adata.write_h5ad(h5ad)
    return h5ad


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------


class TestSpatialInfo:
    def test_runs_and_prints_config(self):
        result = CliRunner().invoke(cli, ["spatial", "info"])
        assert result.exit_code == 0, result.output
        assert "lambda_spatial" in result.output
        assert "random_state" in result.output


# ---------------------------------------------------------------------------
# report (explicit not-implemented)
# ---------------------------------------------------------------------------


class TestSpatialReport:
    def test_report_help(self):
        # The spatial report command is now implemented (publication layer):
        # it builds an HTML report from a results directory.
        result = CliRunner().invoke(cli, ["spatial", "report", "--help"])
        assert result.exit_code == 0
        assert "--results-dir" in result.output

    def test_report_runs_on_results_dir(self, tmp_path):
        import numpy as np
        import pandas as pd

        rdir = tmp_path / "spatial"
        rdir.mkdir()
        pd.DataFrame(np.eye(2), index=["sp0", "sp1"], columns=["A", "B"]).to_csv(
            rdir / "spatial_spot_proportions.tsv", sep="\t")
        result = CliRunner().invoke(
            cli, ["spatial", "report", "--results-dir", str(rdir)])
        assert result.exit_code == 0, result.output
        assert (rdir / "report.html").exists()


# ---------------------------------------------------------------------------
# benchmark
# ---------------------------------------------------------------------------


class TestSpatialBenchmark:
    def test_benchmark_writes_csv(self, tmp_path):
        out = tmp_path / "bench"
        result = CliRunner().invoke(cli, [
            "spatial", "benchmark",
            "--output", str(out),
            "--n-spots", "40", "--n-types", "3", "--n-genes", "30",
            "--scenarios", "basic", "--max-iter", "15", "--seed", "0",
        ])
        assert result.exit_code == 0, result.output
        csv = out / "benchmark_results.csv"
        assert csv.exists()
        frame = pd.read_csv(csv)
        assert {"scenario", "method", "rmse"}.issubset(frame.columns)


# ---------------------------------------------------------------------------
# run (end-to-end)
# ---------------------------------------------------------------------------


class TestSpatialRun:
    def _invoke(self, tmp_path, extra=None):
        ref_dir = _write_reference(tmp_path)
        h5ad = _write_visium(tmp_path)
        out = tmp_path / "out"
        marker_file = tmp_path / "markers.txt"
        marker_file.write_text("\n".join(GENE_NAMES))
        args = [
            "spatial", "run",
            "--visium", str(h5ad),
            "--reference", str(ref_dir),
            "--output", str(out),
            "--marker-genes", str(marker_file),
            "--min-counts", "0", "--min-genes", "0",
            "--max-iter", "20", "--random-state", "0",
        ]
        if extra:
            args += extra
        result = CliRunner().invoke(cli, args)
        return result, out

    def test_run_end_to_end(self, tmp_path):
        result, out = self._invoke(tmp_path)
        assert result.exit_code == 0, result.output
        # Outputs written
        assert (out / "deconv").is_dir()
        assert (out / "qc").is_dir()
        assert (out / "run_metadata.json").exists()
        # Estimate-type message present (never reported as cell counts)
        assert "spot_rna_composition" in result.output

    def test_run_proportions_saved(self, tmp_path):
        from tissueresolve.results import SpatialDeconvResult

        result, out = self._invoke(tmp_path)
        assert result.exit_code == 0, result.output
        loaded = SpatialDeconvResult.load(out / "deconv")
        assert loaded.proportions.shape == (N_SPOTS, N_TYPES)
        np.testing.assert_allclose(
            loaded.proportions.to_numpy().sum(axis=1), 1.0, atol=1e-4
        )

    def test_run_lambda_override_recorded(self, tmp_path):
        import json

        result, out = self._invoke(tmp_path, extra=["--lambda-spatial", "0.0"])
        assert result.exit_code == 0, result.output
        with (out / "run_metadata.json").open() as fh:
            meta = json.load(fh)
        assert meta["lambda_spatial"] == 0.0

    def test_run_missing_visium_errors(self, tmp_path):
        ref_dir = _write_reference(tmp_path)
        result = CliRunner().invoke(cli, [
            "spatial", "run",
            "--visium", str(tmp_path / "nope.h5ad"),
            "--reference", str(ref_dir),
            "--output", str(tmp_path / "out"),
        ])
        assert result.exit_code != 0
