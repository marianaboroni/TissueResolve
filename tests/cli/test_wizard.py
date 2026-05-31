from pathlib import Path

from tissueresolve import cli


def test_wizard_writes_config(tmp_path, monkeypatch):
    answers = iter(["ref.h5ad", "query.h5ad", "auto", "standard"])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    # run wizard
    rc = cli.wizard([])
    assert rc == 0
    assert Path("tissueresolve_config.yaml").exists()
    # cleanup
    Path("tissueresolve_config.yaml").unlink()
