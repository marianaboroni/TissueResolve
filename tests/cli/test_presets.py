from tissueresolve.presets import get_preset


def test_get_preset_standard():
    p = get_preset("standard")
    assert isinstance(p, dict)
    assert "plots" in p


def test_preset_bootstrap_disabled_sets_n_bootstrap_zero():
    """Presets that declare ``bootstrap: False`` must actually disable bootstrap
    (n_bootstrap == 0), not silently fall back to the config default (200)."""
    from tissueresolve.cli import _configure_from_preset

    for name in ("quick", "standard"):
        cfg = _configure_from_preset(get_preset(name))
        assert cfg.bootstrap.n_bootstrap == 0, name

    assert _configure_from_preset(get_preset("publication")).bootstrap.n_bootstrap == 100
    assert _configure_from_preset(get_preset("diagnostic")).bootstrap.n_bootstrap == 200
