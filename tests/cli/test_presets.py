from tissueresolve.presets import get_preset


def test_get_preset_standard():
    p = get_preset("standard")
    assert isinstance(p, dict)
    assert "plots" in p
