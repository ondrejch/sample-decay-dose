import re

from leaky_box_origen import Case1F71DecaySeries as c1


def test_case_positions_filters_and_sorts(monkeypatch):
    fake_index = {
        3: {"case": "1", "time": "30.0"},
        1: {"case": "2", "time": "10.0"},
        2: {"case": "1", "time": "20.0"},
    }
    monkeypatch.setattr(c1, "get_f71_positions_index", lambda _path: fake_index)
    assert c1._case_positions("dummy.f71", case=1) == [(2, 20.0), (3, 30.0)]


def test_origen_decay_deck_contains_requested_volume_and_time():
    deck = c1._origen_decay_deck(
        atom_file="sample_atom_dens.inp",
        out_f71="origen.f71",
        volume_cm3=6500.0,
        decay_days=2.0,
        decay_steps=30,
    )
    assert "volume=6500.0" in deck
    assert "t=[27I 0.01 2.0]" in deck
    assert re.search(r'file="origen\.f71"', deck) is not None
