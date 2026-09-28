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


def test_run_series_tracks_f18_and_titles_from_parameters(monkeypatch, tmp_path):
    fake_index = {1: {"case": "1", "time": "0.0"}, 2: {"case": "1", "time": "86400.0"}}
    monkeypatch.setattr(c1, "get_f71_positions_index", lambda _path: fake_index)
    monkeypatch.setattr(c1, "get_burned_nuclide_atom_dens", lambda _path, _pos: {"f-19": 0.04, "f-18": 1e-9})

    def fake_run_scale(deck):
        open("origen.f71", "w").close()
        return True

    monkeypatch.setattr(c1, "run_scale", fake_run_scale)
    monkeypatch.setattr(c1, "get_burned_nuclide_data",
                        lambda _path, _pos, f71units="becq": {"f-18": 5.0, "f-19": 0.0, "na-24": 7.0})
    captured = {}

    def fake_plot(df, png_path, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(c1, "_plot_activity", fake_plot)
    df, csv_path, _ = c1.run_case1_decay_series(
        f71_path="dummy.f71", case=1, volume_liters=3.25, decay_days=0.5, out_root=str(tmp_path)
    )
    assert df["tracked_nuclide"].tolist() == ["f-18", "f-18"]
    assert df["tracked_activity_bq"].tolist() == [5.0, 5.0]
    assert df["total_activity_bq"].tolist() == [12.0, 12.0]
    assert "f19_activity_bq" not in df.columns
    assert csv_path.is_file()
    assert captured["title"] == "0.5-day decay of 3.25 L sample from each case-1 timestep"
    assert captured["track_nuclide"] == "f-18"
    deck = (tmp_path / "case1_pos0001_t0000000000.0s" / "origen.inp").read_text()
    assert "volume=3250.0" in deck


def test_plot_title_and_label():
    assert c1._plot_title(decay_days=2.0, volume_liters=6.5, case=1) == \
        "2-day decay of 6.5 L sample from each case-1 timestep"
    assert c1._nuclide_label("kr-85m") == "Kr-85m"


def test_cli_track_nuclide_default():
    args = c1._build_parser().parse_args([])
    assert args.track_nuclide == "f-18"
    args = c1._build_parser().parse_args(["--track-nuclide", "na-24"])
    assert args.track_nuclide == "na-24"
