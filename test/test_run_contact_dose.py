""" Mocked-pipeline tests for the concrete-lid contact dose driver """
import importlib.util
import os
import sys

import pytest

_SPEC = importlib.util.spec_from_file_location('run_contact_dose', os.path.join(
    os.path.dirname(__file__), '..', 'concrete_irrad', 'run_contact_dose.py'))
assert _SPEC is not None and _SPEC.loader is not None
driver = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(driver)

from concrete_irrad.concrete_lid_contact_dose import ConcreteLidContactDose


@pytest.fixture()
def mocked_chain(monkeypatch, tmp_path):
    """ Stub out every SCALE-facing step and record the calculations and MPI counts they receive """
    monkeypatch.chdir(tmp_path)
    record: dict = {'calcs': [], 'activation_nmpi': [], 'mavric_nmpi': []}

    def run_activation(self, nmpi=1):
        record['activation_nmpi'].append(nmpi)

    def run_mavric(self, nmpi=1):
        record['mavric_nmpi'].append(nmpi)

    monkeypatch.setattr(ConcreteLidContactDose, 'write_inputs', lambda self: record['calcs'].append(self))
    monkeypatch.setattr(ConcreteLidContactDose, 'run_activation', run_activation)
    monkeypatch.setattr(ConcreteLidContactDose, 'run_mavric', run_mavric)
    monkeypatch.setattr(ConcreteLidContactDose, 'get_responses',
                        lambda self: self.responses.update({'2': {'value': 10.0 * self.decay_days,
                                                                  'stdev': 0.1}}))
    return record


@pytest.fixture()
def custom_mixer(monkeypatch):
    """ Non-default rebar grid, so the tests can tell whether MIXER_PARAMS reached the model """
    params = {'rebar_diameter_cm': 2.0, 'spacing_x_cm': 20.0, 'spacing_y_cm': 25.0}
    for key, value in params.items():
        monkeypatch.setitem(driver.MIXER_PARAMS, key, value)
    return params


def _assert_mixer(calc, params):
    assert calc.mixer.rebar_diameter_cm == params['rebar_diameter_cm']
    assert calc.mixer.spacing_x_cm == params['spacing_x_cm']
    assert calc.mixer.spacing_y_cm == params['spacing_y_cm']
    assert calc.mixed_thickness_cm == params['rebar_diameter_cm']


def test_collect_dose_scales_with_decay_days(mocked_chain):
    calc = ConcreteLidContactDose(decay_days=3.0, irradiation_lib_f33='fake.f33')
    row = driver.collect_dose(calc)
    assert row == {'decay_days': 3.0, 'dose_rem_per_h': pytest.approx(30.0),
                   'stdev_rem_per_h': pytest.approx(0.1)}


def test_collect_dose_forwards_nmpi(mocked_chain):
    driver.collect_dose(ConcreteLidContactDose(irradiation_lib_f33='fake.f33'), nmpi=4)
    assert mocked_chain['activation_nmpi'] == [4]
    assert mocked_chain['mavric_nmpi'] == [4]


def test_scan_decay_times_rows_and_order(mocked_chain):
    rows = driver.scan_decay_times([90.0, 1.0], 'fake.f33')
    assert [r['decay_days'] for r in rows] == [90.0, 1.0]
    assert all(r['dose_rem_per_h'] == pytest.approx(10.0 * r['decay_days']) for r in rows)


def test_scan_decay_times_wires_nmpi_and_mixer(mocked_chain, custom_mixer):
    driver.scan_decay_times([90.0, 1.0], 'fake.f33', nmpi=3)
    assert mocked_chain['activation_nmpi'] == [3, 3]
    assert mocked_chain['mavric_nmpi'] == [3, 3]
    assert len(mocked_chain['calcs']) == 2
    for calc in mocked_chain['calcs']:
        _assert_mixer(calc, custom_mixer)


def test_main_single_run_wires_nmpi_and_mixer(mocked_chain, custom_mixer, monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['run_contact_dose.py', '--lib', 'fake.f33',
                                      '--decay-days', '2', '--nmpi', '6'])
    driver.main()
    assert mocked_chain['activation_nmpi'] == [6]
    assert mocked_chain['mavric_nmpi'] == [6]
    (calc,) = mocked_chain['calcs']
    assert calc.decay_days == 2.0
    assert calc.irradiation_lib_f33 == 'fake.f33'
    _assert_mixer(calc, custom_mixer)


def test_main_scan_ignores_empty_items(mocked_chain, monkeypatch, tmp_path):
    csv_path = os.path.join(str(tmp_path), 'scan.csv')
    monkeypatch.setattr(sys, 'argv', ['run_contact_dose.py', '--lib', 'fake.f33',
                                      '--scan', '1,30,', '--nmpi', '2', '--csv', csv_path])
    driver.main()
    assert [c.decay_days for c in mocked_chain['calcs']] == [1.0, 30.0]
    assert mocked_chain['activation_nmpi'] == [2, 2]
    with open(csv_path) as f:
        assert len(f.read().splitlines()) == 3


def test_main_scan_without_times_exits(mocked_chain, monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['run_contact_dose.py', '--lib', 'fake.f33', '--scan', ' , '])
    with pytest.raises(SystemExit):
        driver.main()
    assert mocked_chain['calcs'] == []


def test_main_preview_uses_mixer_params(custom_mixer, monkeypatch, tmp_path, capsys):
    """ The preview deck carries the edited rebar diameter as the mixed-layer thickness """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['run_contact_dose.py'])
    driver.main()
    out = capsys.readouterr().out
    mixed_volume = driver.PARAMS['length_cm'] * driver.PARAMS['width_cm'] * custom_mixer['rebar_diameter_cm']
    assert f'volume={mixed_volume:.6e}' in out
    assert os.listdir(str(tmp_path)) == []


@pytest.mark.parametrize('text, expected', [('1,30,90', [1.0, 30.0, 90.0]), ('1,30,', [1.0, 30.0]),
                                            (' 5 , ,7', [5.0, 7.0]), ('', [])])
def test_parse_scan(text, expected):
    assert driver.parse_scan(text) == expected


def test_write_csv_roundtrip(tmp_path):
    rows = [{'decay_days': 1.0, 'dose_rem_per_h': 1.5, 'stdev_rem_per_h': 0.2}]
    path = os.path.join(str(tmp_path), 'scan.csv')
    driver.write_csv(rows, path)
    with open(path) as f:
        content = f.read()
    assert content.splitlines()[0] == 'decay_days,dose_rem_per_h,stdev_rem_per_h'
    assert '1.5' in content
