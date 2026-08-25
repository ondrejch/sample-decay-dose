""" Mocked-pipeline tests for the concrete-lid contact dose driver """
import importlib.util
import os

import pytest

_SPEC = importlib.util.spec_from_file_location('run_contact_dose', os.path.join(
    os.path.dirname(__file__), '..', 'concrete_irrad', 'run_contact_dose.py'))
assert _SPEC is not None and _SPEC.loader is not None
driver = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(driver)

from concrete_irrad.concrete_lid_contact_dose import ConcreteLidContactDose


@pytest.fixture()
def mocked_chain(monkeypatch):
    """ Stub out every SCALE-facing step on the prototype class """
    monkeypatch.setattr(ConcreteLidContactDose, 'write_inputs', lambda self: None)
    monkeypatch.setattr(ConcreteLidContactDose, 'run_activation', lambda self: None)
    monkeypatch.setattr(ConcreteLidContactDose, 'run_mavric', lambda self: None)
    monkeypatch.setattr(ConcreteLidContactDose, 'get_responses',
                        lambda self: self.responses.update({'2': {'value': 10.0 * self.decay_days,
                                                                  'stdev': 0.1}}))


def test_collect_dose_scales_with_decay_days(mocked_chain):
    calc = ConcreteLidContactDose(decay_days=3.0, irradiation_lib_f33='fake.f33')
    row = driver.collect_dose(calc)
    assert row == {'decay_days': 3.0, 'dose_rem_per_h': pytest.approx(30.0),
                   'stdev_rem_per_h': pytest.approx(0.1)}


def test_scan_decay_times_rows_and_order(mocked_chain):
    rows = driver.scan_decay_times([90.0, 1.0], 'fake.f33')
    assert [r['decay_days'] for r in rows] == [90.0, 1.0]
    assert all(r['dose_rem_per_h'] == pytest.approx(10.0 * r['decay_days']) for r in rows)


def test_write_csv_roundtrip(tmp_path):
    rows = [{'decay_days': 1.0, 'dose_rem_per_h': 1.5, 'stdev_rem_per_h': 0.2}]
    path = os.path.join(str(tmp_path), 'scan.csv')
    driver.write_csv(rows, path)
    with open(path) as f:
        content = f.read()
    assert content.splitlines()[0] == 'decay_days,dose_rem_per_h,stdev_rem_per_h'
    assert '1.5' in content
