import importlib.util
import json
from pathlib import Path

import pytest

from sample_decay_dose import download_NIST_nuclide_data as nist
from sample_decay_dose.data import ISOTOPIC_DATA
from sample_decay_dose.ValencyMapper import ValencyMapper


def _load_rel_iso_mass() -> dict:
    """Load isotopes.py from its file under a private name.

    test_SampleDose installs a dummy 'sample_decay_dose.isotopes' only when that module is not imported yet.
    Loading it privately keeps the two tests independent of collection order.
    """
    path = Path(__file__).resolve().parent.parent / 'sample_decay_dose' / 'isotopes.py'
    spec = importlib.util.spec_from_file_location('_isotopes_under_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.rel_iso_mass


@pytest.mark.parametrize("estimate_type", ['upper', 'lower', 'doligez'])
def test_valency_fluoride_cations(estimate_type):
    mapper = ValencyMapper(estimate_type)
    assert mapper['Pb'] == 2  # PbF2
    assert mapper['Bi'] == 3  # BiF3
    assert mapper['Ac'] == 3  # AcF3
    assert mapper['U'] == 4
    assert mapper['Xe'] == 0


def test_modern_superheavy_names():
    rel_iso_mass = _load_rel_iso_mass()
    assert rel_iso_mass['mc-291'] == pytest.approx(291.19707)
    assert rel_iso_mass['ts-294'] == pytest.approx(294.21046)
    assert rel_iso_mass['uup-291'] == rel_iso_mass['mc-291']  # legacy alias
    assert 291 in ISOTOPIC_DATA['mc'] and 294 in ISOTOPIC_DATA['ts']
    assert 'uup' not in ISOTOPIC_DATA and 'uus' not in ISOTOPIC_DATA


def test_nist_parser_writes_to_given_path_independent_of_cwd(tmp_path, monkeypatch):
    source = tmp_path / 'listing.txt'
    source.write_text(
        "Atomic Number = 1\nAtomic Symbol = D\nMass Number = 2\nRelative Atomic Mass = 2.01410177812(12)\n"
        "Isotopic Composition = 0.000115(70)\n\n"
        "Atomic Number = 115\nAtomic Symbol = Uup\nMass Number = 291\nRelative Atomic Mass = 291.19707(88#)\n"
        "Isotopic Composition = \n",
        encoding='utf-8')
    work = tmp_path / 'elsewhere'
    work.mkdir()
    monkeypatch.chdir(work)
    out = tmp_path / 'out' / 'isotopic_data.json'
    nist.download_and_parse_nist_data(source, out)
    data = json.loads(out.read_text(encoding='utf-8'))
    assert data == {'h': {'2': {'mass': 2.01410177812, 'abundance': 0.000115}},
                    'mc': {'291': {'mass': 291.19707, 'abundance': 0.0}}}
    assert not (work / 'data').exists()


def test_nist_defaults_point_at_package_data():
    package_data = Path(__file__).resolve().parent.parent / 'sample_decay_dose' / 'data'
    assert nist.DEFAULT_OUTPUT == package_data / 'isotopic_data.json'
    assert nist.DEFAULT_SOURCE == package_data / 'aw.html'


def test_shipped_json_matches_parser():
    """The shipped JSON is exactly what the parser produces from the shipped NIST listing."""
    parsed = nist.parse_nist_text(nist.DEFAULT_SOURCE.read_text(encoding='utf-8'))
    shipped = json.loads(nist.DEFAULT_OUTPUT.read_text(encoding='utf-8'))
    assert json.loads(json.dumps(parsed)) == shipped
