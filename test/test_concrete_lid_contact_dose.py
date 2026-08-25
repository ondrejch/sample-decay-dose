""" Structural and mocked-pipeline tests for the concrete-lid contact dose prototype """
import importlib.util
import os

import pytest

_SPEC = importlib.util.spec_from_file_location('concrete_lid_contact_dose', os.path.join(
    os.path.dirname(__file__), '..', 'concrete_irrad', 'concrete_lid_contact_dose.py'))
assert _SPEC is not None and _SPEC.loader is not None
cld = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cld)


@pytest.fixture()
def calc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return cld.ConcreteLidContactDose(name='lid', irradiation_lib_f33='spectrum.f33')


def test_origen_deck_chains_both_regions(calc):
    deck = calc.origen_deck()
    for case_name in ('case(irr_bottom_slab)', 'case(dec_bottom_slab)',
                      'case(irr_mixed_layer)', 'case(dec_mixed_layer)'):
        assert case_name in deck
    assert deck.count('file="lid_bottom.f71"') == 1
    assert deck.count('file="lid_mixed.f71"') == 1
    assert 'flux=[' in deck and f'{calc.irradiation_steps}R 10000000000.0' in deck
    assert 'end7dec' in deck


def test_origen_deck_requires_spectrum_library():
    with pytest.raises(ValueError):
        cld.ConcreteLidContactDose().origen_deck()


def test_fresh_atom_densities(calc):
    assert calc.fresh_adens['mixed_layer']['co-59'] > 0
    assert calc.fresh_adens['bottom_slab']['o-16'] > 0
    assert 'fe-56' in calc.fresh_adens['mixed_layer']


def test_mavric_deck_structure(calc):
    calc.source_strengths = {'bottom_slab': 1.0e8, 'mixed_layer': 5.0e9}
    deck = calc.mavric_deck()
    assert '<lid_mixed_adens_mavric.inp' in deck
    assert '<lid_bottom_adens_mavric.inp' in deck
    assert 'cuboid 11' in deck and 'cuboid 12' in deck
    assert 'strength=5.000000e+09' in deck and 'strength=1.000000e+08' in deck
    assert deck.count('origensBinaryConcentrationFile') == 2
    assert 'filename="lid_mixed.f71"' in deck and 'filename="lid_bottom.f71"' in deck
    assert 'doseData=9505' in deck
    assert f'position 0.0 0.0 {calc.bottom_thickness_cm + calc.mixed_thickness_cm + 0.1}' in deck


def test_source_strengths_use_gato_and_volume(calc, monkeypatch):
    monkeypatch.setattr(cld, 'get_burned_nuclide_data',
                        lambda f71, pos, units: {'co-60': 2.0, 'eu-152': 3.0})
    strengths = calc._photon_source_strengths()
    vol_mix = calc.length_cm * calc.width_cm * calc.mixed_thickness_cm
    assert strengths['mixed_layer'] == pytest.approx(5.0 * vol_mix)
    assert strengths['bottom_slab'] > 0


def test_write_inputs_creates_case_files(calc):
    calc.write_inputs()
    case = os.path.join(calc.cwd, calc.case_dir)
    assert os.path.isfile(os.path.join(case, calc.ORIGEN_input_file_name))
    for region in calc.regions.values():
        assert os.path.isfile(os.path.join(case, region['adens_file']))


def test_full_pipeline_with_mocked_scale(calc, monkeypatch):
    fake_decay = {'co-60': 1.0e-10, 'eu-154': 5.0e-12}
    runs: list[str] = []

    def fake_run(deck_file, nmpi=1):
        runs.append(os.path.join(os.getcwd(), deck_file))
        return True

    monkeypatch.setattr(cld, 'run_scale_or_raise', fake_run)
    monkeypatch.setattr(cld, 'get_burned_nuclide_atom_dens', lambda f71, pos: dict(fake_decay))
    monkeypatch.setattr(cld, 'get_burned_nuclide_data', lambda f71, pos, units:
                        {'co-60': 1000.0} if 'mixed' in f71 else {'co-60': 100.0})

    calc.run_activation()
    assert set(calc.decayed_adens) == {'bottom_slab', 'mixed_layer'}

    for region in calc.regions.values():
        with open(os.path.join(calc.cwd, calc.case_dir, region['f71']), 'w') as f:
            f.write('fake f71')
    calc.run_mavric()

    case_mavric = os.path.join(calc.cwd, calc.case_dir_mavric)
    assert len(runs) == 2
    assert os.path.isfile(os.path.join(case_mavric, calc.MAVRIC_input_file_name))
    with open(os.path.join(case_mavric, 'lid_mixed_adens_mavric.inp')) as f:
        comp_line = f.readline()
        assert comp_line.startswith('co-60 1 ')

    with open(os.path.join(case_mavric, calc.MAVRIC_out_file_name), 'w') as f:
        f.write('junk\nFinal Tally Results Summary\n\n'
                '   response  2   1.2345e+01  2.500e-02\n')

    calc.get_responses()
    dose = calc.contact_dose
    assert dose['value'] == pytest.approx(12.345)
    assert dose['stdev'] == pytest.approx(0.025)


def test_mavric_before_activation_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = cld.ConcreteLidContactDose()
    with pytest.raises(RuntimeError):
        calc.run_mavric()


def test_missing_response_raises(calc, tmp_path):
    os.mkdir(os.path.join(str(tmp_path), calc.case_dir_mavric))
    with open(os.path.join(str(tmp_path), calc.case_dir_mavric, calc.MAVRIC_out_file_name), 'w') as f:
        f.write('no tally here')
    with pytest.raises((RuntimeError, FileNotFoundError)):
        calc.get_responses()
