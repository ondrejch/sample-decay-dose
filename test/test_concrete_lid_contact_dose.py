""" Structural and mocked-pipeline tests for the concrete-lid contact dose prototype """
import importlib.util
import os
import re
import subprocess
import sys

import pytest

_SPEC = importlib.util.spec_from_file_location('concrete_lid_contact_dose', os.path.join(
    os.path.dirname(__file__), '..', 'concrete_irrad', 'concrete_lid_contact_dose.py'))
assert _SPEC is not None and _SPEC.loader is not None
cld = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cld)


@pytest.fixture()
def calc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with open(os.path.join(str(tmp_path), 'spectrum.f33'), 'wb') as f:
        f.write(b'fake f33 payload')
    return cld.ConcreteLidContactDose(name='lid', irradiation_lib_f33='spectrum.f33')


def _cuboid_bounds(deck: str, body: str) -> list[float]:
    """ Six numeric bounds of a named geometry cuboid (e.g. '11') or of the n-th source cuboid ('src1') """
    if body.startswith('src'):
        sources_block = deck.split('read sources')[1].split('end sources')[0]
        lines = re.findall(r'^\s+cuboid ([-\d.e+ ]+)$', sources_block, flags=re.M)
        return [float(v) for v in lines[int(body[3:]) - 1].split()]
    m = re.search(rf'^\s+cuboid {body} ([-\d.e+ ]+)$', deck, flags=re.M)
    assert m is not None, f"cuboid {body} missing"
    return [float(v) for v in m.group(1).split()]


def _decay_grids(deck: str) -> list[tuple[int, float, float]]:
    """ (interpolated steps, first time, last time) of each decay-case ORIGEN time list """
    grids = []
    for block in re.findall(r'case\(dec_\w+\) \{(.*?)\n\}', deck, flags=re.S):
        m = re.search(r't=\[(\d+)L (\S+) (\S+)\]', block)
        assert m is not None
        grids.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return grids


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
    deck = calc.mavric_deck()
    assert '<lid_mixed_adens_mavric.inp' in deck
    assert '<lid_bottom_adens_mavric.inp' in deck
    assert 'cuboid 11' in deck and 'cuboid 12' in deck
    assert deck.count('origensBinaryConcentrationFile') == 2
    assert 'filename="lid_mixed.f71"' in deck and 'filename="lid_bottom.f71"' in deck
    assert 'doseData=9505' in deck
    assert f'position 0.0 0.0 {calc.bottom_thickness_cm + calc.mixed_thickness_cm + 0.1}' in deck


def test_sources_use_norm_const_of_their_own_region_distribution(calc):
    """ Each src takes its strength from its own region's F71 photon distribution """
    deck = calc.mavric_deck()
    assert 'strength=' not in deck and 'multiplier=' not in deck
    sources = re.findall(r'src (\d+)\n(.*?)end src', deck, flags=re.S)
    assert len(sources) == 2
    z_top = calc.bottom_thickness_cm + calc.mixed_thickness_cm
    expected = {'mixed_layer': ('lid_mixed.f71', z_top, calc.bottom_thickness_cm),
                'bottom_slab': ('lid_bottom.f71', calc.bottom_thickness_cm, 0.0)}
    for src_id, body in sources:
        assert 'useNormConst' in body
        region = re.search(r'title="Decay photons, (\w+)"', body).group(1)
        dist_id = re.search(r'eDistributionID=(\d+)', body).group(1)
        dist = re.search(rf'distribution {dist_id}\n(.*?)end distribution', deck, flags=re.S).group(1)
        f71, z_hi, z_lo = expected[region]
        assert f'filename="{f71}"' in dist
        assert f'parameters {calc.f71_position} 5 end' in dist
        assert _cuboid_bounds(deck, f'src{src_id}')[4:] == pytest.approx([z_hi, z_lo])


def test_cuboids_list_bounds_max_first(calc):
    """ Geometry bodies and source cuboids use +X -X +Y -Y +Z -Z """
    deck = calc.mavric_deck()
    for body in ('11', '12', '99999', 'src1', 'src2'):
        b = _cuboid_bounds(deck, body)
        assert b[0] > b[1] and b[2] > b[3] and b[4] > b[5], (body, b)
    assert _cuboid_bounds(deck, '11') == pytest.approx([50.0, -50.0, 50.0, -50.0, 5.0, 0.0])
    assert _cuboid_bounds(deck, '12') == pytest.approx([50.0, -50.0, 50.0, -50.0, 5.0 + 1.59, 5.0])


def test_non_square_lid_boundary_uses_both_half_widths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = cld.ConcreteLidContactDose(length_cm=100.0, width_cm=300.0, irradiation_lib_f33='x.f33')
    deck = calc.mavric_deck()
    lid = _cuboid_bounds(deck, '12')
    boundary = _cuboid_bounds(deck, '99999')
    assert lid[:4] == pytest.approx([50.0, -50.0, 150.0, -150.0])
    assert boundary[:4] == pytest.approx([60.1, -60.1, 160.1, -160.1])
    assert boundary[4] > lid[4] and boundary[5] < 0.0


def test_decay_time_grid_is_strictly_increasing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for days in (30.0, 0.05, 0.0):
        calc = cld.ConcreteLidContactDose(decay_days=days, irradiation_lib_f33='x.f33')
        grids = _decay_grids(calc.origen_deck())
        assert len(grids) == 2
        for steps, t_first, t_last in grids:
            assert steps == calc.f71_position - 3
            assert 0.0 < t_first < t_last
            if days > 0:
                assert t_last == pytest.approx(days)


def test_zero_decay_days_reads_decay_case_start(tmp_path, monkeypatch):
    """ Zero cooling reads F71 position 1, the end-of-irradiation state of the decay case """
    monkeypatch.chdir(tmp_path)
    zero = cld.ConcreteLidContactDose(decay_days=0.0, irradiation_lib_f33='x.f33')
    assert zero.decayed_f71_position == 1
    assert zero.mavric_deck().count('parameters 1 5 end') == 2
    cooled = cld.ConcreteLidContactDose(decay_days=30.0, irradiation_lib_f33='x.f33')
    assert cooled.decayed_f71_position == cooled.f71_position

    open('x.f33', 'wb').close()
    positions: list[int] = []
    monkeypatch.setattr(cld, 'run_scale_or_raise', lambda deck, nmpi=1, context_dir='': True)
    monkeypatch.setattr(cld, 'get_burned_nuclide_atom_dens',
                        lambda f71, pos: positions.append(pos) or {'na-24': 1e-9})
    zero.run_activation()
    assert positions == [1, 1]


def test_invalid_times_raise():
    with pytest.raises(ValueError):
        cld.ConcreteLidContactDose(decay_days=-1.0)
    with pytest.raises(ValueError):
        cld.ConcreteLidContactDose(irradiation_days=0.0)


def test_write_inputs_creates_case_files(calc):
    calc.write_inputs()
    case = os.path.join(calc.cwd, calc.case_dir)
    assert os.path.isfile(os.path.join(case, calc.ORIGEN_input_file_name))
    for region in calc.regions.values():
        assert os.path.isfile(os.path.join(case, region['adens_file']))


def test_write_inputs_stages_f33_by_basename(tmp_path, monkeypatch):
    """ The F33 lands in the case directory and the deck copies and names it by basename """
    monkeypatch.chdir(tmp_path)
    os.mkdir('libs')
    with open(os.path.join('libs', 'cavity.f33'), 'wb') as f:
        f.write(b'cavity spectrum')
    calc = cld.ConcreteLidContactDose(irradiation_lib_f33=os.path.join('libs', 'cavity.f33'))
    calc.write_inputs()
    case = os.path.join(calc.cwd, calc.case_dir)
    with open(os.path.join(case, 'cavity.f33'), 'rb') as f:
        assert f.read() == b'cavity spectrum'
    with open(os.path.join(case, calc.ORIGEN_input_file_name)) as f:
        deck = f.read()
    shell = deck.split('=origen')[0]
    assert 'cp -r ${INPDIR}/cavity.f33 .' in shell
    assert deck.count('file="cavity.f33"') == 2
    assert 'libs/' not in deck


def test_write_inputs_missing_f33_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calc = cld.ConcreteLidContactDose(irradiation_lib_f33='absent.f33')
    with pytest.raises(FileNotFoundError):
        calc.write_inputs()


def test_full_pipeline_with_mocked_scale(calc, monkeypatch):
    fake_decay = {'co-60': 1.0e-10, 'eu-154': 5.0e-12}
    runs: list[tuple[str, int]] = []
    read_positions: list[int] = []

    def fake_run(deck_file, nmpi=1, context_dir=''):
        runs.append((os.path.join(os.getcwd(), deck_file), nmpi))
        return True

    def fake_adens(f71, pos):
        read_positions.append(pos)
        return dict(fake_decay)

    monkeypatch.setattr(cld, 'run_scale_or_raise', fake_run)
    monkeypatch.setattr(cld, 'get_burned_nuclide_atom_dens', fake_adens)

    calc.run_activation(nmpi=3)
    assert set(calc.decayed_adens) == {'bottom_slab', 'mixed_layer'}
    assert read_positions == [calc.f71_position, calc.f71_position]
    assert os.path.isfile(os.path.join(calc.cwd, calc.case_dir, 'spectrum.f33'))

    for region in calc.regions.values():
        with open(os.path.join(calc.cwd, calc.case_dir, region['f71']), 'w') as f:
            f.write(f"fake f71 {region['f71']}")
    calc.run_mavric(nmpi=5)

    case_mavric = os.path.join(calc.cwd, calc.case_dir_mavric)
    assert [nmpi for _, nmpi in runs] == [3, 5]
    assert os.path.isfile(os.path.join(case_mavric, calc.MAVRIC_input_file_name))
    for region in calc.regions.values():
        with open(os.path.join(case_mavric, region['f71'])) as f:
            assert f.read() == f"fake f71 {region['f71']}"
    with open(os.path.join(case_mavric, 'lid_mixed_adens_mavric.inp')) as f:
        comp_line = f.readline()
        assert comp_line.startswith('co-60 1 ')

    with open(os.path.join(case_mavric, calc.MAVRIC_out_file_name), 'w') as f:
        f.write('junk\n Final Tally Results Summary\n\n'
                ' Photon Point Detector 2.  photon contact dose at lid top center\n'
                '    tally/quantity        value       deviation    uncert   (/min)   1 2 3 4 5 6\n'
                '    total flux          1.09857E+06  9.94570E+03  0.00905  7.37E+03  X - X - X -\n'
                '    response 2          1.23385E+00  8.82139E-03  0.00715  1.18E+04  X - X - X -\n')

    calc.get_responses()
    dose = calc.contact_dose
    assert dose['value'] == pytest.approx(1.23385)
    assert dose['stdev'] == pytest.approx(8.82139e-3)  # absolute rem/h, not the 0.00715 relative column


_SCRIPT = os.path.join(os.path.dirname(__file__), '..', 'concrete_irrad', 'concrete_lid_contact_dose.py')


def _run_script(tmp_path, *args):
    env = dict(os.environ, PYTHONPATH=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    return subprocess.run([sys.executable, _SCRIPT, *args], cwd=str(tmp_path), env=env,
                          capture_output=True, text=True, timeout=120)


def test_script_preview_labels_placeholder_library(tmp_path):
    result = _run_script(tmp_path)
    assert result.returncode == 0, result.stderr
    assert 'cavity_spectrum.f33 (placeholder; pass --lib for real runs)' in result.stdout


def test_script_run_without_library_is_rejected(tmp_path):
    result = _run_script(tmp_path, '--run')
    assert result.returncode == 2
    assert '--run requires --lib' in result.stderr


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
