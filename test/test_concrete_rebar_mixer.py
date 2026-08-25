""" Unit tests for the concrete-rebar mixed-layer homogenizer """
import importlib.util
import math
import os

import pytest

_SPEC = importlib.util.spec_from_file_location('concrete_rebar_mixer', os.path.join(
    os.path.dirname(__file__), '..', 'concrete_irrad', 'concrete_rebar_mixer.py'))
assert _SPEC is not None and _SPEC.loader is not None
crm = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(crm)


def make_default_mixer():
    return crm.ConcreteRebarMixer(rebar_diameter_cm=1.59, spacing_x_cm=30.48, spacing_y_cm=15.24)


def test_steel_volume_fraction_matches_unit_cell_geometry():
    mixer = make_default_mixer()
    d, sx, sy = 1.59, 30.48, 15.24
    expected = math.pi * d * (sx + sy) / (4.0 * sx * sy)
    assert mixer.steel_volume_fraction() == pytest.approx(expected)


def test_mixture_density_is_volume_weighted():
    mixer = make_default_mixer()
    f = mixer.steel_volume_fraction()
    expected = f * crm.REBAR_STEEL_DENSITY + (1.0 - f) * crm.CONCRETE_DENSITY
    assert mixer.mixture_density() == pytest.approx(expected)


def test_atom_densities_iron_conservation():
    """ Iron number density must match a hand calculation from both constituents """
    mixer = make_default_mixer()
    f = mixer.steel_volume_fraction()
    m_fe = sum(d['mass'] * d['abundance'] for d in crm.ConcreteRebarMixer.ISOTOPIC_DATA['fe'].values()
               if d['abundance'] > 0)
    steel_co = mixer.steel_co59_wt_fraction
    concrete_imp_total = sum(mixer.concrete_impurities_wt.values())
    steel_fe = f * 7.85 * 0.98 * (1.0 - steel_co) / m_fe * 6.02214076e23 / 1e24
    concrete_fe = (1.0 - f) * 2.30 * 0.0105 * (1.0 - concrete_imp_total) / m_fe * 6.02214076e23 / 1e24
    adens = mixer.get_atom_densities()
    assert adens['fe-56'] == pytest.approx(0.91754 * (steel_fe + concrete_fe), rel=1e-9)


def test_atom_densities_sum_over_isotopes_matches_elemental():
    mixer = make_default_mixer()
    adens = mixer.get_elemental_atom_densities()
    nuclides = mixer.get_atom_densities()
    for symbol in ('H', 'O', 'Si', 'Ca', 'Fe', 'Mn', 'Co', 'Eu'):
        isotope_sum = sum(v for k, v in nuclides.items() if k.startswith(f"{symbol.lower()}-"))
        assert isotope_sum == pytest.approx(adens[symbol], rel=1e-12)


def test_activation_nuclides_present():
    """ Co-59, Eu-151, and Eu-153 must appear in the mixed-layer atom densities """
    nuclides = make_default_mixer().get_atom_densities()
    assert nuclides['co-59'] > 0
    assert nuclides['eu-151'] > 0 and nuclides['eu-153'] > 0


def test_steel_co59_weight_fraction_is_preserved():
    """ Layer Co equals steel-mass share times steel Co-59 plus concrete trace Co """
    mixer = make_default_mixer()
    f = mixer.steel_volume_fraction()
    steel_term = f * crm.REBAR_STEEL_DENSITY * mixer.steel_co59_wt_fraction
    concrete_term = (1.0 - f) * crm.CONCRETE_DENSITY * mixer.concrete_impurities_wt['co']
    assert mixer.get_weight_fractions()['Co'] == pytest.approx(
        (steel_term + concrete_term) / mixer.mixture_density(), rel=1e-9)


def test_steel_co59_out_of_band_raises():
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=0.005)
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=0.0005)


def test_steel_co59_configurable_within_band():
    mixer = crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=0.004)
    assert mixer.steel_co59_wt_fraction == pytest.approx(0.004)
    assert 'co' in {k.lower() for k in mixer.steel_wt_fractions}


def test_custom_concrete_impurities_passed_through():
    impurities = {'co': 20e-6, 'eu': 2e-6}
    mixer = crm.ConcreteRebarMixer(1.59, 30.48, 15.24, concrete_impurities_wt=impurities)
    assert sum(mixer.get_weight_fractions().values()) == pytest.approx(1.0)
    elemental = mixer.get_elemental_atom_densities()
    assert elemental['Eu'] > 0


def test_weight_fractions_sum_to_one():
    assert sum(make_default_mixer().get_weight_fractions().values()) == pytest.approx(1.0)


def test_zero_spacing_raises():
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 0.0, 15.24)


def test_bad_weight_fractions_raise():
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_wt_fractions={'fe': 2.0})


def test_impossible_geometry_raises():
    mixer = crm.ConcreteRebarMixer(10.0, 11.0, 11.0)
    with pytest.raises(ValueError):
        mixer.steel_volume_fraction()
