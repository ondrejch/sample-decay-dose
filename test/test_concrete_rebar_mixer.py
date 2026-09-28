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
    """ Two bar segments per cell minus one Steinmetz crossing volume 2 d^3 / 3 """
    mixer = make_default_mixer()
    d, sx, sy = 1.59, 30.48, 15.24
    expected = (math.pi * d ** 2 * (sx + sy) / 4.0 - 2.0 * d ** 3 / 3.0) / (sx * sy * d)
    assert mixer.steel_volume_fraction() == pytest.approx(expected, rel=1e-12)


def test_steel_volume_fraction_hand_value():
    """ d = 2 cm, 10 cm square grid: pi/10 - 8/300 = 0.2874925987 """
    mixer = crm.ConcreteRebarMixer(rebar_diameter_cm=2.0, spacing_x_cm=10.0, spacing_y_cm=10.0)
    assert mixer.steel_volume_fraction() == pytest.approx(0.2874925987, rel=1e-9)


def test_steel_volume_fraction_removes_one_crossing_at_defaults():
    """ The crossing overlap is 2 d^2 / (3 sx sy), about 2.95% of the uncorrected fraction """
    d, sx, sy = 1.59, 30.48, 15.24
    uncorrected = math.pi * d * (sx + sy) / (4.0 * sx * sy)
    overlap = 2.0 * d ** 2 / (3.0 * sx * sy)
    assert make_default_mixer().steel_volume_fraction() == pytest.approx(uncorrected - overlap, rel=1e-12)
    assert overlap / uncorrected == pytest.approx(0.0295, abs=5e-4)


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
    concrete_fe = (1.0 - f) * 2.30 * crm.CONCRETE_WT_FRACTIONS['fe'] * (1.0 - concrete_imp_total) \
        / m_fe * 6.02214076e23 / 1e24
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


def test_steel_co59_outside_physical_limits_raises():
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=0.02)
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=-1e-4)


@pytest.mark.parametrize('co_fraction', [93e-6, 151e-6, 0.005])
def test_steel_co59_outside_design_band_accepted_with_note(co_fraction):
    """ Measured carbon-steel Co (93--151 ppm) and values above the band are accepted """
    with pytest.warns(UserWarning, match='design band'):
        mixer = crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=co_fraction)
    assert mixer.steel_co59_wt_fraction == pytest.approx(co_fraction)


def test_steel_co59_configurable_within_band(recwarn):
    mixer = crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=0.004)
    assert mixer.steel_co59_wt_fraction == pytest.approx(0.004)
    assert 'co' in {k.lower() for k in mixer.steel_wt_fractions}
    assert not [w for w in recwarn.list if 'design band' in str(w.message)]


def test_concrete_composition_is_pnnl_15870_material_96():
    """ PNNL-15870 Rev. 1 #96 weight fractions, normalized from their tabulated sum 0.999993 """
    tabulated = {'h': 0.005558, 'o': 0.498457, 'na': 0.017125, 'mg': 0.002565, 'al': 0.045746,
                 'si': 0.315092, 's': 0.001283, 'k': 0.019231, 'ca': 0.082705, 'fe': 0.012231}
    assert set(crm.CONCRETE_WT_FRACTIONS) == set(tabulated)
    for symbol, fraction in tabulated.items():
        assert crm.CONCRETE_WT_FRACTIONS[symbol] == pytest.approx(fraction / 0.999993, rel=1e-12)
    assert crm.CONCRETE_DENSITY == pytest.approx(2.30)


def test_custom_concrete_impurities_passed_through():
    """ Layer Eu weight fraction equals the concrete mass share times the 2 ppm input """
    impurities = {'co': 20e-6, 'eu': 2e-6}
    mixer = crm.ConcreteRebarMixer(1.59, 30.48, 15.24, concrete_impurities_wt=impurities)
    f = mixer.steel_volume_fraction()
    concrete_share = (1.0 - f) * crm.CONCRETE_DENSITY / mixer.mixture_density()
    assert mixer.get_weight_fractions()['Eu'] == pytest.approx(concrete_share * 2e-6, rel=1e-9)
    assert mixer.get_elemental_atom_densities()['Eu'] > 0


def test_atom_densities_reproduce_mixture_density():
    """ Summing nuclide atom densities times nuclide masses recovers the volume-averaged density """
    mixer = make_default_mixer()
    grams_per_cc = 0.0
    for nuclide, density in mixer.get_atom_densities().items():
        symbol, mass_number = nuclide.split('-')
        grams_per_cc += density * 1e24 / 6.02214076e23 \
            * crm.ConcreteRebarMixer.ISOTOPIC_DATA[symbol][int(mass_number)]['mass']
    assert grams_per_cc == pytest.approx(mixer.mixture_density(), rel=1e-9)


def test_zero_spacing_raises():
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 0.0, 15.24)


def test_bad_weight_fractions_raise():
    with pytest.raises(ValueError):
        crm.ConcreteRebarMixer(1.59, 30.48, 15.24, steel_wt_fractions={'fe': 2.0})


def test_impossible_geometry_raises():
    """ A diameter larger than the bar spacing makes parallel bars overlap """
    mixer = crm.ConcreteRebarMixer(10.0, 9.0, 12.0)
    with pytest.raises(ValueError):
        mixer.steel_volume_fraction()


def test_touching_bars_stay_below_unity():
    """ At d = spacing the fraction is pi/2 - 2/3 """
    mixer = crm.ConcreteRebarMixer(10.0, 10.0, 10.0)
    assert mixer.steel_volume_fraction() == pytest.approx(math.pi / 2.0 - 2.0 / 3.0, rel=1e-12)
