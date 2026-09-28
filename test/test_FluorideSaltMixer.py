import pytest
from sample_decay_dose.FluorideSaltMixer import FluorideSalt, FlibeSalt, FlibeUF4Salt, FlibeUF4AdmixtureSalt


def test_fluoridesalt_basic_init():
    salt = FluorideSalt("70%LiF-30%BeF2")
    assert salt.components['LiF'] == pytest.approx(70.0)
    assert salt.components['BeF2'] == pytest.approx(30.0)


def test_fluoridesalt_xxx_composition():
    """Tests the XXX% feature to automatically calculate the remainder."""
    salt = FluorideSalt("12%NaF-XXX%ZrF4")
    assert salt.components['NaF'] == pytest.approx(12.0)
    assert salt.components['ZrF4'] == pytest.approx(88.0)


def test_fluoridesalt_generic_formula_parsing():
    """Tests parsing of a valid chemical formula that isn't in the known salts list."""
    # FeF2 is not in SALT_COEFFICIENTS, but should be parsed via regex
    salt = FluorideSalt("50%LiF-50%FeF2")
    assert salt.components['FeF2'] == pytest.approx(50.0)
    # Ensure it calculated a molar mass for the generic salt
    assert salt.molar_masses['FeF2'] > 0


def test_fluoridesalt_lanthanide_fallback():
    """
    Tests that using a lanthanide not in SALT_COEFFICIENTS (e.g., EuF3)
    triggers the LaF3 molar volume fallback without crashing.
    """
    salt = FluorideSalt("90%LiF-10%EuF3")
    rho = salt.density(1000)
    assert rho > 0
    # Sanity check: approximate density range for molten salts (1.5 - 5 g/cm3)
    assert 1.5 < rho < 5.0


def test_fluoridesalt_elemental_impurity():
    """Tests if the base class correctly handles an elemental impurity with natural abundance."""
    salt = FluorideSalt("100%LiF", impurities_wt={'Al': 1e-5})
    densities = salt.get_atom_densities(900)
    assert 'Al-27' in densities
    assert densities['Al-27'] > 0


def test_fluoridesalt_uf_ratio():
    salt_no_ratio = FluorideSalt("20%UF4-80%LiF")
    salt_with_ratio = FluorideSalt("20%UF4-80%LiF", uf3_to_uf4_ratio=0.1)
    densities_no_ratio = salt_no_ratio.get_atom_densities(950)
    densities_with_ratio = salt_with_ratio.get_atom_densities(950)

    # F-19 should be reduced because some UF4 became UF3 (less Fluorine)
    assert densities_with_ratio['F-19'] < densities_no_ratio['F-19']
    # Uranium density should remain constant
    assert densities_with_ratio['U-238'] == pytest.approx(densities_no_ratio['U-238'])


def test_fluoridesalt_validations():
    """Tests various error conditions."""
    # Sum != 100% (Standard validation)
    # The error message is "Molar percentages must sum to 100, but sum to X"
    with pytest.raises(ValueError, match="Molar percentages must sum to 100"):
        FluorideSalt("60%LiF-50%BeF2")

    # XXX% calculation overflow (Parsing validation)
    # This triggers the specific "exceed 100%" check inside _parse_composition
    with pytest.raises(ValueError, match="exceed 100%"):
        FluorideSalt("60%LiF-50%BeF2-XXX%NaF")

    # Negative percentage / Malformed separator
    # Since '-' is the delimiter, starting with '-' creates an empty split result (''),
    # causing a parsing error before the numeric value is checked.
    with pytest.raises(ValueError, match="Could not parse component"):
        FluorideSalt("-10%LiF-110%BeF2")

    # Total impurity > 100%
    with pytest.raises(ValueError, match="Total weight fraction"):
        FluorideSalt("100%LiF", impurities_wt={'Fe': 1.1})


def test_from_atom_densities_reconstruction():
    """
    Integration test: Create a salt, get densities, and try to reconstruct the salt object.
    This validates the round-trip logic.
    """
    # 1. Create initial salt
    original_salt = FluorideSalt("67%LiF-33%BeF2")
    temp = 900
    densities = original_salt.get_atom_densities(temp)

    # 2. Reconstruct
    # We rely on ValencyMapper working correctly in the environment
    reconstructed_salt = FluorideSalt.from_atom_densities(densities)

    # 3. Compare Components
    # Note: Reconstruction might have small floating point deviations
    assert reconstructed_salt.components['LiF'] == pytest.approx(67.0, abs=0.1)
    assert reconstructed_salt.components['BeF2'] == pytest.approx(33.0, abs=0.1)


# --- Tests for the FlibeSalt child class ---

def test_flibesalt_init():
    salt = FlibeSalt("TestFlibe", temperature=923.0, Li7_enr=0.99, impurities_wt={'Fe-56': 1e-6})
    assert salt.name == "TestFlibe"
    assert salt.temperature == 923.0
    assert salt.components['LiF'] == pytest.approx(100.0 * 2.0 / 3.0, 5e-6)
    assert salt.components['BeF2'] == pytest.approx(100.0 * 1.0 / 3.0, 5e-6)
    assert 'Fe-56' in salt.impurities_wt


# --- Tests for the FlibeUF4Salt child class (Refactored) ---

def test_flibeuf4salt_init():
    salt = FlibeUF4Salt(name="TestFlibeUF4", temperature=950, u_enr_weight={'u-235': 0.2, 'u-238': 0.8},
        u_impurities_weight={'Al': 1e-5}, UF4=0.22, Li7_enr=0.999)
    assert salt.name == "TestFlibeUF4"
    assert salt.components['UF4'] == pytest.approx(22.0)
    assert 'U' in salt.enrichments
    assert 'Al' in salt.impurities_wt  # Check that the key is passed to the parent


def test_flibeuf4salt_direct_admixture():
    """
    Tests the NEW capability of FlibeUF4Salt to handle admixtures directly,
    bypassing the need for the wrapper class.
    """
    salt = FlibeUF4Salt(name="DirectAdmixture", temperature=900, u_enr_weight=None, u_impurities_weight=None, UF4=0.20,
        admixtures={'LuF3': {'mol_frac': 0.05}})
    # Check normalization:
    # UF4 = 20%, LuF3 = 5%. Remaining = 75%.
    # LiF = 75 * 2/3 = 50%
    # BeF2 = 75 * 1/3 = 25%
    assert salt.components['LuF3'] == pytest.approx(5.0)
    assert salt.components['UF4'] == pytest.approx(20.0)
    assert salt.components['LiF'] == pytest.approx(50.0)
    assert salt.components['BeF2'] == pytest.approx(25.0)


def test_flibeuf4salt_impurity_conversion():
    """Verify that impurities relative to U-mass are correctly converted to fractions of total salt mass."""
    li7_enrichment = 0.99995
    salt = FlibeUF4Salt(name="ImpurityTest", temperature=950, u_enr_weight={'u-235': 0.2, 'u-238': 0.8},
        u_impurities_weight={'Al': 275e-6},  # 275 ppm of Al in U
        UF4=0.22, Li7_enr=li7_enrichment)
    # Manually calculate the expected final weight fraction
    # 1. Avg mass of U
    avg_u_mass = 0.2 * 235.043 + 0.8 * 238.050
    # 2. Mass of U per mole of salt
    m_u_per_mole = 0.22 * avg_u_mass
    # 3. Mass of other components per mole of pure salt
    m_lif_per_mole = (1.0 - 0.22) * (2.0 / 3.0) * (li7_enrichment * 7.016 + (1.0 - li7_enrichment) * 6.015 + 18.998)
    m_bef2_per_mole = (1.0 - 0.22) * (1.0 / 3.0) * (9.012 + 2.0 * 18.998)
    m_f_in_uf4_per_mole = 0.22 * 4.0 * 18.998
    # The default UF3/UF4 ratio of 0.05 removes one F per U3+ from the salt mass.
    m_f_removed_by_uf3 = 0.22 * (0.05 / 1.05) * 18.998
    pure_salt_mass = m_u_per_mole + m_lif_per_mole + m_bef2_per_mole + m_f_in_uf4_per_mole - m_f_removed_by_uf3

    # 4. Absolute mass of impurity
    m_al_impurity = 275e-6 * m_u_per_mole

    # 5. Final weight fraction of impurity
    total_mass = pure_salt_mass + m_al_impurity
    expected_frac = m_al_impurity / total_mass

    assert salt.impurities_wt['Al'] == pytest.approx(expected_frac, 5e-5)
    # Check that the atom density is calculated
    densities = salt.get_atom_densities()
    assert 'Al-27' in densities and densities['Al-27'] > 0


def test_flibeuf4salt_uranium_enrichment_conversion():
    salt = FlibeUF4Salt("Test", 950, {'u-235': 0.2, 'u-238': 0.8}, {}, 0.2)
    u_enrich = salt.enrichments['U']
    # Check that mole fractions are not the same as weight fractions
    assert u_enrich[235] != 0.2
    # Heavier isotope should have a lower mole fraction for the same weight fraction
    assert u_enrich[235] > u_enrich[238] * (0.2 / 0.8)


# --- Tests for the FlibeUF4AdmixtureSalt child class (Wrapper) ---

@pytest.fixture
def admixture_salt():
    """This fixture now correctly initializes the admixture salt without error."""
    return FlibeUF4AdmixtureSalt(name="AdmixtureTest", temperature=950, u_enr_weight={'u-235': 0.2, 'u-238': 0.8},
        u_impurities_weight={'Ni-58': 5e-6}, UF4=0.20,
        admixtures={'LuF3': {'mol_frac': 0.01, 'enrichment': {176: 1.0}}, 'YbF3': {'mol_frac': 0.02}}, Li7_enr=0.9999)


def test_admixturesalt_init(admixture_salt):
    assert admixture_salt.name == "AdmixtureTest"
    assert admixture_salt.components['UF4'] == pytest.approx(20.0)
    assert admixture_salt.components['LuF3'] == pytest.approx(1.0)
    assert admixture_salt.components['YbF3'] == pytest.approx(2.0)


def test_admixturesalt_enrichments(admixture_salt):
    enrich = admixture_salt.enrichments
    assert 'U' in enrich
    assert 'Lu' in enrich
    assert enrich['Lu'] == {176: 1.0}


def test_admixturesalt_densities(admixture_salt):
    densities = admixture_salt.get_atom_densities()
    assert 'Lu-176' in densities
    assert 'Lu-175' not in densities
    assert 'Yb-174' in densities
    assert 'U-235' in densities
    assert 'Ni-58' in densities
    assert densities['Ni-58'] > 0


# --- Review 2026-09-27 (A11, E) regression tests ---

N_A = 6.02214076e23


def _mass_density_from_atoms(atom_densities: dict) -> float:
    """Mass density [g/cm3] implied by atom densities [atoms/barn-cm]; isomers use the ground-state mass."""
    total = 0.0
    for name, dens in atom_densities.items():
        symbol, mass = name.split('-')
        total += dens * FluorideSalt.ISOTOPIC_DATA[symbol.lower()][int(mass.rstrip('m'))]['mass']
    return total * 1e24 / N_A


def _decayed_fuel_salt_densities() -> dict:
    salt = FlibeUF4Salt("f", 950, {'u-235': 0.1975, 'u-238': 0.8025}, {'Al': 275e-6}, 0.22,
                        UF3_to_UF4=0.08, Li7_enr=0.9999)
    return salt.get_atom_densities()


def test_lanthanide_density_independent_of_call_order():
    salt = FluorideSalt("90%LiF-10%CeF3")
    salt.density(900)
    assert salt.density(1200) == pytest.approx(FluorideSalt("90%LiF-10%CeF3").density(1200), rel=1e-12)
    assert salt.density(900) > salt.density(1200)


def test_admixture_single_letter_cation_enrichment():
    """'KF' used to become cation 'Kf' through re.IGNORECASE."""
    salt = FlibeUF4Salt("K", 950, None, None, 0.2, admixtures={'KF': {'mol_frac': 0.05, 'enrichment': {39: 0.5, 41: 0.5}}})
    assert salt.enrichments['K'] == {39: 0.5, 41: 0.5}
    densities = salt.get_atom_densities()
    assert densities['K-41'] == pytest.approx(densities['K-39'])
    assert 'K-40' not in densities


def test_zero_abundance_cation_requires_enrichment():
    with pytest.raises(ValueError, match="no natural isotopic abundance"):
        FluorideSalt("95%LiF-5%PmF3")
    salt = FluorideSalt("95%LiF-5%PmF3", enrichments={'Pm': {147: 1.0}})
    assert salt.molar_masses['PmF3'] > 200


def test_enrichment_mole_fraction_form():
    salt = FluorideSalt("78%LiF-22%UF4", enrichments={'U': {235: 0.05, 238: 0.95}})
    assert salt.enrichments['U'] == {235: 0.05, 238: 0.95}
    with pytest.raises(ValueError, match="List every isotope"):
        FluorideSalt("78%LiF-22%UF4", enrichments={'U': {235: 0.05}})


def test_uf3_ratio_is_mass_consistent():
    salt = FluorideSalt("78%LiF-22%UF4", uf3_to_uf4_ratio=0.08)
    densities = salt.get_atom_densities(950)
    assert _mass_density_from_atoms(densities) == pytest.approx(salt.density(950), rel=1e-12)
    # Constant molar volume: the removed fluorine lowers the mass density.
    assert salt.density(950) < FluorideSalt("78%LiF-22%UF4").density(950)
    with pytest.raises(ValueError, match="uf3_to_uf4_ratio"):
        FluorideSalt("78%LiF-22%UF4", uf3_to_uf4_ratio=-0.1)


def test_from_atom_densities_impurity_keeps_uf3_ratio_and_round_trips():
    densities = _decayed_fuel_salt_densities()
    salt = FluorideSalt.from_atom_densities(densities)
    assert salt.uf3_to_uf4_ratio == pytest.approx(0.08, rel=1e-12)
    assert set(salt.components) == {'LiF', 'BeF2', 'UF4'}
    assert 'Al-27' in salt.impurities_wt
    rebuilt = salt.get_atom_densities(950)
    assert set(rebuilt) == set(densities)
    for name, dens in densities.items():
        assert rebuilt[name] == pytest.approx(dens, rel=1e-12)
    assert salt.excluded_species == {}


def test_from_atom_densities_drops_zero_valence_with_report():
    densities = _decayed_fuel_salt_densities()
    densities.update({'xe-135': 1e-6, 'Mo-99': 2e-6})
    with pytest.warns(UserWarning, match="Xe"):
        salt = FluorideSalt.from_atom_densities(densities)
    assert not any(formula.startswith(('Xe', 'Mo')) for formula in salt.components)
    assert not any(key.startswith(('Xe', 'Mo')) for key in salt.impurities_wt)
    total = sum(densities.values())
    assert salt.excluded_species['Xe']['atom_fraction'] == pytest.approx(1e-6 / total)
    assert salt.excluded_species['Mo']['nuclides'] == {'Mo-99': pytest.approx(2e-6 / total)}
    assert salt.uf3_to_uf4_ratio == pytest.approx(0.08, rel=1e-9)
    kept = FluorideSalt.from_atom_densities(densities, drop_zero_valence=False)
    assert 'Mo-99' in kept.impurities_wt and kept.excluded_species == {}


def test_from_atom_densities_trace_anions_kept_major_anion_raises():
    densities = _decayed_fuel_salt_densities()
    densities['I-131'] = 1e-7
    salt = FluorideSalt.from_atom_densities(densities)
    assert 'I-131' in salt.impurities_wt
    assert salt.get_atom_densities(950)['I-131'] == pytest.approx(1e-7, rel=1e-9)
    densities['Cl-35'] = 0.2 * densities['Li-7']
    with pytest.raises(ValueError, match="anion"):
        FluorideSalt.from_atom_densities(densities)


def test_from_atom_densities_keeps_isomers_distinct():
    densities = _decayed_fuel_salt_densities()
    densities.update({'am-242': 2e-9, 'am-242m': 1e-9})
    salt = FluorideSalt.from_atom_densities(densities)
    rebuilt = salt.get_atom_densities(950)
    assert rebuilt['Am-242'] == pytest.approx(2e-9, rel=1e-9)
    assert rebuilt['Am-242m'] == pytest.approx(1e-9, rel=1e-9)
    # As a salt component, the isomer gets its own enrichment key.
    major = FluorideSalt.from_atom_densities(densities, impurity_cutoff=0.0)
    assert major.enrichments['Am'] == {242: pytest.approx(2 / 3), '242m': pytest.approx(1 / 3)}
    assert 'AmF3' in major.components


def test_from_atom_densities_random_salts_sum_to_100():
    """Six-decimal rounding failed the 1e-6 sum check for some random salts."""
    import random
    rng = random.Random(20260927)
    formulas = ['LiF', 'BeF2', 'NaF', 'KF', 'ZrF4', 'ThF4', 'UF4']
    for _ in range(300):
        chosen = rng.sample(formulas, rng.randint(2, len(formulas)))
        weights = [rng.random() for _ in chosen]
        total = sum(weights)
        composition = "-".join(f"{100 * w / total!r}%{f}" for f, w in zip(chosen, weights))
        ratio = rng.uniform(0.0, 0.2) if 'UF4' in chosen else None
        original = FluorideSalt(composition, uf3_to_uf4_ratio=ratio)
        salt = FluorideSalt.from_atom_densities(original.get_atom_densities(900), impurity_cutoff=0.0)
        assert sum(salt.components.values()) == pytest.approx(100.0, abs=1e-9)
        for formula, pct in original.components.items():
            assert salt.components[formula] == pytest.approx(pct, rel=1e-9)
        if ratio:
            assert salt.uf3_to_uf4_ratio == pytest.approx(ratio, rel=1e-6, abs=1e-12)


def test_from_atom_densities_excess_fluorine_kept_as_impurity():
    densities = FluorideSalt("78%LiF-22%UF4").get_atom_densities(950)
    # Remove 1% of the uranium, as fission does; its fluorine stays in the salt.
    for name in [n for n in densities if n.startswith('U-')]:
        densities[name] *= 0.99
    salt = FluorideSalt.from_atom_densities(densities)
    assert salt.uf3_to_uf4_ratio is None
    assert salt.excess_fluorine_atom_fraction > 0
    assert 'F-19' in salt.impurities_wt
    rebuilt = salt.get_atom_densities(950)
    assert rebuilt['F-19'] / rebuilt['Li-7'] == pytest.approx(densities['F-19'] / densities['Li-7'], rel=1e-12)


def test_from_atom_densities_actionable_errors():
    densities = FluorideSalt("78%LiF-22%UF4").get_atom_densities(950)
    with pytest.raises(ValueError, match="Could not parse nuclide"):
        FluorideSalt.from_atom_densities({**densities, 'bogus': 1.0})
    with pytest.raises(ValueError, match="No fluorine"):
        FluorideSalt.from_atom_densities({k: v for k, v in densities.items() if not k.startswith('F-')})
    # Hf has no valence in the maps; above the cutoff it needs an override.
    with_hf = {**densities, 'Hf-180': 0.01 * densities['Li-7']}
    with pytest.raises(ValueError, match="valence_overrides"):
        FluorideSalt.from_atom_densities(with_hf)
    with_hf['F-19'] += 4 * with_hf['Hf-180']
    assert 'HfF4' in FluorideSalt.from_atom_densities(with_hf, valence_overrides={'Hf': 4}).components
    # A fluorine deficit without uranium cannot be expressed as UF3.
    flibe = FluorideSalt("67%LiF-33%BeF2").get_atom_densities(900)
    flibe['F-19'] *= 0.99
    with pytest.raises(ValueError, match="no UF4 component"):
        FluorideSalt.from_atom_densities(flibe)


def test_isomer_impurity_key_is_preserved():
    salt = FluorideSalt("100%LiF", impurities_wt={'am-242m': 1e-6, 'Ac-227': 1e-9})
    assert set(salt.impurities_wt) == {'Am-242m', 'Ac-227'}
    densities = salt.get_atom_densities(900)
    assert densities['Am-242m'] > 0 and 'Am-242' not in densities
    with pytest.raises(ValueError, match="by nuclide"):
        FluorideSalt("100%LiF", impurities_wt={'Pm': 1e-6}).get_atom_densities(900)
