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
    pure_salt_mass = m_u_per_mole + m_lif_per_mole + m_bef2_per_mole + m_f_in_uf4_per_mole

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
