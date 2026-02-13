import re
from collections import defaultdict
from sample_decay_dose.data import ISOTOPIC_DATA as _isotopic_data
from sample_decay_dose.ValencyMapper import ValencyMapper


class FluorideSalt:
    """
    Represents a mixed fluoride salt, calculating its density and isotopic atom densities.
    - Molar masses are calculated from isotopic data.
    - Supports custom isotope enrichments and impurities by weight fraction.
    - Allows 'XXX%' to specify a component's fraction as the remainder to 100%.
    - Can adjust fluorine content based on a specified UF3/UF4 ratio.
    """

    SALT_COEFFICIENTS = {"LiF": {"A": 2.37, "B": 5e-4}, "BeF2": {"A": 1.97, "B": 1.45e-5},
        "UF4": {"A": 7.78, "B": 9.92e-4}, "LaF3": {"A": 5.79, "B": 6.82e-4}, "KF": {"A": 2.64, "B": 6.57e-4},
        "NaF": {"A": 2.70, "B": 5.90e-4}, "SrF2": {"A": 4.78, "B": 7.51e-4}, "ThF4": {"A": 7.11, "B": 7.59e-4},
        "ZrF4": {"A": 5.36, "B": 1.23e-3}}
    AVOGADRO_NUMBER = 6.02214076e23  # atoms/mol
    BARN_CM_CONVERSION = 1e24  # (cm^2 / barn)
    LANTHANIDE_SYMBOLS = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
    ISOTOPIC_DATA = _isotopic_data

    def __init__(self, composition_str: str, enrichments: dict = None, impurities_wt: dict = None,
                 uf3_to_uf4_ratio: float = None):
        """
        Initializes the FluorideSalt object.

        Args:
            composition_str: String defining the molar composition of the salt (e.g., "78%LiF-22%UF4").
            enrichments: Dictionary defining custom isotopic enrichments. Can be mole fractions (e.g., {'U': {235: 0.05}})
                         or special weight fraction keys (e.g., {'Li7_enr': 0.9999}).
            impurities_wt: Dictionary defining impurities by weight fraction of the *total salt mass*.
                           Example: {'Fe-56': 1e-5}.
            uf3_to_uf4_ratio: Molar ratio of U3+ to U4+ to set the fluorine potential.
        """
        # This mapping is the key to robust, case-insensitive salt parsing.
        self._known_salts = {s.lower(): s for s in self.SALT_COEFFICIENTS.keys()}
        for symbol in self.LANTHANIDE_SYMBOLS:
            formula = f"{symbol}F3"
            if formula.lower() not in self._known_salts:
                self._known_salts[formula.lower()] = formula

        # Process enrichments from user-friendly formats (like wt%) to internal mole%
        self.enrichments = self._process_enrichments(enrichments)

        processed_impurities = {}
        if impurities_wt:
            for key, val in impurities_wt.items():
                parts = re.match(r"([a-zA-Z]+)-?(\d*)", key)
                if parts:
                    symbol, mass = parts.groups()
                    new_key = f"{symbol.capitalize()}-{mass}" if mass else symbol.capitalize()
                    processed_impurities[new_key] = val
                else:
                    processed_impurities[key] = val
        self.impurities_wt = processed_impurities

        self.uf3_to_uf4_ratio = uf3_to_uf4_ratio

        self.total_impurity_frac = sum(self.impurities_wt.values())
        if self.total_impurity_frac >= 1.0:
            raise ValueError("Total weight fraction of impurities must be less than 1.")

        self.components = self._parse_composition(composition_str)
        self.elemental_masses = self._calculate_elemental_masses()
        self.molar_masses = {}
        self._validate_and_calculate_masses()
        self._ref_molar_volume = None  # To cache the reference molar volume

    def _process_enrichments(self, user_enrichments: dict) -> dict:
        """
        Processes the user-provided enrichment dictionary to convert special cases
        (like weight fractions) into the internal mole fraction format.
        """
        if user_enrichments is None:
            return {}

        processed = {str(k).capitalize(): v for k, v in user_enrichments.items()}

        # Handle Li-7 weight fraction enrichment
        if 'Li7_enr' in processed:
            w7 = processed.pop('Li7_enr')
            if not (0 <= w7 <= 1):
                raise ValueError("Li-7 weight fraction must be between 0 and 1.")
            m7 = self.ISOTOPIC_DATA['li'][7]['mass']
            m6 = self.ISOTOPIC_DATA['li'][6]['mass']
            w6 = 1.0 - w7
            n7 = w7 / m7
            n6 = w6 / m6
            total_li_moles = n7 + n6
            x7 = n7 / total_li_moles
            x6 = n6 / total_li_moles
            processed['Li'] = {7: x7, 6: x6}

        # Handle Uranium weight fraction enrichment
        if 'U' in processed and isinstance(processed['U'], dict):
            # Check if it's in the user-friendly {'u-235': wt_frac} format
            if any(isinstance(k, str) for k in processed['U'].keys()):
                u_enr_weight = processed.pop('U')
                u_enr_mole = {}
                for key, wt_frac in u_enr_weight.items():
                    match = re.match(r"[Uu]-(\d+)", key)
                    if not match: raise ValueError(f"Invalid key in u_enr_weight: '{key}'")
                    mass_num = int(match.group(1))
                    isotope_mass = self.ISOTOPIC_DATA['u'][mass_num]['mass']
                    u_enr_mole[mass_num] = wt_frac / isotope_mass
                total_moles = sum(u_enr_mole.values())
                if total_moles > 0:
                    processed['U'] = {mass: mol / total_moles for mass, mol in u_enr_mole.items()}

        return processed

    @classmethod
    def from_atom_densities(cls, atom_densities: dict, valency_estimate_type: str = 'upper'):
        """
        Reconstructs the salt's molar composition from a dictionary of nuclide atom densities.

        Args:
            atom_densities (dict): A dictionary mapping nuclide names (e.g., 'Li-7') to their
                                   atom densities in atoms/barn-cm.
            valency_estimate_type (str): The valency map to use ('upper', 'lower', 'doligez'). Defaults to 'upper'.

        Returns:
            A new FluorideSalt instance.
        """
        elemental_densities = defaultdict(float)
        nuclide_breakdown = defaultdict(dict)
        valence_mapper = ValencyMapper(estimate_type=valency_estimate_type)

        for nuclide, density in atom_densities.items():
            match = re.match(r"([A-Z][a-z]*)-(\d+)", nuclide)
            if not match:
                raise ValueError(f"Could not parse nuclide name: '{nuclide}'")
            symbol, mass_num = match.groups()
            mass_num = int(mass_num)

            elemental_densities[symbol] += density
            nuclide_breakdown[symbol][mass_num] = density

        # Reconstruct enrichments (mole fractions) from densities
        enrichments = {}
        for symbol, breakdowns in nuclide_breakdown.items():
            if symbol == 'F':
                continue
            total_density = elemental_densities[symbol]
            if total_density > 0:
                enrichments[symbol] = {mass: dens / total_density for mass, dens in breakdowns.items()}

        salt_components_moles = {}
        fluorine_consumed = 0

        for symbol, density in elemental_densities.items():
            if symbol == 'F':
                continue

            valence = valence_mapper.get(symbol)
            if valence is None or valence < 0:
                raise ValueError(
                    f"Cannot determine salt component for element '{symbol}'. Valence is unknown or negative.")

            formula = f"{symbol}F{valence if valence > 1 else ''}"
            salt_components_moles[formula] = density
            fluorine_consumed += density * valence

        # Check for fluorine discrepancy, which implies a mixed valence state (e.g., U+3/U+4)
        total_fluorine_density = elemental_densities.get('F', 0)
        fluorine_discrepancy = fluorine_consumed - total_fluorine_density
        uf3_to_uf4_ratio = None

        if abs(fluorine_discrepancy) > 1e-9 and 'U' in elemental_densities:
            # The reduction in F is equal to the amount of U+3.
            moles_u3 = fluorine_discrepancy
            moles_u4 = elemental_densities['U'] - moles_u3
            if moles_u4 > 0:
                uf3_to_uf4_ratio = moles_u3 / moles_u4

        # Build the composition string
        total_salt_moles = sum(salt_components_moles.values())
        if total_salt_moles == 0:
            raise ValueError("Cannot reconstruct salt from zero cation densities.")

        composition_str = "-".join(
            [f"{(moles / total_salt_moles) * 100:.6f}%{formula}"
             for formula, moles in salt_components_moles.items()]
        )

        # Impurities are not reconstructed, as their original definition (wt% of U or total) is lost.
        return cls(composition_str, enrichments=enrichments, uf3_to_uf4_ratio=uf3_to_uf4_ratio)

    def _parse_composition(self, s: str) -> dict:
        # This method now uses the _known_salts mapping for robust parsing.
        components = {}
        parts = s.split('-')
        xxx_component_salt = None
        total_specified_percentage = 0.0
        for part in parts:
            match = re.match(r"((?:\d+(?:\.\d+)?)|XXX)%([a-zA-Z0-9]+)", part, re.IGNORECASE)
            if not match:
                raise ValueError(f"Could not parse component: '{part}'")
            percentage_str, salt_input = match.group(1), match.group(2)

            # Use the known salts mapping to find the canonical formula.
            canonical_salt = self._known_salts.get(salt_input.lower())
            # If not a known salt, try to parse it as a generic formula. This allows
            # reconstruction from atom densities where density data may not exist.
            if canonical_salt is None:
                # Attempt to parse as a valid chemical formula (e.g., FeF2)
                parsed_elements = re.findall(r'([A-Z][a-z]?)(\d*)', salt_input)
                if not parsed_elements or "".join(f"{s}{c}" for s, c in parsed_elements) != salt_input:
                    raise ValueError(f"Unsupported or malformed salt component: '{salt_input}'")
                # Use the user's input, but with the first letter capitalized as a convention
                canonical_salt = salt_input[0].upper() + salt_input[1:]

            if percentage_str.upper() == "XXX":
                if xxx_component_salt is not None:
                    raise ValueError("Only one 'XXX%' component allowed.")
                xxx_component_salt = canonical_salt
                components[canonical_salt] = 0
            else:
                percentage = float(percentage_str)
                if percentage < 0:
                    raise ValueError(f"Percentage cannot be negative: {part}")
                components[canonical_salt] = percentage
                total_specified_percentage += percentage
        if xxx_component_salt:
            remainder = 100.0 - total_specified_percentage
            if remainder < -1e-9:
                raise ValueError(f"Specified percentages ({total_specified_percentage}%) exceed 100%.")
            components[xxx_component_salt] = remainder
        return components

    def _calculate_elemental_masses(self) -> dict:
        elemental_masses = {}
        all_elements = set()
        # Scan components, impurities, AND enrichments to find all relevant elements.
        for salt_formula in self.components.keys():
            # The salt_formula is now guaranteed to be in canonical form.
            for symbol, count in re.findall(r'([A-Z][a-z]?)(\d*)', salt_formula):
                all_elements.add(symbol)
        for isotope_name in self.impurities_wt.keys():
            symbol = re.match(r"([A-Z][a-z]*)", isotope_name).group(1)
            all_elements.add(symbol)
        for symbol in self.enrichments.keys():
            all_elements.add(symbol)

        for symbol in all_elements:
            if symbol in elemental_masses:
                continue

            iso_dist = {}
            if symbol in self.enrichments:
                custom_dist = self.enrichments[symbol]
                # The enrichment dict should have integer mass numbers as keys.
                if not all(isinstance(k, int) for k in custom_dist.keys()):
                    raise TypeError(f"Enrichment dictionary for {symbol} has non-integer keys: {custom_dist.keys()}")
                if abs(sum(custom_dist.values()) - 1.0) > 1e-6:
                    raise ValueError(f"Enrichment fractions for {symbol} must sum to 1.")
                iso_dist = custom_dist
            else:
                try:
                    iso_data = self.ISOTOPIC_DATA[symbol.lower()]
                    iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}
                except KeyError:
                    # Check if the element is part of a salt component before raising an error.
                    if any(symbol in f for f in self.components.keys()):
                        raise ValueError(f"Isotopic data not found for element: {symbol}")
                    else:
                        iso_dist = {}

            avg_mass = 0.0
            if not iso_dist and symbol.lower() in self.ISOTOPIC_DATA:
                elemental_masses[symbol] = 0
                continue

            for mass_num, fraction in iso_dist.items():
                try:
                    isotope_mass = self.ISOTOPIC_DATA[symbol.lower()][int(mass_num)]['mass']
                    avg_mass += isotope_mass * fraction
                except KeyError:
                    raise ValueError(f"Isotope {symbol}-{mass_num} not found in data.")
            elemental_masses[symbol] = avg_mass
        return elemental_masses

    def _calculate_molar_mass(self, formula: str) -> float:
        # This method is now simpler as it can assume formulas are canonical.
        total_mass = 0.0
        elements = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        if not elements: raise ValueError(f"Cannot parse formula: {formula}")
        for symbol, count_str in elements:
            count = int(count_str) if count_str else 1
            if symbol not in self.elemental_masses:
                raise ValueError(f"Elemental mass for '{symbol}' not found. Check isotopic data.")
            total_mass += self.elemental_masses[symbol] * count
        return total_mass

    def _validate_and_calculate_masses(self):
        # This method is now simpler as it can assume component keys are canonical.
        total_percentage = sum(self.components.values())
        if abs(total_percentage - 100.0) > 1e-6:
            raise ValueError(f"Molar percentages must sum to 100, but sum to {total_percentage}")
        for salt in self.components:
            # The salt key is already canonical from _parse_composition
            self.molar_masses[salt] = self._calculate_molar_mass(salt)

    def density(self, temperature: float) -> float:
        # This method now uses a molar volume fallback for unknown lanthanide salts.
        mixture_molar_mass = 0.0
        mixture_molar_volume = 0.0
        for salt, percentage in self.components.items():
            molar_fraction = percentage / 100.0
            molar_mass = self.molar_masses[salt]
            molar_volume_pure = 0.0

            if salt in self.SALT_COEFFICIENTS:
                coeffs = self.SALT_COEFFICIENTS[salt]
                density_pure = coeffs['A'] - coeffs['B'] * temperature
                if density_pure <= 0:
                    raise ValueError(f"Density of {salt} is non-positive at {temperature} K.")
                molar_volume_pure = molar_mass / density_pure
            else:
                # Fallback for lanthanides not in SALT_COEFFICIENTS, using LaF3 molar volume.
                if any(salt.startswith(s) for s in self.LANTHANIDE_SYMBOLS):
                    if self._ref_molar_volume is None:  # Calculate and cache if not already done
                        ref_salt = 'LaF3'
                        # Ensure the elemental mass for the reference salt's cation is available.
                        if 'La' not in self.elemental_masses:
                            la_data = self.ISOTOPIC_DATA['la']
                            avg_mass = sum(d['abundance'] * d['mass'] for d in la_data.values())
                            self.elemental_masses['La'] = avg_mass
                        if 'F' not in self.elemental_masses:  # Should always be present, but safe to check
                            f_data = self.ISOTOPIC_DATA['f']
                            avg_mass = sum(d['abundance'] * d['mass'] for d in f_data.values())
                            self.elemental_masses['F'] = avg_mass

                        ref_coeffs = self.SALT_COEFFICIENTS[ref_salt]
                        ref_molar_mass = self.molar_masses.get(ref_salt)
                        if ref_molar_mass is None:
                            ref_molar_mass = self._calculate_molar_mass(ref_salt)
                            self.molar_masses[ref_salt] = ref_molar_mass

                        ref_density = ref_coeffs['A'] - ref_coeffs['B'] * temperature
                        if ref_density <= 0:
                            raise ValueError(f"Reference density for {ref_salt} is non-positive at {temperature} K.")

                        self._ref_molar_volume = ref_molar_mass / ref_density

                    molar_volume_pure = self._ref_molar_volume
                else:
                    raise ValueError(f"Density coefficients not found for non-lanthanide salt: {salt}")

            mixture_molar_mass += molar_fraction * molar_mass
            mixture_molar_volume += molar_fraction * molar_volume_pure

        if mixture_molar_volume <= 0:
            raise ValueError("Mixture molar volume is non-positive.")

        pure_salt_density = mixture_molar_mass / mixture_molar_volume
        final_density = pure_salt_density / (
                    1.0 - self.total_impurity_frac) if self.total_impurity_frac < 1.0 else float('inf')

        return final_density

    def get_atom_densities(self, temperature: float) -> dict:
        # This method is updated to handle elemental impurities correctly.
        bulk_density = self.density(temperature)
        salt_weight_fraction = 1.0 - self.total_impurity_frac
        salt_mass_density = bulk_density * salt_weight_fraction

        mixture_molar_mass = sum((p / 100.0) * self.molar_masses[s] for s, p in self.components.items())
        total_salt_number_density = (salt_mass_density / mixture_molar_mass) * self.AVOGADRO_NUMBER \
            if mixture_molar_mass > 0 else 0

        isotope_densities = defaultdict(float)

        for salt_formula, salt_percentage in self.components.items():
            elements_in_formula = re.findall(r'([A-Z][a-z]?)(\d*)', salt_formula)
            for symbol, count_str in elements_in_formula:
                atoms_per_molecule = int(count_str) if count_str else 1
                element_total_density = total_salt_number_density * (salt_percentage / 100.0) * atoms_per_molecule

                iso_dist = self.enrichments.get(symbol, None)
                if iso_dist is None:
                    iso_data = self.ISOTOPIC_DATA[symbol.lower()]
                    iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}

                for mass_num, fraction in iso_dist.items():
                    if fraction > 0:
                        isotope_name = f"{symbol}-{mass_num}"
                        isotope_densities[isotope_name] += element_total_density * fraction

        if self.uf3_to_uf4_ratio is not None:
            total_uranium_density = sum(v for k, v in isotope_densities.items() if k.startswith('U-'))
            if total_uranium_density > 0:
                ratio = self.uf3_to_uf4_ratio
                frac_u3 = ratio / (1.0 + ratio)
                fluorine_density_reduction = total_uranium_density * frac_u3
                isotope_densities['F-19'] -= fluorine_density_reduction
                if isotope_densities['F-19'] < 0:
                    raise ValueError("Fluorine potential adjustment resulted in negative fluorine density.")

        for impurity_key, weight_frac in self.impurities_wt.items():
            impurity_mass_density = bulk_density * weight_frac
            match = re.fullmatch(r"([A-Z][a-z]*)-(\d+)", impurity_key)
            if match:
                symbol, mass_num_str = match.groups()
                mass_num = int(mass_num_str)
                try:
                    isotope_mass = self.ISOTOPIC_DATA[symbol.lower()][mass_num]['mass']
                    impurity_atom_density = (impurity_mass_density / isotope_mass) * self.AVOGADRO_NUMBER
                    isotope_densities[impurity_key] += impurity_atom_density
                except KeyError:
                    raise ValueError(f"Isotopic data for impurity '{impurity_key}' not found.")
            else:
                symbol = impurity_key
                avg_atomic_mass = self.elemental_masses.get(symbol)
                if avg_atomic_mass is None or avg_atomic_mass == 0:
                    raise ValueError(f"Could not determine average atomic mass for elemental impurity '{symbol}'.")

                element_atom_density = (impurity_mass_density / avg_atomic_mass) * self.AVOGADRO_NUMBER
                iso_data = self.ISOTOPIC_DATA[symbol.lower()]
                for mass_num, data in iso_data.items():
                    if data.get('abundance', 0) > 0:
                        isotope_name = f"{symbol}-{mass_num}"
                        isotope_densities[isotope_name] += element_atom_density * data['abundance']
                    # Handle cases like Ac-227 which has 0 natural abundance but is the only isotope
                    elif len(iso_data) == 1:
                        isotope_name = f"{symbol}-{mass_num}"
                        isotope_densities[isotope_name] += element_atom_density
                        break

        for iso, dens in isotope_densities.items():
            isotope_densities[iso] = dens / self.BARN_CM_CONVERSION

        return dict(sorted(isotope_densities.items()))


class FlibeUF4Salt(FluorideSalt):
    """
    A specialized class for LiF-BeF2-UF4 salts (FLiBe-UF4) with optional admixtures.
    Assumes a 2:1 molar ratio of LiF to BeF2 for the solvent.
    """

    def __init__(self, name: str, temperature: float,
                 u_enr_weight: dict, u_impurities_weight: dict,
                 UF4: float, admixtures: dict = None,
                 UF3_to_UF4: float = 0.05,
                 Li7_enr: float = 0.99995):
        """
        Initializes the FLiBe-UF4 salt.

        Args:
            name (str): A name for the salt material.
            temperature (float): The operating temperature in Kelvin.
            u_enr_weight (dict): Dictionary of Uranium isotope names to their weight fractions.
            u_impurities_weight (dict): Dictionary of impurities relative to the uranium mass.
            UF4 (float): Molar fraction of UF4 in the salt.
            admixtures (dict): Optional dictionary defining extra salts.
                               Example: {'LuF3': {'mol_frac': 0.01, 'enrichment': {176: 1.0}}}
            UF3_to_UF4 (float): Molar ratio of U3+ to U4+.
            Li7_enr (float): Weight fraction enrichment of Li-7.
        """
        self.name = name
        self.temperature = temperature

        # --- Handle Admixtures ---
        total_admixture_frac = 0.0
        if admixtures:
            total_admixture_frac = sum(v.get('mol_frac', 0) for v in admixtures.values())

        if not (0 <= UF4 + total_admixture_frac < 1):
            raise ValueError("Sum of UF4 and admixture molar fractions must be less than 1.")

        flibe_frac = 1.0 - UF4 - total_admixture_frac

        # --- Assemble Composition String ---
        # 2:1 ratio for LiF:BeF2
        components = {
            'UF4': UF4,
            'LiF': flibe_frac * (2.0 / 3.0),
            'BeF2': flibe_frac * (1.0 / 3.0)
        }

        if admixtures:
            for salt, props in admixtures.items():
                components[salt] = props.get('mol_frac', 0)

        composition_str = "-".join([f"{frac * 100:.6f}%{salt}" for salt, frac in components.items() if frac > 0])

        # --- Assemble Enrichments ---
        enrichments = {'Li7_enr': Li7_enr}
        if u_enr_weight:
            enrichments['U'] = u_enr_weight

        if admixtures:
            for salt, props in admixtures.items():
                if 'enrichment' in props:
                    # Extract cation from salt formula (e.g. "LuF3" -> "Lu")
                    cation_match = re.match(r"([A-Z][a-z]?)", salt, re.IGNORECASE)
                    if not cation_match:
                        raise ValueError(f"Could not parse cation from admixture salt: {salt}")
                    cation = cation_match.group(1).capitalize()
                    enrichments[cation] = props['enrichment']

        # --- Calculate Impurity Weight Fractions (Relative to Total Salt Mass) ---
        # We need a temporary object to calculate the mass of Uranium per mole of salt
        temp_salt = FluorideSalt(composition_str, enrichments=enrichments)
        avg_u_mass = temp_salt.elemental_masses.get('U', 0)

        # Calculate total mass of 1 mole of the salt mixture
        total_salt_mass_per_mole = sum(
            (frac / 100.0) * temp_salt.molar_masses[s]
            for s, frac in temp_salt.components.items()
        )

        # Calculate mass of Uranium in 1 mole of salt
        # UF4 component percentage / 100 * Atomic Mass of U
        m_u_per_mole_salt = (temp_salt.components.get('UF4', 0) / 100.0) * avg_u_mass

        final_impurities_wt = {}
        if m_u_per_mole_salt > 0 and u_impurities_weight:
            impurity_masses = {}
            for key, wt_frac_of_u in u_impurities_weight.items():
                impurity_masses[key] = wt_frac_of_u * m_u_per_mole_salt

            total_impurity_mass = sum(impurity_masses.values())
            total_system_mass = total_salt_mass_per_mole + total_impurity_mass

            if total_system_mass > 0:
                final_impurities_wt = {k: m / total_system_mass for k, m in impurity_masses.items()}

        # --- Initialize Base Class ---
        super().__init__(
            composition_str=composition_str,
            enrichments=enrichments,
            impurities_wt=final_impurities_wt,
            uf3_to_uf4_ratio=UF3_to_UF4
        )

    def get_atom_densities(self):
        """Convenience method to get atom densities at the stored temperature."""
        return super().get_atom_densities(self.temperature)

    def get_density(self):
        """Convenience method to get density at the stored temperature."""
        return super().density(self.temperature)


class FlibeSalt(FluorideSalt):
    """
    A specialized class for pure 2LiF-BeF2 eutectic salts (FLiBe).
    Provides a simplified interface for defining enrichment and impurities.
    """

    def __init__(self, name: str, temperature: float,
                 Li7_enr: float = 0.99995, impurities_wt: dict = None):
        """
        Initializes the pure FLiBe salt.

        Args:
            name (str): A name for the salt material.
            temperature (float): The operating temperature in Kelvin.
            Li7_enr (float): Weight fraction enrichment of Li-7.
            impurities_wt (dict, optional): Dictionary of impurity isotope/element names to their
                                            weight fractions of the *total salt*. Defaults to None.
        """
        self.name = name
        self.temperature = temperature
        lif_frac_eutectic = 2.0 / 3.0
        bef2_frac_eutectic = 1.0 / 3.0
        composition_str = f"{lif_frac_eutectic * 100:.4f}%LiF-{bef2_frac_eutectic * 100:.4f}%BeF2"

        enrichments = {'Li7_enr': Li7_enr}

        super().__init__(
            composition_str=composition_str,
            enrichments=enrichments,
            impurities_wt=impurities_wt,
            uf3_to_uf4_ratio=None
        )

    def get_atom_densities(self):
        """Convenience method to get atom densities at the stored temperature."""
        return super().get_atom_densities(self.temperature)

    def get_density(self):
        """Convenience method to get density at the stored temperature."""
        return super().density(self.temperature)


class FlibeUF4AdmixtureSalt(FlibeUF4Salt):
    """
    Extends FlibeUF4Salt to include arbitrary fluoride salt admixtures.
    This class is now a thin wrapper around FlibeUF4Salt, which handles admixtures directly.
    """

    def __init__(self, name: str, temperature: float,
                 u_enr_weight: dict, u_impurities_weight: dict,
                 UF4: float, admixtures: dict, UF3_to_UF4: float = 0.05,
                 Li7_enr: float = 0.99995):
        """
        Initializes a FLiBe-UF4 salt with additional fluoride admixtures.
        Delegates completely to FlibeUF4Salt.
        """
        super().__init__(
            name=name,
            temperature=temperature,
            u_enr_weight=u_enr_weight,
            u_impurities_weight=u_impurities_weight,
            UF4=UF4,
            admixtures=admixtures,
            UF3_to_UF4=UF3_to_UF4,
            Li7_enr=Li7_enr
        )


if __name__ == "__main__":
    try:
        # Example using the FlibeUF4Salt child class with case-insensitive API
        print("\n--- Testing FlibeUF4Salt Child Class (case-insensitive) ---")

        flibe_uf4_salt = FlibeUF4Salt(
            name="FuelSalt_LEU",
            temperature=950,
            u_enr_weight={'u-235': 0.1975, 'u-238': 0.8025},
            u_impurities_weight={'al': 275e-6, 'ac-227': 1e-9},  # lowercase keys
            UF4=0.22,
            UF3_to_UF4=0.08,
            Li7_enr=0.9999
        )

        print(f"Successfully created salt: '{flibe_uf4_salt.name}'")
        print(f"Composition: {flibe_uf4_salt.components}")
        print(f"Impurities (total wt frac): {flibe_uf4_salt.impurities_wt}")

        atom_densities_flibe_uf4 = flibe_uf4_salt.get_atom_densities()

        print(f"\n--- Atom Densities for {flibe_uf4_salt.name} at {flibe_uf4_salt.temperature} K,"
              f" rho = {flibe_uf4_salt.get_density():.5f} g/cm3 ---")
        for iso, dens in atom_densities_flibe_uf4.items():
            print(f"  - {iso:<8}: {dens:.4e}")

        # Example using the FlibeSalt child class
        print("\n--- Testing FlibeSalt Child Class ---")

        pure_flibe_salt = FlibeSalt(
            name="CoolantFlibe",
            temperature=923,
            Li7_enr=0.99995,
            impurities_wt={'fe-56': 10e-6}  # lowercase key
        )

        print(f"\nSuccessfully created salt: '{pure_flibe_salt.name}'")
        atom_densities_flibe = pure_flibe_salt.get_atom_densities()

        print(f"--- Atom Densities for {pure_flibe_salt.name} at {pure_flibe_salt.temperature} K, "
              f"rho = {pure_flibe_salt.get_density():.5f} g/cm3---")
        for iso, dens in atom_densities_flibe.items():
            print(f"  - {iso:<8}: {dens:.4e}")

        # Example using the FlibeUF4AdmixtureSalt class (now wrapper)
        print("\n--- Testing FlibeUF4AdmixtureSalt Child Class (case-insensitive) ---")

        admixture_salt = FlibeUF4AdmixtureSalt(
            name="AdvFlibeSalt",
            temperature=950,
            u_enr_weight={'u-235': 0.1975, 'u-238': 0.8025},
            u_impurities_weight={'al': 275e-6},
            UF4=0.20,
            admixtures={
                'luf3': {'mol_frac': 0.01, 'enrichment': {176: 1.0}},  # Lowercase salt
                'YbF3': {'mol_frac': 0.02}  # Mixed case salt
            },
            UF3_to_UF4=0.08,
            Li7_enr=0.9999
        )

        print(f"\nSuccessfully created salt: '{admixture_salt.name}'")
        print(f"Composition: {admixture_salt.components}")
        print(f"Enrichments: {admixture_salt.enrichments}")
        print(f"Impurities (total wt frac): {admixture_salt.impurities_wt}")

        atom_densities_admix = admixture_salt.get_atom_densities()

        print(f"--- Atom Densities for {admixture_salt.name} at {admixture_salt.temperature} K,"
              f" rho = {admixture_salt.get_density():.5f} g/cm3 ---")
        for iso, dens in atom_densities_admix.items():
            print(f"  - {iso:<8}: {dens:.4e}")

        # --- Test the from_atom_densities method ---
        print("\n--- Testing from_atom_densities Reconstruction ---")
        reconstructed_salt = FluorideSalt.from_atom_densities(atom_densities_admix)
        # print(atom_densities_admix)
        print(f"Original Composition: {admixture_salt.components}")
        print(f"Reconstructed Composition: {reconstructed_salt.components}")
        # Compare the two dictionaries, allowing for small float inaccuracies
        original_norm = {k: round(v, 4) for k, v in admixture_salt.components.items()}
        reconstructed_norm = {k: round(v, 4) for k, v in reconstructed_salt.components.items()}
        if original_norm == reconstructed_norm:
            print("SUCCESS: Reconstructed composition matches original.")
        else:
            print("FAILURE: Reconstructed composition does NOT match original.")

    except (ValueError, KeyError, RuntimeError, TypeError) as e:
        import traceback
        print(f"\nAn error occurred: {e}")
        traceback.print_exc()