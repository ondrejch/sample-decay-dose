import re
from collections import defaultdict
from sample_decay_dose import ISOTOPIC_DATA


class FluorideSalt:
    """
    Represents a mixed fluoride salt, calculating its density and isotopic atom densities.
    - Molar masses are calculated from isotopic data.
    - Supports custom isotope enrichments and impurities by weight fraction.
    - Allows 'XXX%' to specify a component's fraction as the remainder to 100%.
    """

    SALT_COEFFICIENTS = {
        "LiF":   {"A": 2.37, "B": 5e-4},
        "BeF2":  {"A": 1.97, "B": 1.45e-5},
        "UF4":   {"A": 7.78, "B": 9.92e-4},
        "LaF3":  {"A": 5.79, "B": 6.82e-4},
        "KF":    {"A": 2.64, "B": 6.57e-4},
        "NaF":   {"A": 2.70, "B": 5.90e-4},
        "SrF2":  {"A": 4.78, "B": 7.51e-4},
        "ThF4":  {"A": 7.11, "B": 7.59e-4},
        "ZrF4":  {"A": 5.36, "B": 1.23e-3}
    }
    AVOGADRO_NUMBER = 6.02214076e23     # atoms/mol
    BARN_CM_CONVERSION = 1e24           # (cm^2 / barn)

    def __init__(self, composition_str: str, enrichments: dict = None, impurities_wt: dict = None):
        """
        Initializes the FluorideSalt object.

        Args:
            composition_str: String defining the molar composition of the salt (e.g., "78%LiF-22%UF4").
            enrichments: Dictionary defining custom isotopic enrichments (e.g., {'U': {235: 0.05, 238: 0.95}}).
            impurities_wt: Dictionary defining impurities by weight fraction (e.g., {'Fe-56': 1e-5}).
        """
        # Programmatically add lanthanide coefficients
        lanthanide_symbols = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
        for symbol in lanthanide_symbols:
            formula = f"{symbol}F3"
            if formula not in self.SALT_COEFFICIENTS:
                self.SALT_COEFFICIENTS[formula] = self.SALT_COEFFICIENTS["LaF3"]

        self.enrichments = enrichments or {}
        self.impurities_wt = impurities_wt or {}
        self.total_impurity_frac = sum(self.impurities_wt.values())
        if self.total_impurity_frac >= 1.0:
            raise ValueError("Total weight fraction of impurities must be less than 1.")
        
        self.components = self._parse_composition(composition_str)
        self.elemental_masses = self._calculate_elemental_masses()
        self.molar_masses = {}
        self._validate_and_calculate_masses()

    def _parse_composition(self, s: str) -> dict:
        # This function remains the same as before
        components = {}
        parts = s.split('-')
        xxx_component_salt = None
        total_specified_percentage = 0.0
        for part in parts:
            match = re.match(r"((?:\d+(?:\.\d+)?)|XXX)%([A-Za-z0-9]+)", part, re.IGNORECASE)
            if not match:
                raise ValueError(f"Could not parse component: '{part}'")
            percentage_str, salt = match.group(1), match.group(2)
            if percentage_str.upper() == "XXX":
                if xxx_component_salt is not None:
                    raise ValueError("Only one 'XXX%' component allowed.")
                xxx_component_salt = salt
                components[salt] = 0
            else:
                percentage = float(percentage_str)
                if percentage < 0:
                    raise ValueError(f"Percentage cannot be negative: {part}")
                components[salt] = percentage
                total_specified_percentage += percentage
        if xxx_component_salt:
            remainder = 100.0 - total_specified_percentage
            if remainder < -1e-9:
                raise ValueError(f"Specified percentages ({total_specified_percentage}%) exceed 100%.")
            components[xxx_component_salt] = remainder
        return components

    def _calculate_elemental_masses(self) -> dict:
        """Calculate avg atomic mass for each element based on enrichment or natural abundance."""
        elemental_masses = {}
        # Find all unique elements in the salt composition AND impurities
        all_elements = set()
        for salt_formula in self.components.keys():
            for symbol, count in re.findall(r'([A-Z][a-z]*)(\d*)', salt_formula):
                all_elements.add(symbol)
        for isotope_name in self.impurities_wt.keys():
            symbol = re.match(r"([A-Z][a-z]*)", isotope_name).group(1)
            all_elements.add(symbol)


        for symbol in all_elements:
            iso_dist = {}
            if symbol in self.enrichments:
                # Use user-provided enrichment
                custom_dist = self.enrichments[symbol]
                if abs(sum(custom_dist.values()) - 1.0) > 1e-6:
                    raise ValueError(f"Enrichment fractions for {symbol} must sum to 1.")
                iso_dist = custom_dist
            else:
                # Use natural abundance
                try:
                    iso_data = ISOTOPIC_DATA[symbol.lower()]
                    iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}
                except KeyError:
                    # If data is not found, it might be an impurity with no natural abundance defined
                    # which is fine as long as it's not part of the main salt composition.
                    if any(symbol in f for f in self.components.keys()):
                        raise ValueError(f"Isotopic data not found for element: {symbol}")
                    else:
                        continue


            # Calculate weighted average atomic mass
            avg_mass = 0.0
            for mass_num, fraction in iso_dist.items():
                try:
                    isotope_mass = ISOTOPIC_DATA[symbol.lower()][mass_num]['mass']
                    avg_mass += isotope_mass * fraction
                except KeyError:
                    raise ValueError(f"Isotope {symbol}-{mass_num} not found in data.")
            elemental_masses[symbol] = avg_mass
        return elemental_masses

    def _calculate_molar_mass(self, formula: str) -> float:
        """Calculates molar mass of a formula using pre-calculated elemental masses."""
        total_mass = 0.0
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        if not elements: raise ValueError(f"Cannot parse formula: {formula}")
        for symbol, count_str in elements:
            count = int(count_str) if count_str else 1
            if symbol not in self.elemental_masses:
                raise ValueError(f"Elemental mass for '{symbol}' not found. Check isotopic data.")
            total_mass += self.elemental_masses[symbol] * count
        return total_mass

    def _validate_and_calculate_masses(self):
        # This function remains mostly the same
        total_percentage = sum(self.components.values())
        if abs(total_percentage - 100.0) > 1e-6:
            raise ValueError(f"Molar percentages must sum to 100, but sum to {total_percentage}")
        for salt in self.components:
            if salt not in self.SALT_COEFFICIENTS:
                raise ValueError(f"Unsupported salt component: '{salt}'")
            self.molar_masses[salt] = self._calculate_molar_mass(salt)

    def density(self, temperature: float) -> float:
        """
        Calculates the bulk density [g/cm^3] of the salt mixture, including impurities.
        Assumes impurities do not change the volume of the salt.
        """
        mixture_molar_mass = 0.0
        mixture_molar_volume = 0.0
        for salt, percentage in self.components.items():
            molar_fraction = percentage / 100.0
            coeffs = self.SALT_COEFFICIENTS[salt]
            molar_mass = self.molar_masses[salt]
            density_pure = coeffs['A'] - coeffs['B'] * temperature
            if density_pure <= 0:
                raise ValueError(f"Density of {salt} is non-positive at {temperature} K.")
            molar_volume_pure = molar_mass / density_pure
            mixture_molar_mass += molar_fraction * molar_mass
            mixture_molar_volume += molar_fraction * molar_volume_pure
        
        if mixture_molar_volume <= 0:
            raise ValueError("Mixture molar volume is non-positive.")
            
        # Density of the pure salt mixture
        pure_salt_density = mixture_molar_mass / mixture_molar_volume
        
        # Adjust density for the mass of impurities, assuming volume is constant
        final_density = pure_salt_density / (1.0 - self.total_impurity_frac)
        
        return final_density

    def get_atom_densities(self, temperature: float) -> dict:
        """
        Calculates the atom density of each isotope in the mixture, including impurities.

        Args:
            temperature: The temperature in Kelvin.

        Returns:
            A dictionary mapping isotope names (e.g., 'U-235') to their
            atom densities in [atoms / barn-cm].
        """
        bulk_density = self.density(temperature) # g/cm^3
        salt_weight_fraction = 1.0 - self.total_impurity_frac
        salt_mass_density = bulk_density * salt_weight_fraction

        # Calculate overall mixture molar mass (already weighted by molar fractions)
        mixture_molar_mass = sum(
            (p / 100.0) * self.molar_masses[s] for s, p in self.components.items()
        )

        # Total number density of the salt "molecules" [molecules/cm^3]
        total_salt_number_density = (salt_mass_density / mixture_molar_mass) * self.AVOGADRO_NUMBER

        isotope_densities = defaultdict(float)
        # 1. Calculate densities for base salt components
        for salt_formula, salt_percentage in self.components.items():
            elements_in_formula = re.findall(r'([A-Z][a-z]*)(\d*)', salt_formula)
            for symbol, count_str in elements_in_formula:
                atoms_per_molecule = int(count_str) if count_str else 1

                # Get isotopic distribution for this element
                iso_dist = {}
                if symbol in self.enrichments:
                    iso_dist = self.enrichments[symbol]
                else:
                    iso_data = ISOTOPIC_DATA[symbol.lower()]
                    iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}

                for mass_num, fraction in iso_dist.items():
                    if fraction == 0.0:  # Skip zero-fraction nuclides
                        continue
                    isotope_name = f"{symbol}-{mass_num}"
                    # Contribution to atom density from this salt component
                    density_contrib = (total_salt_number_density *    # [molecules/cm^3]
                                       (salt_percentage / 100.0) *    # [unit-less fraction]
                                       atoms_per_molecule *           # [atoms of element / molecule]
                                       fraction)                      # [atoms of isotope / atoms of element]

                    isotope_densities[isotope_name] += density_contrib
        
        # 2. Add densities for impurities
        for isotope_name, weight_frac in self.impurities_wt.items():
            match = re.match(r"([A-Z][a-z]*)-(\d+)", isotope_name)
            if not match:
                raise ValueError(f"Could not parse impurity isotope name: {isotope_name}")
            symbol, mass_num_str = match.groups()
            mass_num = int(mass_num_str)

            try:
                isotope_mass = ISOTOPIC_DATA[symbol.lower()][mass_num]['mass']
            except KeyError:
                raise ValueError(f"Isotopic data for impurity '{isotope_name}' not found.")

            impurity_mass_density = bulk_density * weight_frac
            impurity_atom_density = (impurity_mass_density / isotope_mass) * self.AVOGADRO_NUMBER
            isotope_densities[isotope_name] += impurity_atom_density


        # Convert from atoms/cm^3 to atoms/barn-cm
        for iso, dens in isotope_densities.items():
            isotope_densities[iso] = dens / self.BARN_CM_CONVERSION

        return dict(sorted(isotope_densities.items()))


# Example Usage:
if __name__ == "__main__":
    try:
        # Define a salt with 5% enriched Uranium and some impurities
        enrichment_leu = {'U': {235: 0.05, 238: 0.95}}
        # Impurities are defined as weight fractions (e.g., 10 ppm Fe-56)
        impurities = {'Fe-56': 10e-6, 'Cr-52': 5e-6, 'Ni-58': 2e-6}
        
        leu_salt_str = "78%LiF-XXX%UF4"

        leu_salt = FluorideSalt(leu_salt_str, enrichments=enrichment_leu, impurities_wt=impurities)

        print(f"Original Composition String: '{leu_salt_str}'")
        print(f"With Enrichments: {enrichment_leu}")
        print(f"With Impurities (wt. frac): {impurities}")
        
        print("\nCalculated Component Percentages:")
        for salt, perc in leu_salt.components.items():
            print(f"  - {salt}: {perc:.2f}%")

        print("\nAverage elemental masses (g/mol):")
        for el, mass in sorted(leu_salt.elemental_masses.items()):
            print(f"  - {el}: {mass:.4f}")

        temp_K = 950
        atom_densities = leu_salt.get_atom_densities(temp_K)

        print(f"\n--- Atom Densities at {temp_K} K ---")
        print(f"Bulk Density: {leu_salt.density(temp_K):.4f} g/cm³")
        print("Isotope Densities (atoms/barn-cm):")
        for iso, dens in atom_densities.items():
            print(f"  - {iso:<6}: {dens:.4e}")

    except (ValueError, KeyError, RuntimeError) as e:
        print(f"\nAn error occurred: {e}")
