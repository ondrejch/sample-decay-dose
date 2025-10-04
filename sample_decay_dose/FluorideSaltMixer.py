import re
from collections import defaultdict


class FluorideSalt:
    """
    Represents a mixed fluoride salt, calculating its density and isotopic atom densities.
    - Molar masses are calculated from isotopic data.
    - Supports custom isotope enrichments.
    - Allows 'XXX%' to specify a component's fraction as the remainder to 100%.
    """

    # Source: Mostly from NIST. Abundances are mole fractions.
    ISOTOPIC_DATA = {
        'li': {6: {'mass': 6.015122, 'abundance': 0.0759}, 7: {'mass': 7.016004, 'abundance': 0.9241}},
        'be': {9: {'mass': 9.012182, 'abundance': 1.0}},
        'f': {19: {'mass': 18.998403, 'abundance': 1.0}},
        'na': {23: {'mass': 22.989770, 'abundance': 1.0}},
        'k': {39: {'mass': 38.963707, 'abundance': 0.932581}, 40: {'mass': 39.963998, 'abundance': 0.000117}, 41: {'mass': 40.961826, 'abundance': 0.067302}},
        'sr': {84: {'mass': 83.913425, 'abundance': 0.0056}, 86: {'mass': 85.909260, 'abundance': 0.0986}, 87: {'mass': 86.908877, 'abundance': 0.0700}, 88: {'mass': 87.905612, 'abundance': 0.8258}},
        'zr': {90: {'mass': 89.904704, 'abundance': 0.5145}, 91: {'mass': 90.905644, 'abundance': 0.1122}, 92: {'mass': 91.905040, 'abundance': 0.1715}, 94: {'mass': 93.906315, 'abundance': 0.1738}, 96: {'mass': 95.908273, 'abundance': 0.0280}},
        'th': {232: {'mass': 232.03805, 'abundance': 1.0}},
        'u': {234: {'mass': 234.040952, 'abundance': 0.000055}, 235: {'mass': 235.043929, 'abundance': 0.007204}, 238: {'mass': 238.050788, 'abundance': 0.992741}},
        # Lanthanides
        'la': {138: {'mass': 137.90711, 'abundance': 0.0009}, 139: {'mass': 138.90635, 'abundance': 0.9991}},
        'ce': {136: {'mass': 135.90717, 'abundance': 0.00185}, 138: {'mass': 137.90599, 'abundance': 0.00251}, 140: {'mass': 139.90544, 'abundance': 0.88450}, 142: {'mass': 141.90924, 'abundance': 0.11114}},
        'pr': {141: {'mass': 140.90765, 'abundance': 1.0}},
        'nd': {142: {'mass': 141.90772, 'abundance': 0.272}, 143: {'mass': 142.90981, 'abundance': 0.122}, 144: {'mass': 143.91009, 'abundance': 0.238}, 145: {'mass': 144.91257, 'abundance': 0.083}, 146: {'mass': 145.91312, 'abundance': 0.172}, 148: {'mass': 147.91689, 'abundance': 0.057}, 150: {'mass': 149.92089, 'abundance': 0.056}},
        'pm': {145: {'mass': 144.912749, 'abundance': 1.0}}, # Most stable isotope
        'sm': {144: {'mass': 143.91200, 'abundance': 0.0307}, 147: {'mass': 146.91490, 'abundance': 0.1499}, 148: {'mass': 147.91482, 'abundance': 0.1124}, 149: {'mass': 148.91718, 'abundance': 0.1382}, 150: {'mass': 149.91728, 'abundance': 0.0738}, 152: {'mass': 151.91973, 'abundance': 0.2675}, 154: {'mass': 153.92221, 'abundance': 0.2275}},
        'eu': {151: {'mass': 150.91985, 'abundance': 0.4781}, 153: {'mass': 152.92123, 'abundance': 0.5219}},
        'gd': {152: {'mass': 151.91979, 'abundance': 0.0020}, 154: {'mass': 153.92086, 'abundance': 0.0218}, 155: {'mass': 154.92262, 'abundance': 0.1480}, 156: {'mass': 155.92212, 'abundance': 0.2047}, 157: {'mass': 156.92396, 'abundance': 0.1565}, 158: {'mass': 157.92410, 'abundance': 0.2484}, 160: {'mass': 159.92705, 'abundance': 0.2186}},
        'tb': {159: {'mass': 158.92535, 'abundance': 1.0}},
        'dy': {156: {'mass': 155.92428, 'abundance': 0.00056}, 158: {'mass': 157.92441, 'abundance': 0.00095}, 160: {'mass': 159.92519, 'abundance': 0.02329}, 161: {'mass': 160.92693, 'abundance': 0.18889}, 162: {'mass': 161.92680, 'abundance': 0.25475}, 163: {'mass': 162.92873, 'abundance': 0.24896}, 164: {'mass': 163.92917, 'abundance': 0.28260}},
        'ho': {165: {'mass': 164.93032, 'abundance': 1.0}},
        'er': {162: {'mass': 161.92878, 'abundance': 0.00139}, 164: {'mass': 163.92920, 'abundance': 0.01601}, 166: {'mass': 165.93029, 'abundance': 0.33503}, 167: {'mass': 166.93205, 'abundance': 0.22869}, 168: {'mass': 167.93237, 'abundance': 0.26978}, 170: {'mass': 169.93546, 'abundance': 0.14910}},
        'tm': {169: {'mass': 168.93421, 'abundance': 1.0}},
        'yb': {168: {'mass': 167.93389, 'abundance': 0.0013}, 170: {'mass': 169.93476, 'abundance': 0.0304}, 171: {'mass': 170.93632, 'abundance': 0.1428}, 172: {'mass': 171.93638, 'abundance': 0.2183}, 173: {'mass': 172.93821, 'abundance': 0.1613}, 174: {'mass': 173.93886, 'abundance': 0.3183}, 176: {'mass': 175.94257, 'abundance': 0.1276}},
        'lu': {175: {'mass': 174.94077, 'abundance': 0.9741}, 176: {'mass': 175.94268, 'abundance': 0.0259}},
    }

    SALT_COEFFICIENTS = {
        "LiF": {"A": 2.37, "B": 5e-4}, "BeF2": {"A": 1.97, "B": 1.45e-5},
        "UF4": {"A": 7.78, "B": 9.92e-4}, "LaF3": {"A": 5.79, "B": 6.82e-4},
        "KF": {"A": 2.64, "B": 6.57e-4}, "NaF": {"A": 2.70, "B": 5.90e-4},
        "SrF2": {"A": 4.78, "B": 7.51e-4}, "ThF4": {"A": 7.11, "B": 7.59e-4},
        "ZrF4": {"A": 5.36, "B": 1.23e-3}
    }
    AVOGADRO_NUMBER = 6.02214076e23     # atoms/mol
    BARN_CM_CONVERSION = 1e24           # (cm^2 / barn)

    def __init__(self, composition_str: str, enrichments: dict = None):
        # Programmatically add lanthanide coefficients
        lanthanide_symbols = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
        for symbol in lanthanide_symbols:
            formula = f"{symbol}F3"
            if formula not in self.SALT_COEFFICIENTS:
                self.SALT_COEFFICIENTS[formula] = self.SALT_COEFFICIENTS["LaF3"]

        self.enrichments = enrichments or {}
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
            if not match: raise ValueError(f"Could not parse component: '{part}'")
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
        # Find all unique elements in the salt composition
        all_elements = set()
        for salt_formula in self.components.keys():
            for symbol, count in re.findall(r'([A-Z][a-z]*)(\d*)', salt_formula):
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
                    iso_data = self.ISOTOPIC_DATA[symbol.lower()]
                    iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}
                except KeyError:
                    raise ValueError(f"Isotopic data not found for element: {symbol}")

            # Calculate weighted average atomic mass
            avg_mass = 0.0
            for mass_num, fraction in iso_dist.items():
                try:
                    isotope_mass = self.ISOTOPIC_DATA[symbol.lower()][mass_num]['mass']
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
        """Calculates the bulk density [g/cm^3] of the salt mixture."""
        mixture_molar_mass = 0.0
        mixture_molar_volume = 0.0
        for salt, percentage in self.components.items():
            molar_fraction = percentage / 100.0
            coeffs = self.SALT_COEFFICIENTS[salt]
            molar_mass = self.molar_masses[salt]
            density_pure = coeffs['A'] - coeffs['B'] * temperature
            if density_pure <= 0: raise ValueError(f"Density of {salt} is non-positive at {temperature} K.")
            molar_volume_pure = molar_mass / density_pure
            mixture_molar_mass += molar_fraction * molar_mass
            mixture_molar_volume += molar_fraction * molar_volume_pure
        if mixture_molar_volume <= 0: raise ValueError("Mixture molar volume is non-positive.")
        return mixture_molar_mass / mixture_molar_volume

    def get_atom_densities(self, temperature: float) -> dict:
        """
        Calculates the atom density of each isotope in the mixture.

        Args:
            temperature: The temperature in Kelvin.

        Returns:
            A dictionary mapping isotope names (e.g., 'U-235') to their
            atom densities in [atoms / barn-cm].
        """
        bulk_density = self.density(temperature) # g/cm^3

        # Calculate overall mixture molar mass (already weighted by molar fractions)
        mixture_molar_mass = sum(
            (p / 100.0) * self.molar_masses[s] for s, p in self.components.items()
        )

        # Total number density of the salt "molecules" [molecules/cm^3]
        total_number_density = (bulk_density / mixture_molar_mass) * self.AVOGADRO_NUMBER

        isotope_densities = defaultdict(float)
        for salt_formula, salt_percentage in self.components.items():
            elements_in_formula = re.findall(r'([A-Z][a-z]*)(\d*)', salt_formula)
            for symbol, count_str in elements_in_formula:
                atoms_per_molecule = int(count_str) if count_str else 1

                # Get isotopic distribution for this element
                iso_dist = {}
                if symbol in self.enrichments:
                    iso_dist = self.enrichments[symbol]
                else:
                    iso_data = self.ISOTOPIC_DATA[symbol.lower()]
                    iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}

                for mass_num, fraction in iso_dist.items():
                    isotope_name = f"{symbol}-{mass_num}"
                    # Contribution to atom density from this salt component
                    density_contrib = (total_number_density *         # [molecules/cm^3]
                                       (salt_percentage / 100.0) *    # [unitless fraction]
                                       atoms_per_molecule *           # [atoms of element / molecule]
                                       fraction)                      # [atoms of isotope / atoms of element]

                    isotope_densities[isotope_name] += density_contrib

        # Convert from atoms/cm^3 to atoms/barn-cm
        for iso, dens in isotope_densities.items():
            isotope_densities[iso] = dens / self.BARN_CM_CONVERSION

        return dict(sorted(isotope_densities.items()))


# Example Usage:
if __name__ == "__main__":
    try:
        # Define a salt with 5% enriched Uranium
        enrichment_leu = {'U': {235: 0.05, 238: 0.95}}
        leu_salt_str = "78%LiF-XXX%UF4"

        leu_salt = FluorideSalt(leu_salt_str, enrichments=enrichment_leu)

        print(f"Original Composition String: '{leu_salt_str}'")
        print(f"With Enrichments: {enrichment_leu}")
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

    except (ValueError, KeyError) as e:
        print(f"\nAn error occurred: {e}")
