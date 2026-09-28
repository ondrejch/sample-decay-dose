import re
import warnings
from collections import defaultdict
from decimal import Decimal
from sample_decay_dose.data import ISOTOPIC_DATA as _isotopic_data
from sample_decay_dose.ValencyMapper import ValencyMapper

# Nuclide names such as 'Li-7', 'li7', 'am-242m' or 'Am242m1'. The optional suffix marks an isomeric state.
_NUCLIDE_RE = re.compile(r"([A-Za-z]{1,3})-?(\d+)(m\d*)?", re.IGNORECASE)
# Keys of an isotopic distribution: an int mass number for a ground state, or a string such as '242m' for an isomer.
_ISOTOPE_KEY_RE = re.compile(r"(\d+)(m\d*)?")


def _normalize_isomer_suffix(suffix: str | None) -> str:
    """Return '' for a ground state, 'm' for the first isomer ('m' or 'm1'), and 'm2', 'm3', ... otherwise."""
    if not suffix:
        return ''
    suffix = suffix.lower()
    return 'm' if suffix == 'm1' else suffix


def _parse_nuclide(name: str) -> tuple[str, int, str] | None:
    """Split a nuclide name into (element symbol, mass number, isomer suffix). Returns None if it does not parse."""
    match = _NUCLIDE_RE.fullmatch(str(name).strip())
    if not match:
        return None
    symbol_raw, mass_num, suffix = match.groups()
    return symbol_raw.capitalize(), int(mass_num), _normalize_isomer_suffix(suffix)


def _isotope_key(mass_num: int, suffix: str):
    """Isotopic-distribution key: the int mass number for a ground state, or e.g. '242m' for an isomer."""
    return f"{mass_num}{suffix}" if suffix else mass_num


def _isotope_key_mass_number(key) -> int | None:
    """Mass number of an isotopic-distribution key (int, or a string such as '242m'). None if the key is invalid."""
    if isinstance(key, int) and not isinstance(key, bool):
        return key
    if isinstance(key, str):
        match = _ISOTOPE_KEY_RE.fullmatch(key.strip().lower())
        if match:
            return int(match.group(1))
    return None


def _format_percent(value: float) -> str:
    """Format a percentage in positional notation with full float precision.

    The composition string is parsed back with float(). Using the shortest round-trip representation keeps the
    parsed fractions equal to the computed ones, so their sum stays at 100 to float precision. Rounding every
    component to a fixed number of decimals accumulated up to n * 5e-7 of error with 6 decimals, which failed the
    1e-6 sum check. Exponent notation is avoided because '-' separates components.
    """
    return format(Decimal(repr(float(value))), 'f')


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
            enrichments: Dictionary defining custom isotopic compositions. Three forms are accepted:
                         - Mole fractions keyed by int mass number, listing every isotope of the element so that
                           they sum to 1, e.g. {'U': {235: 0.05, 238: 0.95}}. An isomer uses a string key such
                           as '242m', e.g. {'Am': {241: 0.9, 242: 0.05, '242m': 0.05}}.
                         - Uranium weight fractions keyed by nuclide name, e.g. {'U': {'u-235': 0.05, 'u-238': 0.95}}.
                         - Special weight-fraction keys, e.g. {'Li7_enr': 0.9999, 'U235_enr': 0.1975}.
            impurities_wt: Dictionary defining impurities by weight fraction of the *total salt mass*.
                           Keys are nuclides (e.g. 'Fe-56', 'Am-242m') or elements with natural abundance ('Fe').
                           Impurities are added as bare atoms. They carry no fluorine.
            uf3_to_uf4_ratio: Molar ratio of U3+ to U4+ to set the fluorine potential.
        """
        # This mapping is the key to robust, case-insensitive salt parsing.
        self._known_salts = self._known_salt_map()

        # Process enrichments from user-friendly formats (like wt%) to internal mole%
        self.enrichments = self._process_enrichments(enrichments)

        processed_impurities = {}
        if impurities_wt:
            for key, val in impurities_wt.items():
                parsed = _parse_nuclide(key)
                if parsed:
                    symbol, mass_num, suffix = parsed
                    new_key = f"{symbol}-{mass_num}{suffix}"
                elif re.fullmatch(r"[A-Za-z]{1,3}", str(key)):
                    new_key = str(key).capitalize()
                else:
                    new_key = key
                processed_impurities[new_key] = processed_impurities.get(new_key, 0.0) + val
        self.impurities_wt = processed_impurities

        if uf3_to_uf4_ratio is not None and not (0.0 <= uf3_to_uf4_ratio < float('inf')):
            raise ValueError(f"uf3_to_uf4_ratio must be a finite non-negative number, got {uf3_to_uf4_ratio}.")
        self.uf3_to_uf4_ratio = uf3_to_uf4_ratio

        self.total_impurity_frac = sum(self.impurities_wt.values())
        if self.total_impurity_frac >= 1.0:
            raise ValueError("Total weight fraction of impurities must be less than 1.")

        self.components = self._parse_composition(composition_str)
        self.elemental_masses = self._calculate_elemental_masses()
        self.molar_masses = {}
        self._validate_and_calculate_masses()
        # Species dropped by from_atom_densities(); empty for a salt defined directly.
        self.excluded_species = {}
        self.excess_fluorine_atom_fraction = 0.0

    @classmethod
    def _known_salt_map(cls) -> dict:
        """Map lower-case salt formulas to their canonical spelling (e.g. 'kf' -> 'KF', 'luf3' -> 'LuF3')."""
        known = {s.lower(): s for s in cls.SALT_COEFFICIENTS.keys()}
        for symbol in cls.LANTHANIDE_SYMBOLS:
            formula = f"{symbol}F3"
            known.setdefault(formula.lower(), formula)
        return known

    @classmethod
    def _canonical_formula(cls, formula: str) -> str:
        """Canonical spelling of a known salt formula. Unknown formulas are returned unchanged."""
        return cls._known_salt_map().get(formula.lower(), formula)

    @staticmethod
    def _formula_cation(formula: str) -> str | None:
        """Cation symbol of a canonical formula, e.g. 'KF' -> 'K', 'LuF3' -> 'Lu'. Case-sensitive on purpose."""
        match = re.match(r"[A-Z][a-z]?", formula)
        return match.group(0) if match else None

    def _natural_elemental_mass(self, symbol: str) -> float:
        """Natural-abundance average atomic mass [g/mol] of an element."""
        iso_data = self.ISOTOPIC_DATA[symbol.lower()]
        return sum(d['abundance'] * d['mass'] for d in iso_data.values())

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

        # Handle U-235 weight fraction enrichment (Auto-calculate U-234/U-238)
        if 'U235_enr' in processed:
            self._set_uranium_isotopics_from_enrichment(processed)

        # Handle generic Uranium weight fraction enrichment (Explicit dictionary)
        if 'U' in processed and isinstance(processed['U'], dict):
            # Check if it's in the user-friendly {'u-235': wt_frac} format. Isomer keys such as '235m' are
            # mole-fraction keys, so only keys that start with 'u' select the weight-fraction format.
            if any(isinstance(k, str) and k.strip().lower().startswith('u') for k in processed['U'].keys()):
                u_enr_weight = processed.pop('U')
                u_enr_mole = {}
                for key, wt_frac in u_enr_weight.items():
                    parsed = _parse_nuclide(key) if isinstance(key, str) else None
                    if not parsed or parsed[0] != 'U':
                        raise ValueError(f"Invalid key in u_enr_weight: '{key}'")
                    _, mass_num, suffix = parsed
                    isotope_mass = self.ISOTOPIC_DATA['u'][mass_num]['mass']
                    u_enr_mole[_isotope_key(mass_num, suffix)] = wt_frac / isotope_mass
                total_moles = sum(u_enr_mole.values())
                if total_moles > 0:
                    processed['U'] = {mass: mol / total_moles for mass, mol in u_enr_mole.items()}

        return processed

    def _set_uranium_isotopics_from_enrichment(self, processed_enrichments: dict):
        """
        Calculates U-234, U-235, and U-238 weight fractions based on U-235 enrichment,
        then converts them to mole fractions.
        
        Uses the standard correlation for fresh fuel (typical in SCALE/ASTM):
        w_234 = 0.0089 * w_235
        w_238 = 1.0 - w_235 - w_234
        """
        w235 = processed_enrichments.pop('U235_enr')
        if not (0 <= w235 <= 1):
            raise ValueError("U-235 enrichment must be between 0 and 1.")

        # Correlation for U-234 in fresh enriched uranium
        w234 = 0.0089 * w235
        w238 = 1.0 - w235 - w234
        
        # Handle edge case of very high enrichment where sum might exceed 1.0
        if w238 < 0:
            total_minor = w235 + w234
            w235 /= total_minor
            w234 /= total_minor
            w238 = 0.0

        # Retrieve isotopic masses
        m234 = self.ISOTOPIC_DATA['u'][234]['mass']
        m235 = self.ISOTOPIC_DATA['u'][235]['mass']
        m238 = self.ISOTOPIC_DATA['u'][238]['mass']

        # Convert weight fractions to mole fractions
        n234 = w234 / m234
        n235 = w235 / m235
        n238 = w238 / m238
        
        total_moles = n234 + n235 + n238

        processed_enrichments['U'] = {
            234: n234 / total_moles,
            235: n235 / total_moles,
            238: n238 / total_moles
        }

    @classmethod
    def from_atom_densities(cls, atom_densities: dict, valency_estimate_type: str = 'upper', *,
                            impurity_cutoff: float = 1e-3, drop_zero_valence: bool = True,
                            valence_overrides: dict = None):
        """
        Reconstructs a fluoride salt from nuclide atom densities, e.g. a decayed fuel salt read from an F71 file.

        Each element of the input is assigned to one of four roles:

        1. Salt component. A cation with a positive valence v (from ValencyMapper or ``valence_overrides``)
           whose share of the retained non-fluorine atoms is at least ``impurity_cutoff`` becomes the fluoride
           XF_v. Its isotopic composition goes into ``enrichments``. Isomers keep their own keys, so 'Am-242m'
           stays distinct from 'Am-242'.
        2. Impurity. A retained element below ``impurity_cutoff`` is carried in ``impurities_wt`` nuclide by
           nuclide. Impurities are bare atoms outside the fluorine balance, as in the forward model, so a trace
           metal such as an Al impurity leaves the UF3/UF4 ratio unchanged. Trace anions other than fluorine
           (Cl, Br, I, O, S, ...) are carried the same way. An anion above the cutoff raises ValueError.
        3. Excluded. Elements with valence 0 are dropped when ``drop_zero_valence`` is True. These are the noble
           gases, which are off-gassed, and the noble metals, which plate out. They are listed with their atom
           fraction of the input in ``salt.excluded_species``, and a UserWarning names them.
        4. Fluorine. The fluorine balance against the salt components sets the uranium redox state. A fluorine
           deficit is assigned to U3+ and sets ``uf3_to_uf4_ratio``. A fluorine excess (fission is oxidizing)
           is carried as a fluorine impurity, so no fluorine atoms are lost. Its atom fraction of the input is
           stored in ``salt.excess_fluorine_atom_fraction``.

        A round trip reproduces the retained atom densities: ``from_atom_densities(s.get_atom_densities(T))``
        returns the components, enrichments, impurities and UF3/UF4 ratio of ``s``. A component without density
        coefficients (e.g. CsF above the cutoff) is accepted here, and ``density()`` raises for it.

        Args:
            atom_densities (dict): Nuclide name (e.g. 'Li-7', 'li7', 'am-242m') to atom density [atoms/barn-cm].
            valency_estimate_type (str): The valency map to use ('upper', 'lower', 'doligez'). Defaults to 'upper'.
            impurity_cutoff (float): Share of the retained non-fluorine atoms below which an element is carried
                as an impurity. Defaults to 1e-3 (0.1 mol%). With 0, every cation becomes a salt component.
            drop_zero_valence (bool): Drop valence-0 elements (default). With False, they are kept as impurities.
            valence_overrides (dict): Element symbol to valence, used before the valency map. A valence of 0
                drops the element.

        Returns:
            A new FluorideSalt instance with ``excluded_species`` and ``excess_fluorine_atom_fraction`` set.

        Raises:
            ValueError: for an unparsable nuclide name, a negative density, missing fluorine or cations, an anion
                or an element of unknown valence above the cutoff, or a fluorine deficit that UF3 cannot absorb.
        """
        if impurity_cutoff < 0:
            raise ValueError(f"impurity_cutoff must be non-negative, got {impurity_cutoff}.")
        valence_mapper = ValencyMapper(estimate_type=valency_estimate_type)
        overrides = {str(k).capitalize(): v for k, v in (valence_overrides or {}).items()}

        def valence_of(element: str):
            return overrides[element] if element in overrides else valence_mapper.get(element)

        # Group atom densities by element and isotope key. Isomers get keys such as '242m'.
        nuclides = defaultdict(lambda: defaultdict(float))
        for name, density in atom_densities.items():
            parsed = _parse_nuclide(name)
            if parsed is None:
                raise ValueError(f"Could not parse nuclide name: '{name}'")
            if density < 0:
                raise ValueError(f"Negative atom density for '{name}': {density}")
            if density == 0:
                continue
            symbol, mass_num, suffix = parsed
            nuclides[symbol][_isotope_key(mass_num, suffix)] += density
        elemental = {symbol: sum(isotopes.values()) for symbol, isotopes in nuclides.items()}
        total_atoms = sum(elemental.values())
        fluorine_total = elemental.get('F', 0.0)
        if fluorine_total <= 0:
            raise ValueError("No fluorine in atom_densities, so a fluoride salt cannot be reconstructed.")

        # Valence-0 species are off-gassed (noble gases) or plate out (noble metals) and leave the salt.
        excluded = {}
        retained = []
        for symbol in elemental:
            if symbol == 'F':
                continue
            if drop_zero_valence and valence_of(symbol) == 0:
                excluded[symbol] = {
                    'atom_fraction': elemental[symbol] / total_atoms,
                    'reason': 'valence 0 (noble gas or noble metal)',
                    'nuclides': {f"{symbol}-{key}": dens / total_atoms for key, dens in nuclides[symbol].items()},
                }
            else:
                retained.append(symbol)
        retained_total = sum(elemental[symbol] for symbol in retained)
        if retained_total <= 0:
            raise ValueError("Cannot reconstruct salt from zero cation densities.")

        components = {}  # element symbol -> valence
        impurity_elements = []
        for symbol in retained:
            valence = valence_of(symbol)
            share = elemental[symbol] / retained_total
            if share < impurity_cutoff or valence == 0:
                impurity_elements.append(symbol)
            elif valence is None:
                raise ValueError(
                    f"Element '{symbol}' makes up {share:.3e} of the non-fluorine atoms, which is above "
                    f"impurity_cutoff={impurity_cutoff}, and the '{valency_estimate_type}' valency map has no "
                    f"valence for it. Pass valence_overrides={{'{symbol}': <valence>}}, or 0 to drop it.")
            elif valence < 0:
                raise ValueError(
                    f"Element '{symbol}' (valence {valence}) is an anion making up {share:.3e} of the "
                    f"non-fluorine atoms, which is above impurity_cutoff={impurity_cutoff}. Only fluorine can be "
                    f"a major anion of a fluoride salt. Other anions are carried as trace impurities.")
            else:
                components[symbol] = valence
        if not components:
            raise ValueError("No cation reaches impurity_cutoff, so the salt has no fluoride components.")

        # Fluorine balance. The salt components need sum(n * valence) fluorine atoms with U as U4+.
        fluorine_needed = sum(elemental[symbol] * valence for symbol, valence in components.items())
        deficit = fluorine_needed - fluorine_total
        tolerance = 1e-9 * fluorine_total
        uf3_to_uf4_ratio = None
        excess_fluorine = 0.0
        if deficit > tolerance:
            # Each U3+ carries one fluorine less than U4+, so the deficit equals the U3+ density.
            uranium = elemental['U'] if components.get('U') == 4 else 0.0
            if uranium <= 0:
                raise ValueError(
                    f"The salt components need {deficit:.4e} at/b-cm more fluorine than the input holds, and "
                    f"there is no UF4 component to take up the deficit as UF3. Check valency_estimate_type, "
                    f"valence_overrides, or the input densities.")
            if deficit >= uranium:
                raise ValueError(
                    f"The fluorine deficit ({deficit:.4e} at/b-cm) is at least the uranium density "
                    f"({uranium:.4e} at/b-cm). Absorbing it would need uranium below U3+.")
            uf3_to_uf4_ratio = deficit / (uranium - deficit)
        elif deficit < -tolerance:
            excess_fluorine = -deficit

        def nuclide_mass(element: str, key) -> float:
            # Isomers use the ground-state mass; the excitation energy changes it by less than 1e-5 u.
            try:
                return cls.ISOTOPIC_DATA[element.lower()][_isotope_key_mass_number(key)]['mass']
            except KeyError:
                raise ValueError(f"Isotopic mass not found for {element}-{key}.") from None

        # Weight fractions of the impurities relative to all retained atoms, as the forward model defines them.
        retained_mass = sum(dens * nuclide_mass(symbol, key)
                            for symbol in retained + ['F'] for key, dens in nuclides[symbol].items())
        impurities_wt = {}
        for symbol in impurity_elements:
            for key, dens in nuclides[symbol].items():
                impurities_wt[f"{symbol}-{key}"] = dens * nuclide_mass(symbol, key) / retained_mass
        if excess_fluorine > 0:
            for key, dens in nuclides['F'].items():
                impurities_wt[f"F-{key}"] = (excess_fluorine * dens / fluorine_total
                                             * nuclide_mass('F', key) / retained_mass)

        enrichments = {symbol: {key: dens / elemental[symbol] for key, dens in nuclides[symbol].items()}
                       for symbol in list(components) + ['F']}

        cation_total = sum(elemental[symbol] for symbol in components)
        composition_str = "-".join(
            f"{_format_percent(100.0 * elemental[symbol] / cation_total)}%{symbol}F{valence if valence > 1 else ''}"
            for symbol, valence in components.items())

        salt = cls(composition_str, enrichments=enrichments, impurities_wt=impurities_wt,
                   uf3_to_uf4_ratio=uf3_to_uf4_ratio)
        salt.excluded_species = excluded
        salt.excess_fluorine_atom_fraction = excess_fluorine / total_atoms
        if excluded:
            dropped = ", ".join(f"{symbol} ({info['atom_fraction']:.3e})" for symbol, info in excluded.items())
            warnings.warn(f"from_atom_densities dropped valence-0 species, listed with their atom fraction of the "
                          f"input: {dropped}. Details are in salt.excluded_species.", UserWarning, stacklevel=2)
        return salt

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
        # Elements of the salt components need a defined molar mass. Elements that appear only as impurities or
        # enrichments get one when it is defined. An elemental impurity without a mass raises at use.
        component_elements = set()
        for salt_formula in self.components.keys():
            # The salt_formula is now guaranteed to be in canonical form.
            for symbol, count in re.findall(r'([A-Z][a-z]?)(\d*)', salt_formula):
                component_elements.add(symbol)
        other_elements = set(self.enrichments.keys())
        for isotope_name in self.impurities_wt.keys():
            match = re.match(r"([A-Z][a-z]*)", isotope_name)
            if match:
                other_elements.add(match.group(1))

        for symbol in component_elements | other_elements:
            in_salt = symbol in component_elements
            if symbol in self.enrichments:
                custom_dist = self.enrichments[symbol]
                # Keys are int mass numbers, or strings such as '242m' for isomers.
                mass_numbers = {key: _isotope_key_mass_number(key) for key in custom_dist.keys()}
                if any(mass_num is None for mass_num in mass_numbers.values()):
                    raise TypeError(f"Enrichment dictionary for {symbol} has invalid keys: {list(custom_dist.keys())}."
                                    f" Use int mass numbers, or strings such as '242m' for isomers.")
                total_fraction = sum(custom_dist.values())
                if abs(total_fraction - 1.0) > 1e-6:
                    raise ValueError(f"Enrichment fractions for {symbol} must sum to 1, but sum to {total_fraction:.6g}."
                                     f" List every isotope, e.g. {{'U': {{235: 0.05, 238: 0.95}}}}.")
                avg_mass = 0.0
                for key, fraction in custom_dist.items():
                    try:
                        # Isomers use the ground-state mass; the excitation energy changes it by less than 1e-5 u.
                        avg_mass += self.ISOTOPIC_DATA[symbol.lower()][mass_numbers[key]]['mass'] * fraction
                    except KeyError:
                        raise ValueError(f"Isotope {symbol}-{key} not found in data.") from None
                elemental_masses[symbol] = avg_mass
                continue

            iso_data = self.ISOTOPIC_DATA.get(symbol.lower())
            if iso_data is None:
                if in_salt:
                    raise ValueError(f"Isotopic data not found for element: {symbol}")
                continue
            if sum(data['abundance'] for data in iso_data.values()) <= 0:
                # Elements such as Pm, Tc or the actinides above U have no natural isotopic composition.
                if in_salt:
                    raise ValueError(
                        f"Element '{symbol}' has no natural isotopic abundance, so its molar mass is undefined. "
                        f"Specify its isotopic composition, e.g. enrichments={{'{symbol}': {{<mass number>: 1.0}}}}.")
                continue
            elemental_masses[symbol] = sum(data['abundance'] * data['mass'] for data in iso_data.values())
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

    def _fluorine_removed_per_mole(self) -> float:
        """Moles of fluorine removed per mole of salt formula units by the UF3/UF4 ratio (one F per U3+)."""
        if not self.uf3_to_uf4_ratio:
            return 0.0
        frac_u3 = self.uf3_to_uf4_ratio / (1.0 + self.uf3_to_uf4_ratio)
        uranium_per_mole = 0.0
        for salt, percentage in self.components.items():
            for symbol, count_str in re.findall(r'([A-Z][a-z]?)(\d*)', salt):
                if symbol == 'U':
                    uranium_per_mole += (percentage / 100.0) * (int(count_str) if count_str else 1)
        return uranium_per_mole * frac_u3

    def mixture_molar_mass(self) -> float:
        """Mean molar mass of the pure salt [g per mole of formula units].

        The fluorine removed by the UF3/UF4 ratio is subtracted, so the mass density and the atom densities
        describe the same material.
        """
        molar_mass = sum((p / 100.0) * self.molar_masses[s] for s, p in self.components.items())
        return molar_mass - self._fluorine_removed_per_mole() * self.elemental_masses.get('F', 0.0)

    def _lanthanide_reference_molar_volume(self, temperature: float) -> float:
        """Molar volume [cm3/mol] of natural LaF3 at the given temperature.

        Lanthanide trifluorides without their own density correlation are assigned this molar volume. It is
        evaluated at every call because it depends on temperature. The value used to be cached at the first
        temperature, which made density(T) depend on call order.
        """
        ref_coeffs = self.SALT_COEFFICIENTS['LaF3']
        ref_density = ref_coeffs['A'] - ref_coeffs['B'] * temperature
        if ref_density <= 0:
            raise ValueError(f"Reference density for LaF3 is non-positive at {temperature} K.")
        ref_molar_mass = self._natural_elemental_mass('La') + 3.0 * self._natural_elemental_mass('F')
        return ref_molar_mass / ref_density

    def density(self, temperature: float) -> float:
        """Mass density [g/cm3] of the salt including impurities, at temperature [K].

        Pure-salt molar volumes are mixed ideally. Lanthanide trifluorides without coefficients use the molar
        volume of LaF3. Assumption for the UF3/UF4 ratio: reducing U4+ to U3+ removes fluorine at constant
        molar volume. The UF4 correlation still sets the uranium component's molar volume, and the mass of the
        removed fluorine lowers the density. Number densities of all nuclides other than fluorine are therefore
        independent of uf3_to_uf4_ratio.
        """
        mixture_molar_volume = 0.0
        ref_molar_volume = None
        for salt, percentage in self.components.items():
            molar_fraction = percentage / 100.0
            molar_mass = self.molar_masses[salt]

            if salt in self.SALT_COEFFICIENTS:
                coeffs = self.SALT_COEFFICIENTS[salt]
                density_pure = coeffs['A'] - coeffs['B'] * temperature
                if density_pure <= 0:
                    raise ValueError(f"Density of {salt} is non-positive at {temperature} K.")
                molar_volume_pure = molar_mass / density_pure
            elif self._formula_cation(salt) in self.LANTHANIDE_SYMBOLS:
                # Fallback for lanthanides not in SALT_COEFFICIENTS, using the LaF3 molar volume at this temperature.
                if ref_molar_volume is None:
                    ref_molar_volume = self._lanthanide_reference_molar_volume(temperature)
                molar_volume_pure = ref_molar_volume
            else:
                raise ValueError(f"Density coefficients not found for non-lanthanide salt: {salt}")

            mixture_molar_volume += molar_fraction * molar_volume_pure

        if mixture_molar_volume <= 0:
            raise ValueError("Mixture molar volume is non-positive.")

        pure_salt_density = self.mixture_molar_mass() / mixture_molar_volume
        final_density = pure_salt_density / (
                    1.0 - self.total_impurity_frac) if self.total_impurity_frac < 1.0 else float('inf')

        return final_density

    def get_atom_densities(self, temperature: float) -> dict:
        """Nuclide atom densities [atoms/barn-cm] at temperature [K]. Isomers are named e.g. 'Am-242m'."""
        bulk_density = self.density(temperature)
        salt_weight_fraction = 1.0 - self.total_impurity_frac
        salt_mass_density = bulk_density * salt_weight_fraction

        mixture_molar_mass = self.mixture_molar_mass()
        total_salt_number_density = (salt_mass_density / mixture_molar_mass) * self.AVOGADRO_NUMBER \
            if mixture_molar_mass > 0 else 0

        isotope_densities = defaultdict(float)

        def isotopic_distribution(symbol: str) -> dict:
            iso_dist = self.enrichments.get(symbol, None)
            if iso_dist is None:
                iso_data = self.ISOTOPIC_DATA[symbol.lower()]
                iso_dist = {mass_num: data['abundance'] for mass_num, data in iso_data.items()}
            return iso_dist

        for salt_formula, salt_percentage in self.components.items():
            elements_in_formula = re.findall(r'([A-Z][a-z]?)(\d*)', salt_formula)
            for symbol, count_str in elements_in_formula:
                atoms_per_molecule = int(count_str) if count_str else 1
                element_total_density = total_salt_number_density * (salt_percentage / 100.0) * atoms_per_molecule

                for mass_num, fraction in isotopic_distribution(symbol).items():
                    if fraction > 0:
                        isotope_name = f"{symbol}-{mass_num}"
                        isotope_densities[isotope_name] += element_total_density * fraction

        # One fluorine less per U3+. The mixture molar mass above already excludes this fluorine.
        fluorine_density_reduction = total_salt_number_density * self._fluorine_removed_per_mole()
        if fluorine_density_reduction > 0:
            for mass_num, fraction in isotopic_distribution('F').items():
                if fraction > 0:
                    isotope_name = f"F-{mass_num}"
                    isotope_densities[isotope_name] -= fluorine_density_reduction * fraction
                    if isotope_densities[isotope_name] < 0:
                        raise ValueError("Fluorine potential adjustment resulted in negative fluorine density.")

        for impurity_key, weight_frac in self.impurities_wt.items():
            impurity_mass_density = bulk_density * weight_frac
            parsed = _parse_nuclide(impurity_key)
            if parsed:
                symbol, mass_num, suffix = parsed
                try:
                    # Isomers use the ground-state mass; the excitation energy changes it by less than 1e-5 u.
                    isotope_mass = self.ISOTOPIC_DATA[symbol.lower()][mass_num]['mass']
                except KeyError:
                    raise ValueError(f"Isotopic data for impurity '{impurity_key}' not found.") from None
                impurity_atom_density = (impurity_mass_density / isotope_mass) * self.AVOGADRO_NUMBER
                isotope_densities[f"{symbol}-{mass_num}{suffix}"] += impurity_atom_density
            else:
                symbol = impurity_key
                avg_atomic_mass = self.elemental_masses.get(symbol)
                if avg_atomic_mass is None or avg_atomic_mass == 0:
                    raise ValueError(
                        f"Could not determine average atomic mass for elemental impurity '{symbol}'. An element "
                        f"without natural isotopic abundance must be given by nuclide, e.g. '{symbol}-<mass number>'.")

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

        composition_str = "-".join([f"{_format_percent(frac * 100)}%{salt}"
                                    for salt, frac in components.items() if frac > 0])

        # --- Assemble Enrichments ---
        enrichments = {'Li7_enr': Li7_enr}
        if u_enr_weight:
            enrichments['U'] = u_enr_weight

        if admixtures:
            for salt, props in admixtures.items():
                if 'enrichment' in props:
                    # Extract cation from the canonical salt formula (e.g. "luf3" -> "LuF3" -> "Lu", "KF" -> "K").
                    # The match is case-sensitive: with re.IGNORECASE, "KF" gave the cation "Kf".
                    cation = self._formula_cation(self._canonical_formula(salt))
                    if not cation:
                        raise ValueError(f"Could not parse cation from admixture salt: {salt}")
                    enrichments[cation] = props['enrichment']

        # --- Calculate Impurity Weight Fractions (Relative to Total Salt Mass) ---
        # We need a temporary object to calculate the mass of Uranium per mole of salt
        temp_salt = FluorideSalt(composition_str, enrichments=enrichments, uf3_to_uf4_ratio=UF3_to_UF4)
        avg_u_mass = temp_salt.elemental_masses.get('U', 0)

        # Calculate total mass of 1 mole of the salt mixture, without the fluorine removed by the UF3/UF4 ratio
        total_salt_mass_per_mole = temp_salt.mixture_molar_mass()

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
        composition_str = f"{_format_percent(lif_frac_eutectic * 100)}%LiF-{_format_percent(bef2_frac_eutectic * 100)}%BeF2"

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
