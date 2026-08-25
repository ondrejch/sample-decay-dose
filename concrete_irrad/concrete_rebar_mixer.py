#!/bin/env python3
"""
Homogenized atom densities of a concrete-rebar layer.

Computes the volume-averaged material for a horizontal slice that contains one
mat of reinforcement bars. The slice thickness equals the bar diameter, so the
mat fits inside a single homogenized layer.
"""
import math
from collections import defaultdict

from sample_decay_dose.data import ISOTOPIC_DATA

AVOGADRO_NUMBER: float = 6.02214076e23  # atoms/mol
BARN_CM_CONVERSION: float = 1e24  # cm^2/barn

# ASTM A615 Grade 60 reinforcing steel, nominal elemental weight fractions.
REBAR_STEEL_WT_FRACTIONS: dict = {'fe': 0.98, 'mn': 0.012, 'si': 0.005, 'c': 0.003}
REBAR_STEEL_DENSITY: float = 7.85  # g/cm^3

# Co-59 impurity of the rebar drives the Co-60 source term. Design band is
# 0.1--0.4 wt%; measured carbon steels are typically lower (93--151 ppm,
# ORNL/SPR-2020/1586), so values in this band are conservative.
STEEL_CO59_WT_FRACTION_DEFAULT: float = 0.0025
STEEL_CO59_WT_FRACTION_RANGE: tuple = (0.001, 0.004)

# Ordinary concrete per ANSI/ANS-6.4.3, elemental weight fractions.
CONCRETE_WT_FRACTIONS: dict = {'h': 0.0056, 'o': 0.4985, 'na': 0.0152, 'mg': 0.0012,
    'al': 0.0636, 'si': 0.3041, 'ca': 0.0811, 'k': 0.0193, 's': 0.0009, 'fe': 0.0105}
CONCRETE_DENSITY: float = 2.30  # g/cm^3

# Trace impurities in concrete that dominate long-lived activation (Co-60 from
# Co-59; Eu-152/154 from natural Eu-151/153). Surveyed ordinary concretes span
# roughly Co 0.16-21.9 ppm and Eu 0.05-1.08 ppm (Alhajali et al. 2016, citing
# Suzuki et al. 2001); defaults sit mid-range and are overridable.
CONCRETE_IMPURITIES_WT_DEFAULT: dict = {'co': 10e-6, 'eu': 1e-6}


def _normalize_composition(wt_fractions: dict) -> dict:
    """ Normalize element symbols to 'Xx' form and validate the fractions sum to 1 """
    composition: dict = {}
    for symbol, fraction in wt_fractions.items():
        normalized: str = str(symbol).strip().capitalize()
        if not normalized.isalpha() or len(normalized) > 2:
            raise ValueError(f"Invalid element symbol: '{symbol}'")
        if fraction < 0:
            raise ValueError(f"Negative weight fraction for '{symbol}': {fraction}")
        composition[normalized] = composition.get(normalized, 0.0) + float(fraction)
    total: float = sum(composition.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"Weight fractions must sum to 1.0, got {total}")
    return {symbol: fraction / total for symbol, fraction in sorted(composition.items())}


def _merge_impurities(base_wt_fractions: dict, impurities_wt: dict) -> dict:
    """ Scale base fractions by the impurity remainder and add impurity entries """
    total_impurity: float = sum(float(fraction) for fraction in impurities_wt.values())
    if total_impurity < 0 or total_impurity >= 1.0:
        raise ValueError(f"Total impurity weight fraction must be in [0, 1), got {total_impurity}")
    merged: dict = {symbol: fraction * (1.0 - total_impurity)
                    for symbol, fraction in base_wt_fractions.items()}
    for symbol, fraction in impurities_wt.items():
        merged[symbol] = merged.get(symbol, 0.0) + float(fraction)
    return merged


class ConcreteRebarMixer:
    """
    Homogenizes one concrete-rebar layer into a single material.

    The layer contains bars running along x on a grid with spacing `spacing_x_cm`
    and bars running along y with spacing `spacing_y_cm`, both within a slice of
    height equal to the bar diameter `rebar_diameter_cm`. Steel volume fraction
    follows from the unit cell of area spacing_x * spacing_y and height d:

        f_steel = pi * d * (spacing_x + spacing_y) / (4 * spacing_x * spacing_y)

    Activation-relevant trace impurities are folded into both constituents: the
    rebar carries Co-59 within a configurable design band, and the concrete
    carries Co and Eu at surveyed trace levels. Outputs are nuclide-level atom
    densities in atoms/barn-cm suitable for ORIGEN/MAVRIC material input.
    """

    ISOTOPIC_DATA: dict = ISOTOPIC_DATA

    def __init__(self, rebar_diameter_cm: float, spacing_x_cm: float, spacing_y_cm: float,
                 steel_density_g_cc: float = REBAR_STEEL_DENSITY,
                 steel_wt_fractions: dict | None = None,
                 steel_co59_wt_fraction: float = STEEL_CO59_WT_FRACTION_DEFAULT,
                 concrete_density_g_cc: float = CONCRETE_DENSITY,
                 concrete_wt_fractions: dict | None = None,
                 concrete_impurities_wt: dict | None = None):
        """
        Initializes the mixer geometry and constituent materials.

        Args:
            rebar_diameter_cm: Bar diameter in cm; also the mixed-layer thickness.
            spacing_x_cm: Center-to-center bar spacing along x in cm.
            spacing_y_cm: Center-to-center bar spacing along y in cm.
            steel_density_g_cc: Rebar mass density in g/cm^3.
            steel_wt_fractions: Element -> weight fraction; defaults to REBAR_STEEL_WT_FRACTIONS.
            steel_co59_wt_fraction: Co-59 impurity as weight fraction of the steel; must lie
                                    within STEEL_CO59_WT_FRACTION_RANGE.
            concrete_density_g_cc: Concrete mass density in g/cm^3.
            concrete_wt_fractions: Element -> weight fraction; defaults to CONCRETE_WT_FRACTIONS.
            concrete_impurities_wt: Element -> weight fraction of the total concrete; defaults
                                    to CONCRETE_IMPURITIES_WT_DEFAULT ({'co', 'eu'}).
        """
        for name, value in [('rebar_diameter_cm', rebar_diameter_cm),
                            ('spacing_x_cm', spacing_x_cm), ('spacing_y_cm', spacing_y_cm)]:
            if value <= 0:
                raise ValueError(f"{name} must be positive, got {value}")
        self.rebar_diameter_cm: float = float(rebar_diameter_cm)
        self.spacing_x_cm: float = float(spacing_x_cm)
        self.spacing_y_cm: float = float(spacing_y_cm)
        if steel_density_g_cc <= 0 or concrete_density_g_cc <= 0:
            raise ValueError("Densities must be positive")
        self.steel_density_g_cc: float = float(steel_density_g_cc)
        self.concrete_density_g_cc: float = float(concrete_density_g_cc)
        lo, hi = STEEL_CO59_WT_FRACTION_RANGE
        if not lo <= steel_co59_wt_fraction <= hi:
            raise ValueError(
                f"steel_co59_wt_fraction must be within {STEEL_CO59_WT_FRACTION_RANGE}, "
                f"got {steel_co59_wt_fraction}")
        self.steel_co59_wt_fraction: float = float(steel_co59_wt_fraction)
        self.concrete_impurities_wt: dict = {
            str(symbol).strip().lower(): float(fraction) for symbol, fraction in
            (concrete_impurities_wt if concrete_impurities_wt is not None
             else CONCRETE_IMPURITIES_WT_DEFAULT).items()}
        steel_composition: dict = _merge_impurities(
            steel_wt_fractions if steel_wt_fractions is not None else REBAR_STEEL_WT_FRACTIONS,
            {'co': self.steel_co59_wt_fraction})
        concrete_composition: dict = _merge_impurities(
            concrete_wt_fractions if concrete_wt_fractions is not None else CONCRETE_WT_FRACTIONS,
            self.concrete_impurities_wt)
        self.steel_wt_fractions: dict = _normalize_composition(steel_composition)
        self.concrete_wt_fractions: dict = _normalize_composition(concrete_composition)

    def steel_volume_fraction(self) -> float:
        """ Steel volume fraction of the mixed layer from the unit-cell geometry """
        f_steel: float = math.pi * self.rebar_diameter_cm * (self.spacing_x_cm + self.spacing_y_cm) \
            / (4.0 * self.spacing_x_cm * self.spacing_y_cm)
        if f_steel >= 1.0:
            raise ValueError(f"Rebar volume fraction {f_steel:.4f} >= 1; spacings too small "
                             f"for diameter {self.rebar_diameter_cm} cm")
        return f_steel

    def mixture_density(self) -> float:
        """ Mass density of the homogenized layer in g/cm^3 """
        f_steel: float = self.steel_volume_fraction()
        return f_steel * self.steel_density_g_cc + (1.0 - f_steel) * self.concrete_density_g_cc

    def get_weight_fractions(self) -> dict:
        """ Elemental mass fractions of the whole mixed layer, keyed by element symbol """
        f_steel: float = self.steel_volume_fraction()
        steel_mass: float = f_steel * self.steel_density_g_cc
        concrete_mass: float = (1.0 - f_steel) * self.concrete_density_g_cc
        total_mass: float = steel_mass + concrete_mass
        wt: dict = defaultdict(float)
        for symbol, fraction in self.steel_wt_fractions.items():
            wt[symbol] += steel_mass * fraction
        for symbol, fraction in self.concrete_wt_fractions.items():
            wt[symbol] += concrete_mass * fraction
        return {symbol: mass / total_mass for symbol, mass in sorted(wt.items())}

    def get_elemental_atom_densities(self) -> dict:
        """ Elemental number densities in atoms/barn-cm, keyed by element symbol """
        f_steel: float = self.steel_volume_fraction()
        merged: dict = defaultdict(float)
        for wt_fractions, mass_density in [
                (self.steel_wt_fractions, f_steel * self.steel_density_g_cc),
                (self.concrete_wt_fractions, (1.0 - f_steel) * self.concrete_density_g_cc)]:
            for symbol, density in wt_fractions_to_atom_densities(wt_fractions, mass_density).items():
                merged[symbol] += density
        return dict(sorted(merged.items()))

    def get_atom_densities(self) -> dict:
        """ Nuclide-level atom densities in atoms/barn-cm, keyed like 'fe-56' """
        return expand_to_nuclides(self.get_elemental_atom_densities())


def _elemental_mass(symbol: str) -> float:
    """ Natural-abundance-weighted atomic mass in g/mol """
    iso_data: dict = ISOTOPIC_DATA[symbol.lower()]
    return sum(data['mass'] * data['abundance'] for data in iso_data.values() if data['abundance'] > 0)


def wt_fractions_to_atom_densities(wt_fractions: dict, mass_density_g_cc: float) -> dict:
    """ Elemental number densities [atoms/barn-cm] from elemental weight fractions and mass density """
    densities: dict = {}
    for symbol, fraction in wt_fractions.items():
        try:
            molar_mass: float = _elemental_mass(symbol)
        except KeyError as exc:
            raise ValueError(f"No isotopic data for element '{symbol}'") from exc
        densities[str(symbol).strip().capitalize()] = \
            mass_density_g_cc * fraction / molar_mass * AVOGADRO_NUMBER / BARN_CM_CONVERSION
    return dict(sorted(densities.items()))


def expand_to_nuclides(elemental_densities: dict) -> dict:
    """ Split elemental number densities into natural-abundance nuclides keyed like 'fe-56' """
    nuclide_densities: dict = defaultdict(float)
    for symbol, density in elemental_densities.items():
        for mass_number, data in ISOTOPIC_DATA[symbol.lower()].items():
            if data['abundance'] > 0:
                nuclide_densities[f"{symbol.lower()}-{mass_number}"] += density * data['abundance']
    return dict(sorted(nuclide_densities.items()))


if __name__ == "__main__":
    mixer = ConcreteRebarMixer(rebar_diameter_cm=1.59, spacing_x_cm=30.48, spacing_y_cm=15.24)
    print(f"Steel volume fraction : {mixer.steel_volume_fraction():.5f}")
    print(f"Mixture density       : {mixer.mixture_density():.5f} g/cm^3")
    print(f"Steel Co-59           : {mixer.steel_co59_wt_fraction:.5%} wt")
    print(f"Concrete impurities   : "
          f"{ {k: f'{v:.2e}' for k, v in mixer.concrete_impurities_wt.items()} }")
    print("\nElemental weight fractions:")
    for symbol, fraction in mixer.get_weight_fractions().items():
        print(f"  {symbol:<3}: {fraction:.6f}")
    print("\nNuclide atom densities [atoms/barn-cm]:")
    for nuclide, density in mixer.get_atom_densities().items():
        print(f"  {nuclide:<8}: {density:.6e}")
