# FluorideSalt: physics changes of the 2026-10-21 review fixes

These three changes alter the computed atom densities and densities of
enriched salts. The remaining fixes in the same pass are input validation
and robustness only.

## 1. Density of enriched salts

The pure-component molar volume now uses the natural molar mass, because the
`SALT_COEFFICIENTS` density correlations describe the natural salts.
Isotopic enrichment leaves the molar volume unchanged. An enriched salt keeps
the number densities of the natural salt, and its mass density scales with the
enriched molar mass.

Before the fix the pure molar volume was `M_enriched / rho_natural`, so a Li-7
enriched FLiBe had about 0.16% lower Be and F number densities than natural
FLiBe at the same temperature.

## 2. Isotopic split of elemental impurities

An elemental impurity of an enriched element now follows the enrichment given
for that element. Before the fix it took the enriched mass with the natural
isotopic split. Impurities of non-enriched elements still take the natural
abundance, through the same `isotopic_distribution()` path as the salt
components.

## 3. Generic formula parsing

A salt component that is not in `SALT_COEFFICIENTS` is now accepted as the
fluoride of one element, case-insensensitively: `puf3` resolves to `PuF3`,
`CSF` to `CsF`. Before the fix the generic fallback split the string into
elements, so `PUF3` resolved as P + U + F3, three fictitious components.
Strings that are not a fluoride of one element (two cations, an anion other
than fluorine, an unknown element) raise `ValueError`.

## Validation and robustness (no change of results for previously valid input)

- A uranium weight-fraction dictionary must list every isotope so the weights
  sum to 1. A single-isotope dictionary normalized to 100% of that isotope
  before the fix, which silently produced a U-235 source.
- `U235_enr` together with an explicit `U` entry raises; it overwrote the
  explicit entry before the fix.
- Negative impurity weight fractions raise.
- In `from_atom_densities()`, an anion below `ANION_IMPURITY_FLOOR` (1e-3 of
  the retained non-fluorine atoms) stays a trace impurity even with
  `impurity_cutoff=0`. A larger anion share raises, as before.
- The admixture enrichment dict is copied, so the caller's dict is not
  mutated.

## Deliberately kept

The BeF2 density slope stays at 1.45e-5 g/cm3/K. It is ten times smaller than
the 1.4771e-4 of some compilations, but ideal mixing of the pure-component
molar volumes underpredicts the Janz FLiBe correlation by 2.9% at 800 K, 1.8%
at 900 K and 0.8% at 1000 K with the current slope, and by about twice that
with the larger one. The comparison is documented at `SALT_COEFFICIENTS`. Pass
a corrected `SALT_COEFFICIENTS` entry if a measured value is preferred.

The `Li7_enr` convention is the weight fraction of Li-7 in lithium; vendor
specifications are usually atom percent, and 99.995 at% Li-7 is a weight
fraction of 0.999957.
