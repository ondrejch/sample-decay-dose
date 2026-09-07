#!/usr/bin/env python3
"""Atom density of the human body, derived from elemental composition.

Prints a table of atom densities per element in the units:
  - atoms/cm^3          (number density, N)
  - at/(barn.cm)        (N x 1e-24; the factor that multiplies a
                         microscopic cross section sigma in barns to
                         give a dimensionless interaction probability
                         P = sigma [barn] * N [at/(barn.cm)] * x [cm])

Derivation (per element i):
  mass fraction w_i [g/g]
  -> moles of i per gram:  w_i / M_i          (M_i in g/mol)
  -> atoms of i per gram:  (w_i / M_i) * N_A
  -> number density:       (w_i / M_i) * N_A * rho   [atoms/cm^3]
  -> per barn.cm:          x 1e-24  (since 1 barn.cm = 1e-24 cm^3)

Sanity checks the script prints:
  atom fractions sum to 1; mean mass per atom ~6.5 u, i.e. ~5.1 atoms
  per nucleon, the standard reference-tissue value (ICRU-89 / MedPhy).
"""

# --- constants ----------------------------------------------------------

# Avogadro constant, exact since the 2019 SI redefinition (CODATA 2018,
# N_A = 6.02214076e23 mol^-1, no uncertainty).
N_A = 6.02214076e23

# 1 barn = 1e-24 cm^2 (defined unit in nuclear physics), so
# 1 barn.cm = 1e-24 cm^3. Hence N[at/(barn.cm)] = N[atoms/cm^3] * 1e-24.
BARN_CM2 = 1e-24

# Mean density of the adult human body ~0.98-1.10 g/cm^3 depending on
# composition; 1.05 g/cm^3 is a common reference-tissue value (water
# 1.00, lean muscle ~1.06, bone ~1.8). All results scale linearly with
# this; use ~1.8 for bone.
RHO_BODY = 1.05

# --- input data ----------------------------------------------------------
#
# BODY: element -> (mass fraction [g/g], atomic mass [g/mol])
#
# Mass fractions: ICRU Report 89 "Reference Values of Elemental
# Composition of Human Tissues" (2005), lean adult reference values:
#   O 65%, C 18.5%, H 9.5%, N 3.3%, Ca 1.5%, P 1.0%, K 0.35%,
#   S 0.25%, Cl 0.15%, Na 0.15%, Mg 0.05%  (sums to ~99.75%; the
#   remainder is trace elements and is normalized away here).
# Atomic masses: IUPAC/CIAAW standard atomic weights (2021) in g/mol.
BODY = {
    "O":  (0.650, 15.999),
    "C":  (0.185, 12.011),
    "H":  (0.095, 1.008),
    "N":  (0.033, 14.007),
    "Ca": (0.015, 40.078),
    "P":  (0.010, 30.974),
    "K":  (0.0035, 39.098),
    "S":  (0.0025, 32.065),
    "Cl": (0.0015, 35.453),
    "Na": (0.0015, 22.990),
    "Mg": (0.0005, 24.305),
}


def main():
    # Normalize mass fractions to exactly 1 (covers the ~0.25% of trace
    # elements not listed above).
    wsum = sum(w for w, _ in BODY.values())

    rows = []
    for el, (w, M) in BODY.items():
        w = w / wsum
        # Atoms per gram of tissue: (w_i / M_i) is mol/g, times N_A.
        n_per_g = w / M * N_A
        rows.append([el, w, M, n_per_g])

    n_total_per_g = sum(r[3] for r in rows)
    rows.sort(key=lambda r: r[3], reverse=True)

    header = ["el", "mass %", "atoms/cm^3", "at/(barn.cm)", "atom fraction"]
    widths = (max(len(header[0]), 4), 8, 12, 14, 13)
    w0, w1, w2, w3, w4 = widths

    def line(a):
        # a = [el, mass_fraction, M_g_per_mol, atoms_per_g]
        return " | ".join([
            f"{a[0]:<{w0}}", f"{a[1]*100:.2f}".rjust(w1),
            f"{a[3]*RHO_BODY:.3e}".rjust(w2),
            f"{a[3]*RHO_BODY*BARN_CM2:.4e}".rjust(w3),
            f"{a[3]/n_total_per_g:.4f}".rjust(w4),
        ])

    print("human body atom density (ICRU-89 soft tissue, rho = "
          f"{RHO_BODY} g/cm^3)")
    print(f"1 at/(barn.cm) = {1/BARN_CM2:.0e} atoms/cm^3")
    print()
    print(" | ".join([header[0].ljust(w0), header[1].rjust(w1),
                      header[2].rjust(w2), header[3].rjust(w3),
                      header[4].rjust(w4)]))
    print("-" * 60)
    for r in rows:
        print(line(r))
    print("-" * 60)
    print(line(["TOTAL", 1.0, None, n_total_per_g]))
    print()

    # Mean mass per atom = total mass / total moles = 1 / sum(w_i/M_i).
    # Expected ~6.5 u (i.e. ~0.154 mol atoms per g, ~5.1 atoms per
    # nucleon), matching the standard tissue value used in medical
    # physics (e.g. ICRU-89, MedPhy/Geant4 tissue definitions).
    mol_per_g = sum(w / M for w, M in BODY.values()) / wsum
    mean_u_per_atom = 1.0 / mol_per_g
    print(f"mean mass per atom: {mean_u_per_atom:.2f} u "
          f"({mol_per_g:.3f} atoms/u)")


if __name__ == "__main__":
    main()
