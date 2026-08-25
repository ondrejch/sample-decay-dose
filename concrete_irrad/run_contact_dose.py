#!/bin/env python3
"""
Driver for Erika: contact dose from the irradiated reactor-cavity concrete lid.

Edit the EDIT-ME block below (lid size, rebar grid, irradiation history), make
sure SCALE_BIN points at a SCALE installation, and run:

    PYTHONPATH=. python concrete_irrad/run_contact_dose.py --lib cavity_spectrum.f33

Add --scan 1,30,90,365 to compute the dose at several cooling times; results
land in lid_contact_doses.csv. Without --lib the script only previews decks.
"""
import argparse
import csv
import sys

from concrete_irrad.concrete_lid_contact_dose import ConcreteLidContactDose
from concrete_irrad.concrete_rebar_mixer import ConcreteRebarMixer

# ----------------------------- EDIT ME ------------------------------------
PARAMS: dict = {
    'name': 'lid',
    'length_cm': 100.0,          # lid span along x [cm]
    'width_cm': 100.0,           # lid span along y [cm]
    'bottom_thickness_cm': 6.0,  # concrete between cavity surface and rebar mat [cm]
    'irradiation_flux': 1.0e10,  # cavity neutron flux at the lid [n/cm2/s]
    'irradiation_days': 365.24,  # operating time before shutdown [d]
}
MIXER_PARAMS: dict = {
    'rebar_diameter_cm': 1.59,   # also the mixed-layer thickness [cm]
    'spacing_x_cm': 30.48,
    'spacing_y_cm': 15.24,
}
# ---------------------------------------------------------------------------


def collect_dose(calc: ConcreteLidContactDose) -> dict:
    """ Runs the full chain for one decay time and returns the contact dose row """
    calc.write_inputs()
    calc.run_activation()
    calc.run_mavric()
    calc.get_responses()
    dose: dict = calc.contact_dose
    return {'decay_days': calc.decay_days, 'dose_rem_per_h': dose['value'],
            'stdev_rem_per_h': dose['stdev']}


def scan_decay_times(days_list: list, lib_f33: str, nmpi: int = 1) -> list:
    """ Contact dose at each cooling time; ORIGEN is rerun per point """
    rows: list = []
    for days in days_list:
        calc = ConcreteLidContactDose(decay_days=days, irradiation_lib_f33=lib_f33,
                                      mixer=ConcreteRebarMixer(**MIXER_PARAMS), **PARAMS)
        rows.append(collect_dose(calc))
        print(f"decay {days:>8.2f} d : {rows[-1]['dose_rem_per_h']:.4e} "
              f"+/- {rows[-1]['stdev_rem_per_h']:.2e} rem/h")
    return rows


def write_csv(rows: list, path: str):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Concrete-lid contact dose driver")
    parser.add_argument('--lib', default='', help='Alpha-library F33 with the cavity neutron spectrum')
    parser.add_argument('--scan', default='', help='Comma-separated cooling times in days')
    parser.add_argument('--decay-days', type=float, default=30.0, help='Cooling time for a single run [d]')
    parser.add_argument('--csv', default='lid_contact_doses.csv', help='Output CSV for scans')
    parser.add_argument('--nmpi', type=int, default=1, help='MPI tasks per SCALE run')
    args = parser.parse_args()

    if not args.lib:
        preview = ConcreteLidContactDose(irradiation_lib_f33='cavity_spectrum.f33', **PARAMS)
        print(preview.origen_deck())
        print(preview.mavric_deck())
        print("# Preview only; pass --lib <f33> to execute.", file=sys.stderr)
        return

    if args.scan:
        days: list = [float(d) for d in args.scan.split(',')]
        rows: list = scan_decay_times(days, args.lib, nmpi=args.nmpi)
        write_csv(rows, args.csv)
        print(f"Wrote {args.csv}")
    else:
        calc = ConcreteLidContactDose(decay_days=args.decay_days,
                                      irradiation_lib_f33=args.lib, **PARAMS)
        row: dict = collect_dose(calc)
        print(f"\nContact photon dose after {row['decay_days']} d: "
              f"{row['dose_rem_per_h']:.4e} +/- {row['stdev_rem_per_h']:.2e} rem/h")


if __name__ == "__main__":
    main()
