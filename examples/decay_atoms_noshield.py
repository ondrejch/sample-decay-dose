#!/bin/env python3
"""
Example use case of SampleDose -- calculate decay doses of sample as a function of decay times.
The nuclide inventory is read from a CSV file of <nuclide>, <number of atoms> rows.

Usage:
    python decay_atoms_noshield.py --run-analysis --atoms-csv /tmp/atoms.csv
        runs the ORIGEN/MAVRIC jobs, writes doses_decay_atoms_noshield.json, and prints the doses
    python decay_atoms_noshield.py
        prints the doses stored in doses_decay_atoms_noshield.json
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose, utils
import argparse
import json5
from joblib import Parallel, delayed, cpu_count

n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel

cm2_to_barn: float = 1e24  # 1 cm^2 = 1e24 barn
my_inner_r: float = 37.7825  # IR of 30" schedule 40 pipe
my_thick: float = 0.635  # Thickness of 30" schedule 40 pipe
my_volume: float = 600e3  # 600 liters
my_atoms_file: str = '/tmp/atoms.csv'  # Default nuclide CSV file, override with --atoms-csv
my_doses_file: str = 'doses_decay_atoms_noshield.json'  # Written by run_analysis(), read by print_results()
decay_days_list = [1.0, 2.0, 7.0, 14.0, 30.0]


def print_atoms(my_atom_density: dict):
    """ Debugging """
    tot_atoms: float = 0.0
    tot_atom_density: float = 0.0
    for k, v in my_atom_density.items():
        tot_atom_density += v
        tot_atoms += v * cm2_to_barn * my_volume
    print(f'Total atoms: {tot_atoms}, total atom density {tot_atom_density} atoms / barn-cm')


def get_dose_parallel(decay_days: float, my_atom_density: dict) -> dict:
    origen_decay = SampleDose.OrigenDecayBox(my_atom_density, my_volume)
    origen_decay.set_decay_days(decay_days)  # also names the case directory with the decay time
    origen_decay.SAMPLE_F71_position = 900  # sample decay steps
    origen_decay.write_atom_dens()
    origen_decay.run_decay_sample()
    mavric = SampleDose.DoseEstimatorGenericTank(origen_decay)
    mavric.cyl_r = my_inner_r  # mavric_deck() sets the sample half-height from sample_volume and cyl_r
    # Material composition of additional layers, in dictionaries of atom densities
    mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD]
    # Thicknesses of additional layers [cm]
    mavric.layers_thicknesses = [my_thick, 2.0 * 2.54]
    # Temperatures of additional layers [cm]
    mavric.layers_temperature_K = [873.0, 300.0]
    # Add more planes since the source is large
    # mavric.N_planes_cyl = 20
    # Monaco histories
    mavric.histories_per_batch = 400000
    mavric.batches = 200

    # Run simulation
    mavric.run_mavric()
    mavric.get_responses()
    # Print doses
    print(f'Neutron dose {mavric.responses["1"]["value"]} +- {mavric.responses["1"]["stdev"]}  rem/h')
    print(f'Gamma dose   {mavric.responses["2"]["value"]} +- {mavric.responses["2"]["stdev"]}  rem/h')
    return mavric.responses


def run_analysis(atoms_file: str = my_atoms_file, datafile: str = my_doses_file):
    my_atom_density: dict = utils.read_cvs_atom_dens(atoms_file, my_volume)
    print_atoms(my_atom_density)
    # Parallel MAVRIC jobs
    results = Parallel(n_jobs=n_jobs)(delayed(get_dose_parallel)(d, my_atom_density) for d in decay_days_list)
    print(results)

    with open(datafile, 'w') as file_out:
        json5.dump(results, file_out, indent=4)


def print_results(datafile: str = my_doses_file):
    with open(datafile) as file_in:
        responses = json5.load(file_in)
    # run_analysis() writes one responses dict per entry of decay_days_list, in the same order
    if not isinstance(responses, list) or len(responses) != len(decay_days_list):
        raise ValueError(f'{datafile} must hold the list of {len(decay_days_list)} responses '
                         f'that run_analysis() writes for decay days {decay_days_list}')

    for d, rd in zip(decay_days_list, responses):
        print(f'Decay days {d:4.1f}: neutron dose {rd["1"]["value"]:8.1f} +- {rd["1"]["stdev"]:4.1f}  rem/h')
        print(f'Decay days {d:4.1f}: gamma dose {rd["2"]["value"]:8.1f} +- {rd["2"]["stdev"]:4.1f}  rem/h')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Offgas tank doses versus decay time, with the tank wall and insulation only.'
    )
    parser.add_argument('--run-analysis', action='store_true',
                        help='Run ORIGEN/MAVRIC jobs before printing the results.')
    parser.add_argument('--atoms-csv', default=my_atoms_file,
                        help='Nuclide CSV file (<nuclide>, <number of atoms> per row) for --run-analysis.')
    parser.add_argument('--datafile', default=my_doses_file,
                        help='JSON file with the dose results, written by --run-analysis.')
    args = parser.parse_args()

    if args.run_analysis:
        run_analysis(args.atoms_csv, args.datafile)
    print_results(args.datafile)
