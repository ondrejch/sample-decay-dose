#!/bin/env python3
"""
Example use case of SampleDose -- calculate decay doses of F71 sample as a function of decay times
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5
from joblib import Parallel, delayed, cpu_count

n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel
decay_days: float = 2.0  # 2 days

cm2_to_barn: float = 1e24  # 1 cm^2 = 1e24 barn
my_inner_r: float = 37.7825  # IR of 30" schedule 40 pipe
my_volume: float = 600e3  # 600 liters
my_atoms_file: str = 'atoms.csv'


def print_atoms():
    """ Debugging """
    tot_atoms: float = 0.0
    tot_atom_density: float = 0.0
    for k, v in my_atom_density.items():
        tot_atom_density += v
        tot_atoms += v * cm2_to_barn * my_volume
    print(f'Total atoms: {tot_atoms}, total atom density {tot_atom_density} atoms / barn-cm')


def mavric_process(case: tuple[float, float]) -> dict:
    """ Separating the MAVRIC part into a function for parallel execution """
    steel_cm: float
    concrete_cm: float
    (steel_cm, concrete_cm) = case
    # Calculate dose next to the tank
    mavric = SampleDose.DoseEstimatorGenericTank(origen_decay)
    mavric.cyl_r = my_inner_r
    mavric.sample_h2 = SampleDose.get_cyl_h(my_volume, my_inner_r)
    # Material composition of additional layers, in dictionaries of atom densities
    mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD, SampleDose.ADENS_SS316H_COLD,
                          SampleDose.ADENS_CONCRETE_COLD]
    # Thicknesses of additional layers [cm]
    mavric.layers_thicknesses = [2.54, 2.0 * 2.54, steel_cm, concrete_cm]
    # Temperatures of additional layers [cm]
    mavric.layers_temperature_K = [873.0, 300.0, 300.0, 300.0]
    # Add more planes since the source is large
    mavric.N_planes_cyl = 12
    # Monaco histories
    mavric.histories_per_batch = 200000
    mavric.batches = 40
    # Run simulation
    mavric.run_mavric()
    mavric.get_responses()
    # Print doses
    print(f'Neutron dose {mavric.responses["1"]["value"]} +- {mavric.responses["1"]["stdev"]}  rem/h')
    print(f'Gamma dose   {mavric.responses["2"]["value"]} +- {mavric.responses["2"]["stdev"]}  rem/h')
    _res: dict = {steel_cm: {}}
    _res[steel_cm][concrete_cm] = mavric.responses
    return _res


my_atom_density: dict = SampleDose.read_cvs_atom_dens(my_atoms_file, my_volume)
print_atoms()

origen_decay = SampleDose.OrigenDecayBox(my_atom_density, my_volume)
origen_decay.set_decay_days(decay_days)
origen_decay.SAMPLE_F71_position = 30  # sample decay steps
origen_decay.write_atom_dens()
origen_decay.run_decay_sample()

# Inputs for joblib parallelism have to be iterable
case_inputs: list[tuple[float, float]] = []
d = {}  # This would be better handled with Pandas ..
for steel_shield_thick_in in np.geomspace(0.1, 10, 8):
    s_cm = 2.54 * steel_shield_thick_in
    d[s_cm] = {}
    for concrete_shield_in in np.geomspace(0.1, 10, 8):
        c_cm = 2.54 * concrete_shield_in
        case_inputs.append((s_cm, c_cm))

# Parallel MAVRIC jobs
results = Parallel(n_jobs=n_jobs)(delayed(mavric_process)(case) for case in case_inputs)

print(results)
with open('doses.json', 'w') as file_out:
    json5.dump(results, file_out, indent=4)
