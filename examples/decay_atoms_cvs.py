#!/bin/env python3
"""
Example use case of SampleDose -- calculate decay doses of F71 sample as a function of decay times
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5

cm2_to_barn: float = 1e24  # 1 cm^2 = 1e24 barn
my_inner_r: float = 37.7825  # IR of 30" schedule 40 pipe
my_volume: float = 600e3  # 600 liters
my_atoms_file: str = '/tmp/atom.csv'


def print_atoms():
    """ Debugging """
    tot_atoms: float = 0.0
    tot_atom_density: float = 0.0
    for k, v in my_atom_density.items():
        tot_atom_density += v
        tot_atoms += v * cm2_to_barn * my_volume
    print(f'Total atoms: {tot_atoms}, total atom density {tot_atom_density} atoms / barn-cm')


my_atom_density: dict = SampleDose.read_cvs_atom_dens(my_atoms_file, my_volume)
print_atoms()

origen_decay = SampleDose.OrigenDecayBox(my_atom_density, my_volume)
origen_decay.set_decay_days(1.0)
origen_decay.SAMPLE_F71_position = 30  # sample decay steps
origen_decay.write_atom_dens()
origen_decay.run_decay_sample()

mavric = SampleDose.DoseEstimatorGenericTank(origen_decay)
mavric.cyl_r = my_inner_r
mavric.sample_h2 = SampleDose.get_cyl_h(my_volume, my_inner_r)
# Material composition of additional layers, in dictionaries of atom densities
mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD, SampleDose.ADENS_SS316H_COLD,
                      SampleDose.ADENS_CONCRETE_COLD]

# Thicknesses of additional layers [cm]
mavric.layers_thicknesses = [0.24 * 2.54, 1.0 * 2.54, 0.5*2.54, 2.54]

# Temperatures of additional layers [cm]
mavric.layers_temperature_K = [873.0, 300.0, 300.0, 300.0]

# Add more planes since the source is large
mavric.N_planes_cyl = 12

# Monaco histories
mavric.histories_per_batch = 150000
mavric.batches = 20

# Run simulation
mavric.run_mavric()
mavric.get_responses()

# Print doses
# print(mavric.responses)
print(f'Neutron dose {mavric.responses["1"]["value"]} +- {mavric.responses["1"]["stdev"]}  rem/h')
print(f'Gamma dose   {mavric.responses["2"]["value"]} +- {mavric.responses["2"]["stdev"]}  rem/h')
