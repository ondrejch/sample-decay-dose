#!/bin/env python3
"""
Example use case of SampleDose.DoseEstimatorTank - simple decay doses of F71 sample.
Note that the beta dose is zero if the sample has additional shielding.
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
from sample_decay_dose import SampleDose

sample_mass: float = 1200.0e3  # 1200 kg
decay_days: float = 2.0  # 2 days

# Load and decay the nuclide vector from F71 file
# Set F71 file path and sample mass [g]
origen_triton = SampleDose.OrigenFromTriton('../SCALE_FILE.f71', sample_mass)

# Select F71 file position [seconds]
origen_triton.set_f71_pos(2.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 2 years

# Read burned nuclides
origen_triton.read_burned_material()

# Set hoe long they decay [days]
origen_triton.set_decay_days(decay_days)

# Execute ORIGEN
origen_triton.run_decay_sample()

d = {}  # This would be better handled with Pandas ..
for steel_shield_thick_in in np.geomspace(2, 15, 5):
    s_cm = 2.54 * steel_shield_thick_in
    d[s_cm] = {}
    for concrete_shield_in in np.geomspace(2, 15, 5):
        c_cm = 2.54 * concrete_shield_in

        # Calculate dose next to the tank
        mavric = SampleDose.DoseEstimatorStorageTank(origen_triton)

        # Material composition of additional layers, in dictionaries of atom densities
        mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD, SampleDose.ADENS_SS316H_COLD,
                              SampleDose.ADENS_CONCRETE_COLD]

        # Thicknesses of additional layers [cm]
        mavric.layers_thicknesses = [2.54, 2.0 * 2.54, s_cm, c_cm]

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

        d[s_cm][c_cm] = mavric.responses

with open('doses.json', 'w') as file_out:
    json5.dump(d, file_out, indent=4)
