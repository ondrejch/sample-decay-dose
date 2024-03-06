#!/bin/env python3
"""
Example use case of SampleDose -- calculate decay doses of F71 sample as a function of decay times
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5

F71_2y_burn_file_name: str = '../SCALE_FILE.f71'            # 2 EFPY burn F71 file
F71_60day_burn_file_name: str = '../SCALE_FILE_60days.f71'  # 60 EFPD burn F71 file

r = {}
d = {}
for decay_days in np.geomspace(1./24., 360, 60):
    print(f'***--> Sample decay time {decay_days} days <--***')
    origen_triton = SampleDose.OrigenFromTriton(F71_60day_burn_file_name, 500e3)
    origen_triton.set_f71_pos(30.0 * 24.0 * 60.0 * 60.0)  # 7 days
    origen_triton.read_burned_material()
    origen_triton.set_decay_days(decay_days)
    origen_triton.run_decay_sample()

    mavric = SampleDose.DoseEstimator(origen_triton)
    mavric.det_x = 100.0  # 100 cm = 1 meter away
    mavric.N_planes_box = 10
    mavric.run_mavric()
    mavric.get_responses()

    print(mavric.responses)
    r[decay_days] = mavric.responses
    d[decay_days] = mavric.total_dose

print(r)
print(d)

with open('responses.json', 'w') as fout:
    json5.dump(r, fout, indent=4)

with open('doses.json', 'w') as fout:
    json5.dump(d, fout, indent=4)
