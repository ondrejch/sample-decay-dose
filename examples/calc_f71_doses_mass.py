#!/bin/env python3
"""
Example use case of SampleDose -- calculate doses as a function of sample mass using F71 nuclides
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose

import numpy as np
import json5

r = {}
d = {}

for mass in np.geomspace(1e-3, 1, 30):
    print(f'***--> Sample mass {mass} g <--***')

    origen_triton = SampleDose.OrigenFromTriton('../SCALE_FILE.f71', mass)
    origen_triton.set_f71_pos(365.24 * 24.0 * 60.0 * 60.0)  # 1 year
    origen_triton.read_burned_material()
    origen_triton.run_decay_sample()

    mavric = SampleDose.DoseEstimator(origen_triton)
    mavric.run_mavric()
    mavric.get_responses()
    print(mavric.responses)

    r[mass] = mavric.responses
    d[mass] = mavric.total_dose

print(r)
print(d)

with open('responses.json', 'w') as fout:
    json5.dump(r, fout, indent=4)

with open('doses.json', 'w') as fout:
    json5.dump(d, fout, indent=4)
