#!/bin/env python3
"""
Example use case of DoseF71 -- calculate doses as a function of decay times
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import DoseF71
import numpy as np
import json5

r = {}
d = {}

for decay_days in np.geomspace(1./24., 360, 60):
    print(f'***--> Sample decay time {decay_days} days <--***')
#    de = DoseF71.DoseEstimator('../SCALE_FILE.f71', 0.1)  # 0.1 grams
#    de.set_f71_pos(2.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 2 years
    de = DoseF71.DoseEstimator('../SCALE_FILE_60days.f71', 0.1)  # 0.1 grams
    de.set_f71_pos(4.0 * 24.0 * 60.0 * 60.0)  # 4 days
    de.read_burned_material()
    de.set_decay_time(decay_days)
    de.run_decay_sample()
    de.run_mavric()
    de.get_responses()
    print(de.responses)
    r[decay_days] = de.responses
    d[decay_days] = de.total_dose

print(r)
print(d)

with open('responses.json', 'w') as fout:
    json5.dump(r, fout, indent=4)

with open('doses.json', 'w') as fout:
    json5.dump(d, fout, indent=4)
