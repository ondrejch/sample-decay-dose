#!/bin/env python3
"""
Example use case of DoseF71
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import DoseF71
import numpy as np
import json5

# de = DoseF71.DoseEstimator('../SCALE_FILE.f71', 0.2)  # 0.2 g salt
# de.read_burned_material()
# de.run_decay_sample()
# de.run_mavric()
# de.get_responses()
# print(de.responses)

r = {}
for mass in np.geomspace(1e-3,1,30):
    de = DoseF71.DoseEstimator('../SCALE_FILE.f71', mass)  # 0.2 g salt
    de.read_burned_material()
    de.run_decay_sample()
    de.run_mavric()
    de.get_responses()
    print(de.responses)
    r[mass] = de.responses

print(r)

with open('doses.json','w') as fout:
    json5.dump(r, fout, indent=4)
