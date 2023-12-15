#!/bin/env python3
"""
Example use case of SampleDose - calculate doses as a function of decay times using F33 file
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5

r = {}
d = {}
for decay_days in np.geomspace(1./24., 360, 60):
    irr = SampleDose.OrigenIrradiation('../SCALE_FILE.mix0007.f33', 1.0)
    irr.set_decay_days(decay_days)  # 1 hour
    irr.irradiate_days = 2.0 * 365.24  # 2 years
    irr.irradiate_flux = 1e13  # n/s/cm2
    irr.write_atom_dens()
    irr.run_irradiate_decay_sample()

    mavric = SampleDose.DoseEstimator(irr)
    mavric.run_mavric()
    mavric.get_responses()

    print(mavric.responses)
    # print(mavric.total_dose)

    r[decay_days] = mavric.responses
    d[decay_days] = mavric.total_dose

print(r)
print(d)

with open('responses.json', 'w') as fout:
    json5.dump(r, fout, indent=4)

with open('doses.json', 'w') as fout:
    json5.dump(d, fout, indent=4)
