#!/bin/env python3
"""
Example use case of SampleDose - irradiation of Aluminum
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5

my_al = {'al-27': 6.030726e-02}  # Atom density [at/b-cm] of aluminium at 2.70 g/cm3

r = {}
d = {}
for decay_days in np.geomspace(1. / 24., 360, 60):
    irr = SampleDose.OrigenIrradiation('../SCALE_FILE.mix0001.f33', 50.0)  # 50g
    irr.set_decay_days(decay_days)
    irr.irradiate_days = 30.0  # 30 days
    irr.irradiate_flux = 5e6  # n/s/cm2
    irr.write_atom_dens(my_al)
    irr.run_irradiate_decay_sample()

    mavric = SampleDose.DoseEstimator(irr)
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
