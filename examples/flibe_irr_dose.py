#!/bin/env python3
"""
Example use case of SampleDose - irradiation of FLi7Be
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5

'''atomflibe 1 1.95 3
         3000 2
         4000 1
         9000 4
         1.0 293.0 
         3007 99.995 3006 0.005 end'''

my_flibe = {'3006': 1.383015e-06,
            '3007': 2.371318e-02,
            '4009': 1.185729e-02,
            '9019': 4.742914e-02}
r = {}
d = {}
for decay_days in np.geomspace(1. / 24., 360, 60):
    irr = SampleDose.OrigenIrradiation('../SCALE_FILE.mix0001.f33', 1e3)  # 1 kg
    irr.set_decay_days(decay_days)
    irr.irradiate_days = 30.0  # 30 days
    irr.irradiate_flux = 5e6  # n/s/cm2
    irr.write_atom_dens(my_flibe)
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
