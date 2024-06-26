#!/bin/env python3
"""
Example use case of SampleDose - irradiation of FLi7Be
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5

'''  Atom density from SCALE's mixing table
atomflibeLi7 1 1.95 3
         3000 2
         4000 1
         9000 4
         1.0 293.0
         3007 99.995 3006 0.005 end
atomflibeNat 1 1.95 3
         3000 2
         4000 1
         9000 4
         1.0 293.0
         3007 95.15 3006 4.85 end
         '''

my_flibe_li7 = {'li-6': 1.383015e-06,
                'li-7': 2.371318e-02,
                'be-9': 1.185729e-02,
                'f-19': 4.742914e-02}

my_flibe_nat = {'li-6': 1.332305e-03,
                'li-7': 2.240917e-02,
                'be-9': 1.187073e-02,
                'f-19': 4.748294e-02}

r = {}
d = {}
for decay_days in np.geomspace(1. / 24., 360, 60):
    irr = SampleDose.OrigenIrradiation('../SCALE_FILE.mix0001.f33', 1e3)  # 1 kg
    irr.set_decay_days(decay_days)
    irr.irradiate_days = 30.0  # 30 days
    irr.irradiate_flux = 5e6  # n/s/cm2
    irr.write_atom_dens(my_flibe_nat)
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
