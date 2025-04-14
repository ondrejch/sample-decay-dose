#!/bin/env python3
"""
Irradiation of SS-316 with a specified wt% of cobalt in a steel pipe
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import os
import numpy as np
import json5

cwd: str = os.getcwd()

pipe_or: float = (1.0 / 16.0) * 2.54 / 2.0
pipe_ir: float = 0.0225 * 2.54 / 2.0
pipe_thick: float = pipe_or - pipe_ir
my_mass: float = 2.24222E-06 * 2.64101E-06 * 453.5924

r = {}
d = {}
for decay_days in np.linspace(0, 1, 24):
    irr = SampleDose.OrigenFromTritonMHA('../msrr.f71')
    irr.case_dir = 'run_pipe'
    irr.set_decay_days(decay_days)
    irr.run_decay_sample()

    mavric = SampleDose.MHATank(irr)
    mavric.cyl_r = pipe_ir
    mavric.layers_mats = [SampleDose.ADENS_SS316H_COLD]
    mavric.layers_thicknesses = [pipe_thick]
    mavric.layers_temperature_K = [300.0]
    mavric.det_x = mavric.cyl_r + pipe_thick + 30.0
    mavric.scale_decayed_pipe_material(my_mass)
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
