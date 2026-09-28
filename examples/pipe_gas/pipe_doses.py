#!/bin/env python3
"""
Contact and handling doses from the MHA off-gas nuclides in a 1/16" OD SS-316H tube, as a function of decay time.
OrigenFromTritonMHA loads the MSRR fuel inventory from msrr.f71, retains the elements Te, I, Xe, Br, Kr and H
(processing block of its ORIGEN deck), and decays them.
MHATank scales the decayed gas to 2.24222e-6 lb (about 1.0 mg) and fills the bore of a 10 ft tube with it.
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
my_mass: float = 2.24222E-06 * 453.5924

r = {}
d = {}
for decay_days in np.linspace(0, 1, 24):
    irr = SampleDose.OrigenFromTritonMHA('../msrr.f71')
    irr.set_decay_days(decay_days)  # also names the ORIGEN case directory
    irr.run_decay_sample()

    mavric = SampleDose.MHATank(irr)
    mavric.cyl_r = pipe_ir
    mavric.layers_mats = [SampleDose.ADENS_SS316H_COLD]
    mavric.layers_thicknesses = [pipe_thick]
    mavric.layers_temperature_K = [300.0]
    # MHATank places the contact detector det_standoff_distance (0.1 cm) and the handling detector
    # handling_det_standoff_distance (30 cm) from the tube outer surface
    mavric.scale_decayed_pipe_material(my_mass)
    mavric.run_mavric()
    mavric.get_responses()

    print(mavric.responses)
    # print(mavric.total_dose)

    r[decay_days] = mavric.responses
    d[decay_days] = {'contact': mavric.contact_dose, 'handling': mavric.handling_dose}

print(r)
print(d)

with open('responses.json', 'w') as fout:
    json5.dump(r, fout, indent=4)

with open('doses.json', 'w') as fout:
    json5.dump(d, fout, indent=4)
