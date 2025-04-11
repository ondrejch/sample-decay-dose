#!/bin/env python3
"""
Example use case of SampleDose - irradiation of stellite
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
from sample_decay_dose.SampleDose import extract_flux_values
import numpy as np
import json5
import os
import re

cwd: str = os.getcwd()
stellite_mass: float = float(re.findall(r'_([\d.]+)g', cwd)[0])
irradiation_years: float = float(re.findall(r'dose-([\d.]+)year_', cwd)[0])

stellite_adens: dict = {'c-12': 0.001136285, 'c-13': 1.228975e-05, 'cr-50': 0.001152783, 'cr-52': 0.02223027,
    'cr-53': 0.002520734, 'cr-54': 0.000627464, 'co-59': 0.05488167, 'ni-58': 0.001454611, 'ni-60': 0.0005603136,
    'ni-61': 2.435644e-05, 'ni-62': 7.765899e-05, 'ni-64': 1.977746e-05, 'mo-92': 0.000405459, 'mo-94': 0.0002533775,
    'mo-95': 0.0004364792, 'mo-96': 0.0004578914, 'mo-97': 0.0002624366, 'mo-98': 0.0006640524, 'mo-100': 0.0002654562}

scale_out: str = os.path.expanduser('~/0.02/80-upper-encl-stellite/03-triton-longer/msrr.out')

flux_data = extract_flux_values(scale_out)
stellite_flux: float = flux_data[9000]
print(f'Stellite flux = {stellite_flux} n/cm2/s, mass: {stellite_mass} g, irradiate for {irradiation_years} years')

r = {}
d = {}
for decay_days in np.geomspace(1. / 24., 360, 60):
    irr = SampleDose.OrigenIrradiation('../stellite.f33', stellite_mass)
    irr.set_decay_days(decay_days)
    irr.irradiate_days = irradiation_years * 365.24
    irr.irradiate_flux = stellite_flux  # n/s/cm2
    irr.write_atom_dens(stellite_adens)
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
