#!/bin/env python3
"""
Irradiation of SS-316 with a specified wt% of cobalt in a steel pipe, handling and contact doses
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
from sample_decay_dose.SampleDose import extract_flux_values

import os
import re
import numpy as np
import json5
import argparse

parser = argparse.ArgumentParser(description='Irradiation of SS-316 with a specified wt% of cobalt.')
parser.add_argument('--wtpctCo', required=False, default=1.0, dest='wtpctCo', type=float,
                    help='weight percent of cobalt in SS-316, typically 0.02 to 2.0')
args = parser.parse_args()
w_co = args.wtpctCo * 1e-2  # convert to percent

cwd: str = os.getcwd()
steel_mass: float = float(re.findall(r'_([\d.]+)g', cwd)[0])
irradiation_years: float = float(re.findall(r'dose-([\d.]+)year_', cwd)[0])

scale_out: str = os.path.expanduser('~/0.02/80-upper-encl-stellite/01-triton/msrr.out')
flux_data = extract_flux_values(scale_out)
irradiation_flux: float = flux_data[8140]
print(f'Steel flux = {irradiation_flux} n/cm2/s, mass: {steel_mass} g, irradiate for {irradiation_years} years')


def cobalt_steel(wf_co: float = 0.1e-2) -> dict:
    """ Calculate atom density of SS-316 with cobalt """
    from MSRRpy.mat.material_types import Solid
    ss_composition = [
        ['c', 0.000800],
        ['mn', 0.020000],
        ['p', 0.000450],
        ['s', 0.000300],
        ['si', 0.010000],
        ['cr', 0.170000],
        ['ni', 0.120000],
        ['mo', 0.025000],
        ['fe', 0.653450]
    ]
    my_steel = Solid(name='StainlessSteel', temperature=873.0, composition=ss_composition,
                              composition_mode='weight', alpha=17.2e-6, ref_temperature=20.0 + 273.0, ref_density=8.0)
    my_steel.add_impurity(composition=[['co', 1.0]], composition_mode='weight',
                          fraction=wf_co, fraction_mode='weight')
    return dict(my_steel.mixing_table)


my_SS316_w_cobalt: dict = cobalt_steel(w_co)
pipe1_or: float = 1.315 * 2.54 / 2.0
pipe1_ir: float = 0.957 * 2.54 / 2.0
pipe1_thick: float = pipe1_or - pipe1_ir

r = {}
d = {}
for decay_days in np.geomspace(1. / 24., 30, 5):
    irr = SampleDose.OrigenIrradiation('../steel.f33', steel_mass)
    irr.set_decay_days(decay_days)
    irr.irradiate_days = irradiation_years * 365.24  #  years
    irr.irradiate_flux = irradiation_flux  # n/s/cm2
    irr.write_atom_dens(my_SS316_w_cobalt)
    irr.run_irradiate_decay_sample()

    r[decay_days] = {}
    d[decay_days] = {}
    for pb_shield in np.linspace(1,21,11):
        mavric = SampleDose.HandlingContactDoseEstimatorGenericTank(irr)
        mavric.cyl_r = pipe1_ir
        mavric.sample_h2 = SampleDose.get_cyl_h(mavric.sample_volume, mavric.cyl_r)
        mavric.layers_mats = [SampleDose.ADENS_SS316H_COLD, SampleDose.ADENS_LEAD_COLD]
        mavric.layers_thicknesses = [pipe1_thick, pb_shield]
        mavric.layers_temperature_K = [300.0, 300.0]
        mavric.run_mavric()
        mavric.get_responses()

        print(mavric.responses)

        r[decay_days][pb_shield] = mavric.responses
        d[decay_days][pb_shield] = mavric.total_dose

print(r)
print(d)

with open('responses.json', 'w') as fout:
    json5.dump(r, fout, indent=4)

with open('doses.json', 'w') as fout:
    json5.dump(d, fout, indent=4)
