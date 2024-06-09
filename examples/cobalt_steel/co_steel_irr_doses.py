#!/bin/env python3
"""
Example use case of SampleDose - irradiation of SS-316 with a specified wt% of cobalt
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5
import argparse

parser = argparse.ArgumentParser(description='Irradiation of SS-316 with a specified wt% of cobalt.')
parser.add_argument('--wtpctCo', required=True, dest='wtpctCo', type=float,
                    help='weight percent of cobalt in SS-316, typically 0.02 to 2.0')
args = parser.parse_args()
w_co = args.wtpctCo * 1e-2  # convert to percent


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

r = {}
d = {}
for decay_days in np.geomspace(1. / 24., 360, 60):
    irr = SampleDose.OrigenIrradiation('../SCALE_FILE.mix0007.f33', 1.0)
    irr.set_decay_days(decay_days)
    irr.irradiate_days = 2.0 * 365.24  # 2 years
    irr.irradiate_flux = 1e13  # n/s/cm2
    irr.write_atom_dens(my_SS316_w_cobalt)
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
