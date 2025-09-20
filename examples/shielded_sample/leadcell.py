#!/bin/env python3
"""
Example hotcell is 70 x 70 x 61 cm, close enough to 70x70x70cm
7.5cm of lead shielding
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
import os
from sample_decay_dose.HotCell import HotCellDoses
from sample_decay_dose.SampleDose import ADENS_LEAD_COLD, ADENS_DRYAIR_COLD, OrigenFromTriton
from joblib import Parallel, delayed, cpu_count
n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel

cwd: str = os.getcwd()
hotcell_inner_box2: float = 70.0 / 2.0
hotcell_lead_shield: float = 7.5

r = {}
d = {}
for decay_days in np.linspace(1, 91, 45):
    burned_salt = OrigenFromTriton('./msrr.f71')
    # Get the case number for last burn step with flux > 0
    n_last_fuel_burn: int = [k for k, v in burned_salt.BURNED_MATERIAL_F71_index.items()
                             if v['case'] == '20' and float(v['flux']) > 0][-1]
    burned_salt.BURNED_MATERIAL_F71_position = n_last_fuel_burn
    burned_salt.read_burned_material()
    burned_salt.set_decay_days(decay_days)
    burned_salt.run_decay_sample()

    mavric = HotCellDoses(burned_salt)
    mavric.layers_mats = [ADENS_DRYAIR_COLD, ADENS_LEAD_COLD]
    mavric.layers_thicknesses = [hotcell_inner_box2, hotcell_lead_shield]
    mavric.layers_temperature_K = [300.0, 300.0]
    mavric.N_planes_box = 15
    mavric.N_planes_cyl = 2
    mavric.reuse_adjoint_flux = True
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
