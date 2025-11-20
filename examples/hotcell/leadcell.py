#!/bin/env python3
"""
Example hotcell is 70 x 70 x 61 cm, close enough to 70 x 70 x 70 cm
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
hotcell_inner_box2: float = 70.0 / 2.0  # 70 x 70 x 70 cm
hotcell_lead_shield: float = 7.5        # 7.5 cm of lead shielding
sample_mass: float = 0.7                # 0.1 g sample

r = {}
d = {}
decay_days = np.linspace(1, 128, 128)
#decay_days = np.linspace(1, 128, 2)


def single_run(decay_day: float) -> dict:
    burned_salt = OrigenFromTriton('./msrr.f71', sample_mass)
    # Get the case number for last burn step with flux > 0
    n_last_fuel_burn: int = [k for k, v in burned_salt.BURNED_MATERIAL_F71_index.items()
                             if v['case'] == '20' and float(v['flux']) > 0][-1]
    burned_salt.BURNED_MATERIAL_F71_position = n_last_fuel_burn
    burned_salt.read_burned_material()
    burned_salt.set_decay_days(decay_day)
    burned_salt.run_decay_sample()

    mavric = HotCellDoses(burned_salt)
    mavric.layers_mats = [ADENS_DRYAIR_COLD, ADENS_LEAD_COLD]
    mavric.layers_thicknesses = [hotcell_inner_box2, hotcell_lead_shield]
    mavric.layers_temperature_K = [300.0, 300.0]
    mavric.N_planes_box = 15
    mavric.N_planes_cyl = 2
    mavric.histories_per_batch = int(100000 * np.sqrt(decay_day))
    mavric.reuse_adjoint_flux = True
    mavric.run_mavric()
    mavric.get_responses()

    print(mavric.responses)
    # print(mavric.total_dose)

    r[decay_day] = mavric.responses
    d[decay_day] = mavric.total_dose
    return mavric.responses


def dose_serial(decay_days):
    for decay_day in decay_days:
        single_run(decay_day)


def dose_parallel(decay_days) -> dict:
    results = Parallel(n_jobs=n_jobs)(delayed(single_run)(decay_day) for decay_day in decay_days)
    print(results)
    return results


def main():
#    dose_serial(decay_days)
    mavric_res = dose_parallel(decay_days)
    # Reconstruct r from parallel runs
    r = {decay_days[i]:mavric_res[i] for i in range(len(decay_days))}

    rs = dict(sorted(r.items()))
    print(rs)
    with open('responses.json', 'w') as fout:
        json5.dump(r, fout, indent=4)

#    ds = dict(sorted(d.items()))
#    print(ds)
#    with open('doses.json', 'w') as fout:
#        json5.dump(ds, fout, indent=4)


if __name__ == '__main__':
    main()
