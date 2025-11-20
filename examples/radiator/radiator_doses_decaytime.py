#!/bin/env python3
"""
Example radiator
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
import os
from sample_decay_dose.Radiator import RadiatorBox
from sample_decay_dose.SampleDose import ADENS_SS316H_HOT, ADENS_DRYAIR_COLD, OrigenFromTriton, DAY_IN_SECONDS
from joblib import Parallel, delayed, cpu_count
n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel

cwd: str = os.getcwd()
hotcell_inner_box2: float = 70.0 / 2.0  # 70 x 70 x 70 cm
hotcell_lead_shield: float = 7.5        # 7.5 cm of lead shielding
sample_mass: float = 0.1                # 0.1 g sample

r = {}
d = {}


def single_run() -> dict:
    burned_salt = OrigenFromTriton('./origen.f71', sample_mass)
    burned_salt.set_f71_pos(5.0 * 365.0 * DAY_IN_SECONDS, '1')
    burned_salt.read_burned_material()       # No decay

    radiator = RadiatorBox(burned_salt)
    radiator.layers_mats = [ADENS_SS316H_HOT, ADENS_DRYAIR_COLD]
    radiator.layers_temperature_K = [600.0, 300.0]
    radiator.N_planes_box = 15
    print(f"Radiator inside fluid total volume: {radiator.all_pins_volume} cm3")
    radiator.histories_per_batch = int(100000 * 1)
    radiator.histories_per_batch = int(1000 * 1)
    radiator.reuse_adjoint_flux = False  #True

    radiator.run_mavric()
    radiator.get_responses()
    return radiator.responses

# def dose_serial(decay_days) -> dict:
#     results = {}
#     for decay_day in decay_days:
#         results[decay_day] = single_run(decay_day)
#     return results
#
#
# def dose_parallel(decay_days) -> dict:
#     results = Parallel(n_jobs=n_jobs)(delayed(single_run)(decay_day) for decay_day in decay_days)
#     print(results)
#     return results


def main():
    mavric_res = single_run()
    print(f'Results: {mavric_res}')
    with open('my_response.json', 'w') as fout:
        json5.dump(mavric_res, fout, indent=4)


if __name__ == '__main__':
    main()
