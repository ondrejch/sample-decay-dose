#!/bin/env python3
"""
Example radiator: contact and handling doses from the salt in a tube-bundle radiator, read directly from one
F71 position (no decay step). Writes my_response.json.
Ondrej Chvala <ochvala@utexas.edu>
"""
import json5
from sample_decay_dose.Radiator import RadiatorBox
from sample_decay_dose.SampleDose import ADENS_SS316H_HOT, ADENS_DRYAIR_COLD, OrigenFromTriton, DAY_IN_SECONDS

sample_mass: float = 0.1                # 0.1 g sample; the radiator source is scaled to the tube volume


def single_run() -> dict:
    burned_salt = OrigenFromTriton('./origen.f71', sample_mass)
    # Case 2 starts at the end of the 5-year burn. Its F71 positions carry the photon and neutron spectra
    # that MAVRIC reads; the case-1 position at the same time stores concentrations only.
    burned_salt.set_f71_pos(5.0 * 365.0 * DAY_IN_SECONDS, '2')
    burned_salt.read_burned_material()       # No decay

    radiator = RadiatorBox(burned_salt)
    radiator.layers_mats = [ADENS_SS316H_HOT, ADENS_DRYAIR_COLD]
    radiator.layers_temperature_K = [600.0, 300.0]
    radiator.N_planes_box = 15
    print(f"Radiator inside fluid total volume: {radiator.all_pins_volume} cm3")
    radiator.histories_per_batch = int(100000 * 1)
    radiator.batches = int(50 * 1)
    radiator.reuse_adjoint_flux = False  #True

    radiator.run_mavric()
    radiator.get_responses()
    return radiator.responses


def main():
    mavric_res = single_run()
    print(f'Results: {mavric_res}')
    with open('my_response.json', 'w') as fout:
        json5.dump(mavric_res, fout, indent=4)


if __name__ == '__main__':
    main()
