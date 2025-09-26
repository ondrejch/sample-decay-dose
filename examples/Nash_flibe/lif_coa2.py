#!/bin/env python3
"""
Example use case of SampleDose - irradiation of FLi7Be
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import json5


my_mass: float = 0.01       # grams
my_flux: float = 2e11       # n/cm2/s
irr_time: float = 10.0      # seconds
dec_time: float = 600.0     # seconds
my_adens: dict = {
    'li-7': 0.000232188933773336,
    'li-6': 0.000214502378174925,
    'na-23': 1.76179347295896e-05,
    'ca-40': 5.50092293681612e-10,
    'mg-24': 4.22425489867347e-10,
    'al-27': 2.74002513007875e-10,
    'mg-26': 1.78555877988959e-10,
    'ag-107': 3.81917638712822e-11,
    'ag-109': 3.76233657952604e-11,
    'mg-25': 3.49539722799673e-11,
    'ba-138': 3.46882513576315e-11,
    'rb-85': 3.45855779328456e-11,
    'ca-44': 2.03406951154959e-11,
    'rb-87': 9.08985488937452e-12,
    'sr-88': 7.84372406451334e-12,
    'ba-137': 5.67595814238421e-12,
    'ba-136': 5.4180754462137e-12,
    'ba-135': 3.78860148636177e-12,
    'ca-42': 3.17983976961746e-12,
    'ba-134': 2.81933438448962e-12,
    'ca-48': 1.16590905184715e-12,
    'sr-86': 8.14863429548865e-13,
    'ca-43': 6.77705700370867e-13,
    'sr-87': 5.88271150581067e-13,
    'ba130': 4.81129908278353e-13,
    'ba-132': 5.11320185338408e-14,
    'sr-84': 4.87200347974419e-14,
    'ca-46': 3.8490388071599e-14
}
seconds_to_days: float = 1.0 / (60.0 * 60.0 * 24.0)

irr = SampleDose.OrigenIrradiation('B.f33', my_mass)
irr.irradiate_days = irr_time * seconds_to_days
irr.irradiate_flux = my_flux
irr.set_decay_days(dec_time * seconds_to_days)
irr.write_atom_dens(my_adens)
irr.run_irradiate_decay_sample()

""" Bare sample """
mavric_bare = SampleDose.HandlingContactDoseEstimatorGenericTank(irr)
mavric_bare.cyl_r = SampleDose.get_cyl_r(mavric_bare.sample_volume)
mavric_bare.sample_h2 = mavric_bare.cyl_r
mavric_bare.layers_mats = []
mavric_bare.layers_thicknesses = []
mavric_bare.layers_temperature_K = []
mavric_bare.run_mavric()
mavric_bare.get_responses()
print(mavric_bare.responses)

""" Sample in a 1mm thick steel vial """
mavric_vial = SampleDose.HandlingContactDoseEstimatorGenericTank(irr)
mavric_vial.cyl_r = SampleDose.get_cyl_r(mavric_vial.sample_volume)
mavric_vial.sample_h2 = mavric_vial.cyl_r
mavric_vial.layers_mats = [SampleDose.ADENS_SS316H_COLD]    # Stainless steel 316
mavric_vial.layers_thicknesses = [0.1]                      # 1 mm thick
mavric_vial.layers_temperature_K = [300.0]                  # Room temp
mavric_vial.run_mavric()
mavric_vial.get_responses()
print(mavric_vial.responses)

# Print results
print("** Bare sample **")
print(f"Contact gamma dose: {mavric_bare.responses['2']['value'] * 1e3} ± {mavric_bare.responses['2']['stdev'] * 1e3} mrem/h")
print(f"Handling gamma dose: {mavric_bare.responses['6']['value'] * 1e3} ± {mavric_bare.responses['6']['stdev'] * 1e3} mrem/h")

print("** Sample in a steel vial **")
print(f"Contact gamma dose: {mavric_vial.responses['2']['value'] * 1e3} ± {mavric_vial.responses['2']['stdev'] * 1e3} mrem/h")
print(f"Handling gamma dose: {mavric_vial.responses['6']['value'] * 1e3} ± {mavric_vial.responses['6']['stdev'] * 1e3} mrem/h")



