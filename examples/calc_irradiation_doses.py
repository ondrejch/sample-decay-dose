#!/bin/env python3
"""
Example use case of DoseF71 - simple run
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose

irr = SampleDose.OrigenIrradiation('../SCALE_FILE.mix0007.f33', 1.0)
irr.set_decay_days(1./24.)  # 1 hour
irr.irradiate_days = 2.0 * 365.24  # 2 years
irr.irradiate_flux = 2e12  # n/s/cm2
irr.write_atom_dens()
irr.run_irradiate_decay_sample()

mavric = SampleDose.DoseEstimator(irr)
mavric.run_mavric()
mavric.get_responses()

print(mavric.responses)
print(mavric.total_dose)
