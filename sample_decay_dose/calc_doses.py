#!/bin/env python3
"""
Example use case of DoseF71 - simple run
Ondrej Chvala <ochvala@utexas.edu>
"""

# from sample_decay_dose import DoseF71
#
# de = DoseF71.DoseEstimator('../SCALE_FILE.f71', 1.0)  # 1.0 g sample
# de.read_burned_material()
# de.run_decay_sample()
# de.run_mavric()
# de.get_responses()
# print(de.responses)
# print(de.total_dose)

from sample_decay_dose import SampleDose

origen_triton = SampleDose.OrigenFromTriton('../SCALE_FILE.f71', 1.0)
origen_triton.set_f71_pos(0.5 * 365.24 * 24.0 * 60.0 * 60.0)  # 0.5 years
origen_triton.read_burned_material()
origen_triton.run_decay_sample()

mavric = SampleDose.DoseEstimator(origen_triton)
mavric.run_mavric()
print(mavric.responses)
print(mavric.total_dose)
