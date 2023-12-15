#!/bin/env python3
"""
Example use case of SampleDose - simple decay doses of F71 sample
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose

origen_triton = SampleDose.OrigenFromTriton('../SCALE_FILE.f71', 1.0)
origen_triton.set_f71_pos(0.5 * 365.24 * 24.0 * 60.0 * 60.0)  # 0.5 years
origen_triton.read_burned_material()
origen_triton.run_decay_sample()

mavric = SampleDose.DoseEstimator(origen_triton)
mavric.run_mavric()
mavric.get_responses()
print(mavric.responses)
print(mavric.total_dose)
