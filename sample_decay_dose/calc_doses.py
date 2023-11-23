#!/bin/env python3
"""
Example use case of DoseF71 - simple run
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import DoseF71

de = DoseF71.DoseEstimator('../SCALE_FILE.f71', 1.0)  # 1.0 g sample
de.read_burned_material()
de.run_decay_sample()
de.run_mavric()
de.get_responses()
print(de.responses)
print(de.total_dose)
