#!/bin/env python3
"""
Approximate Th-EIRENE burn in ORIGEN [incomplete]
Ondrej Chvala <ochvala@utexas.edu>
"""
import os
import io
import subprocess
import json5
import numpy as np
from datetime import datetime
from sample_decay_dose import SampleDose


f71_file_name: str = 'ThEIRENE.f71'
f71_position: int = 1
f33_file_name: str = 'ThEIRENE.mix0001.f33'
f33_position: int = 5
fuel_volume: float = 1.06806e+07  # [cm^3]
MTiHM: dict = {' 5.00': 10.4159200656709, '19.75': 10.2647529843048}
fuel_type: dict = {' 5.00': 'LEU+Th', '19.75': 'HALEU+Th'}
flux_per_MW: dict = {' 5.00': 1.3070e+11, '19.75': 1.2441e+11}

SampleDose.ATOM_DENS_MINIMUM = 1e-20
fuel_salt_at_dens_BOC: dict = SampleDose.get_burned_nuclide_atom_dens(f71_file_name, 1)



