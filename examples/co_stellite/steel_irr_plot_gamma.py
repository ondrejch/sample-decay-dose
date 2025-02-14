#!/bin/env python3
"""
Plotting script for calc_doses_mass
Ondrej Chvala <ochvala@utexas.edu>
"""

import os
import re
import numpy as np
import json5
import matplotlib.pyplot as plt
from sample_decay_dose.SampleDose import extract_flux_values

scale_out: str = os.path.expanduser('~/0.02/80-upper-encl-stellite/01-triton/msrr.out')
flux_data = extract_flux_values(scale_out)
stellite_flux: float = flux_data[9000]

cwd: str = os.getcwd()
steel_mass: float = float(re.findall(r'_([\d.]+)g', cwd)[0])
irradiation_years: float = float(re.findall(r'dose-([\d.]+)year_', cwd)[0])
print(f'steel flux: {stellite_flux} n/cm2/s, mass: {steel_mass} g, irradiated for {irradiation_years} years')

# For sample mass dependency, set LABEL='m'.
LABEL = 'irr1'
labels = {'irr1': ['decay time', 'days', f'SS316 {irradiation_years} years irradiation at {stellite_flux:.1e} n/cm2/s, {steel_mass:.1f} g'], }

# particles = {'1': 'Neutron', '2': 'Gamma', '3': 'Beta'}
particles = {'2': 'Gamma'}
data = {'.': 'slategrey',}
#    'co_0.05pct': 'slategrey', 'co_0.10pct': 'cornflowerblue', 'co_0.20pct': 'crimson', 'co_0.40pct': 'orange'}

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses
r = {}

p = '2'  # Gamma only
for d in data.keys():
    with open(os.path.join(d, 'responses.json')) as fin:
        r[d] = json5.load(fin)
        dose[d] = np.array([v[p]['value'] for k, v in r[d].items()], float)
        errd[d] = np.array([v[p]['stdev'] for k, v in r[d].items()], float)

x = np.array(list(r[list(data.keys())[0]].keys()), float)  # x coordinate - sample masses

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title(labels[LABEL][2])
plt.xlabel(f'Sample {labels[LABEL][0]} [{labels[LABEL][1]}]')
plt.ylabel('Dose at 30 cm [rem/h]')

for d in data.keys():
    ytitle = f'SS-316 cylinder, {steel_mass:.1f} g'
    plt.errorbar(x, dose[d], errd[d], ls='none', color=f'{data[d]}', capsize=0.8)
    plt.scatter(x, dose[d], color=f'{data[d]}', s=5, label=ytitle)

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.savefig(f'dose_SS316_{label_file_name}.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_SS316_{label_file_name}-loglog.png', dpi=1000)
