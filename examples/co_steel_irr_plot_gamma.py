#!/bin/env python3
"""
Plotting script for calc_doses_mass
Ondrej Chvala <ochvala@utexas.edu>
"""

import numpy as np
import json5
import matplotlib.pyplot as plt

# For sample mass dependency, set LABEL='m'.
LABEL = 'irr1'
labels = {'irr1': ['decay time', 'days', 'SS-316 coupon, 2 year irradiation at 1e13 n/s/cm2, 1 gram'], }

# particles = {'1': 'Neutron', '2': 'Gamma', '3': 'Beta'}
particles = {'2': 'Gamma'}
data = {'co_0.02pct': 'sandybrown', 'co_0.05pct': 'slategrey', 'co_0.10pct': 'cornflowerblue',
        'co_0.20pct': 'crimson', 'co_0.40pct': 'orange'}

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses
r = {}

p = '2'  # Gamma only
for d in data.keys():
    with open(d + '/responses.json') as fin:
        r[d] = json5.load(fin)
        dose[d] = np.array([v[p]['value'] for k, v in r[d].items()], float)
        errd[d] = np.array([v[p]['stdev'] for k, v in r[d].items()], float)

x = np.array(list(r['co_0.02pct'].keys()), float)  # x coordinate - sample masses

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title(labels[LABEL][2])
plt.xlabel(f'Sample {labels[LABEL][0]} [{labels[LABEL][1]}]')
plt.ylabel('Dose at 30 cm [rem/h]')

for d in data.keys():
    ytitle = d.replace('co_', '').replace('pct', '')
    plt.errorbar(x, dose[d], errd[d], ls='none', color=f'{data[d]}', capsize=0.8)
    plt.scatter(x, dose[d], color=f'{data[d]}', s=5, label=f'SS316 with {ytitle} % Co')

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.savefig(f'dose_Co_{label_file_name}.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_Co_{label_file_name}-loglog.png', dpi=1000)
