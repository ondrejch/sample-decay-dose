#!/bin/env python3
"""
Plotting script for calc_f71_doses_mass.py, calc_f71_doses_decaytime.py and calc_irradiation_doses_decaytime.py
Ondrej Chvala <ochvala@utexas.edu>
"""

import numpy as np
import json5
import matplotlib.pyplot as plt

# For sample mass dependency, set LABEL='m'.
# For decay time dependency, set LABEL='dt'
LABEL = 'irr1'
labels = {'m' : ['mass', 'grams', 'Fuel salt sample, 1 year burn, 30 days decay time'],
          'dt': ['decay time', 'days', 'Fuel salt sample, 1 year burn, 0.1 grams'],
         'dt2': ['decay time', 'days', 'Fuel salt sample, 2 year burn, 0.1 grams'],
        'dt04': ['decay time', 'days', 'Fuel salt sample, 4 day burn, 0.1 grams'],
        'dt12': ['decay time', 'days', 'Fuel salt sample, 12 day burn, 0.1 grams'],
        'dt28': ['decay time', 'days', 'Fuel salt sample, 28 day burn, 0.1 grams'],
        'irr1': ['decay time', 'days', 'SS-316 coupon, 2 year irradiation at 1e13 n/s/cm2, 1 gram'],
          }

# Responses file written by the calculation script that matches each LABEL
datafiles = {'m': 'responses_f71_mass.json',          # calc_f71_doses_mass.py
             'dt': 'responses_f71_decaytime.json',    # calc_f71_doses_decaytime.py
             'dt2': 'responses_f71_decaytime.json',
             'dt04': 'responses_f71_decaytime.json',
             'dt12': 'responses_f71_decaytime.json',
             'dt28': 'responses_f71_decaytime.json',
             'irr1': 'responses_irr_decaytime.json',  # calc_irradiation_doses_decaytime.py
             }

particles = {'1': 'Neutron', '2': 'Gamma', '3': 'Beta'}
colors = {'1': 'sandybrown', '2': 'slategrey', '3': 'cornflowerblue'}

with open(datafiles[LABEL]) as fin:
    r = json5.load(fin)

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses

# Parse dictionary into arrays
x = np.array(list(r.keys()), float)  # x coordinate - sample masses
for p in particles.keys():
    dose[p] = np.array([v[p]['value'] for k, v in r.items()], float)
    errd[p] = np.array([v[p]['stdev'] for k, v in r.items()], float)

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title(labels[LABEL][2])
plt.xlabel(f'Sample {labels[LABEL][0]} [{labels[LABEL][1]}]')
plt.ylabel('Dose 30 cm from the sample centre [rem/h]')
for p in particles.keys():
    if sum(dose[p]) > 1e-6:    # Only plot of there is dose form that particle
        plt.errorbar(x, dose[p], errd[p], ls='none', color=f'{colors[p]}', capsize=1.2)
        plt.scatter(x, dose[p], color=f'{colors[p]}', s=5, label=f'{particles[p]}')

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.savefig(f'dose_{label_file_name}.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_{label_file_name}-loglog.png', dpi=1000)


# neutron_dose    = np.array([v['1']['value'] for k, v in r.items()])
# neutron_dose_e  = np.array([v['1']['stdev'] for k, v in r.items()])
# gamma_dose      = np.array([v['2']['value'] for k, v in r.items()])
# gamma_dose_e    = np.array([v['2']['stdev'] for k, v in r.items()])
# beta_dose       = np.array([v['3']['value'] for k, v in r.items()])
# beta_dose_e     = np.array([v['3']['stdev'] for k, v in r.items()])
