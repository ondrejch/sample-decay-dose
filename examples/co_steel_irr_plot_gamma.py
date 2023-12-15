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
labels = { 'irr1': ['decay time', 'days', 'SS-316 coupon, 2 year irradiation at 1e13 n/s/cm2, 1 gram'], }

#particles = {'1': 'Neutron', '2': 'Gamma', '3': 'Beta'}
particles = { '2': 'Gamma'}
# colors = {'1': 'sandybrown', '2': 'slategrey', '3': 'cornflowerblue', '12' : 'crimson'}
data = { 'co_0.02pct' : 'sandybrown', 'co_0.05pct' : 'slategrey', 'co_0.10pct' : 'cornflowerblue',
           'co_0.20pct' :  'crimson', 'co_0.40pct': 'orange' }

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses
r = {}

p ='2'    # Gamma only
for d in data.keys():
    with open(d+'/responses.json') as fin:
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

"""
p='2'
c1='12'
c2='2'
plt.errorbar(x, dose1[p], errd1[p], ls='none', color=f'{colors[c1]}', capsize=1.2)
plt.scatter(x, dose1[p], color=f'{colors[c1]}', s=5, label=f'{particles[p]}, SS316 with 0.1% Co')
plt.errorbar(x, dose2[p], errd2[p], ls='none', color=f'{colors[c2]}', capsize=1.2)
plt.scatter(x, dose2[p], color=f'{colors[c2]}', s=5, label=f'{particles[p]}, SS316')
"""

for d in data.keys():
    ytitle = d.replace('co_','').replace('pct','')
    plt.errorbar(x, dose[d], errd[d], ls='none', color=f'{data[d]}', capsize=0.8)
    plt.scatter(x, dose[d], color=f'{data[d]}', s=5, label=f'SS316 with {ytitle} % Co')

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.savefig(f'dose_Co_{label_file_name}.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_Co_{label_file_name}-loglog.png', dpi=1000)


# neutron_dose    = np.array([v['1']['value'] for k, v in r.items()])
# neutron_dose_e  = np.array([v['1']['stdev'] for k, v in r.items()])
# gamma_dose      = np.array([v['2']['value'] for k, v in r.items()])
# gamma_dose_e    = np.array([v['2']['stdev'] for k, v in r.items()])
# beta_dose       = np.array([v['3']['value'] for k, v in r.items()])
# beta_dose_e     = np.array([v['3']['stdev'] for k, v in r.items()])
