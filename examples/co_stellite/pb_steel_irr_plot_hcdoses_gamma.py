#!/bin/env python3
"""
Finds minimum lead thickness for shielding
Ondrej Chvala <ochvala@utexas.edu>
"""

import os
import re
import numpy as np
import json5
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve
from sample_decay_dose.SampleDose import extract_flux_values

scale_out: str = os.path.expanduser('~/0.02/80-upper-encl-stellite/01-triton/msrr.out')
flux_data = extract_flux_values(scale_out)
irradiation_flux: float = flux_data[8140]
cwd: str = os.getcwd()
steel_mass: float = float(re.findall(r'_([\d.]+)g', cwd)[0])
irradiation_years: float = float(re.findall(r'dose-([\d.]+)year_', cwd)[0])
print(f'steel flux: {irradiation_flux} n/cm2/s, mass: {steel_mass} g, irradiated for {irradiation_years} years')

LABEL = 'irr2'
labels = {'irr1': ['decay time', 'days',
    f'SS316 in Shed 80 pipe, {irradiation_years} y irradiation, {irradiation_flux:.1e} n/cm2/s, {steel_mass:.1f} g'],
    'irr2': ['lead shielding', 'cm',
        f'SS316 in Shed 80 pipe, {irradiation_years} y irradiation, {irradiation_flux:.1e} n/cm2/s, {steel_mass:.1f} g'], }
particles = {'2': 'Gamma, contact dose', '6': 'Gamma, 30cm handling dose'}
data = {'2': 'slategrey', '6': 'crimson'}

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses
r = {}

# p = ['2', '6']  # Gamma contact/handling
with open(os.path.join(cwd, 'responses.json')) as fin:
    r = json5.load(fin)
decay_days_last: str = list(r.keys())[-1]

pb_thick_list: list = list(r[decay_days_last].keys())
for d in data.keys():
    dose[d] = np.array([r[decay_days_last][pb_thick][d]['value'] for pb_thick in pb_thick_list], float)
    errd[d] = np.array([r[decay_days_last][pb_thick][d]['stdev'] for pb_thick in pb_thick_list], float)

x = np.array(pb_thick_list, float)  # x coordinate - sample masses

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title(labels[LABEL][2])
plt.xlabel(f'Sample {labels[LABEL][0]} [{labels[LABEL][1]}]')
plt.ylabel('Lead shield [cm]')


def exp_f(t: float, a: float, b: float, c: float) -> float:
    """ Exponential fit function """
    return a * np.exp(-b * t + c * t ** 2)


def root_exp_f(t: float, a: float, b: float, c: float, y0: float) -> float:
    """ For root finding """
    return exp_f(t, a, b, c) - y0


max_rem_per_h: float = 80e-3  # maximum rem/h dose
for d in data.keys():
    popt, pcov = curve_fit(exp_f, x, dose[d], p0=[1, 0.5, 1e-5], sigma=errd[d])
    x0: float = fsolve(root_exp_f, x0=5.0, args=(*popt, max_rem_per_h))[0]
    my_title = f'{particles[d]}, Pb thick {max_rem_per_h:.3f} rem/h: {x0:.1f} cm'
    print(my_title)
    plt.errorbar(x, dose[d], errd[d], ls='none', color=f'{data[d]}', capsize=0.8)
    plt.scatter(x, dose[d], color=f'{data[d]}', s=5, label=my_title)
    plt.plot(x, exp_f(x, *popt), ls='dotted', color=f'{data[d]}')

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.yscale('log')
plt.savefig(f'dose_SS316_Pb_{label_file_name}.png', dpi=1000)
plt.xscale('log')
plt.savefig(f'dose_SS316_Pb_{label_file_name}-loglog.png', dpi=1000)
