#!/bin/env python3
"""
Plotting script for pipe_doses.py
Ondrej Chvala <ochvala@utexas.edu>
"""

import os
import numpy as np
import json5
import matplotlib.pyplot as plt
my_mass: float = 2.24222E-06 * 453.5924

cwd: str = os.getcwd()
LABEL = 'pipe'
labels = {'pipe': ['decay time', 'days', f'Offgas pipe, MHA nuclides, {my_mass:.4e} g'] }

# particles = {'1': 'Neutron', '2': 'Gamma', '3': 'Beta'}
particles = {'2': 'Gamma, contact dose', '6': 'Gamma, 30cm handling dose'}
data = {'2': 'slategrey', '6': 'crimson'}
#    'co_0.05pct': 'slategrey', 'co_0.10pct': 'cornflowerblue', 'co_0.20pct': 'crimson', 'co_0.40pct': 'orange'}

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses
r = {}

# p = ['2', '6']  # Gamma contact/handling
for d in data.keys():
    with open(os.path.join(cwd, 'responses.json')) as fin:
        r[d] = json5.load(fin)
        dose[d] = np.array([v[d]['value'] for _, v in r[d].items()], float)
        errd[d] = np.array([v[d]['stdev'] for _, v in r[d].items()], float)

xlist = list(r[list(data.keys())[0]].keys())
x = np.array(xlist, float)  # x coordinate - sample masses
closest_to_1month: str = min(xlist, key=lambda t: abs(float(t) - 30.0))
idx_1month: int = list(xlist).index(closest_to_1month)

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title(labels[LABEL][2])
plt.xlabel(f'Sample {labels[LABEL][0]} [{labels[LABEL][1]}]')
plt.ylabel('Dose [rem/h]')

for d in data.keys():
    my_title = f'{particles[d]}'  # , at {float(closest_to_1month):.1f} days = {dose[d][idx_1month]:.3f} rem/h'
    print(my_title)
    print(dose[d])
    plt.errorbar(x, dose[d], errd[d], ls='none', color=f'{data[d]}', capsize=0.8)
    plt.scatter(x, dose[d], color=f'{data[d]}', s=5, label=my_title)

plt.legend()
plt.tight_layout()
label_file_name = labels[LABEL][0].replace(' ', '_')
plt.savefig(f'dose_pipe_{label_file_name}.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_pipe_{label_file_name}-loglog.png', dpi=1000)


# def merge_doses():
#     import os
#     import json5
#     r: dict = {}
#     my_dirs: list = ['10-run_1d', '11-run_30d', '12-run_1h']
#     for my_dir in my_dirs:
#         with open(os.path.join(my_dir, 'responses.json')) as fin:
#             r.update(json5.load(fin))
#
#     with open('responses.json', 'w') as fout:
#         json5.dump(r, fout, indent=4)
