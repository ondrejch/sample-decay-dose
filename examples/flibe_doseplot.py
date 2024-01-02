#!/bin/env python3
"""
Plotting gamma doses
Ondrej Chvala <ochvala@utexas.edu>
"""

import numpy as np
import json5
import matplotlib.pyplot as plt

with open('responses.json') as fin:
    r = json5.load(fin)

colors = {'1': 'sandybrown', '2': 'slategrey', '3': 'cornflowerblue'}
x = np.array(list(r.keys()), float)  # x coordinate
# Gamma
dose = np.array([v["2"]['value'] for k, v in r.items()], float)
errd = np.array([v["2"]['stdev'] for k, v in r.items()], float)

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title('1kg of Natural FLiBe\nirradiated for 30 days at 5e6 n/cm2/s')
plt.xlabel(f'Sample decay time [days]')
plt.ylabel('Dose at 30 cm [rem/h]')
plt.errorbar(x, dose, errd, ls='none', color=f'{colors["2"]}', capsize=1.2)
plt.scatter(x, dose, color=f'{colors["2"]}', s=5, label=f'gamma')

plt.legend()
plt.tight_layout()
plt.savefig(f'dose_g_NatFLiBe.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_g_NatFLiBe-loglog.png', dpi=1000)

# beta
dose = np.array([v["3"]['value'] for k, v in r.items()], float)
errd = np.array([v["3"]['stdev'] for k, v in r.items()], float)

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title('1kg of Natural FLiBe\nirradiated for 30 days at 5e6 n/cm2/s')
plt.xlabel(f'Sample decay time [days]')
plt.ylabel('Dose at 30 cm [rem/h]')
plt.errorbar(x, dose, errd, ls='none', color=f'{colors["3"]}', capsize=1.2)
plt.scatter(x, dose, color=f'{colors["3"]}', s=5, label=f'beta')

plt.legend()
plt.tight_layout()
plt.savefig(f'dose_b_NatFLiBe.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_b_NatFLiBe-loglog.png', dpi=1000)
plt.close()



