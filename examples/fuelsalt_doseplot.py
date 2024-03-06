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
plt.title('Fuel salt, 500 kg\nirradiated for 8 days at 1 MW')
plt.xlabel(f'Salt decay time [days]')
plt.ylabel('Dose at 100 cm [rem/h]')
plt.errorbar(x, dose, errd, ls='none', color=f'{colors["2"]}', capsize=1.2)
plt.scatter(x, dose, color=f'{colors["2"]}', s=5, label=f'gamma')

plt.legend()
plt.tight_layout()
plt.savefig(f'dose_g_fuelsalt_8MW.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_g_fuelsalt_8MW-loglog.png', dpi=1000)

# beta
dose = np.array([v["1"]['value'] for k, v in r.items()], float)
errd = np.array([v["1"]['stdev'] for k, v in r.items()], float)

# Plots!
plt.close('all')
plt.xscale('linear')
plt.yscale('linear')
plt.grid()
plt.title('Fuel salt, 500 kg\nirradiated for 8 days at 1 MW')
plt.xlabel(f'Salt decay time [days]')
plt.ylabel('Dose at 100 cm [rem/h]')
plt.errorbar(x, dose, errd, ls='none', color=f'{colors["1"]}', capsize=1.2)
plt.scatter(x, dose, color=f'{colors["3"]}', s=5, label=f'neutron')

plt.legend()
plt.tight_layout()
plt.savefig(f'dose_n_fuelsalt_8MW.png', dpi=1000)

plt.xscale('log')
plt.yscale('log')
plt.savefig(f'dose_n_fuelsalt_8MWe-loglog.png', dpi=1000)
plt.close()
