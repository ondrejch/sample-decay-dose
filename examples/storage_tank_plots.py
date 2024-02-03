#!/bin/env python3
"""
Plotting script for storage_tank_doses_scan_parallel.py
Ondrej Chvala <ochvala@utexas.edu>
"""

import os
import numpy as np
import json5
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm, ticker
from numpy import ma

# os.chdir('')
my_data = {
    '1 day': '34-1_decay_days',
    '2 days': '33-2_decay_days',
    '7 days': '35-7_decay_days'
}


def make_plot(title:str, dir: str):
    with open(f'{dir}/doses.json') as fin:
        r = json5.load(fin)

    _steel_cm_list = []
    _concrete_cm_list = []
    for rdict in r:
        for s_cm, _r in rdict.items():
            if s_cm not in _steel_cm_list:
                _steel_cm_list.append(s_cm)
            for c_cm, _d in _r.items():
                if c_cm not in _concrete_cm_list:
                    _concrete_cm_list.append(c_cm)
    # print(len(_steel_cm_list), _steel_cm_list)
    # print(len(_concrete_cm_list), _concrete_cm_list)

    # Setup plotting data
    (x, y) = np.meshgrid(np.array(_steel_cm_list, float), np.array(_concrete_cm_list, float))
    g_mrem_dose = np.zeros((len(_steel_cm_list), len(_concrete_cm_list)))
    g_mrem_stdev = np.zeros((len(_steel_cm_list), len(_concrete_cm_list)))
    for rdict in r:
        for s_cm, _r in rdict.items():
            for c_cm, _d in _r.items():
                i: int = _steel_cm_list.index(s_cm)
                j: int = _concrete_cm_list.index(c_cm)
                g_mrem_dose[i, j] = _d['2']['value'] * 1000.0
                g_mrem_stdev[i, j] = _d['2']['stdev'] * 1000.0
    # print(g_mrem_dose)

    # Plot
    fig, ax = plt.subplots()
    lev_exp = np.arange(np.floor(np.log10(g_mrem_dose.min()) - 1),
                        np.ceil(np.log10(g_mrem_dose.max()) + 1))
    levs = np.power(10, lev_exp)
    # ticker.Locator.major_thresholds = (1, 0.1)
    # ticker.Locator.minor_thresholds = (0.5, 0.25)
    # https://matplotlib.org/stable/users/explain/colors/colormaps.html
    cs = ax.contourf(x, y, g_mrem_dose.T, levs, norm=colors.LogNorm(), locator=ticker.LogLocator(), cmap='jet')
    # cs = ax.contour(x, y, g_dose.T, levs)
    ax.clabel(cs, inline=True, fontsize=10, manual=False, colors=['black'], fmt= "%.0e")
    # cbar = fig.colorbar(cs, format = "%.1e")
    cbar = fig.colorbar(cs, format = "%.05g")
    cbar.set_label('Gamma dose [mrem/h]')
    plt.title(f'Gamma dose from 1.2 t of salt after {title} of decay')
    plt.xlabel('Steel shield thickness [cm]')
    plt.ylabel('Concrete shield thickness [cm]')
    file_name_fig=f'dose_g_storage_tank-{title.replace(" ","_")}_decay.png'
    plt.savefig(file_name_fig, dpi=1000)
    plt.show()


for t, d in my_data.items():
    make_plot(t, d)