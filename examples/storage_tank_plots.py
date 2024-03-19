#!/bin/env python3
"""
Plotting script for storage_tank_doses_scan_parallel.py
Ondrej Chvala <ochvala@utexas.edu>
"""

import os
import numpy as np
import pandas as pd
import json5
import matplotlib.pyplot as plt
from matplotlib import colors, cm, ticker

os.chdir('.')
my_data = {
    '1 day': '34-1_decay_days',
    '2 days': '33-2_decay_days',
    '7 days': '35-7_decay_days',
    '30 days': '36-30_decay_days'
}


def make_plot(title: str, my_dir: str):
    with open(f'{my_dir}/doses.json') as fin:
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
    ax.clabel(cs, inline=True, fontsize=10, manual=False, colors=['black'], fmt="%.0e")
    # cbar = fig.colorbar(cs, format = "%.1e")
    cbar = fig.colorbar(cs, format="%.05g")
    cbar.set_label('Gamma dose [mrem/h]')
    plt.title(f'Gamma dose from 1,200 kg of fuel salt after {title} of decay')
    plt.xlabel('Steel shield thickness [cm]')
    plt.ylabel('Concrete shield thickness [cm]')
    file_name_fig = f'dose_g_storage_tank-{title.replace(" ", "_")}_decay.png'
    plt.savefig(file_name_fig, dpi=1000, bbox_inches='tight', pad_inches=0.1)
    plt.tight_layout()
    # plt.show()

    steel_cm = [f'{float(x):.2f}' for x in _steel_cm_list]
    concrete_cm = [f'{float(x):.2f}' for x in _concrete_cm_list]
    pd_mrem_dose = pd.DataFrame(g_mrem_dose, columns=steel_cm, index=concrete_cm)
    pd_mrem_stdev = pd.DataFrame(g_mrem_stdev, columns=steel_cm, index=concrete_cm)
    return pd_mrem_dose, pd_mrem_stdev


for t, d in my_data.items():
    print(f'Processing {t}, {d}')
    make_plot(t, d)

    writer = pd.ExcelWriter('storage_dose.xlsx')
    for t, d in my_data.items():
        pd_dose, pd_stdev = make_plot(t, d)
        # header = f'Rows = Steel [cm], Columns = Concrete [cm]'
        # pd_dose.columns = pd.MultiIndex.from_product([[header],  pd_dose.columns])
        pd_dose.style.map(lambda v: 'color:#8B0000' if v > 20 else None). \
            map(lambda v: 'font-weight:bold;color:#008000' if 20 > v > 2 else None). \
            to_excel(writer, sheet_name=f'dose (mrem per h), {t}', float_format="%0.1f")
        pd_stdev.to_excel(writer, sheet_name=f'dose ± stdev, {t}', float_format="%0.2f")
    writer.close()
