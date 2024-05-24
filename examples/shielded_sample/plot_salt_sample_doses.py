#!/bin/env python3
"""
Plots for shielded talk
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
import matplotlib.pyplot as plt

low: float = 0.3
hight: float = 0.7
data: dict = {
    '15 minutes': {0.7: '01-fs_5y_0.7g_15min_decay',
                   0.3: '02-fs_5y_0.3g_15min_decay'},
    '30 minutes': {0.7: '11-fs_5y_0.7g_30min_decay',
                   0.3: '12-fs_5y_0.3g_30min_decay'},
    '2 hours': {0.7: '21-fs_5y_0.7g_2h_decay',
                0.3: '22-fs_5y_0.3g_2h_decay'},
}

my_colors: dict = {'hi': '#1995AD', 'lo': '#A1D6E2', 'fill': '#B1B1B2'}


def makefig(decay_time_str: str, d: dict) -> None:
    print(f"Processing decay time {decay_time_str}")
    l: str = d[low]
    h: str = d[hight]
    with open(l + '/doses.json') as f:
        gamma_doses_l = json5.load(f)
    with open(h + '/doses.json') as f:
        gamma_doses_h = json5.load(f)
    _x: list = []
    _y: list = []
    _yerr: list = []
    for dd in gamma_doses_l:
        for t, v in dd.items():
            _x.append(float(t))
            _y.append(float(v['value']))
            _yerr.append(float(v['stdev']))
    xl = np.array(_x)
    yl = np.array(_y)
    ylerr = np.array(_yerr)

    _x: list = []
    _y: list = []
    _yerr: list = []
    for dd in gamma_doses_h:
        for t, v in dd.items():
            _x.append(float(t))
            _y.append(float(v['value']))
            _yerr.append(float(v['stdev']))

    xh = np.array(_x)
    yh = np.array(_y)
    yherr = np.array(_yerr)

    np.testing.assert_array_equal(xh, xl)

    plt.close('all')
    plt.xscale('linear')
    plt.yscale('log')
    plt.grid()
    plt.title(f"Gamma dose from the salt container after {decay_time_str}\n5 EFPY at 1 MWt")
    plt.xlabel(f'SS-316 thickness [cm]')
    plt.ylabel('Dose at 1 cm [rem/h]')
    plt.errorbar(xh, yh, yherr, ls='none', color=my_colors['hi'], capsize=1.2)
    plt.scatter(xh, yh, color=my_colors['hi'], s=5, label='0.7 g sample')
    plt.errorbar(xl, yl, ylerr, ls='none', color=my_colors['lo'], capsize=1.2)
    plt.scatter(xl, yl, color=my_colors['lo'], s=5, label='0.3 g sample')
    plt.fill_between(xl, yl, yh, color=my_colors['fill'], alpha=0.2)

    plt.legend()
    plt.tight_layout()
    time_label: str = time_str.replace(' ', '_')
    label_file_name = f'fs-5y-dec_{decay_time_str}'
    plt.savefig(f'dose_{label_file_name}.png', dpi=1000)
    # plt.show()
    #
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.savefig(f'dose_{label_file_name}-loglog.png', dpi=1000)
    # # plt.show()


for time_str, d in data.items():
    makefig(time_str, d)
