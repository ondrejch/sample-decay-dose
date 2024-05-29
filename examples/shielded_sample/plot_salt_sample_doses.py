#!/bin/env python3
"""
Plots for shielded sample doses.
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
import matplotlib.pyplot as plt

low: float = 0.3
high: float = 0.7
shield: str = 'steel'
distance: float = 30.0
if shield == 'steel' and distance == 1.0:
    data: dict = {
        '15 minutes': {0.7: '01-fs_5y_0.7g_15min_decay',
                       0.3: '02-fs_5y_0.3g_15min_decay'},
        '30 minutes': {0.7: '11-fs_5y_0.7g_30min_decay',
                       0.3: '12-fs_5y_0.3g_30min_decay'},
        '2 hours': {0.7: '21-fs_5y_0.7g_2h_decay',
                    0.3: '22-fs_5y_0.3g_2h_decay'},
        '12 hours': {0.7: '23-fs_5y_0.7g_12h_decay',
                     0.3: '24-fs_5y_0.3g_12h_decay'},
        '24 hours': {0.7: '25-fs_5y_0.7g_24h_decay',
                     0.3: '26-fs_5y_0.3g_24h_decay'},
    }
elif shield == 'lead' and distance == 1.0:
    data: dict = {
        '15 minutes': {0.7: '51-lead-fs_5y_0.7g_15min_decay',
                       0.3: '52-lead-fs_5y_0.3g_15min_decay'},
        '30 minutes': {0.7: '53-lead-fs_5y_0.7g_30min_decay',
                       0.3: '54-lead-fs_5y_0.3g_30min_decay'},
        '2 hours': {0.7: '55-lead-fs_5y_0.7g_2h_decay',
                    0.3: '56-lead-fs_5y_0.3g_2h_decay'},
        '12 hours': {0.7: '57-lead-fs_5y_0.7g_12h_decay',
                     0.3: '58-lead-fs_5y_0.3g_12h_decay'},
        '24 hours': {0.7: '59-lead-fs_5y_0.7g_24h_decay',
                     0.3: '60-lead-fs_5y_0.3g_24h_decay'},
    }
if shield == 'steel' and distance == 30.0:
    data: dict = {
        '15 minutes': {0.7: '71-30cm_5y_0.7g_15min_decay',
                       0.3: '72-30cm_5y_0.3g_15min_decay'},
        '30 minutes': {0.7: '73-30cm_5y_0.7g_30min_decay',
                       0.3: '74-30cm_5y_0.3g_30min_decay'},
        '2 hours': {0.7: '75-30cm_5y_0.7g_2h_decay',
                    0.3: '76-30cm_5y_0.3g_2h_decay'},
        '12 hours': {0.7: '77-30cm_5y_0.7g_12h_decay',
                     0.3: '78-30cm_5y_0.3g_12h_decay'},
        '24 hours': {0.7: '79-30cm_5y_0.7g_24h_decay',
                     0.3: '80-30cm_5y_0.3g_24h_decay'},
    }
elif shield == 'lead' and distance == 30.0:
    data: dict = {
        '15 minutes': {0.7: '81-lead-30cm_5y_0.7g_15min_decay',
                       0.3: '82-lead-30cm_5y_0.3g_15min_decay'},
        '30 minutes': {0.7: '83-lead-30cm_5y_0.7g_30min_decay',
                       0.3: '84-lead-30cm_5y_0.3g_30min_decay'},
        '2 hours': {0.7: '85-lead-30cm_5y_0.7g_2h_decay',
                    0.3: '86-lead-30cm_5y_0.3g_2h_decay'},
        '12 hours': {0.7: '87-lead-30cm_5y_0.7g_12h_decay',
                     0.3: '88-lead-30cm_5y_0.3g_12h_decay'},
        '24 hours': {0.7: '89-lead-30cm_5y_0.7g_24h_decay',
                     0.3: '90-lead-30cm_5y_0.3g_24h_decay'},
    }
burn_years: str = '5'
# data: dict = {
#     '30 minutes': {0.7: '31-fs_1y_0.7g_30min_decay',
#                    0.3: '32-fs_1y_0.3g_30min_decay'},
# }
# burn_years: str = '1'
if shield == 'steel':
    my_colors: dict = {'hi': '#1995AD', 'lo': '#A1D6E2', 'fill': '#B1B1B2'}
elif shield == 'lead':
    my_colors: dict = {'hi': '#962E2A', 'lo': '#E3867D', 'fill': '#CEE6F2'}


def makefig(decay_time_str: str, d: dict) -> None:
    print(f"Processing decay time {decay_time_str}")
    l: str = d[low]
    h: str = d[high]
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
    plt.title(f"Gamma dose from the salt container after {decay_time_str}\n{burn_years} EFPY at 1 MWt")
    if shield == 'steel':
        plt.xlabel(f'SS-316 thickness [cm]')
    elif shield == 'lead':
        plt.xlabel(f'Lead thickness [cm]')
    plt.ylabel(f'Dose at {distance:.1f} cm [rem/h]')
    plt.errorbar(xh, yh, yherr, ls='none', color=my_colors['hi'], capsize=1.2)
    plt.scatter(xh, yh, color=my_colors['hi'], s=5, label='0.7 g sample')
    plt.errorbar(xl, yl, ylerr, ls='none', color=my_colors['lo'], capsize=1.2)
    plt.scatter(xl, yl, color=my_colors['lo'], s=5, label='0.3 g sample')
    plt.fill_between(xl, yl, yh, color=my_colors['fill'], alpha=0.2)

    plt.legend()
    plt.tight_layout()
    time_label: str = decay_time_str.replace(' ', '_')
    label_file_name = f'{distance:.0f}cm_{shield}_fs-{burn_years}y-dec_{time_label}'
    plt.savefig(f'dose_{label_file_name}.png', dpi=1000)
    # plt.show()
    #
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.savefig(f'dose_{label_file_name}-loglog.png', dpi=1000)
    # # plt.show()


for time_str, d in data.items():
    makefig(time_str, d)
