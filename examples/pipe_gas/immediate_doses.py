#!/bin/env python3
"""
Prints dose at t=0
Ondrej Chvala <ochvala@utexas.edu>
"""

import os
import numpy as np
import json5
import matplotlib.pyplot as plt
my_mass: float = 2.24222E-06 * 453.5924

cwd: str = os.getcwd()
particles = {'1': 'Gamma, contact dose', '5': 'Gamma, 30cm handling dose',
    '2': 'Gamma, contact dose', '6': 'Gamma, 30cm handling dose'}

dose = {}  # doses [rem/h]
errd = {}  # stdev of doses
r = {}

with open(os.path.join(cwd, 'responses.json')) as fin:
    r = json5.load(fin)

d0 = r['0.0']
for time in r.keys():
    minutes: float = float(time) *60*24
    print(f'=== time: {minutes:.1f} minutes ===')
    d0 = r[time]
    sep: str = '&'
    print(f'time [minutes]  {sep} \\multicolumn{{2}}{{c|}}{{contact dose}} {sep} \\multicolumn{{2}}{{c}}{{handling dose}} \\\\')
    print(f'{minutes:.1f}    {sep} value {sep} error {sep} value {sep} error \\\\ \\hline')
    print(f'neutron: ', end='')
    for particle in ['1', '5']:
        print(f' {sep} {d0[particle]['value']:.3f} {sep} {d0[particle]['stdev']:.3f} ', end='')
        # print(f' {sep} {d0[particle]['value'] * 1e3:.5f} {sep} {d0[particle]['stdev'] * 1e3:.5f} ', end='')
    print('\\\\')
    print(f'gamma:   ', end='')
    for particle in ['2', '6']:
        print(f' {sep} {d0[particle]['value']:.3f} {sep} {d0[particle]['stdev']:.3f} ', end='')
    print('\\\\ \\hline \\hline')
