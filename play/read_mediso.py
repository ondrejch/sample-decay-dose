#!/bin/env python3
"""
Extract Bq of Ac225 from F71 file
Ondrej Chvala <ochvala@utexas.edu>
"""
import os
import re
import subprocess
import numpy as np
from datetime import datetime

NOW: str = datetime.now().replace(microsecond=0).isoformat()

SCALE_bin_path: str = os.getenv('SCALE_BIN', '/opt/scale6.3.1/bin/')
ATOM_DENS_MINIMUM: float = 1e-60
MAVRIC_NG_XSLIB: str = 'v7.1-28n19g'


def get_f71_positions_index(f71file: str) -> dict:
    """ Read info of SCALE's F71 file """
    output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=info", f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    f71_idx = {}
    skip_data = ['pos', '(-)']  # sip records starting with this
    end_data = 'state definition present'  # stop reading after reaching this
    for line in output:
        data = line.split()
        if data[0].strip() in skip_data:
            continue
        if line.count(end_data) > 0:
            break
        # print(data)
        f71_idx[int(data[0])] = {'time': data[1], 'power': data[2], 'flux': data[3], 'fluence': data[4],
            'energy': data[5], 'initialhm': data[6], 'libpos': data[7], 'case': data[8], 'step': data[9],
            'DCGNAB': data[10]}
    return f71_idx


def get_last_position_for_case(f71file: str, case: int = 1) -> int:
    iii = get_f71_positions_index(f71file)
    maxi = max([i for i, rec in iii.items() if rec['case'] == str(case)])
    return maxi


def get_burned_material_atom_dens(f71file: str, position: int) -> dict:
    """ Read atom density of nuclides from SCALE's F71 file """
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'",
            f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    densities = {}  # densities[nuclide] = (density at position of f71 file)
    skip = ["case", "step", "time", "power", "flux", "volume"]
    regexp = re.compile(r"(?P<elem>[a-zA-Z]+)(?P<num>\d+)(?P<meta>m)?")
    for line in output:
        data = line.split(',')
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            dummy = re.search(regexp, data[0].strip())
            elem = dummy.group("elem").lower()  # convert to all lower cases
            num = int(dummy.group("num"))  # to cut off leading zeros
            if dummy.group("meta"):
                nuclide = elem + "-" + str(num) + "m"  # for metastable isotopes
            else:
                nuclide = elem + "-" + str(num)
            if float(data[position]) > ATOM_DENS_MINIMUM:
                # The [x] here is what causes the code to return only the densities at position x of the f71 file
                densities[nuclide] = float(data[position])
    sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
    return sorted_densities


def get_nuclide_Bq(f71file: str, position: int, my_nuclide: str) -> float:
    my_nuclide = my_nuclide.lower()
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=becq", "-idform='{:Ee}{:AAA}{:m}'",
            f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    # bqs = {}  # densities[nuclide] = (density at position of f71 file)
    skip = ["case", "step", "time", "power", "flux", "volume"]
    regexp = re.compile(r"(?P<elem>[a-zA-Z]+)(?P<num>\d+)(?P<meta>m)?")
    for line in output:
        data = line.split(',')
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            dummy = re.search(regexp, data[0].strip())
            elem = dummy.group("elem").lower()  # convert to all lower cases
            num = int(dummy.group("num"))  # to cut off leading zeros
            if dummy.group("meta"):
                nuclide = elem + "-" + str(num) + "m"  # for metastable isotopes
            else:
                nuclide = elem + "-" + str(num)
            if float(data[position]) > ATOM_DENS_MINIMUM:
                # The [x] here is what causes the code to return only the densities at position x of the f71 file
                # bqs[nuclide] = float(data[position])
                if nuclide == my_nuclide:
                    # print(f"{my_nuclide}: {data[position]}")
                    return float(data[position])
    # sorted_densities = {k: v for k, v in sorted(bqs.items(), key=lambda item: -item[1])}
    return -1.0


if __name__ == "test":
    f71_file_name: str = 'ThEIRENE.f71'
    ipos: int = get_last_position_for_case(f71_file_name, 1)
    ac_225_Bq: float = get_nuclide_Bq(f71_file_name, ipos, 'Ac-225')
    print(f'Ac-225: {ac_225_Bq:.3e} Bq')  # print(f'{ac_225_Bq:.3e}')

if __name__ == "__main__":
    f71_file_name: str = 'ThEIRENE.f71'
    MTiHM: dict = {' 5.00': 10.4159200656709, '19.75': 10.2647529843048}
    # runs: dict = {
    #     '01-1_year_120_steps': {
    #         'power': 400.0, 'enr': '19.75'
    #     },
    #     '02-gasFP_removal-1_year_120_steps': {
    #         'power': 400.0, 'enr': '19.75'
    #     },
    #     '11-5pct-1_year_120_steps': {
    #         'power': 400.0, 'enr': ' 5.00'
    #     },
    #     '12-5ctp_gasFP_removal-1_year_120_steps': {
    #         'power': 400.0, 'enr': ' 5.00'
    #     },
    #     '20-20MWth/05.00pct_1y': {
    #         'power': 20.0, 'enr': ' 5.00'
    #     },
    #     '20-20MWth/05.00pct_1y_gas_removal': {
    #         'power': 20.0, 'enr': ' 5.00'
    #     },
    #     '20-20MWth/19.75pct_1y': {
    #         'power': 20.0, 'enr': '19.75'
    #     },
    #     '20-20MWth/19.75pct_1y_gas_removal': {
    #         'power': 20.0, 'enr': '19.75'
    #     },
    #     '21-01MWth/05.00pct_1y': {
    #         'power': 1.0, 'enr': ' 5.00'
    #     },
    #     '21-01MWth/05.00pct_1y_gas_removal': {
    #         'power': 1.0, 'enr': ' 5.00'
    #     },
    #     '21-01MWth/19.75pct_1y': {
    #         'power': 1.0, 'enr': '19.75'
    #     },
    #     '21-01MWth/19.75pct_1y_gas_removal':{
    #         'power': 1.0, 'enr': '19.75'
    #     },
    # }

    runs: dict = {
        '21-01MWth/05.00pct_1y': {'power':  1.0, 'enr': ' 5.00'},
        '21-01MWth/19.75pct_1y': {'power':  1.0, 'enr': '19.75'},
        '20-20MWth/05.00pct_1y': {'power': 20.0, 'enr': ' 5.00'},
        '20-20MWth/19.75pct_1y': {'power': 20.0, 'enr': '19.75'},
        '11-5pct-1_year_120_steps': {'power': 400.0, 'enr': ' 5.00'},
        '01-1_year_120_steps': {'power': 400.0, 'enr': '19.75'},
    }
    my_nuclide: str = 'Ac-225'
    print(f'{my_nuclide} production')
    for my_path, my_run in runs.items():
        my_f71_file_name = f'{my_path}/{f71_file_name}'
        thermal_power: float = my_run['power']
        enrichment: str = my_run['enr']
        ipos: int = get_last_position_for_case(my_f71_file_name, 1)
        # print(f'ipos = {ipos}')
        ac_225_Bq: float = get_nuclide_Bq(my_f71_file_name, ipos, my_nuclide)
        ac_225_Bq_pery: float = ac_225_Bq * MTiHM[enrichment]
        ac_225_Bq_pery_perMW: float = ac_225_Bq_pery / thermal_power
        doses: float = ac_225_Bq_pery_perMW / 15e3
        print(
            f'Power {thermal_power:5.1f} Uenr {enrichment:5s}   {ac_225_Bq_pery_perMW:.3e} Bq/y/MWth  doses {doses:3.0f} /y/MWth')

"""This is 400 MWth, and 10.265 MTiHM. 
We get about 2e8 Bq of Ac225 and Bi213 each from 1MTiHM in a year of burn, 
so about (2e8 * 10 / 400 = ) 5e6 Bq per 1MWth-year, if I calculate that correctly.
A medical dose is about 15 kBq, so about 330 doses per MWth per year. """
