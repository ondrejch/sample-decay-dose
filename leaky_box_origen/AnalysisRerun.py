#!/bin/env python3
"""
Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

Load analysis results from JSONS and re-generates the spreadsheet.
"""
import numpy as np
import pandas as pd
import json5
from leaky_box_origen.LeakyBox import get_dataframe, PCTperDAY

volume: float = 498663.8595795657
box_A_leak_rate: float = PCTperDAY
box_B_leak_rate: float = PCTperDAY * 0.1

is_xe_136_testing: bool = False  # Use 1 Xe-136 at / b-sm instead of read composition
is_xe_135_testing: bool = True  # Use 1 Xe-135 at / b-sm instead of read composition
if is_xe_135_testing and is_xe_136_testing:
    raise ValueError("Tests are exclusive!")
dtB: float = 0.0  # time offset box B
dtC: float = 43200./2.0  # time offset Bbx C


def main():
    with open('boxA.json5', 'r') as f:  # Save to json
        box_A_adens = json5.load(f)
    with open('boxB.json5', 'r') as f:  # Save to json
        box_B_adens = json5.load(f)
    with open('boxC.json5', 'r') as f:
        box_C_adens = json5.load(f)

    pd_A: pd.DataFrame = get_dataframe(box_A_adens)
    pd_B: pd.DataFrame = get_dataframe(box_B_adens)
    pd_C: pd.DataFrame = get_dataframe(box_C_adens)

    print(pd_B.columns)
    print(pd_C.columns)
    if is_xe_136_testing:  # Add analytical calculations
        pd_B['total'] = pd_B['xe-136'] * volume
        pd_C['total'] = pd_C['xe-136'] * volume
        """
N_i^A(t) & = N_i^A e^{ - \epsilon_i^A t} 
N_i^B(t) & = N_i^A \epsilon_i^A \frac{e^{-\epsilon_i^B t} - e^{-\epsilon_i^A t}} {\epsilon_i^A - \epsilon_i^B}
N_i^C(t) & = N_i^A \frac{ -\epsilon_i^B e^{-\epsilon_i^A t} + \epsilon_i^A (e^{-\epsilon_i^B t} - 1) + \epsilon_i^B}
{\epsilon_i^B - \epsilon_i^A}
        """
        pd_A['analytic'] = (1.0 * np.exp(-box_A_leak_rate * pd_A['time [s]']))
        pd_B['analytic'] = (1.0 *
                            box_A_leak_rate *
                            (np.exp(-box_A_leak_rate * (pd_B['time [s]'] + dtB)) -
                             np.exp(-box_B_leak_rate * (pd_B['time [s]'] + dtB)))
                            / (box_B_leak_rate - box_A_leak_rate)
                            )
        pd_C['analytic'] = (1.0 *
                            (-box_A_leak_rate * np.exp(-box_B_leak_rate * (pd_C['time [s]']+dtC)) +
                             box_B_leak_rate * (np.exp(-box_A_leak_rate * (pd_C['time [s]']+dtC)) - 1.0) +
                             box_A_leak_rate)
                            / (box_A_leak_rate - box_B_leak_rate)
                            )
        pd_A['diff [%]'] = 100.0 * (pd_A['xe-136'] - pd_A['analytic']) / pd_A['analytic']
        pd_B['diff [%]'] = 100.0 * (pd_B['total'] - pd_B['analytic']) / pd_B['analytic']
        pd_C['diff [%]'] = 100.0 * (pd_C['total'] - pd_C['analytic']) / pd_C['analytic']

    if is_xe_135_testing:  # Add analytical calculations
        lambda_xe_135: float = 2.106574217602 * 1e-5  # Xe-135 decay constant [1/s]
        pd_B['total'] = pd_B['xe-135'] * volume
        pd_C['total'] = pd_C['xe-135'] * volume
        """
N_i^A(t) & = N_i^A e^{ - (\lambda_i + \epsilon_i^A) t}  \\
\\
N_i^B(t) & = N_i^A \epsilon_i^A \frac{e^{-(\epsilon_i^A + \lambda_i)  t} - e^{-(\epsilon_i^B + \lambda_i) t}}
{\epsilon_i^B - \epsilon_i^A} \\
% B(t) = -(a N (e^(t (-(a + l))) - e^(t (-(b + l)))))/(a - b)
\\
N_i^C(t) & = N_i^A e^{ - \lambda_i  t} 
\frac{ -\epsilon_i^B e^{-\epsilon_i^A t} + \epsilon_i^A (e^{-\epsilon_i^B t} - 1) + \epsilon_i^B}
{\epsilon_i^B - \epsilon_i^A}
        """
        pd_A['analytic'] = (1.0 * np.exp(-(box_A_leak_rate+lambda_xe_135) * pd_A['time [s]']))
        pd_B['analytic'] = (1.0 *
                            box_A_leak_rate *
                            (np.exp(-(box_A_leak_rate+lambda_xe_135) * (pd_B['time [s]'] + dtB)) - np.exp(
                                -(box_B_leak_rate+lambda_xe_135) * (pd_B['time [s]'] + dtB)))
                            / (box_B_leak_rate - box_A_leak_rate)
                            )
        pd_C['analytic'] = (1.0 * np.exp(-lambda_xe_135 * (pd_C['time [s]'] + dtC)) *
                            (-box_B_leak_rate * np.exp(-box_A_leak_rate * (pd_C['time [s]'] + dtC)) +
                             box_A_leak_rate * (
                                         np.exp(-box_B_leak_rate * (pd_C['time [s]'] + dtC)) - 1.0) + box_B_leak_rate)
                            / (box_B_leak_rate - box_A_leak_rate)
                            )
        pd_A['diff [%]'] = 100.0 * (pd_A['xe-135'] - pd_A['analytic']) / pd_A['analytic']
        pd_B['diff [%]'] = 100.0 * (pd_B['total'] - pd_B['analytic']) / pd_B['analytic']
        pd_C['diff [%]'] = 100.0 * (pd_C['total'] - pd_C['analytic']) / pd_C['analytic']

    writer = pd.ExcelWriter('leaky_boxes1.xlsx')

    pd_A.to_excel(writer, sheet_name='box A')
    pd_B.to_excel(writer, sheet_name='box B')
    pd_C.to_excel(writer, sheet_name='box C')
    writer.close()


if __name__ == "__main__":
    # pass
    main()
