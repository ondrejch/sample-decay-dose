#!/bin/env python3
"""
Example use case of SampleDose -- calculate decay doses of sample as a function of decay times
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import numpy as np
import pandas as pd
import json5
from joblib import Parallel, delayed, cpu_count

n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel
decay_days: float = 2.0  # 2 days

cm2_to_barn: float = 1e24  # 1 cm^2 = 1e24 barn
my_inner_r: float = 37.7825  # IR of 30" schedule 40 pipe
my_volume: float = 600e3  # 600 liters
my_atoms_file: str = '/tmp/atoms.csv'


def print_atoms():
    """ Debugging """
    tot_atoms: float = 0.0
    tot_atom_density: float = 0.0
    for k, v in my_atom_density.items():
        tot_atom_density += v
        tot_atoms += v * cm2_to_barn * my_volume
    print(f'Total atoms: {tot_atoms}, total atom density {tot_atom_density} atoms / barn-cm')


def mavric_process(case: tuple[float, float]) -> dict:
    """ Separating the MAVRIC part into a function for parallel execution """
    steel_cm: float
    concrete_cm: float
    (steel_cm, concrete_cm) = case
    # Calculate dose next to the tank
    mavric = SampleDose.DoseEstimatorGenericTank(origen_decay)
    mavric.cyl_r = my_inner_r
    mavric.sample_h2 = SampleDose.get_cyl_h(my_volume, my_inner_r)
    # Material composition of additional layers, in dictionaries of atom densities
    mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD, SampleDose.ADENS_SS316H_COLD,
                          SampleDose.ADENS_CONCRETE_COLD]
    # Thicknesses of additional layers [cm]
    mavric.layers_thicknesses = [2.54, 2.0 * 2.54, steel_cm, concrete_cm]
    # Temperatures of additional layers [cm]
    mavric.layers_temperature_K = [873.0, 300.0, 300.0, 300.0]
    # Add more planes since the source is large
    mavric.N_planes_cyl = 20
    # Monaco histories
    mavric.histories_per_batch = 400000
    mavric.batches = min(150, 5 + int(steel_cm * concrete_cm / 4))

    # Run simulation
    mavric.run_mavric()
    mavric.get_responses()
    # Print doses
    print(f'Neutron dose {mavric.responses["1"]["value"]} +- {mavric.responses["1"]["stdev"]}  rem/h')
    print(f'Gamma dose   {mavric.responses["2"]["value"]} +- {mavric.responses["2"]["stdev"]}  rem/h')
    _res: dict = {steel_cm: {}}
    _res[steel_cm][concrete_cm] = mavric.responses
    return _res


my_atom_density: dict = SampleDose.read_cvs_atom_dens(my_atoms_file, my_volume)
print_atoms()
origen_decay = SampleDose.OrigenDecayBox(my_atom_density, my_volume)  # This needs to be global scope


def run_2y_decay_only():
    origen_decay.set_decay_days(2.0 * 365.24)
    origen_decay.SAMPLE_F71_position = 300  # sample decay steps
    origen_decay.write_atom_dens()
    origen_decay.run_decay_sample()


def run_analysis():
    origen_decay.set_decay_days(decay_days)
    origen_decay.SAMPLE_F71_position = 30  # sample decay steps
    origen_decay.write_atom_dens()
    origen_decay.run_decay_sample()

    # Inputs for joblib parallelism have to be iterable
    case_inputs: list[tuple[float, float]] = []
    d = {}  # This would be better handled with Pandas ..
    for steel_shield_thick_in in np.geomspace(5, 18, 30):
        s_cm = 2.54 * steel_shield_thick_in
        d[s_cm] = {}
        for concrete_shield_in in np.geomspace(5, 18, 30):
            c_cm = 2.54 * concrete_shield_in
            case_inputs.append((s_cm, c_cm))

    # Parallel MAVRIC jobs
    results = Parallel(n_jobs=n_jobs)(delayed(mavric_process)(case) for case in case_inputs)
    print(results)

    with open('doses.json', 'w') as file_out:
        json5.dump(results, file_out, indent=4)


def make_plot(title: str, my_dir: str):
    import matplotlib.pyplot as plt
    from matplotlib import colors, cm, ticker

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
    print(len(_steel_cm_list), _steel_cm_list)
    print(len(_concrete_cm_list), _concrete_cm_list)

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
    # cs = ax.contourf(x, y, g_mrem_dose.T,  norm=colors.LogNorm(),  cmap='jet')
    # cs = ax.contour(x, y, g_mrem_dose.T, levs)
    ax.clabel(cs, inline=True, fontsize=10, manual=False, colors=['black'], fmt="%.0e")
    # cbar = fig.colorbar(cs, format = "%.1e")
    cbar = fig.colorbar(cs, format="%.05g")
    cbar.set_label('Gamma dose [mrem/h]')
    plt.title(f'Gamma dose from offgas tank after {title} of decay')
    plt.xlabel('Steel shield thickness [cm]')
    plt.ylabel('Concrete shield thickness [cm]')
    file_name_fig = f'dose_g_offgas_tank-{title.replace(" ", "_")}_decay.png'
    plt.savefig(file_name_fig, dpi=1000, bbox_inches='tight', pad_inches=0.1)
    plt.tight_layout()
    # plt.show()

    steel_cm = [f'{float(x):.2f}' for x in _steel_cm_list]
    concrete_cm = [f'{float(x):.2f}' for x in _concrete_cm_list]
    pd_mrem_dose = pd.DataFrame(g_mrem_dose, columns=steel_cm, index=concrete_cm)
    pd_mrem_stdev = pd.DataFrame(g_mrem_stdev, columns=steel_cm, index=concrete_cm)
    return pd_mrem_dose, pd_mrem_stdev


def plot(datafile='doses.json'):
    import os
    os.chdir('/home/o/MSRR-local/92-offgas_dose')
    my_data = {
        '1 day': '01-1day_decay',
        '2 days': '02-2days_decay',
        '7 days': '03-7days_decay',
        '14 days': '04-14days_decay',
        '30 days': '05-30days_decay'
    }
    writer = pd.ExcelWriter('offgas_dose.xlsx')
    for t, d in my_data.items():
        pd_dose, pd_stdev = make_plot(t, d)
        # header = f'Rows = Steel [cm], Columns = Concrete [cm]'
        # pd_dose.columns = pd.MultiIndex.from_product([[header],  pd_dose.columns])
        pd_dose.style.map(lambda v: 'color:#8B0000' if v > 20 else None). \
            map(lambda v: 'font-weight:bold;color:#008000' if 20 > v > 2 else None). \
            to_excel(writer, sheet_name=f'dose (mrem per h), {t}', float_format="%0.1f")
        pd_stdev.to_excel(writer, sheet_name=f'dose ± stdev, {t}', float_format="%0.2f")
    writer.close()


if __name__ == "__main__":
    # run_2y_decay_only()
    # run_analysis()
    plot()
