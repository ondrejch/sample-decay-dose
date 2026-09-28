#!/bin/env python3
"""
Example use case of SampleDose -- gamma dose from an offgas tank as a function of steel and concrete
shield thicknesses. The nuclide inventory is read from a CSV file of <nuclide>, <number of atoms> rows.

Usage:
    python decay_atoms_cvs.py --run-analysis --decay-days 2 --atoms-csv /tmp/atoms.csv
        runs the ORIGEN decay and the parallel MAVRIC scan, and writes doses_decay_atoms_cvs.json
    python decay_atoms_cvs.py --run-2y-decay-only --atoms-csv /tmp/atoms.csv
        runs only a 2-year ORIGEN decay of the inventory
    python decay_atoms_cvs.py
        plots the scans stored in the per-decay-time result folders listed in plot()
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose, utils
import argparse
import numpy as np
import pandas as pd
import json5
from joblib import Parallel, delayed, cpu_count

n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel
decay_days: float = 2.0  # 2 days

cm2_to_barn: float = 1e24  # 1 cm^2 = 1e24 barn
my_inner_r: float = 37.7825  # IR of 30" schedule 40 pipe
my_thick: float = 0.635  # Thickness of 30" schedule 40 pipe
my_volume: float = 600e3  # 600 liters
my_atoms_file: str = '/tmp/atoms.csv'  # Default nuclide CSV file, override with --atoms-csv
my_doses_file: str = 'doses_decay_atoms_cvs.json'  # Written by run_analysis(), read by make_plot()


def print_atoms(my_atom_density: dict):
    """ Debugging """
    tot_atoms: float = 0.0
    tot_atom_density: float = 0.0
    for k, v in my_atom_density.items():
        tot_atom_density += v
        tot_atoms += v * cm2_to_barn * my_volume
    print(f'Total atoms: {tot_atoms}, total atom density {tot_atom_density} atoms / barn-cm')


def load_origen_decay(atoms_file: str = my_atoms_file) -> SampleDose.OrigenDecayBox:
    """ Reads the nuclide CSV file and returns the ORIGEN decay box of the sample """
    my_atom_density: dict = utils.read_cvs_atom_dens(atoms_file, my_volume)
    print_atoms(my_atom_density)
    return SampleDose.OrigenDecayBox(my_atom_density, my_volume)


def mavric_process(_origen_decay: SampleDose.OrigenDecayBox, case: tuple[float, float]) -> dict:
    """ Separating the MAVRIC part into a function for parallel execution """
    steel_cm: float
    concrete_cm: float
    (steel_cm, concrete_cm) = case
    # Calculate dose next to the tank
    mavric = SampleDose.DoseEstimatorGenericTank(_origen_decay)
    mavric.cyl_r = my_inner_r  # mavric_deck() sets the sample half-height from sample_volume and cyl_r
    # Material composition of additional layers, in dictionaries of atom densities
    mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD, SampleDose.ADENS_SS316H_COLD,
                          SampleDose.ADENS_CONCRETE_COLD]
    # Thicknesses of additional layers [cm]
    mavric.layers_thicknesses = [my_thick, 2.0 * 2.54, steel_cm, concrete_cm]
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


def run_2y_decay_only(atoms_file: str = my_atoms_file):
    origen_decay = load_origen_decay(atoms_file)
    origen_decay.set_decay_days(2.0 * 365.24)
    origen_decay.SAMPLE_F71_position = 900  # sample decay steps
    origen_decay.write_atom_dens()
    origen_decay.run_decay_sample()


def run_analysis(atoms_file: str = my_atoms_file, my_decay_days: float = decay_days):
    origen_decay = load_origen_decay(atoms_file)
    origen_decay.set_decay_days(my_decay_days)
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
    results = Parallel(n_jobs=n_jobs)(delayed(mavric_process)(origen_decay, case) for case in case_inputs)
    print(results)

    with open(my_doses_file, 'w') as file_out:
        json5.dump(results, file_out, indent=4)


def make_plot(title: str, my_dir: str):
    import matplotlib.pyplot as plt
    from matplotlib import colors, cm, ticker

    with open(f'{my_dir}/{my_doses_file}') as fin:
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
    # Rows of g_mrem_* follow _steel_cm_list (index i), columns follow _concrete_cm_list (index j)
    pd_mrem_dose = pd.DataFrame(g_mrem_dose, index=pd.Index(steel_cm, name='steel [cm]'),
                                columns=pd.Index(concrete_cm, name='concrete [cm]'))
    pd_mrem_stdev = pd.DataFrame(g_mrem_stdev, index=pd.Index(steel_cm, name='steel [cm]'),
                                 columns=pd.Index(concrete_cm, name='concrete [cm]'))
    return pd_mrem_dose, pd_mrem_stdev


def plot(base_dir='.'):
    from pathlib import Path
    base_path = Path(base_dir).expanduser().resolve()
    my_data = {
        '1 day': '41-1day_decay',
        '2 days': '42-2days_decay',
        '7 days': '43-7days_decay',
        '14 days': '44-14days_decay',
        '30 days': '45-30days_decay'
    }
    writer = pd.ExcelWriter(base_path / 'offgas_dose.xlsx')
    for t, d in my_data.items():
        pd_dose, pd_stdev = make_plot(t, str(base_path / d))
        # header = f'Rows = Steel [cm], Columns = Concrete [cm]'
        # pd_dose.columns = pd.MultiIndex.from_product([[header],  pd_dose.columns])
        pd_dose.style.map(lambda v: 'color:#8B0000' if v > 20 else None). \
            map(lambda v: 'font-weight:bold;color:#008000' if 20 > v > 2 else None). \
            to_excel(writer, sheet_name=f'dose (mrem per h), {t}', float_format="%0.1f")
        pd_stdev.to_excel(writer, sheet_name=f'dose ± stdev, {t}', float_format="%0.2f")
    writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Offgas tank dose scan over steel and concrete shield thicknesses, and plotting. '
                    'Without --run-analysis or --run-2y-decay-only, only the plots are made.'
    )
    parser.add_argument('--run-analysis', action='store_true',
                        help=f'Run the ORIGEN decay and the MAVRIC scan in the current directory, write {my_doses_file}.')
    parser.add_argument('--run-2y-decay-only', action='store_true',
                        help='Run only a 2-year ORIGEN decay of the inventory.')
    parser.add_argument('--atoms-csv', default=my_atoms_file,
                        help='Nuclide CSV file (<nuclide>, <number of atoms> per row) for the ORIGEN runs.')
    parser.add_argument('--decay-days', default=decay_days, type=float,
                        help='Decay time [days] for --run-analysis.')
    parser.add_argument('--base-dir', default='.',
                        help='Directory that holds the per-decay-time result folders listed in plot().')
    args = parser.parse_args()

    if args.run_2y_decay_only:
        run_2y_decay_only(args.atoms_csv)
    if args.run_analysis:
        run_analysis(args.atoms_csv, args.decay_days)
    if not (args.run_analysis or args.run_2y_decay_only):
        plot(args.base_dir)
