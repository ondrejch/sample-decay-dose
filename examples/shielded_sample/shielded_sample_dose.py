#!/bin/env python3
"""
Example use case of SampleDose.DoseEstimatorSquareTank
using parallel execution of the MAVRIC cases.
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
from sample_decay_dose import SampleDose
from joblib import Parallel, delayed, cpu_count

n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel
sample_mass: float = 0.7  # [g]
decay_days: float = 15.0 / (24.0 * 60.0)  # 15 minutes

ss_thick_max: float = 30  # [cm]
ss_thick_steps: int = 64


def my_process(ss_thickness: float) -> dict:
    # Calculate dose next to the tank
    mavric = SampleDose.DoseEstimatorSquareTank(origen_triton)
    mavric.case_dir += f'{ss_thickness:.3f}'
    # Material composition of additional layers, in dictionaries of atom densities
    mavric.layers_mats = [SampleDose.ADENS_HELIUM_COLD, SampleDose.ADENS_SS316H_COLD]
    # Thicknesses of additional layers [cm]
    mavric.layers_thicknesses = [0.1, ss_thickness]
    # Temperatures of additional layers [cm]
    mavric.layers_temperature_K = [300.0, 300.0]
    # Add more planes since the source is large
    # mavric.N_planes_cyl = 12
    # Monaco histories
    mavric.histories_per_batch = 50000
    mavric.batches = 20
    # Run simulation
    mavric.run_mavric()
    mavric.get_responses()
    # Print doses
    print(f'Neutron dose {mavric.responses["1"]["value"]} +- {mavric.responses["1"]["stdev"]}  rem/h')
    print(f'Gamma dose   {mavric.responses["2"]["value"]} +- {mavric.responses["2"]["stdev"]}  rem/h')
    _res: dict = {ss_thickness: mavric.responses["2"]}  # Gamma dose
    return _res


def run_analysis():
    # Load and decay the nuclide vector from F71 file
    # Set F71 file path and sample mass [g]
    global origen_triton
    origen_triton = SampleDose.OrigenFromTriton('../SCALE_FILE.f71', sample_mass)
    # Select F71 file position [seconds]
    print(f'ORIGEN set {decay_days} days')
    origen_triton.set_f71_pos(5.0 * 365.0 * 24.0 * 60.0 * 60.0)  # 5 year burn
    # Read burned nuclides
    origen_triton.read_burned_material()
    # Set hoe long they decay [days]
    origen_triton.set_decay_days(decay_days)
    # Execute ORIGEN
    origen_triton.run_decay_sample()

    # Parallel jobs
    shielding_thicknesses = np.linspace(0.1, ss_thick_max, ss_thick_steps)
    print(shielding_thicknesses)
    results = Parallel(n_jobs=n_jobs)(delayed(my_process)(ss_t) for ss_t in shielding_thicknesses)
    print(results)
    with open('doses.json', 'w') as file_out:
        json5.dump(results, file_out, indent=4)


def plot(datafile='doses.json'):
    import matplotlib.pyplot as plt
    do_plots: bool = False
    with open(datafile) as f:
        gamma_doses = json5.load(f)
    _x: list = []
    _y: list = []
    _yerr: list = []
    for dd in gamma_doses:
        for t, v in dd.items():
            _x.append(float(t))
            _y.append(float(v['value']))
            _yerr.append(float(v['stdev']))

    x = np.array(_x)
    y = np.array(_y)
    yerr = np.array(_yerr)

    # Plots!
    if do_plots:
        plt.close('all')
        plt.xscale('linear')
        plt.yscale('linear')
        plt.grid()
        plt.title(f"Gamma dose from the salt container, {sample_mass:.2f} g, 5 EFPY at 1 MWt")
        plt.xlabel(f'SS-316 thickness [cm]')
        plt.ylabel('Dose at 1 cm [rem/h]')
        plt.errorbar(x, y, yerr, ls='none', color='slategrey', capsize=1.2)
        plt.scatter(x, y, color='slategrey', s=5, label='Gamma')

        plt.legend()
        plt.tight_layout()
        label_file_name = f'fs_{sample_mass:.2f}g_5y'
        plt.savefig(f'dose_{label_file_name}.png', dpi=1000)
        # plt.show()

        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(f'dose_{label_file_name}-loglog.png', dpi=1000)
        # plt.show()


if __name__ == "__main__":
    origen_triton = None
    run_analysis()
    # plot()
    # plot('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/40-decaydays_1200kg/doses.json')
