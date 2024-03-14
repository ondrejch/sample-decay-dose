#!/bin/env python3
"""
Example use case of SampleDose.DoseEstimatorTank - simple decay doses of F71 sample,
using parallel execution of the MAVRIC cases.
Note that the beta dose is zero if the sample has additional shielding.
Ondrej Chvala <ochvala@utexas.edu>
"""
import numpy as np
import json5
from sample_decay_dose import SampleDose
from joblib import Parallel, delayed, cpu_count

n_jobs: int = cpu_count()  # How many MAVRIC cases to run in parallel
sample_mass: float = 1200.0e3  # 1200 kg
max_decay_time_years: float = 4.0
decay_steps: int = 144


def my_process(decay_d: float) -> dict:
    # Load and decay the nuclide vector from F71 file
    # Set F71 file path and sample mass [g]
    origen_triton = SampleDose.OrigenFromTriton('../SCALE_FILE.f71', sample_mass)
    # Select F71 file position [seconds]
    print(f'ORIGEN set {decay_d} days')
    origen_triton.set_f71_pos(2.0 * 365.0 * 24.0 * 60.0 * 60.0)  # 2 year burn
    # Read burned nuclides
    origen_triton.read_burned_material()
    # Set hoe long they decay [days]
    origen_triton.set_decay_days(decay_d)
    # Execute ORIGEN
    origen_triton.run_decay_sample()

    # Calculate dose next to the tank
    mavric = SampleDose.DoseEstimatorStorageTank(origen_triton)
    # Material composition of additional layers, in dictionaries of atom densities
    mavric.layers_mats = [SampleDose.ADENS_SS316H_HOT, SampleDose.ADENS_KAOWOOL_COLD]
    # Thicknesses of additional layers [cm]
    mavric.layers_thicknesses = [2.54, 2.0 * 2.54]
    # Temperatures of additional layers [cm]
    mavric.layers_temperature_K = [873.0, 300.0]
    # Add more planes since the source is large
    mavric.N_planes_cyl = 12
    # Monaco histories
    mavric.histories_per_batch = 200000
    mavric.batches = 40
    # Run simulation
    mavric.run_mavric()
    mavric.get_responses()
    # Print doses
    print(f'Neutron dose {mavric.responses["1"]["value"]} +- {mavric.responses["1"]["stdev"]}  rem/h')
    print(f'Gamma dose   {mavric.responses["2"]["value"]} +- {mavric.responses["2"]["stdev"]}  rem/h')
    _res: dict = {decay_d: mavric.responses["2"]}  # Gamma dose
    return _res


def run_analysis():
    # Parallel jobs
    decay_timeline = np.linspace(1.0, max_decay_time_years * 365.24, decay_steps)
    results = Parallel(n_jobs=n_jobs)(delayed(my_process)(decay_day) for decay_day in decay_timeline)
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

    # Integrate
    dose2y: float = 0
    for i in range(len(x) - 1):
        trapezoid = (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2.0  # Trapezoidal rule
        dose2y += trapezoid * 24.0  # [rem/h] -> [rem], integrating over days
    print(f'Gamma dose over {max_decay_time_years:.1f} years of decay: {dose2y/1e6:.1f} Mrem')
    print(f'Gamma dose at {max_decay_time_years:.1f} years of decay: {y[-1]:.1f} ± {yerr[-1]:.1f} rem/h, '
          f'or {y[-1]*24.0*365.24/1e6:.1f} Mrem/year')

    # Plots!
    if do_plots:
        plt.close('all')
        plt.xscale('linear')
        plt.yscale('linear')
        plt.grid()
        plt.title("Gamma dose from unshielded container, 2 EFPY at 1 MWt")
        plt.xlabel(f'decay time [days]')
        plt.ylabel('Dose at 1 cm [rem/h]')
        plt.errorbar(x, y, yerr, ls='none', color='slategrey', capsize=1.2)
        plt.scatter(x, y, color='slategrey', s=5, label='Gamma')

        plt.legend()
        plt.tight_layout()
        label_file_name = 'fs_days'
        plt.savefig(f'dose_{label_file_name}.png', dpi=1000)
        # plt.show()

        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(f'dose_{label_file_name}-loglog.png', dpi=1000)
        # plt.show()


if __name__ == "__main__":
    # run_analysis()
    # plot()
    plot('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/40-decaydays_1200kg/doses.json')
