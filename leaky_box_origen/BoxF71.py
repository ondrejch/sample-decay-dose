#!/bin/env python3
"""
Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

100% per day  =  1/24/60/60 per second    = 0.0000115740740741
   1% per day  =  1e-2/24/60/60 per second = 0.000000115740740741
 0.1% per day  =  1e-3/24/60/60 per second = 0.0000000115740740741
"""
import copy
import os
import re
import numpy as np
from bisect import bisect_left

from sample_decay_dose.SampleDose import NOW, get_rho_from_atom_density, nicely_print_atom_dens, atom_dens_for_origen, \
    get_f71_positions_index, get_burned_material_total_mass_dens, get_burned_material_atom_dens, run_scale

PCTperDAY: float = 0.000000115740740741


class Origen:
    """ ORIGEN handling parent class """

    def __init__(self):
        self.debug: int = 3  # Debugging flag
        self.cwd: str = os.getcwd()  # Current running fir
        self.ORIGEN_input_file_name: str = 'origen.inp'
        self.decayed_atom_dens: dict = {}  # Atom density of the decayed sample
        self.case_dir: str = ''
        self.ATOM_DENS_file_name_Origen: str = 'my_sample_atom_dens_origen.inp'
        self.F71_file_name: str = self.ORIGEN_input_file_name.replace('inp', 'f71')
        self.DECAY_days: float = 30.0  # Sample decay time [days]
        self.DECAY_steps: int = 30  # Number of steps for ORIGEN decay
        self.salt_weight: float = np.NaN  # Mass of the sample [g]
        self.salt_density: float = np.NaN  # Mass density of the sample [g/cm3]
        self.salt_volume: float = np.NaN  # Sample volume [cm3]
        self.salt_atom_dens: dict = {}

    def set_decay_days(self, decay_days: float = 30.0):
        """ Use this to change decay time, as it also updates the case directory """
        self.DECAY_days = decay_days
        self.case_dir: str = f'run_{NOW}_{self.salt_weight:.5}_g-{decay_days:.5}_days'  # Directory to run the case


class DecaySalt(Origen):
    """ Decays salt in ORIGEN """

    def __init__(self, _f71: str = './SCALE_FILE.f71', _mass: float = 500e3):
        Origen.__init__(self)
        self.ORIGEN_input_file_name: str = 'salt-decay.inp'
        self.BURNED_MATERIAL_F71_file_name: str = _f71  # Burned core F71 file from TRITON
        self.BURNED_MATERIAL_F71_index: dict = get_f71_positions_index(self.BURNED_MATERIAL_F71_file_name)
        self.BURNED_MATERIAL_F71_position: int = 16
        self.salt_weight: float = _mass  # Mass of the sample [g]
        self.nuclide_removal_rates: (dict, None) = None
        self.case_dir = 'fuelsalt'
        self.ORIGEN_input_file_name = 'core-decay.inp'

    def set_f71_pos(self, t: float = 5184000.0, case: str = '1'):
        """ Returns closest position in the F71 file for a case """
        pos_times = [(k, float(v['time'])) for k, v in self.BURNED_MATERIAL_F71_index.items() if v['case'] == case]
        times = [x[1] for x in pos_times]
        t_min = min(times)
        t_max = max(times)
        if t < t_min:
            print(f"Error: Time {t} seconds is less than {t_min} s, the minimum time in records.")
        #    return None
        if t > t_max:
            print(f"Error: Time {t} seconds is longer than {t_max} s, the maximum time in records.")
        #    return None
        self.BURNED_MATERIAL_F71_position = bisect_left(times, t)

    def read_burned_material(self):
        """ Reads atom density and rho from F71 file """
        if self.debug > 0:
            print(f'ORIGEN: reading nuclides from {self.BURNED_MATERIAL_F71_file_name}, '
                  f'position {self.BURNED_MATERIAL_F71_position}')
        self.salt_density = get_burned_material_total_mass_dens(self.BURNED_MATERIAL_F71_file_name,
                                                                self.BURNED_MATERIAL_F71_position)
        self.salt_volume = self.salt_weight / self.salt_density
        self.salt_atom_dens = get_burned_material_atom_dens(self.BURNED_MATERIAL_F71_file_name,
                                                            self.BURNED_MATERIAL_F71_position)
        if self.debug > 2:
            # print(list(self.salt_atom_dens.items())[:25])
            print(f'Salt density {self.salt_density} g/cm3, volume {self.salt_volume} cm3')
            nicely_print_atom_dens(self.salt_atom_dens)

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        with open(self.ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
            f.write(atom_dens_for_origen(self.salt_atom_dens))

        with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.origen_deck())

        if self.debug > 0:
            print(f'ORIGEN: decaying sample for {self.DECAY_days} days')
            print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

        self.decayed_atom_dens = get_burned_material_atom_dens(self.F71_file_name,
                                                               self.DECAY_steps)
        os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def origen_deck(self) -> str:
        """ Sample decay Origen deck """
        removal: str = ''
        if self.nuclide_removal_rates:
            removal += 'processing {\n'
            for k, v in self.nuclide_removal_rates.items():
                removal += f'       removal {{ rate={v} ele=[{k}] }}\n'
            removal += '    }'

        time_interp_steps = self.DECAY_steps - 3
        if time_interp_steps < 1:
            raise ValueError("Too few time steps")

        origen_output = f'''=shell
cp -r ${{INPDIR}}/{self.ATOM_DENS_file_name_Origen} .
end

=origen
' {NOW} 
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
case {{
    gamma=yes
    neutron=yes
    beta=yes
    lib {{ % decay only library
        file="end7dec"
    }}
    mat {{
        iso [
<{self.ATOM_DENS_file_name_Origen}
        ]
        units=ATOMS-PER-BARN-CM
        volume={self.salt_volume}
    }}
    time {{
        units=DAYS
        t=[{time_interp_steps}I 0.0001 {self.DECAY_days}]
        start=0
    }}
{removal} 
    save {{
        file="{self.F71_file_name}"
    }}
}}
end

=opus
data='{self.F71_file_name}'
title='Neutrons'
typarams=nspectrum
units=intensity
time=days
'tmin={self.DECAY_days}
'tmax={self.DECAY_days}
end

=opus
data='{self.F71_file_name}'
title='Gamma'
typarams=gspectrum
units=intensity
time=days
'tmin={self.DECAY_days}
'tmax={self.DECAY_days}
end

=opus
data='{self.F71_file_name}'
title='Beta'
typarams=bspectrum
units=intensity
time=days
'tmin={self.DECAY_days}
'tmax={self.DECAY_days}
end
'''
        return origen_output


class DiffLeak:
    """ Calculate difference between ORIGEN case with and without removals """

    def __init__(self):
        self.decay_base: (None, DecaySalt) = None
        self.decay_leaks: (None, DecaySalt) = None
        self.leaked: dict = {}

        self.nuclide_removal_rates: dict = {
            'H': PCTperDAY,
            'He': PCTperDAY,
            'O': PCTperDAY,
            'N': PCTperDAY,
            'I': PCTperDAY,
            'Ar': PCTperDAY,
            'Kr': PCTperDAY,
            'Xe': PCTperDAY
        }
        self.nuclide_removal_rates: dict = {'Kr': PCTperDAY,
                                            'Xe': PCTperDAY}

    def setup_cases(self, base_case: DecaySalt):
        """ Setup cases based on a base_case """
        self.decay_base = base_case
        self.decay_leaks = copy.deepcopy(base_case)
        self.decay_leaks.case_dir = self.decay_base.case_dir + '_leak'
        self.decay_leaks.nuclide_removal_rates = self.nuclide_removal_rates

    def run_cases(self):
        """ Run cases and get results """
        self.decay_base.run_decay_sample()
        self.decay_leaks.run_decay_sample()

    def get_diff_rate(self):
        """ Calculate difference between sealed and leaky box, and get equivalent ingress rate """
        f71_base_file: str = self.decay_base.case_dir + '/' + self.decay_base.F71_file_name
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        times: dict = get_f71_positions_index(f71_base_file)
        for i, vals in times.items():  # Time-dependent difference between sealed and leaky box
            t = float(vals['time'])
            self.leaked[i]: dict = {}
            self.leaked[i]['time'] = t
            adens_base: dict = get_burned_material_atom_dens(f71_base_file, i)
            adens_leak: dict = get_burned_material_atom_dens(f71_leak_file, i)
            # self.adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base
            #                     if (adens_base[key] - adens_leak.get(key, 0)) > 0}
            adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base}
            self.leaked[i]['adens'] = adens_diff
            # rate1 is calculated from leakage rate and the original concentration
            rate1: dict = {key: adens_base[key] * PCTperDAY for key in adens_base
                           if adens_base[key] > 0 and  # if nonzero
                           next((True for s in [k.lower() for k in self.nuclide_removal_rates.keys()] if s in key),
                                False)  # only for elements that leak
                           }
            self.leaked[i]['feed_rate1'] = rate1
            print(i, t, rate1)

        prev_t: float = -1.0
        prev_adens: dict = {}
        for i, vals in self.leaked.items():  # Calculate leakage, turn it into feed rate
            if i == 1:
                prev_t = vals['time']
                prev_adens = self.leaked[i]['adens']
                continue
            dt: float = vals['time'] - prev_t
            adens: dict = self.leaked[i]['adens']
            # rate2 is calculated using time derivative of concentration
            rate2: dict = {key: (adens[key] - prev_adens[key]) / dt for key in adens  # leak rate
                           if (adens[key] - prev_adens[key]) > 0 and  # if nonzero
                           next((True for s in [k.lower() for k in self.nuclide_removal_rates.keys()] if s in key),
                                False)  # only for elements that leak
                           }
            self.leaked[i]['feed_rate2'] = rate2
            print(i, vals['time'], dt, rate2)

            prev_t = vals['time']
            prev_adens = adens


def main():
    """ Calculate which nuclides leak """
    origen_triton = DecaySalt('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/SCALE_FILE.f71', 500e3)
    origen_triton.set_f71_pos(5.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 5 years
    origen_triton.read_burned_material()
    origen_triton.DECAY_days = 0.1
    origen_triton.DECAY_steps = 10

    dl = DiffLeak()
    dl.setup_cases(origen_triton)
    dl.run_cases()
    dl.get_diff_rate()


if __name__ == "__main__":
    # pass
    main()
