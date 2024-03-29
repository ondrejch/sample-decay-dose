#!/bin/env python3
"""
OBSOLETE broken development version! Use LeakyBox.py

Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

Box A --> Box B --> outside

100% per day  =  1/24/60/60 per second    = 0.0000115740740741
   1% per day  =  1e-2/24/60/60 per second = 0.000000115740740741
 0.1% per day  =  1e-3/24/60/60 per second = 0.0000000115740740741
"""
import copy
import os
import re
import numpy as np
import pandas as pd
from bisect import bisect_left
import json5
from sample_decay_dose.SampleDose import NOW, nicely_print_atom_dens, atom_dens_for_origen, \
    get_f71_positions_index, get_burned_material_total_mass_dens, get_burned_material_atom_dens, run_scale

PCTperDAY: float = 0.000000115740740741


class Origen:
    """ ORIGEN handling parent class """

    def __init__(self):
        self.debug: int = 3  # Debugging flag
        self.cwd: str = os.getcwd()  # Current running fir
        self.ORIGEN_input_file_name: str = 'origen.inp'
        self.atom_dens: dict = {}  # Atom density of the sample
        self.case_dir: str = ''
        self.ATOM_DENS_file_name_Origen: str = 'my_sample_atom_dens_origen.inp'
        self.F71_file_name: str = self.ORIGEN_input_file_name.replace('inp', 'f71')
        self.DECAY_days: float = 30.0  # Sample decay time [days]
        self.DECAY_steps: int = 30  # Number of steps for ORIGEN decay
        self.weight: float = np.NaN  # Mass of the sample [g]
        self.density: float = np.NaN  # Mass density of the sample [g/cm3]
        self.volume: float = np.NaN  # Sample volume [cm3]
        self.final_atom_dens: dict = {}
        self.skip_calculation: bool = False  # Skips actual SCALE calculation, used for re-runs.


class DecayBoxA(Origen):
    """ Decays salt from F71 file in ORIGEN, no feed only removal, meant for box A decay """

    def __init__(self, _f71: str = './SCALE_FILE.f71', _mass: float = 500e3):
        Origen.__init__(self)
        self.BURNED_MATERIAL_F71_file_name: str = _f71  # Burned core F71 file from TRITON
        self.BURNED_MATERIAL_F71_index: dict = get_f71_positions_index(self.BURNED_MATERIAL_F71_file_name)
        self.BURNED_MATERIAL_F71_position: int = 16
        self.weight: float = _mass  # Mass of the sample [g]
        self.nuclide_removal_rates: (dict, None) = None
        self.case_dir = '_box_A'
        self.ORIGEN_input_file_name = 'core-decay.inp'

    def set_f71_pos(self, t: float = 5184000.0, case: str = '1'):
        """ Returns closest position in the F71 file for a case """
        pos_times = [(k, float(v['time'])) for k, v in self.BURNED_MATERIAL_F71_index.items() if v['case'] == case]
        times = [x[1] for x in pos_times]
        t_min = min(times)
        t_max = max(times)
        if t < t_min:
            print(f"Error: Time {t} seconds is less than {t_min} s, the minimum time in records.")
            pos_number: int = 0
        elif t > t_max:
            print(f"Error: Time {t} seconds is longer than {t_max} s, the maximum time in records.")
            pos_number: int = len(times) - 1
        else:
            pos_number: int = bisect_left(times, t)
        print(f'--> Closest F71 position found at slot {pos_number}, {times[pos_number]} seconds')
        self.BURNED_MATERIAL_F71_position = pos_number

    def read_burned_material(self):
        """ Reads atom density and rho from F71 file """
        if self.debug > 0:
            print(f'ORIGEN: reading nuclides from {self.BURNED_MATERIAL_F71_file_name}, '
                  f'position {self.BURNED_MATERIAL_F71_position}')
        self.density = get_burned_material_total_mass_dens(self.BURNED_MATERIAL_F71_file_name,
                                                           self.BURNED_MATERIAL_F71_position)
        self.volume = self.weight / self.density
        self.atom_dens = get_burned_material_atom_dens(self.BURNED_MATERIAL_F71_file_name,
                                                       self.BURNED_MATERIAL_F71_position)
        if self.debug > 2:
            # print(list(self.salt_atom_dens.items())[:25])
            print(f'Salt density {self.density} g/cm3, volume {self.volume} cm3')
            nicely_print_atom_dens(self.atom_dens)

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        if not self.skip_calculation:
            with open(self.ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
                f.write(atom_dens_for_origen(self.atom_dens))

            with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
                f.write(self.origen_deck())

            if self.debug > 0:
                print(f'ORIGEN: decaying sample for {self.DECAY_days} days')
                print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
            run_scale(self.ORIGEN_input_file_name)
        else:
            if not os.path.isfile(self.F71_file_name):
                error_text: str = f'Skip SCALE flag set, output file {self.case_dir}/{self.F71_file_name} is not found!'
                raise ValueError(error_text)

        self.final_atom_dens = get_burned_material_atom_dens(self.F71_file_name,
                                                             self.DECAY_steps)
        os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.final_atom_dens)

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
        volume={self.volume}
    }}
    time {{
        units=DAYS
        t=[{time_interp_steps}I 0.00001 {self.DECAY_days}]
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


class DecayBoxB(Origen):
    """ Decays material from density list in ORIGEN using 1 time step, feed and removal, box B calculation """

    def __init__(self, _atom_dens: dict, _volume: float):
        Origen.__init__(self)
        self.atom_dens = _atom_dens
        self.volume: float = _volume  # Mass of the sample [g]
        self.nuclide_feed_rates: (dict, None) = None
        self.nuclide_removal_rates: (dict, None) = None
        self.case_dir = '_box_B'
        self.ORIGEN_input_file_name = 'origen.inp'
        self.DECAY_start_seconds: float = 0.0
        self.DECAY_end_seconds: float = 1000.0

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        if not self.skip_calculation:
            if len(self.atom_dens) > 0:
                with open(self.ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
                    f.write(atom_dens_for_origen(self.atom_dens))

            with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
                f.write(self.origen_deck())

            if self.debug > 0:
                print(f'ORIGEN: decaying sample for {self.DECAY_end_seconds - self.DECAY_start_seconds} seconds')
                print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
            run_scale(self.ORIGEN_input_file_name)
        else:
            if not os.path.isfile(self.F71_file_name):
                error_text: str = f'Skip SCALE flag set, output file {self.case_dir}/{self.F71_file_name} is not found!'
                raise ValueError(error_text)

        self.final_atom_dens = get_burned_material_atom_dens(self.F71_file_name,
                                                             self.DECAY_steps)
        os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.final_atom_dens)

    def origen_deck(self) -> str:
        """ Sample decay Origen deck """
        if len(self.atom_dens) > 0:
            iso: str = f'''        iso [
<{self.ATOM_DENS_file_name_Origen}
        ]\n'''
        else:
            iso: str = '        iso=0\n'

        removal: str = ''
        if self.nuclide_removal_rates:
            removal += 'processing {\n'
            for k, v in self.nuclide_removal_rates.items():
                removal += f'       removal {{ rate={v} ele=[{k}] }}\n'
            removal += '    }'

        feed: str = ''
        if self.nuclide_feed_rates:
            feed += 'feed=['
            for k, v in self.nuclide_feed_rates.items():
                feed += f'{k}={v} '
            feed += ']'

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
{iso}
        units=ATOMS-PER-BARN-CM
        volume={self.volume}
        {feed}
    }}
    time {{
        units=SECONDS
        t=[{time_interp_steps}I {self.DECAY_start_seconds} {self.DECAY_end_seconds}]
        start={self.DECAY_start_seconds}
    }}
{removal} 
    save {{
        file="{self.F71_file_name}"
    }}
}}
end
'''
        return origen_output


class LeakyBox:
    """ Handles two ORIGEN calculations, one with leakage one without.
    Calculates nuclide leak rate from box N, which is a corresponding feed rate into box N+!
    """

    def __init__(self):
        self.nuc_leak = None
        self.decay_leaks: (None, DecayBoxA) = None
        self.leak_rates: dict = {}  # Leak rate calculated as \epsilon_i^A N_i^A (t)
        self.removal_rate: float = PCTperDAY  # \epsilon_i^A
        self.nuclide_removal_rates: dict = {}
        self.adens: dict = {}  # box_A_adens

    def setup_cases(self, origen_case: DecayBoxA):
        """ Setup cases based on a base_case """
        self.nuclide_removal_rates = {
            'h': self.removal_rate,
            'he': self.removal_rate,
            'o': self.removal_rate,
            'n': self.removal_rate,
            'i': self.removal_rate,
            'ar': self.removal_rate,
            'kr': self.removal_rate,
            'xe': self.removal_rate
        }
        # self.nuclide_removal_rates = {'Kr': self.removal_rate,
        #                               'Xe': self.removal_rate}
        self.decay_leaks = origen_case
        self.decay_leaks.case_dir += '_leak'
        self.decay_leaks.nuclide_removal_rates = self.nuclide_removal_rates
        self.nuc_leak: list = [s.lower() for s in self.nuclide_removal_rates.keys()]

    def run_case(self):
        """ Run cases and get results """
        self.decay_leaks.run_decay_sample()

    def _parallel_leak_calc(self, i, vals) -> [int, float, dict, dict]:
        t = float(vals['time'])
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        adens_tot: dict = get_burned_material_atom_dens(f71_leak_file, i)  # All nuclides
        step_leak_rate: dict = {k: adens_tot[k] for k in adens_tot.keys() if re.sub('-.*', '', k) in self.nuc_leak}
        for k in step_leak_rate.keys():
            step_leak_rate[k] *= self.nuclide_removal_rates[re.sub('-.*', '', k)]  # * self.decay_leaks.volume
            print("LEAK_RATE: ", k, step_leak_rate[k], adens_tot[k])
        return i, t, adens_tot, step_leak_rate

    def get_leak_rate_parallel(self):
        """ Calculate Box B ingress rate = Box A leak rate"""
        from joblib import Parallel, delayed, cpu_count
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        times: dict = get_f71_positions_index(f71_leak_file)
        res_list = Parallel(n_jobs=cpu_count())(delayed(self._parallel_leak_calc)(ii, vv) for ii, vv in times.items())
        for k, res in enumerate(res_list):
            (i, t, adens_tot, step_leak_rate) = res
            self.leak_rates[i]: dict = {}
            self.leak_rates[i]['time'] = t
            self.leak_rates[i]['rate'] = step_leak_rate
            self.adens[i]: dict = {}
            self.adens[i]['time'] = t
            self.adens[i]['adens'] = adens_tot

    def get_leak_rate(self):
        """ Calculate Box B ingress rate = Box A leak rate"""
        # nuc_leak: list = [s.lower() for s in self.nuclide_removal_rates.keys()]
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        times: dict = get_f71_positions_index(f71_leak_file)
        for i, vals in times.items():  # Time-dependent difference between sealed and leaky box
            t = float(vals['time'])
            self.leak_rates[i]: dict = {}
            self.leak_rates[i]['time'] = t
            adens_tot: dict = get_burned_material_atom_dens(f71_leak_file, i)  # All nuclides
            step_leak_rate: dict = {k: adens_tot[k] for k in adens_tot.keys() if re.sub('-.*', '', k) in self.nuc_leak}
            for k in step_leak_rate.keys():
                step_leak_rate[k] *= self.nuclide_removal_rates[re.sub('-.*', '', k)]  # * self.decay_leaks.volume
                print("LEAK_RATE: ", k, step_leak_rate[k], adens_tot[k])
            self.leak_rates[i]['rate'] = step_leak_rate
            self.adens[i]: dict = {}
            self.adens[i]['time'] = t
            self.adens[i]['adens'] = adens_tot


class DiffLeak:
    """ Calculate difference between ORIGEN case with and without removals """

    def __init__(self):
        self.decay_base: (None, DecayBoxA, DecayBoxB) = None
        self.decay_leaks: (None, DecayBoxA, DecayBoxB) = None
        self.leaked: dict = {}
        self.removal_rate: float = PCTperDAY
        self.nuclide_removal_rates: dict = {}

    def setup_cases(self, base_case: (DecayBoxA, DecayBoxB)):
        """ Setup cases based on a base_case """
        self.nuclide_removal_rates = {
            'h': self.removal_rate,
            'he': self.removal_rate,
            'o': self.removal_rate,
            'n': self.removal_rate,
            'i': self.removal_rate,
            'ar': self.removal_rate,
            'kr': self.removal_rate,
            'xe': self.removal_rate
        }
        # self.nuclide_removal_rates = {'Kr': self.removal_rate,
        #                               'Xe': self.removal_rate}
        self.decay_base = base_case
        self.decay_leaks = copy.deepcopy(base_case)
        self.decay_leaks.case_dir = self.decay_base.case_dir + '_leak'
        self.decay_leaks.nuclide_removal_rates = self.nuclide_removal_rates

    def run_cases(self):
        """ Run cases and get results """
        self.decay_base.run_decay_sample()
        self.decay_leaks.run_decay_sample()

    def _parallel_leak_calc(self, i, vals) -> dict:
        f71_base_file: str = self.decay_base.case_dir + '/' + self.decay_base.F71_file_name
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        t = float(vals['time'])
        print(f"Parallel leak: reading {i}, {t} s")
        adens_base: dict = get_burned_material_atom_dens(f71_base_file, i)
        adens_leak: dict = get_burned_material_atom_dens(f71_leak_file, i)
        adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base
                            if (adens_base[key] - adens_leak.get(key, 0)) > 0}
        return adens_diff

    def get_diff_parallel(self):
        from joblib import Parallel, delayed, cpu_count
        """ Calculate difference between sealed and leaky box in parallel """
        f71_base_file: str = self.decay_base.case_dir + '/' + self.decay_base.F71_file_name
        times: dict = get_f71_positions_index(f71_base_file)
        for k, v in times.items():
            self.leaked[k]: dict = {}
            self.leaked[k]['time'] = float(v['time'])
        res_list = Parallel(n_jobs=cpu_count())(delayed(self._parallel_leak_calc)(k, v) for k, v in times.items())
        for i, res in enumerate(res_list):
            print(f'Leak diff: {i+1}, {res}')
            self.leaked[i + 1]['adens'] = res

    def get_diff(self):
        """ Calculate difference between sealed and leaky box """
        f71_base_file: str = self.decay_base.case_dir + '/' + self.decay_base.F71_file_name
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        times: dict = get_f71_positions_index(f71_base_file)
        for i, vals in times.items():  # Time-dependent difference between sealed and leaky box
            t = float(vals['time'])
            self.leaked[i]: dict = {}
            self.leaked[i]['time'] = t
            adens_base: dict = get_burned_material_atom_dens(f71_base_file, i)
            adens_leak: dict = get_burned_material_atom_dens(f71_leak_file, i)
            adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base
                                if (adens_base[key] - adens_leak.get(key, 0)) > 0}
            # adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base}
            self.leaked[i]['adens'] = adens_diff
            print(i, t, adens_diff)

    def get_rate(self):
        """ Calculate difference between sealed and leaky box, and get equivalent ingress rate
            [obsolete] """
        f71_base_file: str = self.decay_base.case_dir + '/' + self.decay_base.F71_file_name
        f71_leak_file: str = self.decay_leaks.case_dir + '/' + self.decay_leaks.F71_file_name
        times: dict = get_f71_positions_index(f71_base_file)
        for i, vals in times.items():  # Time-dependent difference between sealed and leaky box
            t = float(vals['time'])
            self.leaked[i]: dict = {}
            self.leaked[i]['time'] = t
            adens_base: dict = get_burned_material_atom_dens(f71_base_file, i)
            adens_leak: dict = get_burned_material_atom_dens(f71_leak_file, i)
            adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base
                                if (adens_base[key] - adens_leak.get(key, 0)) > 0}
            # adens_diff: dict = {key: adens_base[key] - adens_leak.get(key, 0) for key in adens_base}
            self.leaked[i]['adens'] = adens_diff
            print(i, t, adens_diff)
    #     """ rate1 is calculated from leakage rate and the original concentration """
    #     rate1: dict = {key: adens_base[key] * self.removal_rate for key in adens_base
    #                    if adens_base[key] > 0 and  # if nonzero
    #                    next((True for s in [k.lower() for k in self.nuclide_removal_rates.keys()] if s in key),
    #                         False)  # only for elements that leak
    #                    }
    #     self.leaked[i]['feed_rate1'] = rate1
    #     print(i, t, rate1)
    #
    # prev_t: float = -1.0
    # prev_adens: dict = {}
    # for i, vals in self.leaked.items():  # Calculate leakage, turn it into feed rate
    #     if i == 1:
    #         prev_t = vals['time']
    #         prev_adens = self.leaked[i]['adens']
    #         continue
    #     dt: float = vals['time'] - prev_t
    #     adens: dict = self.leaked[i]['adens']
    #     """ rate2 is calculated using time derivative of concentration """
    #     rate2: dict = {key: (adens[key] - prev_adens[key]) / dt for key in adens  # leak rate
    #                    if (adens[key] - prev_adens[key]) > 0 and  # if nonzero
    #                    next((True for s in [k.lower() for k in self.nuclide_removal_rates.keys()] if s in key),
    #                         False)  # only for elements that leak
    #                    }
    #     self.leaked[i]['feed_rate2'] = rate2
    #     print(i, vals['time'], dt, rate2)
    #     prev_t = vals['time']
    #     prev_adens = adens


def get_dataframe(leaky_box_adens: dict[dict]) -> pd.DataFrame:
    """ Convert """
    _list = []
    for k, v in leaky_box_adens.items():
        a = v['adens']
        a['time [s]'] = float(v['time'])
        a['time [d]'] = float(v['time']) / float(24 * 60 * 60)
        _list.append(a)
    return pd.DataFrame(_list)  # .set_index('time [s]', append=True)


is_xe_136_testing: bool = False  # Use 1 Xe-136 at / b-sm instead of read composition
is_xe_135_testing: bool = True  # Use 1 Xe-135 at / b-sm instead of read composition
if is_xe_135_testing and is_xe_136_testing:
    raise ValueError("Tests are exclusive!")

do_skip_calcs: bool = False     # If True, skips SCALE calculations & re-reads existing data


def main():
    """ Calculate which nuclides leak """
    box_A_leak_rate: float = PCTperDAY
    box_B_leak_rate: float = PCTperDAY * 0.1
    # Baseline decay calculation setup
    origen_triton_box_A = DecayBoxA('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/SCALE_FILE.f71', 1200e3)
    if do_skip_calcs:
        origen_triton_box_A.skip_calculation = True
    origen_triton_box_A.set_f71_pos(5.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 5 years
    origen_triton_box_A.read_burned_material()
    origen_triton_box_A.DECAY_days = 12
    origen_triton_box_A.DECAY_steps = 12 * 2 + 2
    if is_xe_136_testing:
        origen_triton_box_A.atom_dens = {'xe-136': 1.0}
    if is_xe_135_testing:
        origen_triton_box_A.atom_dens = {'xe-135': 1.0}
    volume: float = 1.0  # origen_triton_box_A.volume
    print(volume)

    box_A = LeakyBox()
    box_A.removal_rate = box_A_leak_rate
    box_A.setup_cases(origen_triton_box_A)
    box_A.run_case()
    box_A.get_leak_rate_parallel()
    with open('boxA.json5', 'w') as f:  # Save to json
        json5.dump(box_A.adens, f, indent=4)
    with open('boxA_leak_rate.json5', 'w') as f:  # Save to json
        json5.dump(box_A.leak_rates, f, indent=4)

    origen_box_B: dict = {}  # ORIGEN runs that describe decay of box_B nuclides without leak
    box_B: dict = {}
    box_B_adens: dict = {}  # Summary of box_B atomic density as a function of decay time
    box_B_adens_current: dict = {}  # Content of box B
    dl_C: dict = {}  # Leaking into box_C calculated from box_B_origen runs
    origen_box_C: dict = {}
    box_C_adens: dict = {}  # Summary of box_C atomic density as a function of decay time
    box_C_adens_last: dict = {}  # Content of box C
    prev_time_seconds: float = 0.0

    for k, v in box_A.leak_rates.items():
        print(f'Box B {k}: {v}')
        if k == 1:  # Skip the first time step, there is nothing in box_B
            prev_time_seconds = v['time']
            continue
        origen_box_B[k] = DecayBoxB(box_B_adens_current, volume)
        if do_skip_calcs:
            origen_box_B[k].skip_calculation = True
        origen_box_B[k].nuclide_feed_rates = box_A.leak_rates[k]['rate']
        origen_box_B[k].DECAY_steps = 4
        origen_box_B[k].DECAY_start_seconds = prev_time_seconds
        origen_box_B[k].DECAY_end_seconds = v['time']
        origen_box_B[k].case_dir = f'_box_B_{k:04d}'

        box_B[k] = DiffLeak()
        box_B[k].removal_rate = box_B_leak_rate
        box_B[k].setup_cases(origen_box_B[k])
        box_B[k].run_cases()
        box_B[k].get_diff_parallel()

        box_B_adens_current = box_B[k].decay_leaks.final_atom_dens
        box_B_adens[k] = {}
        box_B_adens[k]['time'] = v['time']
        box_B_adens[k]['adens'] = box_B_adens_current

        # Now feed box C
        leaked_adens = box_B[k].leaked[origen_box_B[k].DECAY_steps]['adens']
        print(f'box_C_adens_last {box_C_adens_last}')
        print(f'leaked_adens {k} {leaked_adens}')
        box_C_adens_new = {iso: box_C_adens_last.get(iso, 0) + leaked_adens.get(iso, 0)
                           for iso in set(box_C_adens_last) | set(leaked_adens)}
        origen_box_C[k] = DecayBoxB(box_C_adens_new, volume)
        origen_box_C[k].DECAY_steps = 4
        origen_box_C[k].DECAY_start_seconds = prev_time_seconds
        origen_box_C[k].DECAY_end_seconds = v['time']
        origen_box_C[k].case_dir = f'_box_C_{k:04d}'
        origen_box_C[k].run_decay_sample()

        box_C_adens_last = origen_box_C[k].final_atom_dens
        box_C_adens[k] = {}
        box_C_adens[k]['time'] = v['time']
        box_C_adens[k]['adens'] = box_C_adens_last

        # dl_C[k] = DiffLeak()
        # dl_C[k].removal_rate = box_B_leak_rate  # Removal rate from box B to box_C
        # dl_C[k].setup_cases(origen_box_B[k])
        # dl_C[k].run_cases()
        #
        # box_B_adens_current = dl_C[k].decay_leaks.final_atom_dens
        # box_B_adens[k] = {}
        # box_B_adens[k]['time'] = v['time']
        # box_B_adens[k]['adens'] = box_B_adens_current
        # # dl_C[k].get_diff()
        # dl_C[k].get_diff_parallel()
        # box_C_adens[k] = {}
        # box_C_adens[k]['time'] = v['time']
        # box_C_adens[k]['adens'] = dl_C[k].leaked[origen_box_B[k].DECAY_steps]['adens']
        prev_time_seconds = v['time']
        print(k, v, box_C_adens[k])

    with open('boxB.json5', 'w') as f:  # Save to json
        json5.dump(box_B_adens, f, indent=4)
    with open('boxC.json5', 'w') as f:
        json5.dump(box_C_adens, f, indent=4)

    pd_A: pd.DataFrame = get_dataframe(box_A.adens)
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
                            (np.exp(-box_A_leak_rate * pd_B['time [s]']) - np.exp(-box_B_leak_rate * pd_B['time [s]']))
                            / (box_B_leak_rate - box_A_leak_rate)
                            )
        pd_C['analytic'] = (1.0 *
                            (-box_A_leak_rate * np.exp(-box_B_leak_rate * pd_C['time [s]']) +
                             box_B_leak_rate * (np.exp(-box_A_leak_rate * pd_C['time [s]']) - 1.0) + box_A_leak_rate)
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
N_i^B(t) & = N_i^A \epsilon_i^A 
\frac{e^{-(\epsilon_i^A + \lambda_i)  t} - e^{-(\epsilon_i^B + \lambda_i) t}}
{\epsilon_i^B - \epsilon_i^A} \\
% B(t) = -(a N (e^(t (-(a + l))) - e^(t (-(b + l)))))/(a - b)
\\
N_i^C(t) & = N_i^A e^{ - \lambda_i  t} 
\frac{ -\epsilon_i^B e^{-\epsilon_i^A t} + \epsilon_i^A (e^{-\epsilon_i^B t} - 1) + \epsilon_i^B}
{\epsilon_i^B - \epsilon_i^A}
        """
        pd_A['analytic'] = (1.0 * np.exp(-(box_A_leak_rate + lambda_xe_135) * pd_A['time [s]']))
        pd_B['analytic'] = (1.0 *
                            box_A_leak_rate *
                            (np.exp(-(box_A_leak_rate + lambda_xe_135) * pd_B['time [s]']) - np.exp(
                                -(box_B_leak_rate + lambda_xe_135) * pd_B['time [s]']))
                            / (box_B_leak_rate - box_A_leak_rate)
                            )
        pd_C['analytic'] = (1.0 * np.exp(-lambda_xe_135 * pd_C['time [s]']) *
                            (-box_B_leak_rate * np.exp(-box_A_leak_rate * pd_C['time [s]']) +
                             box_A_leak_rate * (
                                     np.exp(-box_B_leak_rate * pd_C['time [s]']) - 1.0) + box_B_leak_rate)
                            / (box_B_leak_rate - box_A_leak_rate)
                            )
        pd_A['diff [%]'] = 100.0 * (pd_A['xe-135'] - pd_A['analytic']) / pd_A['analytic']
        pd_B['diff [%]'] = 100.0 * (pd_B['total'] - pd_B['analytic']) / pd_B['analytic']
        pd_C['diff [%]'] = 100.0 * (pd_C['total'] - pd_C['analytic']) / pd_C['analytic']

    writer = pd.ExcelWriter('leaky_boxes.xlsx')

    pd_A.to_excel(writer, sheet_name='box A')
    pd_B.to_excel(writer, sheet_name='box B')
    pd_C.to_excel(writer, sheet_name='box C')
    writer.close()


def main_old():
    """ Calculate which nuclides leak using concentrations in box B """
    # Baseline decay calculation setup
    origen_triton = DecayBoxA('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/SCALE_FILE.f71', 1200e3)
    origen_triton.set_f71_pos(5.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 5 years
    origen_triton.read_burned_material()
    origen_triton.DECAY_days = 12
    origen_triton.DECAY_steps = 120
    if is_xe_136_testing:
        origen_triton.atom_dens = {'xe-136': 1.0}
    if is_xe_135_testing:
        origen_triton.atom_dens = {'xe-135': 1.0}
    volume: float = origen_triton.volume
    print(volume)

    # Calculate nuclide atom density in box B using the (1-2) difference between (1) box_A decay with no leakage
    # and (2) box_A decay with leakage
    dl_B = DiffLeak()
    dl_B.setup_cases(origen_triton)
    dl_B.run_cases()
    # dl_B.get_diff()
    dl_B.get_diff_parallel()

    # Calculate concentrations of nuclides leaking from box_B as the difference between box_B decays without leakage
    # minus box_B decay with leakage.
    # Since the concentrations differ at each timestep, each timestep is its own ORIGEN run.

    box_B_origen: dict = {}  # ORIGEN runs that describe decay of box_B nuclides
    dl_C: dict = {}  # Leaking into box_C calculated from box_B_origen runs
    leaked_adens: dict = {}  # Atomic density that leaked into box_C in the current decay time step
    box_C_adens: dict = {}  # Summary of box_C atomic density as a function of decay time
    for k, v in dl_B.leaked.items():
        print(k, v)
        if k == 1:  # Skip the first time step, there is nothing in box_B
            continue
        box_B_adens_orig = v['adens']
        # Remove nuclides that leaked in the previous step
        box_B_adens_current: dict = {key: box_B_adens_orig[key] - leaked_adens.get(key, 0) for key in box_B_adens_orig}

        box_B_origen[k] = DecayBoxB(box_B_adens_current, volume)
        box_B_origen[k].DECAY_steps = 5
        box_B_origen[k].DECAY_end_seconds = v['time']
        box_B_origen[k].case_dir = f'box_B_{k:03d}'

        dl_C[k] = DiffLeak()
        dl_C[k].removal_rate = 0.1 * PCTperDAY  # Removal rate to box_C
        dl_C[k].setup_cases(box_B_origen[k])
        dl_C[k].run_cases()
        # dl_C[k].get_diff_rate()
        dl_C[k].get_diff_parallel()
        leaked_adens = dl_C[k].leaked[box_B_origen[k].DECAY_steps]['adens']
        box_C_adens[k] = {}
        box_C_adens[k]['time'] = v['time']
        box_C_adens[k]['adens'] = leaked_adens
        print(k, v, box_C_adens[k])

    with open('boxB.json5', 'w') as f:  # Save to json
        json5.dump(dl_B.leaked, f, indent=4)
    with open('boxC.json5', 'w') as f:
        json5.dump(box_C_adens, f, indent=4)

    pd_B = get_dataframe(dl_B.leaked)
    pd_C = get_dataframe(box_C_adens)
    writer = pd.ExcelWriter('leaky_boxes.xlsx')

    pd_B.to_excel(writer, 'box B')
    pd_C.to_excel(writer, 'box C')
    writer.close()


"""
df=pd.read_excel("./leaky_boxes.xlsx")
a=0.000000115740740741
b=0.0000000115740740741
# https://www.wolframalpha.com/input?i=dM%2Fdt+%3D+N*a+-+dU%2Fdt+%3B+dU%2Fdt+%3D+M*b%2C+M%280%29+%3D+0%2C+U%280%29+%3D+0
df['anaB'] = a * (1.0 - np.exp(-b*df['time'])) /b
df['anaC'] = a * np.exp(-b*df['time']) * ( np.exp(b*df['time']) * (b*df['time'] -1.0) +1.0) /b
"""

if __name__ == "__main__":
    # pass
    main()
