#!/bin/env python3
"""
Legacy leaky-box variant using SCALE/Origen. Kept in sync with LeakyBox.py.

Explanation of differences vs LeakyBox.py:
- LeakyBox.py is the maintained implementation and recommended entry point.
- BoxF71.py historically experimented with a DiffLeak (sealed vs leaky difference) approach.
- The refactored main path now mirrors LeakyBox.py logic (interval-averaged A->B feed,
  time-weighted B->C feed, analytic comparisons, plots), but `main_old()` preserves the
  older DiffLeak workflow for reference.

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
from sample_decay_dose.SampleDose import NOW
from sample_decay_dose.utils import nicely_print_atom_dens, get_f71_positions_index, get_burned_nuclide_atom_dens, \
    get_burned_material_total_mass_dens, run_scale, atom_dens_for_origen

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
        self.weight: float = np.nan  # Mass of the sample [g]
        self.density: float = np.nan  # Mass density of the sample [g/cm3]
        self.volume: float = np.nan  # Sample volume [cm3]
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
        self.atom_dens = get_burned_nuclide_atom_dens(self.BURNED_MATERIAL_F71_file_name,
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

        self.final_atom_dens = get_burned_nuclide_atom_dens(self.F71_file_name,
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
        t=[{time_interp_steps}I 0.01 {self.DECAY_days}]
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
            else:
                # Ensure the file exists; the ORIGEN deck always copies it in the shell preamble.
                open(self.ATOM_DENS_file_name_Origen, 'w').close()

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

        self.final_atom_dens = get_burned_nuclide_atom_dens(self.F71_file_name,
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
    """ Handles the ORIGEN calculations with leakage.
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
        adens_tot: dict = get_burned_nuclide_atom_dens(f71_leak_file, i)  # All nuclides
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
            adens_tot: dict = get_burned_nuclide_atom_dens(f71_leak_file, i)  # All nuclides
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
        adens_base: dict = get_burned_nuclide_atom_dens(f71_base_file, i)
        adens_leak: dict = get_burned_nuclide_atom_dens(f71_leak_file, i)
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
            adens_base: dict = get_burned_nuclide_atom_dens(f71_base_file, i)
            adens_leak: dict = get_burned_nuclide_atom_dens(f71_leak_file, i)
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
            adens_base: dict = get_burned_nuclide_atom_dens(f71_base_file, i)
            adens_leak: dict = get_burned_nuclide_atom_dens(f71_leak_file, i)
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


def _sorted_steps_by_time(step_dict: dict) -> list[tuple[int, dict]]:
    # Ensure deterministic chronological order for interval math.
    return sorted(step_dict.items(), key=lambda kv: float(kv[1]['time']))


def _average_rates(prev_rates: dict, cur_rates: dict, volume_ratio: float = 1.0) -> dict:
    # Trapezoidal (endpoint) average of feed rates over a single interval.
    avg_rates: dict = {}
    if prev_rates is None:
        prev_rates = {}
    if cur_rates is None:
        cur_rates = {}
    for iso in set(prev_rates.keys()) | set(cur_rates.keys()):
        avg_rates[iso] = 0.5 * (prev_rates.get(iso, 0.0) + cur_rates.get(iso, 0.0)) * volume_ratio
    return avg_rates


def _time_average_rates(leak_rates: dict, volume_ratio: float = 1.0) -> dict:
    # Time-weighted average of rates across all sub-intervals in a run.
    steps = _sorted_steps_by_time(leak_rates)
    if len(steps) < 2:
        return {}
    total_dt: float = 0.0
    accum: dict = {}
    for (_, prev_v), (_, cur_v) in zip(steps, steps[1:]):
        t0 = float(prev_v['time'])
        t1 = float(cur_v['time'])
        dt = t1 - t0
        if dt <= 0.0:
            continue
        total_dt += dt
        prev_rates = prev_v.get('rate', {})
        cur_rates = cur_v.get('rate', {})
        for iso in set(prev_rates.keys()) | set(cur_rates.keys()):
            r0 = prev_rates.get(iso, 0.0)
            r1 = cur_rates.get(iso, 0.0)
            accum[iso] = accum.get(iso, 0.0) + 0.5 * (r0 + r1) * dt
    if total_dt <= 0.0:
        return {}
    return {iso: (val / total_dt) * volume_ratio for iso, val in accum.items()}


do_skip_calcs: bool = False  # If True, skips SCALE calculations & re-reads existing data
apply_volume_scaling: bool = False  # If True, scales feed by source_volume / dest_volume


def _out_name(prefix: str | None, stem: str, ext: str) -> str:
    if prefix:
        return f"{stem}_{prefix}{ext}"
    return f"{stem}{ext}"


def _add_analytic_columns(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, isotope: str,
                          n0_density: float, n0_total: float, eps_a: float, eps_b: float,
                          lambda_decay: float | None) -> None:
    if lambda_decay is None:
        lambda_decay = 0.0

    pd_A['analytic'] = (n0_density * np.exp(-(eps_a + lambda_decay) * pd_A['time [s]']))
    pd_B['analytic'] = (n0_total *
                        eps_a *
                        (np.exp(-(eps_a + lambda_decay) * pd_B['time [s]']) -
                         np.exp(-(eps_b + lambda_decay) * pd_B['time [s]']))
                        / (eps_b - eps_a)
                        )
    pd_C['analytic'] = (n0_total * np.exp(-lambda_decay * pd_C['time [s]']) *
                        (-eps_b * np.exp(-eps_a * pd_C['time [s]']) +
                         eps_a * (np.exp(-eps_b * pd_C['time [s]']) - 1.0) + eps_b)
                        / (eps_b - eps_a)
                        )

    pd_A['diff [%]'] = 100.0 * (pd_A[isotope] - pd_A['analytic']) / pd_A['analytic']
    pd_B['diff [%]'] = 100.0 * (pd_B['total'] - pd_B['analytic']) / pd_B['analytic']
    pd_C['diff [%]'] = 100.0 * (pd_C['total'] - pd_C['analytic']) / pd_C['analytic']


def plot_results(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, isotope: str, out_prefix: str | None,
                 logy: bool = True) -> str:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 10))

    # Box A (per-volume densities)
    axes[0].plot(pd_A['time [d]'], pd_A[isotope], label='ORIGEN', color='tab:blue')
    if 'analytic' in pd_A.columns:
        axes[0].plot(pd_A['time [d]'], pd_A['analytic'], '--', label='Analytic', color='tab:orange')
    axes[0].set_title(f'Box A ({isotope})')
    axes[0].set_ylabel('atoms/b-cm')

    # Box B (total atoms)
    axes[1].plot(pd_B['time [d]'], pd_B['total'], label='ORIGEN', color='tab:blue')
    if 'analytic' in pd_B.columns:
        axes[1].plot(pd_B['time [d]'], pd_B['analytic'], '--', label='Analytic', color='tab:orange')
    axes[1].set_title(f'Box B ({isotope})')
    axes[1].set_ylabel('total atoms')

    # Box C (total atoms)
    axes[2].plot(pd_C['time [d]'], pd_C['total'], label='ORIGEN', color='tab:blue')
    if 'analytic' in pd_C.columns:
        axes[2].plot(pd_C['time [d]'], pd_C['analytic'], '--', label='Analytic', color='tab:orange')
    axes[2].set_title(f'Box C ({isotope})')
    axes[2].set_ylabel('total atoms')
    axes[2].set_xlabel('time [days]')

    if logy:
        for ax in axes:
            ax.set_yscale('log')

    for ax in axes:
        ax.grid(True, which='both', linestyle=':', alpha=0.5)
        ax.legend()

    fig.tight_layout()
    fname = _out_name(out_prefix, 'leaky_boxes', '.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def run_test(isotope: str, lambda_decay: float | None, apply_volume_scaling: bool,
             do_skip_calcs: bool, out_prefix: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run a single isotope test and return dataframes."""
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
    origen_triton_box_A.atom_dens = {isotope: 1.0}

    volume: float = 1.0  # origen_triton_box_A.volume
    print(volume)

    box_A = LeakyBox()
    box_A.removal_rate = box_A_leak_rate
    box_A.setup_cases(origen_triton_box_A)
    box_A.run_case()
    box_A.get_leak_rate_parallel()

    with open(_out_name(out_prefix, 'boxA', '.json5'), 'w') as f:  # Save to json
        json5.dump(box_A.adens, f, indent=4)
    with open(_out_name(out_prefix, 'boxA_leak_rate', '.json5'), 'w') as f:  # Save to json
        json5.dump(box_A.leak_rates, f, indent=4)

    origen_box_B: dict = {}  # ORIGEN runs that describe decay of box_B nuclides without leak
    box_B: dict = {}
    box_B_adens: dict = {}  # Summary of box_B atomic density as a function of decay time
    box_B_adens_current: dict = {}  # Content of box B
    origen_box_C: dict = {}
    box_C_adens: dict = {}  # Summary of box_C atomic density as a function of decay time
    box_C_adens_current: dict = {}  # Content of box C
    box_B_n_steps: int = 8
    box_C_n_steps: int = 8

    # Iterate by time, not by dict insertion order, to keep interval math correct.
    box_A_steps = _sorted_steps_by_time(box_A.leak_rates)
    for idx in range(1, len(box_A_steps)):
        prev_k, prev_v = box_A_steps[idx - 1]
        k, v = box_A_steps[idx]
        start_t = float(prev_v['time'])
        end_t = float(v['time'])
        print(f'Box B {k}: {v}')
        if end_t <= start_t:
            continue
        origen_box_B[k] = DecayBoxB(box_B_adens_current, volume)
        if do_skip_calcs:
            origen_box_B[k].skip_calculation = True
        volume_ratio = 1.0
        if apply_volume_scaling:
            source_vol = box_A.decay_leaks.volume
            dest_vol = origen_box_B[k].volume
            volume_ratio = source_vol / dest_vol if dest_vol else 1.0
        # Feed into B is the interval-average leak rate from A (per unit volume).
        origen_box_B[k].nuclide_feed_rates = _average_rates(prev_v['rate'], v['rate'], volume_ratio)
        origen_box_B[k].DECAY_steps = box_B_n_steps
        origen_box_B[k].DECAY_start_seconds = start_t
        origen_box_B[k].DECAY_end_seconds = end_t
        origen_box_B[k].case_dir = f'_box_B_{k:04d}'

        box_B[k] = LeakyBox()
        box_B[k].removal_rate = box_B_leak_rate
        box_B[k].setup_cases(origen_box_B[k])
        box_B[k].run_case()
        box_B[k].get_leak_rate_parallel()

        box_B_adens_current = origen_box_B[k].final_atom_dens
        box_B_adens[k] = {}
        box_B_adens[k]['time'] = v['time']
        box_B_adens[k]['adens'] = box_B_adens_current

    # Now feed box C
    for idx in range(1, len(box_A_steps)):
        prev_k, prev_v = box_A_steps[idx - 1]
        k, v = box_A_steps[idx]
        start_t = float(prev_v['time'])
        end_t = float(v['time'])
        print(f'Box C {k}: {v}')
        if end_t <= start_t:
            continue

        origen_box_C[k] = DecayBoxB(box_C_adens_current, volume)
        if do_skip_calcs:
            origen_box_C[k].skip_calculation = True
        volume_ratio = 1.0
        if apply_volume_scaling:
            source_vol = box_B[k].decay_leaks.volume
            dest_vol = origen_box_C[k].volume
            volume_ratio = source_vol / dest_vol if dest_vol else 1.0
        # Feed into C is the time-weighted average of B's leak rates across the whole sub-run.
        origen_box_C[k].nuclide_feed_rates = _time_average_rates(box_B[k].leak_rates, volume_ratio)

        origen_box_C[k].DECAY_steps = box_C_n_steps
        origen_box_C[k].DECAY_start_seconds = start_t
        origen_box_C[k].DECAY_end_seconds = end_t
        origen_box_C[k].case_dir = f'_box_C_{k:04d}'
        origen_box_C[k].run_decay_sample()

        box_C_adens_current = origen_box_C[k].final_atom_dens
        box_C_adens[k] = {}
        box_C_adens[k]['time'] = v['time']
        box_C_adens[k]['adens'] = box_C_adens_current

    with open(_out_name(out_prefix, 'boxB', '.json5'), 'w') as f:  # Save to json
        json5.dump(box_B_adens, f, indent=4)
    with open(_out_name(out_prefix, 'boxC', '.json5'), 'w') as f:
        json5.dump(box_C_adens, f, indent=4)

    pd_A: pd.DataFrame = get_dataframe(box_A.adens)
    pd_B: pd.DataFrame = get_dataframe(box_B_adens)
    pd_C: pd.DataFrame = get_dataframe(box_C_adens)

    pd_B['total'] = pd_B[isotope] * volume
    pd_C['total'] = pd_C[isotope] * volume

    n0_density = float(origen_triton_box_A.atom_dens.get(isotope, 0.0))
    n0_total = n0_density * (origen_triton_box_A.volume if apply_volume_scaling else volume)

    _add_analytic_columns(pd_A, pd_B, pd_C, isotope, n0_density, n0_total,
                          box_A_leak_rate, box_B_leak_rate, lambda_decay)

    writer = pd.ExcelWriter(_out_name(out_prefix, 'leaky_boxes', '.xlsx'))
    pd_A.to_excel(writer, sheet_name='box A')
    pd_B.to_excel(writer, sheet_name='box B')
    pd_C.to_excel(writer, sheet_name='box C')
    writer.close()

    plot_results(pd_A, pd_B, pd_C, isotope, out_prefix)
    return pd_A, pd_B, pd_C


def run_all_tests(apply_volume_scaling: bool, do_skip_calcs: bool):
    tests = {
        'xe-136': None,
        'xe-135': 2.106574217602 * 1e-5,  # Xe-135 decay constant [1/s]
    }
    for iso, lam in tests.items():
        prefix = iso.replace('-', '')
        run_test(iso, lam, apply_volume_scaling, do_skip_calcs, out_prefix=prefix)


def main():
    """Default entrypoint: run Xe-136 test."""
    run_test('xe-136', None, apply_volume_scaling, do_skip_calcs, out_prefix='xe136')


def main_old(isotope: str = 'xe-135'):
    """Legacy method using DiffLeak; kept for reference."""
    # Baseline decay calculation setup
    origen_triton = DecayBoxA('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/SCALE_FILE.f71', 1200e3)
    origen_triton.set_f71_pos(5.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 5 years
    origen_triton.read_burned_material()
    origen_triton.DECAY_days = 12
    origen_triton.DECAY_steps = 120
    origen_triton.atom_dens = {isotope: 1.0}
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
