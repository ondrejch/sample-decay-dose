#!/bin/env python3
"""
Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

Box A --> Box B --> outside

100% per day  =  1/24/60/60 per second    = 0.0000115740740741
   1% per day  =  1e-2/24/60/60 per second = 0.000000115740740741
 0.1% per day  =  1e-3/24/60/60 per second = 0.0000000115740740741
"""
import os
import re
import numpy as np
import pandas as pd
from bisect import bisect_left
import json5
import math
from sample_decay_dose.SampleDose import NOW
from sample_decay_dose.utils import nicely_print_atom_dens, get_f71_positions_index, get_burned_nuclide_atom_dens, \
    get_burned_nuclide_data, get_burned_material_total_mass_dens, run_scale, atom_dens_for_origen

PCTperDAY: float = 0.000000115740740741
DEFAULT_F71_PATH: str = os.getenv(
    'LEAKYBOX_F71_PATH',
    os.path.expanduser('~/0.03/20-burn-MHA/mha-4.5-a4/msrr.f71'),
)
DEFAULT_F71_CASE: str = os.getenv('LEAKYBOX_F71_CASE', '20')


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
        pos_times = sorted(
            (k, float(v['time'])) for k, v in self.BURNED_MATERIAL_F71_index.items() if v['case'] == case
        )
        if not pos_times:
            raise ValueError(f"Case '{case}' not found in F71 index for {self.BURNED_MATERIAL_F71_file_name}")
        times = [x[1] for x in pos_times]
        t_min = min(times)
        t_max = max(times)
        if t < t_min:
            print(f"Error: Time {t} seconds is less than {t_min} s, the minimum time in records.")
            pos_idx: int = 0
        elif t > t_max:
            print(f"Error: Time {t} seconds is longer than {t_max} s, the maximum time in records.")
            pos_idx: int = len(times) - 1
        else:
            pos_idx: int = bisect_left(times, t)
        pos_number: int = pos_times[pos_idx][0]
        print(f'--> Closest F71 position found at slot {pos_number}, {times[pos_idx]} seconds')
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
            'br': self.removal_rate,
            'ar': self.removal_rate,
            'kr': self.removal_rate,
            'xe': self.removal_rate,
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


def get_dataframe(leaky_box_adens: dict[dict]) -> pd.DataFrame:
    """ Convert """
    _list = []
    for k, v in leaky_box_adens.items():
        a = dict(v['adens'])
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


def _write_excel(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, out_prefix: str | None) -> None:
    # Persist A/B/C time series to a workbook for quick inspection.
    writer = pd.ExcelWriter(_out_name(out_prefix, 'leaky_boxes', '.xlsx'))
    pd_A.to_excel(writer, sheet_name='box A')
    pd_B.to_excel(writer, sheet_name='box B')
    pd_C.to_excel(writer, sheet_name='box C')
    writer.close()


def _write_dataframes_json(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, out_prefix: str | None) -> str:
    # JSON snapshot of the three dataframes (column -> list).
    payload = {
        'box_A': pd_A.to_dict(orient='list'),
        'box_B': pd_B.to_dict(orient='list'),
        'box_C': pd_C.to_dict(orient='list'),
    }
    fname = _out_name(out_prefix, 'leaky_boxes', '.json')
    with open(fname, 'w') as f:
        json5.dump(payload, f, indent=2)
    return fname


def _add_analytic_columns(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, isotope: str,
                          n0_density: float, n0_total: float, eps_a: float, eps_b: float,
                          lambda_decay: float | None) -> None:
    # Fill analytic solutions and relative differences for the test isotope.
    if lambda_decay is None:
        lambda_decay = 0.0

    t_a = pd_A['time [s]']
    t_b = pd_B['time [s]']
    t_c = pd_C['time [s]']
    pd_A['analytic'] = n0_density * np.exp(-(eps_a + lambda_decay) * t_a)

    if math.isclose(eps_a, eps_b, rel_tol=0.0, abs_tol=1e-18):
        eps = eps_a
        pd_B['analytic'] = n0_total * eps * t_b * np.exp(-(eps + lambda_decay) * t_b)
        pd_C['analytic'] = n0_total * np.exp(-lambda_decay * t_c) * (1.0 - (1.0 + eps * t_c) * np.exp(-eps * t_c))
    else:
        pd_B['analytic'] = (n0_total *
                            eps_a *
                            (np.exp(-(eps_a + lambda_decay) * t_b) -
                             np.exp(-(eps_b + lambda_decay) * t_b))
                            / (eps_b - eps_a)
                            )
        pd_C['analytic'] = (n0_total * np.exp(-lambda_decay * t_c) *
                            (-eps_b * np.exp(-eps_a * t_c) +
                             eps_a * (np.exp(-eps_b * t_c) - 1.0) + eps_b)
                            / (eps_b - eps_a)
                            )

    def _safe_pct_diff(obs: pd.Series, ana: pd.Series) -> pd.Series:
        denom = np.where(np.abs(ana.to_numpy()) > 0.0, ana.to_numpy(), np.nan)
        return pd.Series(100.0 * (obs.to_numpy() - ana.to_numpy()) / denom, index=ana.index)

    pd_A['diff [%]'] = _safe_pct_diff(pd_A[isotope], pd_A['analytic'])
    pd_B['diff [%]'] = _safe_pct_diff(pd_B['total'], pd_B['analytic'])
    pd_C['diff [%]'] = _safe_pct_diff(pd_C['total'], pd_C['analytic'])


def _line_style(n_points: int, marker: str = 'o', hollow: bool = False) -> dict:
    # Small point markers on all series; downsample markers only for very dense lines.
    markevery = max(1, n_points // 200)
    style = {
        'linewidth': 1.4,
        'marker': marker,
        'markersize': 2.4,
        'markeredgewidth': 0.6,
        'markevery': markevery,
    }
    if hollow:
        style['markerfacecolor'] = 'none'
    return style


def plot_results(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, isotope: str, out_prefix: str | None,
                 logy: bool = True) -> str:
    # Plot analytic vs ORIGEN for A (density) and B/C (totals).
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 10))

    # Box A (per-volume densities)
    axes[0].plot(pd_A['time [d]'], pd_A[isotope], label='ORIGEN', color='tab:blue',
                 **_line_style(len(pd_A), marker='o'))
    if 'analytic' in pd_A.columns:
        axes[0].plot(pd_A['time [d]'], pd_A['analytic'], '--', label='Analytic', color='tab:orange',
                     **_line_style(len(pd_A), marker='s', hollow=True))
    axes[0].set_title(f'Box A ({isotope})')
    axes[0].set_ylabel('atoms/b-cm')

    # Box B (total atoms)
    axes[1].plot(pd_B['time [d]'], pd_B['total'], label='ORIGEN', color='tab:blue',
                 **_line_style(len(pd_B), marker='o'))
    if 'analytic' in pd_B.columns:
        axes[1].plot(pd_B['time [d]'], pd_B['analytic'], '--', label='Analytic', color='tab:orange',
                     **_line_style(len(pd_B), marker='s', hollow=True))
    axes[1].set_title(f'Box B ({isotope})')
    axes[1].set_ylabel('total atoms')

    # Box C (total atoms)
    axes[2].plot(pd_C['time [d]'], pd_C['total'], label='ORIGEN', color='tab:blue',
                 **_line_style(len(pd_C), marker='o'))
    if 'analytic' in pd_C.columns:
        axes[2].plot(pd_C['time [d]'], pd_C['analytic'], '--', label='Analytic', color='tab:orange',
                     **_line_style(len(pd_C), marker='s', hollow=True))
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


def _plot_isotopes(pd_A: pd.DataFrame, pd_B: pd.DataFrame, pd_C: pd.DataFrame, isotopes: list[str],
                   volume: float, out_prefix: str | None) -> None:
    # Generate ORIGEN-only plots for selected isotopes from F71 runs.
    for iso in isotopes:
        if iso not in pd_A.columns or iso not in pd_B.columns or iso not in pd_C.columns:
            continue
        pd_A_plot = pd_A.copy()
        pd_B_plot = pd_B.copy()
        pd_C_plot = pd_C.copy()
        pd_B_plot['total'] = pd_B_plot[iso] * volume
        pd_C_plot['total'] = pd_C_plot[iso] * volume
        plot_prefix = f"{out_prefix}_{iso.replace('-', '')}" if out_prefix else iso.replace('-', '')
        plot_results(pd_A_plot, pd_B_plot, pd_C_plot, iso, plot_prefix)


def _write_box_c_timeseries_json(pd_C: pd.DataFrame, volume_cm3: float, out_prefix: str | None,
                                 source_volume_cm3: float | None = None) -> str:
    # Save Box C atom densities vs time with volume metadata.
    payload: dict = {
        'volume_cm3': float(volume_cm3),
        'time_s': [float(x) for x in pd_C['time [s]'].tolist()],
        'time_d': [float(x) for x in pd_C['time [d]'].tolist()],
        'atom_densities': {},
    }
    if source_volume_cm3 is not None:
        payload['source_volume_cm3'] = float(source_volume_cm3)
    for col in pd_C.columns:
        if col in ('time [s]', 'time [d]'):
            continue
        payload['atom_densities'][col] = [float(x) for x in pd_C[col].tolist()]

    fname = _out_name(out_prefix, 'boxC_timeseries', '.json')
    with open(fname, 'w') as f:
        json5.dump(payload, f, indent=2)
    return fname


def _find_case_f71(case_dirs: list[str], f71_name: str = 'origen.f71') -> str | None:
    # Locate an ORIGEN F71 file for a given case directory (supports legacy _leak suffix).
    roots = [os.getcwd(), os.path.join(os.getcwd(), 'leaky_box_origen')]
    for case_dir in case_dirs:
        for root in roots:
            candidate = os.path.join(root, case_dir, f71_name)
            if os.path.isfile(candidate):
                return candidate
    return None


def _total_activity_from_f71(f71_path: str, position: int = -1) -> float:
    # Sum nuclide activities (Bq) at a given position from an ORIGEN F71 file.
    activity = get_burned_nuclide_data(f71_path, position, f71units='becq')
    return float(sum(activity.values()))


def _activity_timeseries_from_box_json(box_json_path: str, case_prefix: str,
                                       position: int = -1) -> pd.DataFrame:
    # Build total activity vs time from per-case F71 outputs referenced by a box json5 file.
    with open(box_json_path, 'r') as f:
        box_data = json5.load(f)
    steps = sorted(box_data.items(), key=lambda kv: float(kv[1]['time']))
    rows: list[dict] = []
    for key, rec in steps:
        try:
            k_int = int(key)
        except ValueError:
            k_int = int(float(key))
        time_s = float(rec['time'])
        case_dirs = [
            f'_{case_prefix}_{k_int:04d}',
            f'_{case_prefix}_{k_int:04d}_leak',
        ]
        f71_path = _find_case_f71(case_dirs)
        if f71_path is None:
            print(f"Warning: no F71 found for {case_prefix} step {k_int:04d}")
            total_bq = 0.0
        else:
            total_bq = _total_activity_from_f71(f71_path, position)
        rows.append({
            'time [s]': time_s,
            'time [d]': time_s / float(24 * 60 * 60),
            'activity [Bq]': total_bq,
        })
    return pd.DataFrame(rows)


def _write_activity_csv(pd_B_act: pd.DataFrame, pd_C_act: pd.DataFrame, out_prefix: str | None) -> str:
    # Save per-box and combined activity tables.
    fname_b = _out_name(out_prefix, 'boxB_activity', '.csv')
    fname_c = _out_name(out_prefix, 'boxC_activity', '.csv')
    pd_B_act.to_csv(fname_b, index=False)
    pd_C_act.to_csv(fname_c, index=False)

    merged = pd.merge_asof(
        pd_B_act.sort_values('time [s]'),
        pd_C_act.sort_values('time [s]'),
        on='time [s]',
        suffixes=('_B', '_C'),
    )
    merged.rename(columns={
        'time [d]_B': 'time [d]',
        'activity [Bq]_B': 'boxB_activity [Bq]',
        'activity [Bq]_C': 'boxC_activity [Bq]',
    }, inplace=True)
    fname = _out_name(out_prefix, 'leaky_boxes_activity', '.csv')
    merged.to_csv(fname, index=False)
    return fname


def _activity_timeseries_per_nuclide_from_box_json(box_json_path: str, case_prefix: str,
                                                   position: int = -1) -> pd.DataFrame:
    # Build per-nuclide activity (Bq) vs time from per-case F71 outputs.
    with open(box_json_path, 'r') as f:
        box_data = json5.load(f)
    steps = sorted(box_data.items(), key=lambda kv: float(kv[1]['time']))
    rows: list[dict] = []
    for key, rec in steps:
        try:
            k_int = int(key)
        except ValueError:
            k_int = int(float(key))
        time_s = float(rec['time'])
        case_dirs = [
            f'_{case_prefix}_{k_int:04d}',
            f'_{case_prefix}_{k_int:04d}_leak',
        ]
        f71_path = _find_case_f71(case_dirs)
        activity = {}
        if f71_path is None:
            print(f"Warning: no F71 found for {case_prefix} step {k_int:04d}")
        else:
            activity = get_burned_nuclide_data(f71_path, position, f71units='becq')
        row = {'time [s]': time_s, 'time [d]': time_s / float(24 * 60 * 60)}
        row.update(activity)
        rows.append(row)
    return pd.DataFrame(rows)


def _load_dcf_csv(path: str) -> dict[str, float]:
    # Load dose coefficients (DCF) from CSV with columns: nuclide, dcf_sv_bq or dcf_rem_bq.
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if 'nuclide' not in cols:
        raise ValueError("DCF CSV must include 'nuclide' column")
    if 'dcf_sv_bq' in cols:
        dcf_col = cols['dcf_sv_bq']
        scale = 1.0
    elif 'dcf_rem_bq' in cols:
        dcf_col = cols['dcf_rem_bq']
        scale = 0.01  # rem -> Sv
    else:
        raise ValueError("DCF CSV must include 'dcf_sv_bq' or 'dcf_rem_bq'")
    dcf_map: dict[str, float] = {}
    max_dcf = 1e-2  # Sv/Bq sanity cap to discard OCR or parsing artifacts.
    for _, row in df.iterrows():
        nuclide = str(row[cols['nuclide']]).strip().lower()
        try:
            val = float(row[dcf_col]) * scale
        except (TypeError, ValueError):
            continue
        if not math.isfinite(val) or val <= 0.0 or val > max_dcf:
            continue
        dcf_map[nuclide] = val
    return dcf_map


def _load_dcf_immersion_csv(path: str) -> dict[str, float]:
    # Load immersion (cloudshine) dose coefficients from CSV with columns:
    # nuclide, dcf_sv_per_bq_m3_day (Sv/day per Bq/m^3).
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if 'nuclide' not in cols or 'dcf_sv_per_bq_m3_day' not in cols:
        raise ValueError("Immersion DCF CSV must include 'nuclide' and 'dcf_sv_per_bq_m3_day' columns")
    dcf_map: dict[str, float] = {}
    max_dcf = 1e-2  # Sv/day per (Bq/m^3) sanity cap.
    for _, row in df.iterrows():
        nuclide = str(row[cols['nuclide']]).strip().lower()
        try:
            val = float(row[cols['dcf_sv_per_bq_m3_day']])
        except (TypeError, ValueError):
            continue
        if not math.isfinite(val) or val <= 0.0 or val > max_dcf:
            continue
        dcf_map[nuclide] = val
    return dcf_map


def _chi_q_at_time(time_s: float, chi_q_schedule: float | list[tuple[float, float]]) -> float:
    # Resolve chi/Q at a given time. Schedule entries are (t_end_s, chi_q).
    if isinstance(chi_q_schedule, (int, float)):
        return float(chi_q_schedule)
    for t_end_s, chi_q in chi_q_schedule:
        if time_s <= t_end_s:
            return float(chi_q)
    return float(chi_q_schedule[-1][1])


def chi_q_schedule_rg145_400m() -> list[tuple[float, float]]:
    # RG 1.145 (PAVAN) 95th-percentile ground-level chi/Q at 400 m.
    # Values: 0-2h, 0-8h, 8-24h, 1-4d, 4-30d, annual.
    # Convert to piecewise: 0-2h, 2-8h, 8-24h, 1-4d, 4-30d, >30d.
    chi_0_2 = 3.35e-3
    chi_0_8 = 1.91e-3
    chi_8_24 = 1.49e-3
    chi_1_4d = 8.75e-4
    chi_4_30d = 4.08e-4
    chi_annual = 1.60e-4
    # Back-calculate 2-8h average from 0-8h and 0-2h.
    chi_2_8 = (chi_0_8 * 8.0 - chi_0_2 * 2.0) / 6.0
    return [
        (2.0 * 3600.0, chi_0_2),
        (8.0 * 3600.0, chi_2_8),
        (24.0 * 3600.0, chi_8_24),
        (4.0 * 24.0 * 3600.0, chi_1_4d),
        (30.0 * 24.0 * 3600.0, chi_4_30d),
        (365.0 * 24.0 * 3600.0, chi_annual),
    ]


def compute_inhalation_dose_timeseries(activity_df: pd.DataFrame, dcf_map: dict[str, float],
                                       chi_q_s_m3: float, breathing_rate_m3_s: float) -> pd.DataFrame:
    """Compute inhalation dose from a nuclide activity time series.

    Assumes activity_df columns (excluding time) represent release rates in Bq/s for each nuclide.
    Dose rate (Sv/s) = chi/Q [s/m^3] * breathing_rate [m^3/s] * sum_i( Q_i * DCF_i [Sv/Bq] ).
    """
    if 'time [s]' not in activity_df.columns:
        raise ValueError("activity_df must include 'time [s]' column")
    nuclide_cols = [c for c in activity_df.columns if c not in ('time [s]', 'time [d]')]
    dose_rate = []
    for _, row in activity_df.iterrows():
        total = 0.0
        for nuclide in nuclide_cols:
            if nuclide in dcf_map:
                try:
                    total += float(row[nuclide]) * dcf_map[nuclide]
                except (TypeError, ValueError):
                    continue
        dose_rate.append(chi_q_s_m3 * breathing_rate_m3_s * total)
    out = pd.DataFrame({
        'time [s]': activity_df['time [s]'],
        'time [d]': activity_df['time [s]'] / float(24 * 60 * 60),
        'dose_rate [Sv/s]': dose_rate,
    })
    # Integrate cumulative dose (Sv) using trapezoidal rule.
    dose = [0.0]
    for i in range(1, len(out)):
        dt = float(out.loc[i, 'time [s]'] - out.loc[i - 1, 'time [s]'])
        avg_rate = 0.5 * (out.loc[i, 'dose_rate [Sv/s]'] + out.loc[i - 1, 'dose_rate [Sv/s]'])
        dose.append(dose[-1] + avg_rate * dt)
    out['dose [Sv]'] = dose
    out['dose [rem]'] = out['dose [Sv]'] * 100.0
    return out


def compute_dose_timeseries(activity_df: pd.DataFrame, dcf_inhal_map: dict[str, float],
                            chi_q_schedule: float | list[tuple[float, float]],
                            breathing_rate_m3_s: float,
                            dcf_immersion_map: dict[str, float] | None = None) -> pd.DataFrame:
    """Compute inhalation + immersion (cloudshine) dose from a nuclide activity time series.

    Assumes activity_df columns (excluding time) represent release rates in Bq/s for each nuclide.
    Inhalation: dose_rate = chi/Q * breathing_rate * sum_i( Q_i * DCF_i [Sv/Bq] ).
    Immersion: dose_rate = chi/Q * sum_i( Q_i * DCF_imm_i [Sv/day per (Bq/m^3)] ) / 86400.
    """
    if 'time [s]' not in activity_df.columns:
        raise ValueError("activity_df must include 'time [s]' column")
    nuclide_cols = [c for c in activity_df.columns if c not in ('time [s]', 'time [d]')]
    dcf_immersion_map = dcf_immersion_map or {}

    dose_rate_inh = []
    dose_rate_imm = []
    chi_q_vals = []
    for _, row in activity_df.iterrows():
        time_s = float(row['time [s]'])
        chi_q = _chi_q_at_time(time_s, chi_q_schedule)
        chi_q_vals.append(chi_q)
        inhal_sum = 0.0
        imm_sum = 0.0
        for nuclide in nuclide_cols:
            try:
                val = float(row[nuclide])
            except (TypeError, ValueError):
                continue
            if math.isnan(val):
                continue
            if nuclide in dcf_inhal_map:
                inhal_sum += val * dcf_inhal_map[nuclide]
            if nuclide in dcf_immersion_map:
                imm_sum += val * dcf_immersion_map[nuclide]
        dose_rate_inh.append(chi_q * breathing_rate_m3_s * inhal_sum)
        dose_rate_imm.append(chi_q * (imm_sum / 86400.0))

    out = pd.DataFrame({
        'time [s]': activity_df['time [s]'],
        'time [d]': activity_df['time [s]'] / float(24 * 60 * 60),
        'chi/Q [s/m^3]': chi_q_vals,
        'dose_rate_inhalation [Sv/s]': dose_rate_inh,
        'dose_rate_immersion [Sv/s]': dose_rate_imm,
    })
    out['dose_rate [Sv/s]'] = out['dose_rate_inhalation [Sv/s]'] + out['dose_rate_immersion [Sv/s]']

    # Integrate cumulative doses (Sv) using trapezoidal rule.
    def integrate(col: str) -> list[float]:
        vals = [0.0]
        for i in range(1, len(out)):
            dt = float(out.loc[i, 'time [s]'] - out.loc[i - 1, 'time [s]'])
            avg_rate = 0.5 * (out.loc[i, col] + out.loc[i - 1, col])
            vals.append(vals[-1] + avg_rate * dt)
        return vals

    out['dose_inhalation [Sv]'] = integrate('dose_rate_inhalation [Sv/s]')
    out['dose_immersion [Sv]'] = integrate('dose_rate_immersion [Sv/s]')
    out['dose [Sv]'] = integrate('dose_rate [Sv/s]')
    out['dose [rem]'] = out['dose [Sv]'] * 100.0
    return out


def _inventory_to_release_rate(activity_df: pd.DataFrame, removal_rate_s: float) -> pd.DataFrame:
    """Convert nuclide inventory activity (Bq) into release rate (Bq/s)."""
    if removal_rate_s <= 0.0:
        raise ValueError("removal_rate_s must be positive when converting inventory to release rate")
    out = activity_df.copy()
    nuclide_cols = [c for c in out.columns if c not in ('time [s]', 'time [d]')]
    for col in nuclide_cols:
        out[col] = out[col] * removal_rate_s
    return out


def plot_dose_timeseries(pd_dose: pd.DataFrame, out_prefix: str | None, logy: bool = True) -> str:
    # Plot cumulative dose vs time.
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    if 'dose_inhalation [Sv]' in pd_dose.columns:
        ax.plot(pd_dose['time [d]'], pd_dose['dose_inhalation [Sv]'],
                label='Inhalation [Sv]', color='tab:blue', **_line_style(len(pd_dose), marker='o'))
    if 'dose_immersion [Sv]' in pd_dose.columns:
        ax.plot(pd_dose['time [d]'], pd_dose['dose_immersion [Sv]'],
                label='Immersion [Sv]', color='tab:orange', **_line_style(len(pd_dose), marker='s', hollow=True))
    ax.plot(pd_dose['time [d]'], pd_dose['dose [Sv]'], label='Total [Sv]', color='tab:red',
            **_line_style(len(pd_dose), marker='^'))
    ax.set_xlabel('time [days]')
    ax.set_ylabel('dose [Sv]')
    ax.set_title('Cumulative Dose vs Time')
    if logy:
        ax.set_yscale('log')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fname = _out_name(out_prefix, 'leaky_boxes_dose', '.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def plot_dose_comparison(csv_a: str, csv_b: str, label_a: str, label_b: str,
                         out_prefix: str | None, title: str,
                         logy: bool = True) -> str:
    # Plot cumulative dose comparison between two datasets.
    import matplotlib.pyplot as plt

    pd_a = pd.read_csv(csv_a)
    pd_b = pd.read_csv(csv_b)

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(pd_a['time [d]'], pd_a['dose [Sv]'], label=label_a, color='tab:blue',
            **_line_style(len(pd_a), marker='o'))
    ax.plot(pd_b['time [d]'], pd_b['dose [Sv]'], label=label_b, color='tab:orange',
            **_line_style(len(pd_b), marker='s', hollow=True))
    ax.set_xlabel('time [days]')
    ax.set_ylabel('dose [Sv]')
    ax.set_title(title)
    if logy:
        ax.set_yscale('log')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fname = _out_name(out_prefix, 'leaky_boxes_dose_compare', '.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def plot_dose_rate_comparison(csv_a: str, csv_b: str, label_a: str, label_b: str,
                              out_prefix: str | None, title: str,
                              logy: bool = True) -> str:
    # Plot dose rate comparison between two datasets.
    import matplotlib.pyplot as plt

    pd_a = pd.read_csv(csv_a)
    pd_b = pd.read_csv(csv_b)

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(pd_a['time [d]'], pd_a['dose_rate [Sv/s]'], label=label_a, color='tab:blue',
            **_line_style(len(pd_a), marker='o'))
    ax.plot(pd_b['time [d]'], pd_b['dose_rate [Sv/s]'], label=label_b, color='tab:orange',
            **_line_style(len(pd_b), marker='s', hollow=True))
    ax.set_xlabel('time [days]')
    ax.set_ylabel('dose rate [Sv/s]')
    ax.set_title(title)
    if logy:
        ax.set_yscale('log')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fname = _out_name(out_prefix, 'leaky_boxes_dose_rate_compare', '.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def compute_dose_from_box(box_json_path: str, case_prefix: str, dcf_inhal_csv: str,
                          chi_q_schedule: float | list[tuple[float, float]],
                          breathing_rate_m3_s: float, dcf_immersion_csv: str | None = None,
                          out_prefix: str | None = None,
                          activity_representation: str = 'inventory_bq',
                          removal_rate_s: float | None = None) -> pd.DataFrame:
    """Compute and save inhalation+immersion dose from per-case F71 activities.

    Args:
        activity_representation: `inventory_bq` or `release_rate_bq_s`.
            - `inventory_bq`: activity in each nuclide column is Bq and is converted using `removal_rate_s`.
            - `release_rate_bq_s`: activity in each nuclide column is already Bq/s.
        removal_rate_s: Required when `activity_representation='inventory_bq'`.
    """
    dcf_inhal = _load_dcf_csv(dcf_inhal_csv)
    dcf_imm = _load_dcf_immersion_csv(dcf_immersion_csv) if dcf_immersion_csv else {}
    activity_df = _activity_timeseries_per_nuclide_from_box_json(box_json_path, case_prefix)
    if activity_representation == 'inventory_bq':
        if removal_rate_s is None:
            raise ValueError("removal_rate_s must be provided when activity_representation='inventory_bq'")
        activity_df = _inventory_to_release_rate(activity_df, removal_rate_s)
    elif activity_representation != 'release_rate_bq_s':
        raise ValueError("activity_representation must be 'inventory_bq' or 'release_rate_bq_s'")
    pd_dose = compute_dose_timeseries(activity_df, dcf_inhal, chi_q_schedule,
                                      breathing_rate_m3_s, dcf_imm)
    fname = _out_name(out_prefix, 'leaky_boxes_dose', '.csv')
    pd_dose.to_csv(fname, index=False)
    plot_dose_timeseries(pd_dose, out_prefix)
    return pd_dose

def plot_total_activity(pd_B_act: pd.DataFrame, pd_C_act: pd.DataFrame,
                        out_prefix: str | None, logy: bool = True) -> str:
    # Plot total activity in Box B and Box C vs time.
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(pd_B_act['time [d]'], pd_B_act['activity [Bq]'], label='Box B', color='tab:blue',
            **_line_style(len(pd_B_act), marker='o'))
    ax.plot(pd_C_act['time [d]'], pd_C_act['activity [Bq]'], label='Box C', color='tab:orange',
            **_line_style(len(pd_C_act), marker='s', hollow=True))
    ax.set_xlabel('time [days]')
    ax.set_ylabel('activity [Bq]')
    ax.set_title('Total Activity vs Time')
    if logy:
        ax.set_yscale('log')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fname = _out_name(out_prefix, 'leaky_boxes_activity', '.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def compute_inhalation_dose_from_box(box_json_path: str, case_prefix: str, dcf_csv: str,
                                     chi_q_s_m3: float, breathing_rate_m3_s: float,
                                     out_prefix: str | None = None,
                                     activity_representation: str = 'inventory_bq',
                                     removal_rate_s: float | None = None) -> pd.DataFrame:
    """Compute and save inhalation dose from per-case F71 activities.

    This assumes activity columns represent release rates in Bq/s. For Box C, this is a
    conservative bounding assumption if Box C is treated as an instantaneous release.
    """
    dcf_map = _load_dcf_csv(dcf_csv)
    activity_df = _activity_timeseries_per_nuclide_from_box_json(box_json_path, case_prefix)
    if activity_representation == 'inventory_bq':
        if removal_rate_s is None:
            raise ValueError("removal_rate_s must be provided when activity_representation='inventory_bq'")
        activity_df = _inventory_to_release_rate(activity_df, removal_rate_s)
    elif activity_representation != 'release_rate_bq_s':
        raise ValueError("activity_representation must be 'inventory_bq' or 'release_rate_bq_s'")
    pd_dose = compute_inhalation_dose_timeseries(activity_df, dcf_map, chi_q_s_m3, breathing_rate_m3_s)
    fname = _out_name(out_prefix, 'leaky_boxes_dose', '.csv')
    pd_dose.to_csv(fname, index=False)
    plot_dose_timeseries(pd_dose, out_prefix)
    return pd_dose

def _get_case_max_time(index: dict, case: str) -> float:
    # Max time for a given case in an F71 index.
    case_str = str(case)
    times = [float(v['time']) for v in index.values() if str(v.get('case')) == case_str]
    if not times:
        raise ValueError(f"No times found for case '{case}' in F71 index")
    return max(times)


def _get_case_time_from_end(index: dict, case: str, offset_from_end: int) -> float:
    # Time selection by offset from the end (e.g., pre-last = 1).
    case_str = str(case)
    times = sorted(float(v['time']) for v in index.values() if str(v.get('case')) == case_str)
    if len(times) <= offset_from_end:
        raise ValueError(f"Not enough times for case '{case}' to select offset {offset_from_end} from end")
    return times[-1 - offset_from_end]


def _setup_box_a(f71_path: str, f71_case: str, f71_time_seconds: float | None,
                 sample_mass_g: float, do_skip_calcs: bool,
                 decay_days: float, steps_per_day: float) -> DecayBoxA:
    # Configure Box A from F71 (case/time selection + initial decay settings).
    box_a = DecayBoxA(f71_path, sample_mass_g)
    if do_skip_calcs:
        box_a.skip_calculation = True
    if f71_time_seconds is None:
        # Use the pre-last position by default (second-to-last time in the case).
        f71_time_seconds = _get_case_time_from_end(box_a.BURNED_MATERIAL_F71_index, f71_case, 1)
    box_a.set_f71_pos(f71_time_seconds, case=f71_case)
    box_a.read_burned_material()
    box_a.DECAY_days = decay_days
    box_a.DECAY_steps = int(decay_days * steps_per_day) + 2
    return box_a


def _run_simulation(box_a: DecayBoxA, apply_volume_scaling: bool, do_skip_calcs: bool,
                    out_prefix: str | None, box_B_n_steps: int = 8, box_C_n_steps: int = 8
                    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float, float, float]:
    # Core leaky-box simulation (A->B->C) shared by tests and F71 runs.
    box_A_leak_rate: float = PCTperDAY
    box_B_leak_rate: float = PCTperDAY * 0.1

    volume: float = 1.0  # box_a.volume
    print(volume)

    box_A = LeakyBox()
    box_A.removal_rate = box_A_leak_rate
    box_A.setup_cases(box_a)
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

    return pd_A, pd_B, pd_C, volume, box_A_leak_rate, box_B_leak_rate


def run_test(isotope: str, lambda_decay: float | None, apply_volume_scaling: bool,
             do_skip_calcs: bool, out_prefix: str | None = None,
             f71_path: str = DEFAULT_F71_PATH, f71_case: str = DEFAULT_F71_CASE,
             f71_time_seconds: float | None = None,
             sample_mass_g: float = 1200e3, decay_days: float = 12.0,
             steps_per_day: float = 2.0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run a single isotope test and return dataframes."""
    # Baseline decay calculation setup (forces a single-isotope inventory).
    origen_triton_box_A = _setup_box_a(
        f71_path, f71_case, f71_time_seconds, sample_mass_g, do_skip_calcs, decay_days, steps_per_day
    )
    origen_triton_box_A.atom_dens = {isotope: 1.0}
    pd_A, pd_B, pd_C, volume, box_A_leak_rate, box_B_leak_rate = _run_simulation(
        origen_triton_box_A, apply_volume_scaling, do_skip_calcs, out_prefix
    )

    pd_B['total'] = pd_B[isotope] * volume
    pd_C['total'] = pd_C[isotope] * volume

    n0_density = float(origen_triton_box_A.atom_dens.get(isotope, 0.0))
    n0_total = n0_density * (origen_triton_box_A.volume if apply_volume_scaling else volume)

    _add_analytic_columns(pd_A, pd_B, pd_C, isotope, n0_density, n0_total,
                          box_A_leak_rate, box_B_leak_rate, lambda_decay)

    _write_excel(pd_A, pd_B, pd_C, out_prefix)
    _write_dataframes_json(pd_A, pd_B, pd_C, out_prefix)
    plot_results(pd_A, pd_B, pd_C, isotope, out_prefix)
    return pd_A, pd_B, pd_C


def run_from_f71(apply_volume_scaling: bool, do_skip_calcs: bool, out_prefix: str | None = None,
                 f71_path: str = DEFAULT_F71_PATH, f71_case: str = DEFAULT_F71_CASE,
                 f71_time_seconds: float | None = None, sample_mass_g: float = 1200e3,
                 plot_isotopes: list[str] | None = None, decay_days: float = 12.0,
                 steps_per_day: float = 2.0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run a simulation using the actual F71 nuclide composition (no analytic comparison)."""
    # Uses full F71 composition; analytic columns are not added.
    origen_triton_box_A = _setup_box_a(
        f71_path, f71_case, f71_time_seconds, sample_mass_g, do_skip_calcs, decay_days, steps_per_day
    )
    pd_A, pd_B, pd_C, volume, _, _ = _run_simulation(
        origen_triton_box_A, apply_volume_scaling, do_skip_calcs, out_prefix
    )

    _write_excel(pd_A, pd_B, pd_C, out_prefix)
    _write_dataframes_json(pd_A, pd_B, pd_C, out_prefix)
    _write_box_c_timeseries_json(pd_C, volume, out_prefix, source_volume_cm3=origen_triton_box_A.volume)

    # Total activity (Bq) for Box B and Box C from per-case F71 outputs.
    try:
        boxB_json = _out_name(out_prefix, 'boxB', '.json5')
        boxC_json = _out_name(out_prefix, 'boxC', '.json5')
        pd_B_act = _activity_timeseries_from_box_json(boxB_json, 'box_B')
        pd_C_act = _activity_timeseries_from_box_json(boxC_json, 'box_C')
        _write_activity_csv(pd_B_act, pd_C_act, out_prefix)
        plot_total_activity(pd_B_act, pd_C_act, out_prefix)
    except (FileNotFoundError, RuntimeError, ValueError, OSError) as exc:
        print(f"Warning: activity plot generation failed: {exc}")

    if plot_isotopes:
        _plot_isotopes(pd_A, pd_B, pd_C, plot_isotopes, volume, out_prefix)

    return pd_A, pd_B, pd_C


def run_all_tests(apply_volume_scaling: bool, do_skip_calcs: bool,
                  decay_days: float = 12.0, steps_per_day: float = 2.0):
    # Convenience wrapper for Xe-136 and Xe-135 test cases.
    tests = {
        'xe-136': None,
        'xe-135': 2.106574217602 * 1e-5,  # Xe-135 decay constant [1/s]
    }
    for iso, lam in tests.items():
        prefix = iso.replace('-', '')
        run_test(iso, lam, apply_volume_scaling, do_skip_calcs, out_prefix=prefix,
                 decay_days=decay_days, steps_per_day=steps_per_day)


def main():
    """Default entrypoint: run Xe-136 test."""
    run_test('xe-136', None, apply_volume_scaling, do_skip_calcs, out_prefix='xe136')


if __name__ == "__main__":
    # pass
    main()
