#!/bin/env python3
"""
Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

Box A --> Box B --> Box C (outside)

Box A leaks into box B, and box B leaks into box C. Box C accumulates the release to the
environment and has no removal. Only the elements in LEAKED_ELEMENTS are removed from a box, all at
one rate by default; a box can also be given one removal rate per element (see _element_removal_rates).

100% per day  =  1/24/60/60 per second    = 0.0000115740740741
   1% per day  =  1e-2/24/60/60 per second = 0.000000115740740741
 0.1% per day  =  1e-3/24/60/60 per second = 0.0000000115740740741
"""
import numbers
import os
import re
import warnings
from collections.abc import Mapping
from contextlib import contextmanager
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from bisect import bisect_left
import json5
import math
from sample_decay_dose.SampleDose import NOW
from sample_decay_dose.utils import nicely_print_atom_dens, get_f71_positions_index, get_burned_nuclide_atom_dens, \
    get_burned_nuclide_data, get_burned_material_total_mass_dens, get_last_position_for_case, run_scale, \
    atom_dens_for_origen

PCTperDAY: float = 0.000000115740740741
DEFAULT_F71_PATH: str = os.getenv(
    'LEAKYBOX_F71_PATH',
    os.path.expanduser('~/0.03/20-burn-MHA/mha-4.5-a4/msrr.f71'),
)
DEFAULT_F71_CASE: str = os.getenv('LEAKYBOX_F71_CASE', '20')
LEAKY_BOX_ORIGEN_DIR: Path = Path(__file__).resolve().parent
LEAKY_BOX_DATA_DIR: Path = LEAKY_BOX_ORIGEN_DIR / 'data'
# Elements that leak out of a box. Each one gets an ORIGEN `removal { rate=... ele=[x] }` entry;
# all other elements (including decay daughters such as Rb, Cs, Sr, Ba) are retained in the box.
LEAKED_ELEMENTS: tuple[str, ...] = ('h', 'he', 'o', 'n', 'i', 'br', 'ar', 'kr', 'xe')
# Fail-closed default: raise when more than 5% of the release rate has no DCF.
# Pass max_missing_dcf_fraction=None to warn instead of raising.
DEFAULT_MAX_MISSING_DCF_FRACTION: float = 0.05
# ICRP-72 inhalation table for paper runs: maximum adult DCF across types F/M/S per nuclide
# (journal method), the upper bound over particulate forms. The Type-F table
# (dcf_icrp72_inhalation_adult.csv) is lower for insoluble forms. Gases and vapours (Table A.3) sit
# outside it: elemental iodine vapour is 2.7 times the Type F value for I-131.
DEFAULT_ICRP72_INHALATION_CSV: str = 'dcf_icrp72_inhalation_adult_max.csv'
# ICRP-72 Table A.3 adult coefficients for elemental iodine vapour (I-120 to I-135), the iodine-vapour sensitivity.
ICRP72_IODINE_VAPOUR_CSV: str = 'dcf_icrp72_iodine_vapour_adult.csv'
# ICRP-72 adult immersion (cloudshine) table, the default for the noble-gas dose.
DEFAULT_ICRP72_IMMERSION_CSV: str = 'dcf_icrp72_immersion_adult.csv'
# Adult breathing rate of the journal method: 22.2 m^3/d (ICRP-66 adult male daily average, the ICRP-89
# reference value) in m^3/s.
BREATHING_RATE_M3_S: float = 22.2 / 86400.0


def _element_removal_rates(removal_rate: float | Mapping[str, float],
                           default_elements: tuple[str, ...] = LEAKED_ELEMENTS) -> dict[str, float]:
    """ Normalize a box leak rate into an {element: rate [1/s]} dict.

    A scalar applies to every element in default_elements, which is the journal method. A mapping
    carries one rate per element, e.g. the PSAR MHA enclosure leak rates for Ar/Kr/Xe at 0.05 %/d
    and the other released elements at 0.005 %/d. Keys are lower-cased to match nuclide column
    names; values must be finite and non-negative (zero retains that element in the box).
    """
    if isinstance(removal_rate, numbers.Real):
        rates: dict[str, float] = {ele: float(removal_rate) for ele in default_elements}
    elif isinstance(removal_rate, Mapping) and removal_rate:
        rates = {str(ele).strip().lower(): float(rate) for ele, rate in removal_rate.items()}
        if len(rates) != len(removal_rate):
            raise ValueError(f"removal_rate keys collide after lower-casing: {sorted(removal_rate)}")
    else:
        raise ValueError("removal_rate must be a non-negative number or a non-empty {element: rate [1/s]} mapping")
    for ele, rate in rates.items():
        if not math.isfinite(rate) or rate < 0.0:
            raise ValueError(f"removal rate of '{ele}' must be a finite non-negative number per second, got {rate!r}")
    return rates


def _run_removal_rates(box_A_leak_rate: float | Mapping[str, float], box_B_leak_rate: float | Mapping[str, float],
                       box_C_leak_rate: float | Mapping[str, float]
                       ) -> tuple[dict[str, float], dict[str, float], dict[str, float], list[str]]:
    """ Positive {element: rate [1/s]} removal rates that ORIGEN applies to boxes A, B and C, and the run's
    released elements.

    A scalar rate applies to every element the run releases. That set is the elements with a positive rate in
    the per-element rates of any box, or LEAKED_ELEMENTS when every box takes one rate (the journal method).
    A scalar box among per-element boxes therefore removes the same elements as its neighbours, so an element
    that leaks into a box also leaks out of it. An empty dict leaves the box without removal (box C as a
    collector).
    """
    given = (box_A_leak_rate, box_B_leak_rate, box_C_leak_rate)
    run_elements: tuple[str, ...] = tuple(dict.fromkeys(
        ele for rate in given if isinstance(rate, Mapping)
        for ele, value in _element_removal_rates(rate).items() if value > 0.0)) or LEAKED_ELEMENTS
    rates_A, rates_B, rates_C = ({ele: value for ele, value in _element_removal_rates(rate, run_elements).items()
                                  if value > 0.0} for rate in given)
    return rates_A, rates_B, rates_C, list(dict.fromkeys(rates_A | rates_B | rates_C))


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
        self.case_path: str = ''  # Absolute path of case_dir, set when run_decay_sample() runs

    def _enter_case_dir(self) -> str:
        """ Create case_dir under the current cwd, chdir into it, and return the caller's cwd """
        prev_cwd: str = os.getcwd()
        self.case_path = os.path.abspath(self.case_dir)
        os.makedirs(self.case_path, exist_ok=True)
        os.chdir(self.case_path)
        return prev_cwd


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

    def set_f71_pos(self, t: float = 5184000.0, case: str | int = '1'):
        """ Sets the F71 position closest in time to t for a case (ties go to the earlier slot) """
        case_str: str = str(case)
        pos_times = sorted(
            ((k, float(v['time'])) for k, v in self.BURNED_MATERIAL_F71_index.items()
             if str(v['case']) == case_str),
            key=lambda x: (x[1], x[0])
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
            # bisect_left gives the first slot at or after t; the previous slot may be closer.
            pos_idx: int = bisect_left(times, t)
            if pos_idx > 0 and (t - times[pos_idx - 1]) <= (times[pos_idx] - t):
                pos_idx -= 1
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
        prev_cwd: str = self._enter_case_dir()
        try:
            if not self.skip_calculation:
                with open(self.ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
                    f.write(atom_dens_for_origen(self.atom_dens))

                with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
                    f.write(self.origen_deck())

                if self.debug > 0:
                    print(f'ORIGEN: decaying sample for {self.DECAY_days} days')
                    print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
                if not run_scale(self.ORIGEN_input_file_name):
                    raise RuntimeError(f"SCALE run failed for {self.case_dir}/{self.ORIGEN_input_file_name}")
            else:
                if not os.path.isfile(self.F71_file_name):
                    error_text: str = f'Skip SCALE flag set, output file {self.case_dir}/{self.F71_file_name} is not found!'
                    raise ValueError(error_text)

            self.final_atom_dens = get_burned_nuclide_atom_dens(self.F71_file_name,
                                                                self.DECAY_steps)
        finally:
            os.chdir(prev_cwd)  # restore the caller's cwd, not the construction-time one
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
        self.volume: float = _volume  # Sample volume [cm3]
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
        prev_cwd: str = self._enter_case_dir()
        try:
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
                if not run_scale(self.ORIGEN_input_file_name):
                    raise RuntimeError(f"SCALE run failed for {self.case_dir}/{self.ORIGEN_input_file_name}")
            else:
                if not os.path.isfile(self.F71_file_name):
                    error_text: str = f'Skip SCALE flag set, output file {self.case_dir}/{self.F71_file_name} is not found!'
                    raise ValueError(error_text)

            self.final_atom_dens = get_burned_nuclide_atom_dens(self.F71_file_name,
                                                                self.DECAY_steps)
        finally:
            os.chdir(prev_cwd)  # restore the caller's cwd, not the construction-time one
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

        # DECAY_steps F71 positions: the start state plus DECAY_steps - 1 evenly spaced times after it. The list
        # begins one step after start, because ORIGEN stores the start state itself as position 1.
        time_interp_steps = self.DECAY_steps - 3
        if time_interp_steps < 1:
            raise ValueError("Too few time steps")
        first_time = self.DECAY_start_seconds + (self.DECAY_end_seconds - self.DECAY_start_seconds) / (self.DECAY_steps - 1)

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
        t=[{time_interp_steps}I {first_time} {self.DECAY_end_seconds}]
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
        self.debug: int = 0  # Debugging flag; > 2 prints every per-nuclide leak rate
        self.nuc_leak = None
        self.decay_leaks: (None, DecayBoxA) = None
        self.leak_rates: dict = {}  # Leak rate calculated as \epsilon_i^A N_i^A (t)
        self.removal_rate: float | dict[str, float] = PCTperDAY  # \epsilon_i^A, one rate or one per element
        self.nuclide_removal_rates: dict = {}
        self.adens: dict = {}  # box_A_adens

    def setup_cases(self, origen_case: DecayBoxA):
        """ Setup cases based on a base_case """
        self.nuclide_removal_rates = _element_removal_rates(self.removal_rate)
        self.decay_leaks = origen_case
        if not self.decay_leaks.case_dir.endswith('_leak'):  # Idempotent on repeated calls
            self.decay_leaks.case_dir += '_leak'
        self.decay_leaks.nuclide_removal_rates = self.nuclide_removal_rates
        self.nuc_leak: list = [s.lower() for s in self.nuclide_removal_rates.keys()]

    def run_case(self):
        """ Run cases and get results """
        self.decay_leaks.run_decay_sample()

    def _leak_f71_path(self) -> str:
        """ Absolute path of the leak-case F71 file.
        Resolved in the parent process: reused joblib workers keep the cwd they started with,
        so relative paths inside workers can point to a previous run directory.
        """
        case_path: str = self.decay_leaks.case_path or os.path.abspath(self.decay_leaks.case_dir)
        return os.path.join(case_path, self.decay_leaks.F71_file_name)

    def _parallel_leak_calc(self, i, vals, f71_leak_file: str | None = None) -> [int, float, dict, dict]:
        t = float(vals['time'])
        if f71_leak_file is None:
            f71_leak_file = self._leak_f71_path()
        adens_tot: dict = get_burned_nuclide_atom_dens(f71_leak_file, i)  # All nuclides
        step_leak_rate: dict = {k: adens_tot[k] for k in adens_tot.keys() if re.sub('-.*', '', k) in self.nuc_leak}
        for k in step_leak_rate.keys():
            step_leak_rate[k] *= self.nuclide_removal_rates[re.sub('-.*', '', k)]  # * self.decay_leaks.volume
            if self.debug > 2:
                print("LEAK_RATE: ", k, step_leak_rate[k], adens_tot[k])
        return i, t, adens_tot, step_leak_rate

    def get_leak_rate_parallel(self):
        """ Calculate Box B ingress rate = Box A leak rate"""
        from joblib import Parallel, delayed, cpu_count
        f71_leak_file: str = self._leak_f71_path()
        times: dict = get_f71_positions_index(f71_leak_file)
        res_list = Parallel(n_jobs=cpu_count())(delayed(self._parallel_leak_calc)(ii, vv, f71_leak_file)
                                                for ii, vv in times.items())
        for k, res in enumerate(res_list):
            (i, t, adens_tot, step_leak_rate) = res
            self.leak_rates[i]: dict = {}
            self.leak_rates[i]['time'] = t
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


def _resolve_run_dir(run_dir: str | None = None, run_date: str | None = None) -> str:
    # Default run artifacts location: leaky_box_origen/run_YYYY-MM-DD.
    if run_dir:
        resolved = Path(run_dir).expanduser()
        if not resolved.is_absolute():
            resolved = (Path.cwd() / resolved).resolve()
        else:
            resolved = resolved.resolve()
    else:
        if run_date is None:
            run_date = date.today().isoformat()
        resolved = LEAKY_BOX_ORIGEN_DIR / f'run_{run_date}'
    resolved.mkdir(parents=True, exist_ok=True)
    return str(resolved)


def _resolve_input_file(path: str, fallback_dir: Path | None = None) -> str:
    # Resolve inputs from cwd first, then optionally from a project data directory by basename.
    resolved = Path(path).expanduser()
    if resolved.exists():
        return str(resolved.resolve())
    if not resolved.is_absolute() and fallback_dir is not None:
        candidate = (fallback_dir / resolved.name).resolve()
        if candidate.exists():
            return str(candidate)
    return str(resolved)


def _resolve_run_artifact(path: str) -> str:
    # Resolve artifacts from cwd first, then today's run directory, then the newest run directory.
    resolved = Path(path).expanduser()
    if resolved.exists():
        return str(resolved.resolve())
    if resolved.is_absolute():
        return str(resolved)

    run_candidates = [LEAKY_BOX_ORIGEN_DIR / f'run_{date.today().isoformat()}']
    run_candidates.extend(sorted(LEAKY_BOX_ORIGEN_DIR.glob('run_*'), reverse=True))
    for run_dir in run_candidates:
        if run_dir.is_dir():
            candidate = (run_dir / resolved).resolve()
            if candidate.exists():
                return str(candidate)
    return str(resolved)


@contextmanager
def _in_directory(path: str):
    # Execute a block in a target directory and restore cwd on exit.
    prev_cwd = os.getcwd()
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


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


def _find_case_f71(case_dirs: list[str], f71_name: str = 'origen.f71',
                   search_roots: list[str] | None = None) -> str | None:
    # Locate an ORIGEN F71 file for a given case directory (supports legacy _leak suffix).
    roots = search_roots or [os.getcwd(), os.path.join(os.getcwd(), 'leaky_box_origen')]
    for case_dir in case_dirs:
        for root in roots:
            candidate = os.path.join(root, case_dir, f71_name)
            if os.path.isfile(candidate):
                return candidate
    return None


def _step_activity_from_f71(f71_path: str, position: int | None = None) -> dict:
    # Nuclide activities (Bq) at a position of a box B/C step F71. None reads the end of the step: each step deck
    # runs one ORIGEN case, and its last position holds the final_atom_dens state recorded in the box json.
    if position is None:
        position = get_last_position_for_case(f71_path)
    return get_burned_nuclide_data(f71_path, position, f71units='becq')


def _total_activity_from_f71(f71_path: str, position: int | None = None) -> float:
    # Sum nuclide activities (Bq) at a given position from an ORIGEN F71 file.
    return float(sum(_step_activity_from_f71(f71_path, position).values()))


def _case_dir_name(box: str, run_prefix: str | None, step: int | None = None) -> str:
    # ORIGEN case directory name for a box ('box_A', 'box_B', 'box_C').
    # Runs with an output prefix (e.g. the xe135 and xe136 tests) get their own case directories,
    # so one run cannot overwrite the F71 files of another run in the same run directory.
    name = f'_{box}_{run_prefix}' if run_prefix else f'_{box}'
    if step is not None:
        name += f'_{step:04d}'
    return name


def _run_prefix_from_box_json(box_json_path: str) -> str | None:
    # Infer the run prefix from a box json5 name written by _out_name, e.g. boxB_xe135.json5 -> xe135.
    name = Path(box_json_path).name
    stem = name[:-len('.json5')] if name.endswith('.json5') else Path(box_json_path).stem
    if '_' in stem:
        return stem.split('_', 1)[1] or None
    return None


def _step_case_dirs(case_prefix: str, step: int, run_prefix: str | None) -> list[str]:
    # Candidate case directories for one step: prefixed names first, then legacy unprefixed names.
    case_dirs: list[str] = []
    if run_prefix:
        base = _case_dir_name(case_prefix, run_prefix, step)
        case_dirs += [base, f'{base}_leak']
    base = _case_dir_name(case_prefix, None, step)
    case_dirs += [base, f'{base}_leak']
    return case_dirs


def _activity_timeseries_from_box_json(box_json_path: str, case_prefix: str,
                                       position: int | None = None, run_prefix: str | None = None) -> pd.DataFrame:
    # Build total activity vs time from per-case F71 outputs referenced by a box json5 file.
    # run_prefix defaults to the prefix encoded in the json name (boxB_xe135.json5 -> xe135).
    box_root = os.path.dirname(os.path.abspath(box_json_path))
    if run_prefix is None:
        run_prefix = _run_prefix_from_box_json(box_json_path)
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
        case_dirs = _step_case_dirs(case_prefix, k_int, run_prefix)
        f71_path = _find_case_f71(case_dirs, search_roots=[box_root])
        if f71_path is None:
            warnings.warn(f"No F71 found for {case_prefix} step {k_int:04d} (searched {case_dirs} under {box_root}); "
                          f"its activity is NaN", stacklevel=2)
            total_bq = float('nan')
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
        pd_C_act.sort_values('time [s]').drop(columns=['time [d]'], errors='ignore'),
        on='time [s]',
        suffixes=('_B', '_C'),
    )
    merged.rename(columns={
        'activity [Bq]_B': 'boxB_activity [Bq]',
        'activity [Bq]_C': 'boxC_activity [Bq]',
    }, inplace=True)
    fname = _out_name(out_prefix, 'leaky_boxes_activity', '.csv')
    merged.to_csv(fname, index=False)
    return fname


def _activity_timeseries_per_nuclide_from_box_json(box_json_path: str, case_prefix: str,
                                                   position: int | None = None,
                                                   run_prefix: str | None = None) -> pd.DataFrame:
    # Build per-nuclide activity (Bq) vs time from per-case F71 outputs.
    # run_prefix defaults to the prefix encoded in the json name (boxB_xe135.json5 -> xe135).
    box_root = os.path.dirname(os.path.abspath(box_json_path))
    if run_prefix is None:
        run_prefix = _run_prefix_from_box_json(box_json_path)
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
        case_dirs = _step_case_dirs(case_prefix, k_int, run_prefix)
        f71_path = _find_case_f71(case_dirs, search_roots=[box_root])
        if f71_path is None:  # a zero row here would integrate a hole into the cumulative dose
            raise FileNotFoundError(f"No F71 found for {case_prefix} step {k_int:04d}; searched {case_dirs} under "
                                    f"{box_root}. Rerun the leaky-box calculation or pass run_prefix.")
        activity = _step_activity_from_f71(f71_path, position)
        row = {'time [s]': time_s, 'time [d]': time_s / float(24 * 60 * 60)}
        row.update(activity)
        rows.append(row)
    return pd.DataFrame(rows)


def _dcf_missing_hint(resolved: Path) -> str:
    # Recovery hint for a missing DCF table.
    if 'fgr11' in resolved.name:
        return (f"DCF CSV not found: {resolved}. The FGR-11 tables ship with the repository under "
                f"{LEAKY_BOX_DATA_DIR}/; restore them from version control "
                f"(git checkout -- leaky_box_origen/data/) or reinstall the package. Regenerating them with "
                f"'python -m leaky_box_origen.extract_fgr11_dcf --pdf <path>' needs a local copy of the "
                f"EPA FGR-11 PDF, which is not distributed with the repository "
                f"(and pdftotext, pdftoppm and tesseract).")
    if 'fgr12' in resolved.name:
        return (f"DCF CSV not found: {resolved}. The FGR-12 table ships with the repository under "
                f"{LEAKY_BOX_DATA_DIR}/; restore it from version control "
                f"(git checkout -- leaky_box_origen/data/) or regenerate it with "
                f"'python -m leaky_box_origen.extract_fgr12_dcf --pdf <path>' from the EPA PDF "
                f"402-r-93-081.pdf (needs pdftotext).")
    if 'site-packages' in str(LEAKY_BOX_DATA_DIR):
        return (f"DCF CSV not found: {resolved}. Expected under {LEAKY_BOX_DATA_DIR}/; "
                f"reinstall the sample_decay_dose package to restore bundled data.")
    return (f"DCF CSV not found: {resolved}. Expected under {LEAKY_BOX_DATA_DIR}/; "
            f"restore it from version control (git checkout -- leaky_box_origen/data/).")


# Plausibility window for tabulated DCF values. The upper cap drops OCR or parsing artifacts
# (no inhalation coefficient approaches 1e-2 Sv/Bq). The lower bound rejects OCR exponent garbage
# such as 1e-97; the smallest published values are ~1e-15 (ICRP-72 elemental tritium,
# Sv/Bq) and ~1e-15 Sv/day per Bq/m^3 (Ar-37 immersion).
MAX_DCF: float = 1e-2
MIN_DCF: float = 1e-20


def _load_dcf_csv(path: str) -> dict[str, float]:
    # Load dose coefficients (DCF) from CSV with columns: nuclide, dcf_sv_bq or dcf_rem_bq.
    resolved = Path(_resolve_input_file(path, LEAKY_BOX_DATA_DIR))
    if not resolved.is_file():
        raise FileNotFoundError(_dcf_missing_hint(resolved))
    df = pd.read_csv(resolved)
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
    max_dcf = MAX_DCF  # Sv/Bq sanity cap to discard OCR or parsing artifacts.
    min_dcf = MIN_DCF  # Sv/Bq lower bound to discard garbled exponents (e.g. 1e-97).
    for _, row in df.iterrows():
        nuclide = str(row[cols['nuclide']]).strip().lower()
        try:
            val = float(row[dcf_col]) * scale
        except (TypeError, ValueError):
            continue
        if not math.isfinite(val) or val < min_dcf or val > max_dcf:
            continue
        dcf_map[nuclide] = val
    return dcf_map


def _load_dcf_immersion_csv(path: str) -> dict[str, float]:
    # Load immersion (cloudshine) dose coefficients from CSV with columns:
    # nuclide, dcf_sv_per_bq_m3_day (Sv/day per Bq/m^3).
    resolved = Path(_resolve_input_file(path, LEAKY_BOX_DATA_DIR))
    if not resolved.is_file():
        raise FileNotFoundError(_dcf_missing_hint(resolved))
    df = pd.read_csv(resolved)
    cols = {c.lower(): c for c in df.columns}
    if 'nuclide' not in cols or 'dcf_sv_per_bq_m3_day' not in cols:
        raise ValueError("Immersion DCF CSV must include 'nuclide' and 'dcf_sv_per_bq_m3_day' columns")
    dcf_map: dict[str, float] = {}
    max_dcf = MAX_DCF  # Sv/day per (Bq/m^3) sanity cap.
    min_dcf = MIN_DCF  # Sv/day per (Bq/m^3) lower bound to discard garbled exponents.
    for _, row in df.iterrows():
        nuclide = str(row[cols['nuclide']]).strip().lower()
        try:
            val = float(row[cols['dcf_sv_per_bq_m3_day']])
        except (TypeError, ValueError):
            continue
        if not math.isfinite(val) or val < min_dcf or val > max_dcf:
            continue
        dcf_map[nuclide] = val
    return dcf_map


def _chi_q_at_time(time_s: float, chi_q_schedule: float | list[tuple[float, float]]) -> float:
    # Resolve chi/Q at a given time. Schedule entries are (t_end_s, chi_q).
    # A scalar (any real number, including numpy int/float scalars) is a constant chi/Q.
    if isinstance(chi_q_schedule, numbers.Real):
        return float(chi_q_schedule)
    for t_end_s, chi_q in chi_q_schedule:
        if time_s <= t_end_s:
            return float(chi_q)
    return float(chi_q_schedule[-1][1])


def _chi_q_boundaries(chi_q_schedule: float | list[tuple[float, float]]) -> list[float]:
    # Times at which a piecewise-constant chi/Q schedule changes value; none for a scalar.
    if isinstance(chi_q_schedule, numbers.Real):
        return []
    return sorted(float(t_end) for t_end, _ in chi_q_schedule)


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


def compute_dose_timeseries(activity_df: pd.DataFrame, dcf_inhal_map: dict[str, float],
                            chi_q_schedule: float | list[tuple[float, float]],
                            breathing_rate_m3_s: float,
                            dcf_immersion_map: dict[str, float] | None = None,
                            max_missing_dcf_fraction: float | None = DEFAULT_MAX_MISSING_DCF_FRACTION) -> pd.DataFrame:
    """Compute inhalation + immersion (cloudshine) dose from a nuclide activity time series.

    Assumes activity_df columns (excluding time) represent release rates in Bq/s for each nuclide.
    Inhalation: dose_rate [Sv/s] = chi/Q [s/m^3] * breathing_rate [m^3/s] * sum_i( Q_i [Bq/s] * DCF_i [Sv/Bq] ).
    Immersion: dose_rate [Sv/s] = chi/Q * sum_i( Q_i * DCF_imm_i [Sv/day per (Bq/m^3)] ) / 86400.

    Rows are sorted by time, so any index works. The cumulative doses integrate chi/Q(t) times the
    DCF-weighted release sum, with the sum interpolated linearly between rows and every interval split at the
    chi/Q schedule boundaries. The result is exact for a linear release under a stepwise chi/Q, also when a
    boundary (2 h, 8 h, 24 h, ...) falls between two ORIGEN time points. The dose-rate columns hold the
    instantaneous rates at the row times.
    Column 'activity_fraction_without_dcf [-]' is the fraction of the positive release rate at each
    step carried by nuclides that have neither an inhalation nor an immersion DCF; that dose is
    missing from the totals. A warning names the largest such nuclides whenever the fraction is
    nonzero, and max_missing_dcf_fraction (default 0.05) raises ValueError when it is exceeded.
    Pass None to warn only.
    """
    if 'time [s]' not in activity_df.columns:
        raise ValueError("activity_df must include 'time [s]' column")
    activity_df = activity_df.sort_values('time [s]', kind='stable').reset_index(drop=True)
    nuclide_cols = [c for c in activity_df.columns if c not in ('time [s]', 'time [d]')]
    dcf_immersion_map = dcf_immersion_map or {}

    # Release rates [Bq/s]; non-numeric and NaN entries contribute nothing.
    vals = np.zeros((len(activity_df), len(nuclide_cols)))
    for j, nuclide in enumerate(nuclide_cols):
        col = pd.to_numeric(activity_df[nuclide], errors='coerce').to_numpy(dtype=float)
        vals[:, j] = np.where(np.isfinite(col), col, 0.0)
    dcf_inh_vec = np.array([dcf_inhal_map.get(n, 0.0) for n in nuclide_cols], dtype=float)
    dcf_imm_vec = np.array([dcf_immersion_map.get(n, 0.0) for n in nuclide_cols], dtype=float)
    no_dcf = np.array([n not in dcf_inhal_map and n not in dcf_immersion_map for n in nuclide_cols], dtype=bool)

    time_s = activity_df['time [s]'].astype(float).to_numpy()
    chi_q_vals = np.array([_chi_q_at_time(t, chi_q_schedule) for t in time_s], dtype=float)
    inhal_sum = vals @ dcf_inh_vec if nuclide_cols else np.zeros(len(activity_df))
    imm_sum = vals @ dcf_imm_vec if nuclide_cols else np.zeros(len(activity_df))

    # Fraction of the (positive) release rate without any DCF at each time step.
    positive = np.clip(vals, 0.0, None)
    total_rate = positive.sum(axis=1)
    missing_rate = positive[:, no_dcf].sum(axis=1)
    missing_frac = np.divide(missing_rate, total_rate, out=np.zeros_like(total_rate), where=total_rate > 0.0)
    max_missing = float(missing_frac.max()) if len(missing_frac) else 0.0
    if max_missing > 0.0:
        with np.errstate(divide='ignore', invalid='ignore'):
            share = np.where(total_rate[:, None] > 0.0, positive / total_rate[:, None], 0.0)
        worst = share.max(axis=0)
        missing = sorted(((worst[j], n) for j, n in enumerate(nuclide_cols) if no_dcf[j] and worst[j] > 0.0),
                         reverse=True)
        top = ', '.join(f"{n} ({w:.2e})" for w, n in missing[:5])
        msg = (f"Release without a DCF reaches {max_missing:.3e} of the total release rate; "
               f"largest nuclides without DCF (max step fraction): {top}")
        if max_missing_dcf_fraction is not None and max_missing > max_missing_dcf_fraction:
            raise ValueError(f"{msg} (limit {max_missing_dcf_fraction:.3e})")
        warnings.warn(msg, stacklevel=2)

    out = pd.DataFrame({
        'time [s]': time_s,
        'time [d]': time_s / float(24 * 60 * 60),
        'chi/Q [s/m^3]': chi_q_vals,
        'dose_rate_inhalation [Sv/s]': chi_q_vals * breathing_rate_m3_s * inhal_sum,
        'dose_rate_immersion [Sv/s]': chi_q_vals * (imm_sum / 86400.0),
        'activity_fraction_without_dcf [-]': missing_frac,
    })
    out['dose_rate [Sv/s]'] = out['dose_rate_inhalation [Sv/s]'] + out['dose_rate_immersion [Sv/s]']

    # Cumulative doses [Sv]: integral of chi/Q(t) * factor * S(t), S the DCF-weighted release sum, linear between
    # rows. Each row interval is split at the schedule boundaries, where chi/Q is constant on every piece.
    boundaries = _chi_q_boundaries(chi_q_schedule)

    def integrate(weighted_sum: np.ndarray, factor: float) -> np.ndarray:
        cumulative = np.zeros(len(time_s))
        for i in range(1, len(time_s)):
            t0, t1 = time_s[i - 1], time_s[i]
            cumulative[i] = cumulative[i - 1]
            if not t1 > t0:
                continue
            s0, s1 = weighted_sum[i - 1], weighted_sum[i]
            cuts = [t0] + [b for b in boundaries if t0 < b < t1] + [t1]
            for a, b in zip(cuts, cuts[1:]):
                s_a = s0 + (s1 - s0) * (a - t0) / (t1 - t0)
                s_b = s0 + (s1 - s0) * (b - t0) / (t1 - t0)
                chi = _chi_q_at_time(0.5 * (a + b), chi_q_schedule)
                cumulative[i] += chi * factor * 0.5 * (s_a + s_b) * (b - a)
        return cumulative

    out['dose_inhalation [Sv]'] = integrate(inhal_sum, breathing_rate_m3_s)
    out['dose_immersion [Sv]'] = integrate(imm_sum, 1.0 / 86400.0)
    out['dose [Sv]'] = out['dose_inhalation [Sv]'] + out['dose_immersion [Sv]']
    out['dose [rem]'] = out['dose [Sv]'] * 100.0
    return out


def _with_empty_start(activity_df: pd.DataFrame, case_prefix: str) -> pd.DataFrame:
    # Boxes B and C start empty at t = 0 and their json holds end-of-step states only, so a zero row at t = 0
    # closes the first interval of the dose integral. Box A, and frames that already start at t = 0, pass through.
    if (not case_prefix.lower().startswith(('box_b', 'box_c')) or activity_df.empty
            or float(activity_df['time [s]'].min()) <= 0.0):
        return activity_df
    start = pd.DataFrame([{c: 0.0 for c in activity_df.columns}])
    return pd.concat([start, activity_df], ignore_index=True)


def _element_of(nuclide: str) -> str:
    # 'kr-85m' -> 'kr'
    return re.sub('-.*', '', str(nuclide)).strip().lower()


def _inventory_to_release_rate(activity_df: pd.DataFrame, removal_rate_s: float | Mapping[str, float],
                               released_elements: list[str] | tuple[str, ...] | None = None) -> pd.DataFrame:
    """Convert nuclide inventory activity (Bq) into release rate (Bq/s).

    removal_rate_s is the box leak constant [1/s]: a scalar applies one rate to every released
    element, and a mapping carries one rate per element (PSAR MHA style). Only released_elements
    leak: the default is LEAKED_ELEMENTS for a scalar and the mapping keys for a mapping, and a
    scalar with an explicit list releases exactly that list. All other nuclides are retained in the
    box and get zero release. Every applied rate must be positive.
    """
    rates: dict[str, float] = _element_removal_rates(removal_rate_s)
    if released_elements is not None:
        released = {str(e).strip().lower() for e in released_elements}
        if isinstance(removal_rate_s, numbers.Real):
            rates = {ele: float(removal_rate_s) for ele in released}
        else:
            rates = {ele: rate for ele, rate in rates.items() if ele in released}
        if not rates:
            raise ValueError(f"released_elements {sorted(released)} shares no element with the removal rates")
    for ele, rate in rates.items():
        if not rate > 0.0:
            raise ValueError(f"removal_rate_s must be positive when converting inventory to release rate "
                             f"('{ele}' is {rate})")
    out = activity_df.copy()
    nuclide_cols = [c for c in out.columns if c not in ('time [s]', 'time [d]')]
    for col in nuclide_cols:
        out[col] = out[col] * rates.get(_element_of(col), 0.0)
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


def _removal_rate_from_run(box_json_path: str, case_prefix: str,
                           run_prefix: str | None) -> float | dict[str, float]:
    # The removal rate ORIGEN applied to this box, from the rates file _run_simulation writes next to the box json.
    # Older journal runs record a scalar. Per-element rates come back as a dict of the positive entries, because a
    # zero rate leaves that element in the box.
    rates_path = os.path.join(os.path.dirname(os.path.abspath(box_json_path)),
                              _out_name(run_prefix, 'leaky_boxes_rates', '.json5'))

    def no_rate() -> ValueError:
        return ValueError(f"{rates_path} has no usable removal rate for {case_prefix}")

    if not os.path.isfile(rates_path):
        raise ValueError(f"removal_rate_s is needed for activity_representation='inventory_bq', and the run did not "
                         f"leave {rates_path}. Pass removal_rate_s explicitly (the ORIGEN removal rate of "
                         f"{case_prefix}, PCTperDAY for box A and 0.1 * PCTperDAY for box B in run_test).")
    with open(rates_path, 'r') as f:
        rates = json5.load(f)
    try:
        applied: float | dict[str, float] = rates[case_prefix]['removal_rate_s']
    except (KeyError, TypeError):
        raise no_rate() from None
    if isinstance(applied, Mapping):
        try:
            applied = {str(ele).strip().lower(): float(rate) for ele, rate in applied.items()}
        except (TypeError, ValueError):
            raise no_rate() from None
        applied = {ele: rate for ele, rate in applied.items() if rate > 0.0}
    else:
        try:
            applied = float(applied)
        except (TypeError, ValueError):
            raise no_rate() from None
    if not (applied if isinstance(applied, Mapping) else applied > 0.0):
        raise ValueError(f"{case_prefix} has no removal in the ORIGEN run ({rates_path}), so its inventory is not "
                         f"released at a rate of its own. Pass removal_rate_s explicitly to use rate * A(t) as a "
                         f"release-rate surrogate, or use case_prefix='box_B' for the environmental release.")
    return applied


def _run_box_leaks(box_json_path: str, case_prefix: str, run_prefix: str | None) -> bool:
    # True when the run gave this box a positive removal rate of its own, i.e. it is a leaking
    # compartment rather than a collector, so rate * A(t) is its physical release rather than a surrogate.
    try:
        _removal_rate_from_run(box_json_path, case_prefix, run_prefix)
    except ValueError:
        return False
    return True


def compute_dose_from_box(box_json_path: str, case_prefix: str,
                          dcf_inhal_csv: str = DEFAULT_ICRP72_INHALATION_CSV,
                          chi_q_schedule: float | list[tuple[float, float]] | None = None,
                          breathing_rate_m3_s: float = BREATHING_RATE_M3_S,
                          dcf_immersion_csv: str | None = DEFAULT_ICRP72_IMMERSION_CSV,
                          out_prefix: str | None = None,
                          activity_representation: str = 'inventory_bq',
                          removal_rate_s: float | dict[str, float] | None = None,
                          released_elements: list[str] | tuple[str, ...] | None = None,
                          activity_scale: float = 1.0,
                          run_prefix: str | None = None,
                          max_missing_dcf_fraction: float | None = DEFAULT_MAX_MISSING_DCF_FRACTION) -> pd.DataFrame:
    """Compute and save inhalation+immersion dose from per-case F71 activities.

    The default inhalation table is the conservative ICRP-72 maximum across types F/M/S
    (DEFAULT_ICRP72_INHALATION_CSV), matching the journal method for unknown chemical forms.
    Pass dcf_icrp72_inhalation_adult.csv explicitly for the Type-F (fast-clearance) table.
    The defaults of the journal method: chi/Q from chi_q_schedule_rg145_400m() (RG 1.145, 400 m), the breathing
    rate BREATHING_RATE_M3_S (22.2 m^3/d) and the ICRP-72 immersion table for the noble gases.

    The release to the environment is the leak out of box B, removal_rate_s(box B) * A_B(t), so the
    environmental dose comes from case_prefix='box_B' with the default 'inventory_bq'. Box C collects that
    release and, in the journal method, has no removal, so its inventory is a cumulative amount. The journal
    method uses removal_rate_s * A_C(t) as a release-rate surrogate; box C with 'inventory_bq' then needs an
    explicit removal_rate_s and warns that the result is a surrogate. Integrating dA_C/dt = eps_B A_B - lambda A_C
    over a run long enough for A_C to decay gives int A_C dt = (eps_B / lambda) int A_B dt, so the time-integrated
    surrogate is eps_B / lambda of the box-B release (1.16 % for I-131 at 0.1 %/day), and a chi/Q schedule that
    weights early times lowers the dose ratio further. The instantaneous ratio A_C / A_B grows with time. A run that
    gives box C a positive removal rate of its own (a leaking third compartment, see `_run_simulation`) releases
    ε_C * A_C(t) physically, and that case_prefix carries no surrogate warning.

    Args:
        chi_q_schedule: A constant chi/Q [s/m^3] or a list of (t_end_s, chi_q) steps; None (default) selects
            chi_q_schedule_rg145_400m().
        dcf_immersion_csv: Immersion table; None drops the immersion dose of the noble gases.
        activity_representation: `inventory_bq` or `release_rate_bq_s`.
            - `inventory_bq`: activity in each nuclide column is Bq and is converted using `removal_rate_s`.
            - `release_rate_bq_s`: activity in each nuclide column is already Bq/s.
        removal_rate_s: Release constant [1/s] of the box for `inventory_bq`, either one rate for every released
            element or an {element: rate} mapping for per-element rates. None reads the rate ORIGEN applied
            to the box from `leaky_boxes_rates_<prefix>.json5`, which `_run_simulation` writes next to the box json.
        released_elements: Elements that leak at `removal_rate_s` (default LEAKED_ELEMENTS with a scalar rate, the
            `removal_rate_s` mapping keys with a per-element rate). Other nuclides are retained and give no release.
        activity_scale: Multiplier applied to the F71 activities before the dose calculation.
            The default leaky-box run (apply_volume_scaling=False) models boxes B and C with a volume
            of 1 cm^3 fed by the leak from 1 cm^3 of box A, so B/C activities and the resulting doses are
            per cm^3 of box A. Pass the box-A volume in cm^3 to get the whole-inventory dose.
            The value is written to the 'activity_scale [-]' output column.
        run_prefix: Case-directory prefix of the run (default: inferred from the json name,
            e.g. boxB_xe135.json5 -> xe135), with fallback to unprefixed legacy directories.
        max_missing_dcf_fraction: Raise ValueError when the release fraction without any DCF exceeds it
            (default DEFAULT_MAX_MISSING_DCF_FRACTION = 0.05; None warns only).
            All dose and dose-rate columns scale with activity_scale, so with the default
            activity_scale = 1.0 they are per cm^3 of box A; check 'activity_scale [-]'.
    """
    box_json_path = _resolve_run_artifact(box_json_path)
    if not os.path.isfile(box_json_path):
        raise FileNotFoundError(f"Box JSON not found: {box_json_path}")
    if chi_q_schedule is None:
        chi_q_schedule = chi_q_schedule_rg145_400m()
    if activity_representation == 'inventory_bq' and case_prefix.lower().startswith('box_c'):
        if not _run_box_leaks(box_json_path, case_prefix, run_prefix or _run_prefix_from_box_json(box_json_path)):
            warnings.warn("Box C has no removal, so removal_rate_s * A_C(t) is a release-rate surrogate (journal "
                          "method). Integrated over the run it is removal_rate_s / lambda of the physical release "
                          "for a nuclide with decay constant lambda. The environmental release is the box-B "
                          "leak: case_prefix='box_B'.", stacklevel=2)
    if activity_representation == 'inventory_bq' and removal_rate_s is None:
        removal_rate_s = _removal_rate_from_run(box_json_path, case_prefix,
                                                run_prefix or _run_prefix_from_box_json(box_json_path))
    dcf_inhal_csv = _resolve_input_file(dcf_inhal_csv, LEAKY_BOX_DATA_DIR)
    dcf_immersion_csv = _resolve_input_file(dcf_immersion_csv, LEAKY_BOX_DATA_DIR) if dcf_immersion_csv else None
    with _in_directory(os.path.dirname(os.path.abspath(box_json_path))):
        dcf_inhal = _load_dcf_csv(dcf_inhal_csv)
        dcf_imm = _load_dcf_immersion_csv(dcf_immersion_csv) if dcf_immersion_csv else {}
        activity_df = _with_empty_start(
            _activity_timeseries_per_nuclide_from_box_json(box_json_path, case_prefix, run_prefix=run_prefix),
            case_prefix)
        if activity_representation == 'inventory_bq':
            activity_df = _inventory_to_release_rate(activity_df, removal_rate_s, released_elements)
        elif activity_representation != 'release_rate_bq_s':
            raise ValueError("activity_representation must be 'inventory_bq' or 'release_rate_bq_s'")
        nuclide_cols = [c for c in activity_df.columns if c not in ('time [s]', 'time [d]')]
        activity_df[nuclide_cols] = activity_df[nuclide_cols] * float(activity_scale)
        pd_dose = compute_dose_timeseries(activity_df, dcf_inhal, chi_q_schedule,
                                          breathing_rate_m3_s, dcf_imm,
                                          max_missing_dcf_fraction=max_missing_dcf_fraction)
        pd_dose['activity_scale [-]'] = float(activity_scale)
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
                    out_prefix: str | None, box_B_n_steps: int = 8, box_C_n_steps: int = 8,
                    box_A_leak_rate: float | dict[str, float] = PCTperDAY,
                    box_B_leak_rate: float | dict[str, float] = PCTperDAY * 0.1,
                    box_C_leak_rate: float | dict[str, float] = 0.0
                    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    # Core leaky-box simulation (A->B->C) shared by tests and F71 runs. A leak rate is either one
    # removal constant for every released element or an {element: rate} mapping (PSAR MHA style);
    # _run_removal_rates resolves both to the per-element rates ORIGEN applies.
    # Box C collects the Box B leak and has no removal by default; a positive box_C_leak_rate turns it
    # into a real third compartment that leaks on its own, as the PSAR MHA pit is read to do.
    box_A_rates, box_B_rates, box_C_rates, released_elements = _run_removal_rates(
        box_A_leak_rate, box_B_leak_rate, box_C_leak_rate)

    # Boxes B and C are modelled as 1 cm^3 volumes. Box A atom densities are per cm^3, so each B/C
    # feed is the leak from 1 cm^3 of box A (unless apply_volume_scaling rescales it). B/C
    # inventories, activities and doses derived from them are therefore per cm^3 of box A.
    volume: float = 1.0

    # Case directories carry the output prefix so that runs sharing a run directory
    # (e.g. the xe135 and xe136 tests) do not overwrite each other's F71 files.
    box_a.case_dir = _case_dir_name('box_A', out_prefix)
    box_A = LeakyBox()
    box_A.removal_rate = box_A_rates or 0.0
    box_A.setup_cases(box_a)
    box_A.run_case()
    box_A.get_leak_rate_parallel()

    with open(_out_name(out_prefix, 'boxA', '.json5'), 'w') as f:  # Save to json
        json5.dump(box_A.adens, f, indent=4)
    with open(_out_name(out_prefix, 'boxA_leak_rate', '.json5'), 'w') as f:  # Save to json
        json5.dump(box_A.leak_rates, f, indent=4)
    with open(_out_name(out_prefix, 'leaky_boxes_rates', '.json5'), 'w') as f:  # removal rates for the dose step
        # The per-element rates ORIGEN applied, so the dose step releases only what each box lost. A box
        # without removal records 0.0.
        json5.dump({'box_A': {'removal_rate_s': box_A_rates or 0.0}, 'box_B': {'removal_rate_s': box_B_rates or 0.0},
                    'box_C': {'removal_rate_s': box_C_rates or 0.0}, 'volume_cm3': volume,
                    'released_elements': released_elements}, f, indent=4)

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
        origen_box_B[k].case_dir = _case_dir_name('box_B', out_prefix, k)

        box_B[k] = LeakyBox()
        box_B[k].removal_rate = box_B_rates or 0.0
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
        if box_C_rates:  # a leaking third compartment: ORIGEN removes these elements from box C too
            origen_box_C[k].nuclide_removal_rates = box_C_rates
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
        origen_box_C[k].case_dir = _case_dir_name('box_C', out_prefix, k)
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

    return pd_A, pd_B, pd_C, volume


def run_test(isotope: str, lambda_decay: float | None, apply_volume_scaling: bool,
             do_skip_calcs: bool, out_prefix: str | None = None,
             f71_path: str = DEFAULT_F71_PATH, f71_case: str = DEFAULT_F71_CASE,
             f71_time_seconds: float | None = None,
             sample_mass_g: float = 1200e3, decay_days: float = 12.0,
             steps_per_day: float = 2.0,
             run_dir: str | None = None, run_date: str | None = None
             ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run a single isotope test and return dataframes."""
    # Baseline decay calculation setup (forces a single-isotope inventory).
    f71_path = os.path.abspath(os.path.expanduser(f71_path))
    with _in_directory(_resolve_run_dir(run_dir=run_dir, run_date=run_date)):
        origen_triton_box_A = _setup_box_a(
            f71_path, f71_case, f71_time_seconds, sample_mass_g, do_skip_calcs, decay_days, steps_per_day
        )
        origen_triton_box_A.atom_dens = {isotope: 1.0}
        # The analytic solution takes one rate per box, so the test runs the scalar journal rates.
        box_A_leak_rate, box_B_leak_rate = PCTperDAY, PCTperDAY * 0.1
        pd_A, pd_B, pd_C, volume = _run_simulation(
            origen_triton_box_A, apply_volume_scaling, do_skip_calcs, out_prefix,
            box_A_leak_rate=box_A_leak_rate, box_B_leak_rate=box_B_leak_rate
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
                 steps_per_day: float = 2.0, run_dir: str | None = None,
                 run_date: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run a simulation using the actual F71 nuclide composition (no analytic comparison)."""
    # Uses full F71 composition; analytic columns are not added.
    f71_path = os.path.abspath(os.path.expanduser(f71_path))
    with _in_directory(_resolve_run_dir(run_dir=run_dir, run_date=run_date)):
        origen_triton_box_A = _setup_box_a(
            f71_path, f71_case, f71_time_seconds, sample_mass_g, do_skip_calcs, decay_days, steps_per_day
        )
        pd_A, pd_B, pd_C, volume = _run_simulation(
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
                  decay_days: float = 12.0, steps_per_day: float = 2.0,
                  run_dir: str | None = None, run_date: str | None = None):
    # Convenience wrapper for Xe-136 and Xe-135 test cases.
    tests = {
        'xe-136': None,
        'xe-135': 2.106574217602 * 1e-5,  # Xe-135 decay constant [1/s]
    }
    for iso, lam in tests.items():
        prefix = iso.replace('-', '')
        run_test(iso, lam, apply_volume_scaling, do_skip_calcs, out_prefix=prefix,
                 decay_days=decay_days, steps_per_day=steps_per_day,
                 run_dir=run_dir, run_date=run_date)


def main():
    """Default entrypoint: run Xe-136 test."""
    run_test('xe-136', None, apply_volume_scaling, do_skip_calcs, out_prefix='xe136')


if __name__ == "__main__":
    main()
