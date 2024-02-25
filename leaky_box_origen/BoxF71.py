#!/bin/env python3
"""
Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

Box A --> Box B --> outside

100% per day  =  1/24/60/60 per second    = 0.0000115740740741
   1% per day  =  1e-2/24/60/60 per second = 0.000000115740740741
 0.1% per day  =  1e-3/24/60/60 per second = 0.0000000115740740741
"""
import copy
import os
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


class DecaySalt(Origen):
    """ Decays salt from F71 file in ORIGEN, ment for box A decay """

    def __init__(self, _f71: str = './SCALE_FILE.f71', _mass: float = 500e3):
        Origen.__init__(self)
        self.BURNED_MATERIAL_F71_file_name: str = _f71  # Burned core F71 file from TRITON
        self.BURNED_MATERIAL_F71_index: dict = get_f71_positions_index(self.BURNED_MATERIAL_F71_file_name)
        self.BURNED_MATERIAL_F71_position: int = 16
        self.weight: float = _mass  # Mass of the sample [g]
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

        with open(self.ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
            f.write(atom_dens_for_origen(self.atom_dens))

        with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.origen_deck())

        if self.debug > 0:
            print(f'ORIGEN: decaying sample for {self.DECAY_days} days')
            print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

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


class DecayBox(Origen):
    """ Decays material from density list in ORIGEN using 1 time step, box B calculation """

    def __init__(self, _atom_dens: dict, _volume: float):
        Origen.__init__(self)
        self.atom_dens = _atom_dens
        self.volume: float = _volume  # Mass of the sample [g]
        self.nuclide_removal_rates: (dict, None) = None
        self.case_dir = 'box_B'
        self.ORIGEN_input_file_name = 'origen.inp'
        self.DECAY_seconds: float = 1000.0  # Time in seconds is saved from the first step (box A)

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        with open(self.ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
            f.write(atom_dens_for_origen(self.atom_dens))

        with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.origen_deck())

        if self.debug > 0:
            print(f'ORIGEN: decaying sample for {self.DECAY_seconds} seconds')
            print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

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
        units=SECONDS
        t=[{time_interp_steps}I 0.0001 {self.DECAY_seconds}]
        start=0
    }}
{removal} 
    save {{
        file="{self.F71_file_name}"
    }}
}}
end
'''
        return origen_output


class DiffLeak:
    """ Calculate difference between ORIGEN case with and without removals """

    def __init__(self):
        self.decay_base: (None, DecaySalt) = None
        self.decay_leaks: (None, DecaySalt) = None
        self.leaked: dict = {}
        self.removal_rate: float = PCTperDAY
        self.nuclide_removal_rates: dict = {}

    def setup_cases(self, base_case: DecaySalt):
        """ Setup cases based on a base_case """
        self.nuclide_removal_rates = {
            'H': self.removal_rate,
            'He': self.removal_rate,
            'O': self.removal_rate,
            'N': self.removal_rate,
            'I': self.removal_rate,
            'Ar': self.removal_rate,
            'Kr': self.removal_rate,
            'Xe': self.removal_rate
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
        print(f"Reading {i}, {t} s")
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
        res_list = Parallel(n_jobs=cpu_count())(delayed(self._parallel_leak_calc)(k, v) for k,v in times.items())
        for i, res in enumerate(res_list):
            print(i, res)
            self.leaked[i+1]['adens'] = res

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
        a['time'] = float(v['time'])
        _list.append(a)
    return pd.DataFrame(_list).set_index('time')


def main():
    """ Calculate which nuclides leak """
    # Baseline decay calculation setup
    origen_triton = DecaySalt('/home/o/MSRR-local/53-Ko1-cr2half/10-burn/33-SalstDose/SCALE_FILE.f71', 1200e3)
    origen_triton.set_f71_pos(5.0 * 365.24 * 24.0 * 60.0 * 60.0)  # 5 years
    origen_triton.read_burned_material()
    origen_triton.DECAY_days = 12
    origen_triton.DECAY_steps = 12
    volume: float = origen_triton.volume
    print(volume)

    # Calculate nuclide atom density in box B using the (1-2) difference between (1) box_A decay with no leakage
    # and (2) box_A decay with leakage
    dl_B = DiffLeak()
    dl_B.setup_cases(origen_triton)
    dl_B.run_cases()
    #dl_B.get_diff_rate()
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

        box_B_origen[k] = DecayBox(box_B_adens_current, volume)
        box_B_origen[k].DECAY_steps = 5
        box_B_origen[k].DECAY_seconds = v['time']
        box_B_origen[k].case_dir = f'box_B_{k:03d}'

        dl_C[k] = DiffLeak()
        dl_C[k].removal_rate = 0.1 * PCTperDAY  # Removal rate to box_C
        dl_C[k].setup_cases(box_B_origen[k])
        dl_C[k].run_cases()
        #dl_C[k].get_diff_rate()
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

if __name__ == "__main__":
    # pass
    main()
