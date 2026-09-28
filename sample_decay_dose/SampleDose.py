#!/bin/env python3
"""
Handling dose [rem/h] in SCALE
Ondrej Chvala <ochvala@utexas.edu>
"""
import os
import shutil
import re
import hashlib
import warnings
from datetime import datetime
import numpy as np
from sample_decay_dose.read_opus import integrate_opus
from sample_decay_dose.utils import nicely_print_atom_dens, get_rho_from_atom_density, scale_adens, \
    get_f71_positions_index, get_last_position_for_case, get_burned_nuclide_atom_dens, get_burned_nuclide_data, \
    get_burned_material_total_mass_dens, get_F33_num_sets, get_cyl_r, get_cyl_r_4_1, get_fill_height_4_1, get_cyl_h, \
    run_scale, atom_dens_for_origen, atom_dens_for_mavric

NOW: str = datetime.now().replace(microsecond=0).isoformat()
MAVRIC_NG_XSLIB: str = 'v7.1-28n19g'
DAY_IN_SECONDS: float = 24.0 * 60.0 * 60.0
# Zero-day decay: ORIGEN needs a positive time step after start=0, so the deck writes one step of this length.
# The zero-day state is the decay-case start, F71 position 1, which ORIGEN saves with its emission spectra.
ZERO_DAY_STEP_DAYS: float = 1.0e-6
MAVRIC_TALLY_SUMMARY: str = 'Final Tally Results Summary'
_NUMBER: str = r'[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?'
_POINT_DETECTOR_RE = re.compile(r'Point Detector\s+(\d+)\.\s*(.*?)\s*$')
_RESPONSE_LINE_RE = re.compile(rf'^\s*response\s+(\d+)\s+({_NUMBER})(?:\s+({_NUMBER}))?')


def run_scale_or_raise(deck_file: str, nmpi: int = 1, context_dir: str = ''):
    """ Run a SCALE deck via utils.run_scale() and raise RuntimeError on failure """
    if not run_scale(deck_file, nmpi):
        out_file: str = os.path.join(context_dir, deck_file.replace('.inp', '.out'))
        raise RuntimeError(f"SCALE run failed for '{context_dir}/{deck_file}', see output file: {out_file}")


def decay_time_grid(decay_days: float, n_positions: int, t_first_days: float = 1.0e-4) -> str:
    """ ORIGEN `t=[...]` list [days] for a decay case with start=0.
    Position 1 of the saved F71 is the start state. A positive decay time gives n_positions positions:
    t_first, n_positions - 3 log-spaced points, and decay_days at position n_positions.
    t_first is capped at decay_days / 1e3, so the grid stays strictly increasing for short decays.
    A zero decay time writes one ZERO_DAY_STEP_DAYS step, and the zero-day state is read at position 1.
    """
    if decay_days < 0.0:
        raise ValueError(f"Decay time cannot be negative: {decay_days} days")
    if decay_days == 0.0:
        return f't=[{ZERO_DAY_STEP_DAYS}]'
    time_interp_steps: int = n_positions - 3
    if time_interp_steps < 1:
        raise ValueError("Too few time steps")
    t_first: float = min(t_first_days, decay_days / 1e3)
    if not 0.0 < t_first < decay_days:
        raise ValueError(f"Decay time {decay_days} days is too short for a log-spaced ORIGEN time grid")
    return f't=[{time_interp_steps}L {t_first} {decay_days}]'


def decayed_f71_position(decay_days: float, n_positions: int) -> int:
    """ F71 position of the decayed state for a grid from decay_time_grid() """
    return 1 if decay_days == 0.0 else n_positions


def combine_dose_responses(responses: dict, keys: list[str], correlated: tuple[str, ...] = ()) -> dict:
    """ Sum of dose responses [rem/h] and its 1-sigma uncertainty.
    The point-detector tallies are independent, so their sigmas add in quadrature.
    The beta response is beta_over_gamma times the photon response. The (photon, beta) pair listed in
    `correlated` is therefore fully correlated, and its sigmas add linearly before the quadrature sum.
    """
    missing: list[str] = [k for k in keys if k not in responses]
    if missing:
        raise RuntimeError(f"Missing dose responses {missing}, run get_responses() first")
    value: float = sum(responses[k]['value'] for k in keys)
    correlated_sigma: float = sum(responses[k]['stdev'] for k in keys if k in correlated)
    variance: float = correlated_sigma ** 2 + sum(responses[k]['stdev'] ** 2 for k in keys if k not in correlated)
    return {'value': value, 'stdev': float(np.sqrt(variance))}


def parse_point_detector_responses(mavric_output: str) -> dict:
    """ Point-detector results from the MAVRIC 'Final Tally Results Summary' block, keyed by detector ID.
    Each detector block starts with a line like "Neutron Point Detector 5.  neutron detector" and holds a line
    "response 1   <value>  <stdev>  <rel. unc.> ..." with the response ID, value, and absolute 1-sigma.
    Returns {detector_id: {'particle', 'pid', 'value', 'stdev'}}. The particle is the word before "detector"
    in the tally title.
    """
    i_summary: int = mavric_output.find(MAVRIC_TALLY_SUMMARY)
    text: str = mavric_output[i_summary:] if i_summary >= 0 else mavric_output
    results: dict = {}
    detector_id: str | None = None
    particle: str = ''
    for line in text.splitlines():
        m_det = _POINT_DETECTOR_RE.search(line)
        if m_det:
            detector_id = m_det.group(1)
            m_particle = re.search(r'(\w+)\s+detector', m_det.group(2))
            particle = m_particle.group(1) if m_particle else line.split()[0].lower()
            continue
        if detector_id is None:
            continue
        m_resp = _RESPONSE_LINE_RE.match(line)
        if m_resp:
            results[detector_id] = {'particle': particle, 'pid': m_resp.group(1), 'value': float(m_resp.group(2)),
                                    'stdev': 0.0 if m_resp.group(3) is None else float(m_resp.group(3))}
            detector_id = None
    return results

# Predefined lists of atom densities for convenience. Check and make yours :)
# https://www.sandmeyersteel.com/316H.html
ADENS_SS316H_HOT: dict = {'c': 0.0002384, 'n-14': 0.000254609, 'n-15': 9.30164e-07, 'al-27': 5.30544e-05,
    'si-28': 1.56475e-05, 'si-29': 7.94907e-07, 'si-30': 5.24622e-07, 'p-31': 6.93324e-05, 's-32': 4.23788e-05,
    's-33': 3.34604e-07, 's-34': 1.89609e-06, 's-36': 4.46139e-09, 'ti-46': 3.28985e-06, 'ti-47': 2.96684e-06,
    'ti-48': 2.93973e-05, 'ti-49': 2.15734e-06, 'ti-50': 2.06562e-06, 'cr-50': 0.000677912, 'cr-52': 0.0130729,
    'cr-53': 0.00148236, 'cr-54': 0.00036899, 'mn-55': 1.73977e-05, 'fe-54': 0.00390032, 'fe-56': 0.0612267,
    'fe-57': 0.00141399, 'fe-58': 0.000188176, 'ni-58': 6.64309e-05, 'ni-60': 2.55891e-05, 'ni-61': 1.11234e-06,
    'ni-62': 3.54662e-06, 'ni-64': 9.03221e-07, 'mo-92': 0.00018364, 'mo-94': 0.00011476, 'mo-95': 0.00019769,
    'mo-96': 0.000207388, 'mo-97': 0.000118863, 'mo-98': 0.000300762, 'mo-100': 0.00012023}
ADENS_SS316H_COLD: dict = {'c': 0.000311167, 'si-28': 0.00153402, 'si-29': 7.79293e-05, 'si-30': 5.14317e-05,
    'p-31': 6.78722e-05, 's-32': 4.15187e-05, 's-33': 3.27814e-07, 's-34': 1.85761e-06, 's-36': 4.37085e-09,
    'cr-50': 0.000663653, 'cr-52': 0.0127979, 'cr-53': 0.00145118, 'cr-54': 0.000361229, 'mn-55': 0.00170071,
    'fe-54': 0.0031951, 'fe-56': 0.0501563, 'fe-57': 0.00115833, 'fe-58': 0.000154152, 'co-59': 7.93501e-05,
    'ni-58': 0.00650228, 'ni-60': 0.00250467, 'ni-61': 0.000108876, 'ni-62': 0.000347145, 'ni-64': 8.84075e-05,
    'mo-92': 0.000179806, 'mo-94': 0.000112364, 'mo-95': 0.000193563, 'mo-96': 0.000203058, 'mo-97': 0.000116381,
    'mo-98': 0.000294483, 'mo-100': 0.00011772}
ADENS_HDPE_COLD: dict = {'c-12': 3.992647e-02, 'c-13': 4.318339e-04, 'h-poly': 8.071660e-02}
ADENS_HELIUM_COLD: dict = {'he-3': 2.680585e-11, 'he-4': 2.693156e-05}
ADENS_HELIUM_HOT: dict = {'he-3': 8.27507e-12, 'he-4': 8.27506e-06}
ADENS_CONCRETE_COLD: dict = {'h-1': 0.01373939, 'o-16': 0.04606872, 'na-23': 0.001747024, 'al-27': 0.001745235,
    'si-28': 0.01532717, 'si-29': 0.0007786319, 'si-30': 0.0005138805, 'ca-40': 0.001474023, 'ca-42': 9.837866e-06,
    'ca-43': 2.052723e-06, 'ca-44': 3.171838e-05, 'ca-46': 6.082143e-08, 'ca-48': 2.843401e-06, 'fe-54': 2.029158e-05,
    'fe-56': 0.0003185345, 'fe-57': 7.356349e-06, 'fe-58': 9.789951e-07, 'co-59': 2.350275e-07, 'eu-151': 4.357685e-08,
    'eu-153': 4.756903e-08}
ADENS_KAOWOOL_COLD: dict = {'b-10': 4.39982e-07, 'b-11': 1.77098e-06, 'o-16': 0.00300423, 'o-17': 1.14439e-06,
    'o-18': 6.17368e-06, 'al-27': 0.000850507, 'si-28': 0.000770832, 'si-29': 3.91588e-05, 'si-30': 2.5844e-05,
    'ca-40': 1.66602e-06, 'ca-42': 1.11193e-08, 'ca-43': 2.32009e-09, 'ca-44': 3.58497e-08, 'ca-46': 6.87435e-11,
    'ca-48': 3.21376e-09, 'ti-46': 1.69203e-06, 'ti-47': 1.5259e-06, 'ti-48': 1.51195e-05, 'ti-49': 1.10956e-06,
    'ti-50': 1.06239e-06, 'fe-54': 7.05374e-07, 'fe-56': 1.10729e-05, 'fe-57': 2.55721e-07, 'fe-58': 3.40317e-08}
ADENS_LEAD_COLD: dict = {'pb-204': 4.615515e-04, 'pb-206': 7.945277e-03, 'pb-207': 7.285918e-03, 'pb-208': 1.727521e-02}
ADENS_DRYAIR_COLD: dict = {'c-12': 7.499998e-09, 'c-13': 8.111795e-11, 'n-14': 3.932965e-05, 'n-15': 1.436829e-07,
    'o-16': 1.057930e-05, 'o-17': 4.029883e-09, 'o-18': 2.174041e-08}


class Origen:
    """ ORIGEN handling parent class """

    def __init__(self):
        self.debug: int = 3  # Debugging flag
        self.cwd: str = os.getcwd()  # Current running fir
        self.ORIGEN_input_file_name: str = 'origen.inp'
        self.decayed_atom_dens: dict = {}  # Atom density of the decayed sample
        self.case_dir: str = ''
        self.SAMPLE_ATOM_DENS_file_name_Origen: str = 'my_sample_atom_dens_origen.inp'
        self.SAMPLE_F71_file_name: str = self.ORIGEN_input_file_name.replace('inp', 'f71')
        self.SAMPLE_F71_position: int = 12  # Sample decay steps; the last one holds the decayed state
        self.SAMPLE_DECAY_days: float = 30.0  # Sample decay time [days]
        self.sample_weight: float = np.nan  # Mass of the sample [g]
        self.sample_density: float = np.nan  # Mass density of the sample [g/cm3]
        self.sample_volume: float = np.nan  # Sample volume [cm3]

    @property
    def decayed_F71_position(self) -> int:
        """ Position of the decayed sample on the saved F71: SAMPLE_F71_position, or 1 for a zero-day decay """
        return decayed_f71_position(self.SAMPLE_DECAY_days, self.SAMPLE_F71_position)

    def decay_time_line(self, t_first_days: float = 1.0e-4) -> str:
        """ ORIGEN decay time list for SAMPLE_DECAY_days, see decay_time_grid() """
        return decay_time_grid(self.SAMPLE_DECAY_days, self.SAMPLE_F71_position, t_first_days)

    def opus_blocks(self) -> str:
        """ OPUS neutron, gamma, and beta spectra of the decayed state. The .plt files are numbered in this
        order (000...000 neutrons, 000...001 gamma, 000...002 beta), as get_neutron_integral() and
        get_beta_to_gamma() expect. The spectra are selected by F71 position, not by time. """
        out: str = ''
        for title, typarams in (('Neutrons', 'nspectrum'), ('Gamma', 'gspectrum'), ('Beta', 'bspectrum')):
            out += (f"\n=opus\n"
                    f"data='{self.SAMPLE_F71_file_name}'\n"
                    f"title='{title}'\n"
                    f"typarams={typarams}\n"
                    f"units=intensity\n"
                    f"npos={self.decayed_F71_position} end\n"
                    f"end\n")
        return out

    def set_decay_days(self, decay_days: float = 30.0):
        """ Use this to change decay time, as it also updates the case directory """
        self.SAMPLE_DECAY_days = decay_days
        self.case_dir: str = f'run_{NOW}_{self.sample_weight:.5}_g-{decay_days:.5}_days'  # Directory to run the case

    def get_beta_to_gamma(self) -> float:
        """ Calculates beta / gamma dose ratio as a ratio of respective spectral integrals """
        my_inp = self.ORIGEN_input_file_name  # temp variable to make the code PEP-8 compliant...
        gamma_spectrum_file = self.cwd + '/' + self.case_dir + '/' + my_inp.replace(".inp", ".000000000000000001.plt")
        beta_spectrum_file = self.cwd + '/' + self.case_dir + '/' + my_inp.replace(".inp", ".000000000000000002.plt")
        if not os.path.isfile(gamma_spectrum_file):
            raise FileNotFoundError("Expected OPUS file file:" + gamma_spectrum_file)
        if not os.path.isfile(beta_spectrum_file):
            raise FileNotFoundError("Expected OPUS file file:" + beta_spectrum_file)
        gamma_spectrum_integral = integrate_opus(gamma_spectrum_file)
        beta_spectrum_integral = integrate_opus(beta_spectrum_file)
        if gamma_spectrum_integral > 0:
            return beta_spectrum_integral / gamma_spectrum_integral
        else:
            return 0.0

    def get_neutron_integral(self) -> float:
        """ Calculates spectral integral of neutrons to see if neutrons shoudl be transported by MAVRIC """
        my_inp = self.ORIGEN_input_file_name  # temp variable to make the code PEP-8 compliant...
        neutron_spectrum_file = self.cwd + '/' + self.case_dir + '/' + my_inp.replace(".inp", ".000000000000000000.plt")
        if not os.path.isfile(neutron_spectrum_file):
            raise FileNotFoundError("Expected OPUS file file:" + neutron_spectrum_file)
        return integrate_opus(neutron_spectrum_file)


class OrigenFromTriton(Origen):
    """ Decays material generated by SCALE/TRITON sequence """

    def __init__(self, _f71: str = './SCALE_FILE.f71', _mass: float = 0.1):
        Origen.__init__(self)
        self.BURNED_MATERIAL_F71_file_name: str = _f71  # Burned core F71 file from TRITON
        self.BURNED_MATERIAL_F71_index: dict = get_f71_positions_index(self.BURNED_MATERIAL_F71_file_name)
        self.BURNED_MATERIAL_F71_position: int = 16
        self.burned_atom_dens: dict = {}  # Atom density of the burned material from F71 file
        self.sample_weight: float = _mass  # Mass of the sample [g]
        self.case_dir: str = f'run_{_mass:.5}_g'  # Directory to run the case

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
        else:  # closest time; a tie goes to the earlier time, then to the lower position
            pos_idx: int = min(range(len(times)), key=lambda i: (abs(times[i] - t), times[i]))
        pos_number: int = pos_times[pos_idx][0]
        print(f'--> Closest F71 position found at slot {pos_number}, {times[pos_idx]} seconds, '
              f'{times[pos_idx] / (60.0 * 60.0 * 24)} days.')

        self.BURNED_MATERIAL_F71_position = pos_number

    def read_burned_material(self):
        """ Reads atom density and rho from F71 file """
        if self.debug > 0:
            print(f'ORIGEN: reading nuclides from {self.BURNED_MATERIAL_F71_file_name}, '
                  f'position {self.BURNED_MATERIAL_F71_position}')
        self.sample_density = get_burned_material_total_mass_dens(self.BURNED_MATERIAL_F71_file_name,
                                                                  self.BURNED_MATERIAL_F71_position)
        self.sample_volume = self.sample_weight / self.sample_density
        self.burned_atom_dens = get_burned_nuclide_atom_dens(self.BURNED_MATERIAL_F71_file_name,
                                                             self.BURNED_MATERIAL_F71_position)
        if self.debug > 2:
            # print(list(self.burned_atom_dens.items())[:25])
            print(f'Sample density {self.sample_density} g/cm3, volume {self.sample_volume} cm3')
            nicely_print_atom_dens(self.burned_atom_dens)

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
                f.write(atom_dens_for_origen(self.burned_atom_dens))

            with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
                f.write(self.origen_deck())

            if self.debug > 0:
                print(f'ORIGEN: decaying sample for {self.SAMPLE_DECAY_days} days')
                print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
            run_scale_or_raise(self.ORIGEN_input_file_name, context_dir=self.case_dir)

            self.decayed_atom_dens = get_burned_nuclide_atom_dens(self.SAMPLE_F71_file_name, self.decayed_F71_position)
        finally:
            os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def origen_deck(self) -> str:
        """ Sample decay Origen deck """
        my_line_time: str = self.decay_time_line(1.0e-4)
        origen_output: str = f'''
=shell
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_Origen} .
end

=origen
' OrigenFromTriton {NOW} 
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
<{self.SAMPLE_ATOM_DENS_file_name_Origen}
        ]
        units=ATOMS-PER-BARN-CM
        volume={self.sample_volume}
    }}
    time {{
        units=DAYS
        {my_line_time}
        start=0
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end
'''
        origen_output += self.opus_blocks()
        return origen_output


class OrigenFromTritonMHA(OrigenFromTriton):
    """ Decays material generated by SCALE/TRITON sequence """

    def __init__(self, _f71: str = '../msrr.f71', _MTiHM: float = 0.40647261972527676, _f71_case: int = 20):
        super().__init__(_f71, _MTiHM * 1e6)
        self.BURNED_MATERIAL_F71_file_name: str = _f71  # Burned core F71 file from TRITON
        self.BURNED_MATERIAL_F71_index: dict = get_f71_positions_index(self.BURNED_MATERIAL_F71_file_name)
        self.BURNED_MATERIAL_F71_position: int = get_last_position_for_case(self.BURNED_MATERIAL_F71_file_name,
                                                                            _f71_case) - 1  # last step is just decay
        self.burned_atom_dens: dict = {}  # Atom density of the burned material from F71 file
        self.MTiHM: float = _MTiHM
        self.case_dir: str = f'run'  # Directory to run the case

    def origen_deck(self) -> str:
        """ Sample decay Origen deck """
        if self.debug > 0:
            print(
                f'[OrigenFromTritonMHA] Reading {self.BURNED_MATERIAL_F71_file_name} position {self.BURNED_MATERIAL_F71_position}, scaling {self.MTiHM}')
        # A zero-day deck reads the decay-case start (position 1). ORIGEN saves it after the "retained"
        # processing, so it holds the retained elements only.
        my_line_time: str = self.decay_time_line(1.0e-4)

        basename_burned_f71: str = os.path.basename(self.BURNED_MATERIAL_F71_file_name)
        origen_output: str = f'''
=shell
cp -r ${{INPDIR}}/../{self.BURNED_MATERIAL_F71_file_name} .
end

=origen
' OrigenFromTritonMHA {NOW} 
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
        load {{ 
            file="{basename_burned_f71}" pos={self.BURNED_MATERIAL_F71_position} 
        }}
    }}
    processing {{
        retained=[Te={self.MTiHM}, I={self.MTiHM}, Xe={self.MTiHM}, Br={self.MTiHM}, Kr={self.MTiHM}, H={self.MTiHM}]
    }}
    time {{
        units=DAYS
        {my_line_time}
        start=0
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end
'''
        origen_output += self.opus_blocks()
        return origen_output


class OrigenIrradiation(Origen):
    """ Irradiate and decay sample in Origen. Uses F33 file and total flux to irradiate a sample.
    F33 does not provide time index, by default the last one is taken.
    If the index needs to be set, use a corresponding F71 file.
    """

    def __init__(self, _f33: str = './SCALE_FILE.mix0007.f33', _mass: float = 0.1):
        Origen.__init__(self)
        self.F33_file_name: str = _f33  # Burned core F71 file from TRITON
        self.F33_position: int = get_F33_num_sets(self.F33_file_name)
        self.sample_weight: float = _mass  # Mass of the sample [g]
        self.case_dir: str = f'irr_{_mass:.5}_g'  # Directory to run the case
        self.irradiate_flux: float = 1.0e12  # Irradiation total flux [n/cm2/s]
        self.irradiate_days: float = 365.24  # Sample irradiation time [days]
        self.irradiate_steps: int = 30  # How many irradiation steps
        self.irradiate_F71_file_name: str = 'irradiate.f71'  # Saves the sample irradiation

    def write_atom_dens(self, atom_dens: dict = None):
        """ Writes atom density of the fresh material to be irradiated """
        if atom_dens is None:
            atom_dens = ADENS_SS316H_HOT
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            self.sample_density = get_rho_from_atom_density(atom_dens)
            with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write at-dens input for Origen irradiation
                f.write(atom_dens_for_origen(atom_dens))
        finally:
            os.chdir(self.cwd)

    def run_irradiate_decay_sample(self):
        """  Writes Origen input file, runs Origen to irradiate and decay it, and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir + '/' + self.SAMPLE_ATOM_DENS_file_name_Origen):
            raise FileNotFoundError("Write atom density for Origen first")
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
                f.write(self.origen_deck())

            if self.debug > 0:
                print(f'ORIGEN: burning sample for {self.irradiate_days} days at {self.irradiate_flux} n/cm2/s, '
                      f'then decaying for {self.SAMPLE_DECAY_days} days')
                print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
            run_scale_or_raise(self.ORIGEN_input_file_name, context_dir=self.case_dir)

            self.decayed_atom_dens = get_burned_nuclide_atom_dens(self.SAMPLE_F71_file_name, self.decayed_F71_position)
        finally:
            os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def read_irradiated_material_density(self):
        """ Reads rho [g/cm3] at the end of irradiation, the last position of the irradiation F71 in case_dir """
        irr_f71: str = os.path.join(self.cwd, self.case_dir, self.irradiate_F71_file_name)
        if not os.path.isfile(irr_f71):
            raise FileNotFoundError(f"Expected irradiation F71 file: {irr_f71}, run run_irradiate_decay_sample() first")
        last_position: int = max(get_f71_positions_index(irr_f71))
        self.sample_density = get_burned_material_total_mass_dens(irr_f71, last_position)
        self.sample_volume = self.sample_weight / self.sample_density

    def origen_deck(self) -> str:
        """ Sample irradiation and decay Origen deck """
        self.sample_volume = self.sample_weight / self.sample_density
        if not self.irradiate_days > 0.0:
            raise ValueError(f"Irradiation time must be positive: {self.irradiate_days} days")
        if self.irradiate_steps < 3:
            raise ValueError("Too few irradiation steps")
        # The self.F33_file_name potentially includes path to that file.
        # It gets copied to temp directory, where it is just F33_file
        F33_path, F33_file = os.path.split(self.F33_file_name)
        irr_days_min: float = min(0.0001, self.irradiate_days/1e3)
        my_line_time: str = self.decay_time_line(1.0e-4)
        origen_output: str = f'''
=shell
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_Origen} .
cp -r {self.cwd}/{self.F33_file_name} .
end

=origen
' OrigenIrradiation {NOW} 
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
case(irrad) {{
    lib {{ 
        file="{F33_file}" 
        pos={self.F33_position} 
    }}
    mat{{
        iso=[
<{self.SAMPLE_ATOM_DENS_file_name_Origen}
]
        units=ATOMS-PER-BARN-CM
        volume={self.sample_volume}
    }}
    time {{
        units=DAYS
        start=0
        t=[{self.irradiate_steps - 2}L {irr_days_min} {self.irradiate_days}]
    }}
    flux=[{self.irradiate_steps}R {self.irradiate_flux}]
    save {{
        file="{self.irradiate_F71_file_name}"
    }}
}}
case(decay) {{
    gamma=yes
    neutron=yes
    beta=yes
    lib {{
        file="end7dec"
    }}    
    time {{
        units=DAYS
        start=0
        {my_line_time}
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end
'''
        origen_output += self.opus_blocks()
        return origen_output


class OrigenDecayBox(Origen):
    """ Origen decay from a simple dict of atom density and volume [cm] """

    def __init__(self, _adens: (None, dict) = None, _vol: float = 0):
        Origen.__init__(self)
        self.SAMPLE_ATOM_DENSITY: (None, dict) = _adens
        self.sample_volume = _vol
        self.case_dir: str = self._decaybox_case_dir()  # Directory to run the case

    def _decaybox_case_dir(self) -> str:
        """ Case directory named by volume [cm3], decay time [days], and a hash of the composition """
        adens: dict = self.SAMPLE_ATOM_DENSITY or {}
        composition_hash: str = hashlib.sha1(repr(sorted(adens.items())).encode()).hexdigest()[:8]
        return f'_decaybox_{self.sample_volume:.5}_cm3-{self.SAMPLE_DECAY_days:.5}_days-{composition_hash}'

    def set_decay_days(self, decay_days: float = 30.0):
        """ Use this to change decay time, as it also updates the case directory """
        self.SAMPLE_DECAY_days = decay_days
        self.case_dir = self._decaybox_case_dir()

    def write_atom_dens(self):
        """ Writes atom density of the sample to decay """
        if self.SAMPLE_ATOM_DENSITY is None or self.sample_volume == 0:
            raise ValueError(f'Adens is None or Volume is zero: {self.SAMPLE_ATOM_DENSITY}, {self.sample_volume}')
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            self.sample_density = get_rho_from_atom_density(self.SAMPLE_ATOM_DENSITY)
            with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write at-dens input for Origen irradiation
                f.write(atom_dens_for_origen(self.SAMPLE_ATOM_DENSITY))
        finally:
            os.chdir(self.cwd)

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it, and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir + '/' + self.SAMPLE_ATOM_DENS_file_name_Origen):
            raise FileNotFoundError("Write atom density for Origen first")
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
                f.write(self.origen_deck())

            if self.debug > 0:
                print(f'ORIGEN: decaying for {self.SAMPLE_DECAY_days} days')
                print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
            run_scale_or_raise(self.ORIGEN_input_file_name, context_dir=self.case_dir)

            if self.debug > 1:
                print(self.SAMPLE_F71_file_name, self.decayed_F71_position)
            self.decayed_atom_dens = get_burned_nuclide_atom_dens(self.SAMPLE_F71_file_name, self.decayed_F71_position)
        finally:
            os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def origen_deck(self) -> str:
        """ Sample decay Origen deck """
        my_line_time: str = self.decay_time_line(1.0e-3)
        origen_output: str = f'''
=shell
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_Origen} .
end

=origen
' OrigenDecayBox {NOW} 
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
case(decay) {{
    gamma=yes
    neutron=yes
    beta=yes
    mat{{
        iso=[
<{self.SAMPLE_ATOM_DENS_file_name_Origen}
]
        units=ATOMS-PER-BARN-CM
        volume={self.sample_volume}
    }}
    lib {{
        file="end7dec"
    }}    
    time {{
        units=DAYS
        start=0
        {my_line_time}
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end
'''
        origen_output += self.opus_blocks()
        return origen_output


class DoseEstimator:
    """ MAVRIC calculation of rem/h doses from the decayed sample.
    The base geometry is a bare square cylinder of the sample in void. The point detector sits at det_x,
    measured from the sample centre (30 cm by default), so the distance from the sample surface is det_x - cyl_r.
    Responses: '1' neutron and '2' photon point-detector doses [rem/h]; '3' beta dose estimate, see get_responses().
    """

    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        # Manual-detector flags must exist before any det_x assignment is observed.
        object.__setattr__(self, '_det_x_manual', False)
        object.__setattr__(self, '_handling_det_x_manual', False)
        object.__setattr__(self, '_recomputing_detectors', False)
        self.debug: int = 3  # Debugging flag
        self.MAVRIC_input_file_name: str = 'my_dose.inp'
        self.SAMPLE_ATOM_DENS_file_name_MAVRIC: str = 'my_sample_atom_dens_mavric.inp'
        self.responses: dict = {}  # Dose responses 1: neutron, 2: gamma, 3: beta
        self.det_x: float = 30.0  # Detector distance from the sample centre [cm]
        self.N_planes_box: int = 5  # Planes per box
        self.N_planes_cyl: int = 8  # Planes per cylinder
        self.histories_per_batch: int = 100000  # Monaco hist per batch
        self.batches: int = 10  # Monaco number of batches in total
        self.box_a: float = np.nan
        self.cyl_r: float = np.nan
        self.sample_temperature_K: float = 873.0  # Sample temperature [K]
        self.decayed_atom_dens: dict = {}  # Atom density of the decayed sample
        self.beta_over_gamma: float | None = None  # Beta over gamma spectral ratio; None if unknown
        self.neutron_intensity: float | None = None  # Integral of neutron spectra [n/s]; None if unknown
        if _o is not None:
            self.sample_weight: float = _o.sample_weight  # Mass of the sample [g]
            self.sample_density: float = _o.sample_density  # Mass density of the sample [g/cm3]
            self.sample_volume: float = getattr(_o, 'sample_volume', np.nan)  # Sample volume [cm3], if exists
            self.DECAYED_SAMPLE_F71_file_name: str = _o.SAMPLE_F71_file_name
            self.DECAYED_SAMPLE_F71_position: int = decayed_f71_position(_o.SAMPLE_DECAY_days,
                                                                         _o.SAMPLE_F71_position)
            self.DECAYED_SAMPLE_days: float = _o.SAMPLE_DECAY_days  # Sample decay time [days]
            self.decayed_atom_dens: dict = _o.decayed_atom_dens  # Atom density of the decayed sample
            if _o.decayed_atom_dens:  # if the ORIGEN case was run...
                self.beta_over_gamma = _o.get_beta_to_gamma()  # Beta over gamma spectral ratio
                self.neutron_intensity = _o.get_neutron_integral()  # Integral of neutron spectra
            else:  # the OPUS spectra exist if the ORIGEN case was run earlier; otherwise the values stay unknown
                self.beta_over_gamma = self._spectral_integral_or_none(_o.get_beta_to_gamma)
                self.neutron_intensity = self._spectral_integral_or_none(_o.get_neutron_integral)
                if self.beta_over_gamma is not None or self.neutron_intensity is not None:
                    if self._opus_spectra_are_stale(_o):
                        warnings.warn(f"OPUS spectra in {_o.case_dir} predate {_o.SAMPLE_F71_file_name}; "
                                      f"they may describe a different decay time. Re-run the ORIGEN decay "
                                      f"before creating the dose estimator, or set beta_over_gamma / "
                                      f"neutron_intensity explicitly. Treating both as unknown.",
                                      stacklevel=2)
                        self.beta_over_gamma = None
                        self.neutron_intensity = None
                    else:
                        warnings.warn(f"OPUS spectra read from {_o.case_dir} on disk; they describe the decay "
                                      f"stored there. Re-run the ORIGEN decay after changing SAMPLE_DECAY_days.",
                                      stacklevel=2)
            self.ORIGEN_dir: str = _o.case_dir  # Directory to run the case
            self.case_dir: str = self.ORIGEN_dir + '_MAVRIC'
            self.cwd: str = _o.cwd  # Current running directory
        # __init__ assignments above are defaults, not manual overrides.
        object.__setattr__(self, '_det_x_manual', False)
        object.__setattr__(self, '_handling_det_x_manual', False)

    @staticmethod
    def _spectral_integral_or_none(read_integral) -> float | None:
        """ Spectral integral from the ORIGEN case OPUS files, or None when those files do not exist """
        try:
            return read_integral()
        except FileNotFoundError:
            return None

    @staticmethod
    def _opus_spectra_are_stale(_o) -> bool:
        """ True when on-disk OPUS .plt spectra predate the ORIGEN F71, so they may not match
        the current decay. OPUS runs after ORIGEN in the same deck, so fresh spectra are at
        least as new as the F71. A 2 s tolerance covers filesystem timestamp granularity. """
        try:
            case_dir = os.path.join(_o.cwd, _o.case_dir)
            f71_path = os.path.join(case_dir, _o.SAMPLE_F71_file_name)
            base = _o.ORIGEN_input_file_name.replace('.inp', '')
            plt_names = (base + '.000000000000000000.plt', base + '.000000000000000001.plt',
                         base + '.000000000000000002.plt')
            f71_mtime = os.path.getmtime(f71_path)
            for name in plt_names:
                path = os.path.join(case_dir, name)
                if os.path.isfile(path) and os.path.getmtime(path) < f71_mtime - 2.0:
                    return True
        except (AttributeError, OSError):
            return False
        return False

    def __setattr__(self, name, value):
        # Track a user-assigned detector position so tank decks can warn instead of
        # silently discarding it. Internal recomputation sets _recomputing_detectors.
        if name in ('det_x', 'handling_det_x') and not getattr(self, '_recomputing_detectors', False):
            if '_det_x_manual' in self.__dict__:
                flag = '_det_x_manual' if name == 'det_x' else '_handling_det_x_manual'
                object.__setattr__(self, flag, True)
        object.__setattr__(self, name, value)

    def _set_computed_detectors(self, det_x: float, handling_det_x: float | None = None):
        """ Store mavric_deck() detector positions, warning when a manual det_x is discarded.

        Tank estimators place detectors from det_standoff_distance /
        handling_det_standoff_distance. A direct det_x assignment cannot survive
        the deck build, so it warns and points at the standoff attribute.
        """
        if getattr(self, '_det_x_manual', False):
            try:
                differs = abs(float(self.det_x) - float(det_x)) > 1e-9
            except (TypeError, ValueError):
                differs = True
            if differs:
                warnings.warn(f"det_x={self.det_x} is discarded by mavric_deck(); "
                              f"set det_standoff_distance instead (computed det_x={det_x:.5g})",
                              stacklevel=3)
        if handling_det_x is not None and getattr(self, '_handling_det_x_manual', False):
            try:
                h_differs = abs(float(self.handling_det_x) - float(handling_det_x)) > 1e-9
            except (TypeError, ValueError, AttributeError):
                h_differs = True
            if h_differs:
                warnings.warn(f"handling_det_x={getattr(self, 'handling_det_x', None)} is discarded by "
                              f"mavric_deck(); set handling_det_standoff_distance instead "
                              f"(computed handling_det_x={handling_det_x:.5g})",
                              stacklevel=3)
        object.__setattr__(self, '_recomputing_detectors', True)
        try:
            self.det_x = det_x
            if handling_det_x is not None:
                self.handling_det_x = handling_det_x
        finally:
            object.__setattr__(self, '_recomputing_detectors', False)
            object.__setattr__(self, '_det_x_manual', False)
            object.__setattr__(self, '_handling_det_x_manual', False)

    @property
    def MAVRIC_out_file_name(self) -> str:
        return self.MAVRIC_input_file_name.replace('inp', 'out')

    @property
    def beta_applies(self) -> bool:
        """ True when the beta estimate (beta_over_gamma times the photon dose) applies to the detector.
        The estimate holds for a bare sample only. Monaco does not transport electrons, so behind any layer
        the beta response is reported as zero. The base geometry is a bare sample. """
        return True

    def _include_neutron_source(self, note: bool = False) -> bool:
        """ Neutron source and distribution go into the deck unless the ORIGEN neutron spectrum is known to be empty.
        An unknown intensity (None) includes the source. """
        if self.neutron_intensity is None:
            if note and self.debug > 0:
                print('Note: ORIGEN neutron intensity unknown (no OPUS spectra), including the F71 neutron source')
            return True
        return self.neutron_intensity > 0.0

    def _beta_response(self, gamma_response: dict) -> dict:
        """ Beta dose estimate from the photon response. Zero when the beta estimate does not apply.
        The value and stdev are beta_over_gamma times the photon values, so beta is fully correlated with gamma. """
        if not self.beta_applies:
            return {'value': 0.0, 'stdev': 0.0}
        if self.beta_over_gamma is None:
            raise RuntimeError('The beta/gamma ratio is unknown, because the ORIGEN OPUS spectra were not found. '
                               'Run the ORIGEN decay in this process before creating the dose estimator, '
                               'or set beta_over_gamma explicitly.')
        return {'value': self.beta_over_gamma * gamma_response['value'],
                'stdev': self.beta_over_gamma * gamma_response['stdev']}

    def run_mavric(self):
        """ Writes Mavric inputs and runs the case """
        if not os.path.isfile(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name):
            raise FileNotFoundError(
                "Expected decayed sample F71 file: \n" + self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name)
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            shutil.copy2(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name,
                         self.cwd + '/' + self.case_dir)
            os.chdir(self.cwd + '/' + self.case_dir)

            with open(self.SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
                f.write(atom_dens_for_mavric(self.decayed_atom_dens, 1, self.sample_temperature_K))

            with open(self.MAVRIC_input_file_name, 'w') as f:  # write MAVRIC input deck
                f.write(self.mavric_deck())

            if self.debug > 0:
                print(f"MAVRIC: running case {self.case_dir}/{self.MAVRIC_input_file_name}")
            run_scale_or_raise(self.MAVRIC_input_file_name, context_dir=self.case_dir)
        finally:
            os.chdir(self.cwd)

    def get_responses(self):
        """ Reads over the MAVRIC output and returns responses for rem/h doses """
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name):
            raise FileNotFoundError(
                "Expected decayed sample MAVRIC output file: \n" + self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)
        try:
            tally_sep: str = 'Final Tally Results Summary'
            is_in_tally: bool = False
            response_re = re.compile(
                r'\bresponse\s+(\d+)\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)?'
            )

            self.responses = {}
            with open(self.MAVRIC_out_file_name, 'r') as f:
                for line in f.read().splitlines():
                    if tally_sep in line:
                        is_in_tally = True
                        continue
                    if not is_in_tally:
                        continue
                    m = response_re.search(line)
                    if not m:
                        continue
                    rid: str = m.group(1)
                    value: float = float(m.group(2))
                    stdev: float = 0.0 if m.group(3) in (None, '') else float(m.group(3))
                    self.responses[rid] = {'value': value, 'stdev': stdev}

            for rid, particle in (('1', 'neutron'), ('2', 'photon')):
                if rid not in self.responses:
                    raise RuntimeError(f'Failed to parse {particle} response {rid} from {self.MAVRIC_out_file_name}')
            self.responses['3'] = self._beta_response(self.responses['2'])
        finally:
            os.chdir(self.cwd)
        if self.debug > 3:
            print(self.responses)

    def print_response(self):
        """ Prints dose responses """
        if self.responses:
            r1: dict = self.responses['1']
            r2: dict = self.responses['2']
            r3: dict = self.responses['3']
            print(self.sample_weight, r1['value'], r1['stdev'], r2['value'], r2['stdev'], r3['value'], r3['stdev'])

    @property
    def total_dose(self) -> dict:
        """ Total dose [rem/h], the sum of all responses, with its 1-sigma uncertainty.
        Beta '3' is correlated with photon '2', see combine_dose_responses(). """
        if not self.responses:
            raise RuntimeError("No dose responses, run get_responses() first")
        return combine_dose_responses(self.responses, list(self.responses), correlated=('2', '3'))

    def _deck_head(self, title_line: str, comp_pre: str = '', comp_extra: str = '') -> str:
        """ Common deck preamble: shell copies, header, parameters, composition """
        adjoint_flux_file: str = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        return (f'\n'
                f'=shell\n'
                f'cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .\n'
                f'cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .\n'
                f"'cp -r ${{INPDIR}}/{adjoint_flux_file} .\n"
                f'end\n'
                f'\n'
                f'=mavric parm=(   )\n'
                f'{NOW} {title_line}\n'
                f'{MAVRIC_NG_XSLIB}\n'
                f'\n'
                f'read parameters\n'
                f'    randomSeed=0000000100000001\n'
                f'    ceLibrary="ce_v7.1_endf.xml"\n'
                f'    neutrons  photons\n'
                f'    fissionMult=1  secondaryMult=1\n'
                f'    perBatch={self.histories_per_batch} batches={self.batches}\n'
                f'end parameters\n'
                f'\n'
                f'read comp\n'
                f'{comp_pre}'
                f'<{self.SAMPLE_ATOM_DENS_file_name_MAVRIC}\n'
                f'{comp_extra}'
                f'end comp\n')

    def _deck_distributions(self) -> str:
        """ ORIGEN spectral distributions; neutron distribution unless the neutron intensity is known to be zero """
        out: str = ''
        if self._include_neutron_source(note=True):
            out += (f'\n'
                    f'    distribution 1\n'
                    f'        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, neutrons"\n'
                    f'        special="origensBinaryConcentrationFile"\n'
                    f'        parameters {self.DECAYED_SAMPLE_F71_position} 1 end\n'
                    f'        filename="{self.DECAYED_SAMPLE_F71_file_name}"\n'
                    f'    end distribution')
        out += (f'\n'
                f'    distribution 2\n'
                f'        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, photons"\n'
                f'        special="origensBinaryConcentrationFile"\n'
                f'        parameters {self.DECAYED_SAMPLE_F71_position} 5 end\n'
                f'        filename="{self.DECAYED_SAMPLE_F71_file_name}"\n'
                f'    end distribution\n')
        return out

    def _deck_tail(self, neutron_source_cylinder: str, photon_source_cylinder: str,
                   multiplier: float | None = None, handling_detector: bool = False) -> str:
        """ Common deck tail: sources, importance map, tallies, closing """
        neutron_src: str = ''
        if self._include_neutron_source():
            mult_n: str = f'\n        multiplier={multiplier}' if multiplier is not None else ''
            neutron_src = (f'\n'
                           f'    src 1\n'
                           f'        title="Sample neutrons"\n'
                           f'        neutron\n'
                           f'        useNormConst\n'
                           f'        cylinder {neutron_source_cylinder}\n'
                           f'        eDistributionID=1'
                           f'{mult_n}\n'
                           f'    end src')
        mult_p: str = f'\n        multiplier={multiplier}        ' if multiplier is not None else ''
        out: str = (f'\nread sources'
                    f'{neutron_src}'
                    f'\n'
                    f'    src 2\n'
                    f'        title="Sample photons"\n'
                    f'        photon\n'
                    f'        useNormConst\n'
                    f'        cylinder {photon_source_cylinder}\n'
                    f'        eDistributionID=2'
                    f'{mult_p}\n'
                    f'    end src\n'
                    f'end sources\n'
                    f'\n'
                    f'read importanceMap\n'
                    f'   gridGeometryID=1\n'
                    f"'   adjointFluxes=\"{self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')}\"\n"
                    f'   adjointSource 1\n'
                    f'        locationID=1\n'
                    f'        responseID=1\n'
                    f'   end adjointSource\n'
                    f'   adjointSource 2\n'
                    f'        locationID=1\n'
                    f'        responseID=2\n'
                    f'   end adjointSource\n')
        if handling_detector:
            out += (f'   adjointSource 5\n'
                    f'        locationID=2\n'
                    f'        responseID=1\n'
                    f'   end adjointSource\n'
                    f'   adjointSource 6\n'
                    f'        locationID=2\n'
                    f'        responseID=2\n'
                    f'   end adjointSource\n')
        out += (f'end importanceMap\n'
                f'\n'
                f'read tallies\n'
                f'    pointDetector 1\n'
                f'        title="neutron detector"\n'
                f'        neutron\n'
                f'        locationID=1\n'
                f'        responseID=1\n'
                f'    end pointDetector\n'
                f'    pointDetector 2\n'
                f'        title="photon detector"\n'
                f'        photon\n'
                f'        locationID=1\n'
                f'        responseID=2\n'
                f'    end pointDetector\n')
        if handling_detector:
            out += (f'    pointDetector 5\n'
                    f'        title="neutron detector"\n'
                    f'        neutron\n'
                    f'        locationID=2\n'
                    f'        responseID=1\n'
                    f'    end pointDetector\n'
                    f'    pointDetector 6\n'
                    f'        title="photon detector"\n'
                    f'        photon\n'
                    f'        locationID=2\n'
                    f'        responseID=2\n'
                    f'    end pointDetector\n'
                    f'\n'
                    f'end tallies\n'
                    f'\n'
                    f'end data\n'
                    f'end\n')
        else:
            out += (f'end tallies\n'
                    f'\n'
                    f'end data\n'
                    f'end\n')
        return out

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file. The detector is at det_x from the sample centre. """
        self.cyl_r = get_cyl_r(self.sample_volume)
        self.box_a: float = self.det_x + 10.0  # Problem box distance [cm]
        if not self.det_x > self.cyl_r:
            raise ValueError(f'Detector at det_x={self.det_x} cm from the sample centre lies inside the sample '
                             f'(radius and half-height {self.cyl_r:.5g} cm). Increase det_x.')
        if not self.box_a > self.cyl_r:
            raise ValueError(f'Problem boundary at {self.box_a} cm cuts the sample (radius {self.cyl_r:.5g} cm)')
        mavric_output: str = self._deck_head(
            title_line=f'Sample dose, {self.sample_weight} g, at x={self.det_x} cm')
        mavric_output += f'''
read geometry
global unit 1
    cylinder 1 {self.cyl_r} 2p {self.cyl_r}
    cuboid 99  6p {self.box_a}
    media 1 1 1
    media 0 1 99 -1
boundary 99
end geometry

read definitions
     location 1
        position {self.det_x} 0 0
    end location
    response 1
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
'''

        mavric_output += self._deck_distributions()
        mavric_output += f'''
    gridGeometry 1
        title="Grid over the problem"
        xLinear {self.N_planes_box} -{self.box_a} {self.box_a}
        yLinear {self.N_planes_box} -{self.box_a} {self.box_a}
        zLinear {self.N_planes_box} -{self.box_a} {self.box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
    end gridGeometry
end definitions
'''
        mavric_output += self._deck_tail(
            neutron_source_cylinder=f'{self.cyl_r} {self.cyl_r} -{self.cyl_r}',
            photon_source_cylinder=f'{self.cyl_r} {self.cyl_r} -{self.cyl_r}')
        return mavric_output


class DoseEstimatorSquareTank(DoseEstimator):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a square tank made of materials.
    layers_thicknesses, layers_mats, and layers_temperature_K describe the shielding layers from the inside out
    and must have equal lengths. Empty lists give a bare sample. The detector is det_standoff_distance
    from the outer surface of the last layer; mavric_deck() recomputes det_x from it.
    """

    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        super().__init__(_o)  # Init DoseEstimator
        self.cyl_r: (None, float) = None  # tank cylinder inner radius = sample radius
        self.box_a: float = 10.0  # Problem box distance offset [cm]
        self.det_standoff_distance: float = 1.0  # Detector is 1 cm off the tank
        self.planes_xy_around_det: float = 0.5  # Additional +-dX/dY planes around detector for CADIS
        self.layers_thicknesses: list[float] = [2.54, 2.0 * 2.54, 3.0 * 2.54]
        self.layers_mats: list[dict] = [ADENS_SS316H_HOT, ADENS_HDPE_COLD, ADENS_SS316H_COLD]
        self.layers_temperature_K: list[float] = [873.0, 300.0, 300.0]
        self._case_dir_suffix: str = ''  # layer suffix that run_mavric() appended to case_dir

    @property
    def beta_applies(self) -> bool:
        """ The beta estimate applies to a bare sample only, i.e., when there are no layers """
        return not self.layers_mats

    @property
    def outer_body(self) -> int:
        """ KENO-VI body ID of the outermost cylinder: the sample is body 1, layer k is body k + 2 """
        return len(self.layers_mats) + 1

    def _validate_layers(self):
        """ The layer lists are usually assigned after __init__, so they are checked before each deck """
        n_layers: int = len(self.layers_mats)
        if len(self.layers_thicknesses) != n_layers or len(self.layers_temperature_K) != n_layers:
            raise ValueError(f"layers_thicknesses ({len(self.layers_thicknesses)}), layers_mats ({n_layers}), and "
                             f"layers_temperature_K ({len(self.layers_temperature_K)}) need the same length.")
        if any(not t > 0.0 for t in self.layers_thicknesses):
            raise ValueError(f"Layer thicknesses must be positive: {self.layers_thicknesses}")

    def _layer_suffix(self) -> str:
        """ Layer thicknesses in the case directory name, unique IDs for parallel runs """
        return "".join([f'_{s:.3f}' for s in self.layers_thicknesses])

    def _layer_case_dir(self) -> str:
        """ case_dir with the current layer suffix, replacing the suffix of an earlier run_mavric() call """
        base: str = self.case_dir
        if self._case_dir_suffix and base.endswith(self._case_dir_suffix):
            base = base[:-len(self._case_dir_suffix)]
        return base + self._layer_suffix()

    def run_mavric(self, nmpi: int = 1):
        """ Writes Mavric inputs and runs the case """
        self._validate_layers()
        self.case_dir = self._layer_case_dir()
        self._case_dir_suffix = self._layer_suffix()
        if not os.path.isfile(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name):
            raise FileNotFoundError(
                "Expected decayed sample F71 file: \n" + self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name)
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            shutil.copy2(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name,
                         self.cwd + '/' + self.case_dir)
            os.chdir(self.cwd + '/' + self.case_dir)

            with open(self.SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
                f.write(atom_dens_for_mavric(self.decayed_atom_dens, 1, self.sample_temperature_K))
                for k in range(len(self.layers_mats)):
                    f.write(atom_dens_for_mavric(self.layers_mats[k], k + 10, self.layers_temperature_K[k]))

            with open(self.MAVRIC_input_file_name, 'w') as f:  # write MAVRICinput deck
                f.write(self.mavric_deck())

            if self.debug > 0:
                print(f"MAVRIC: running case {self.case_dir}/{self.MAVRIC_input_file_name}")
            run_scale_or_raise(self.MAVRIC_input_file_name, nmpi, context_dir=self.case_dir)
        finally:
            os.chdir(self.cwd)

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self._validate_layers()
        self.cyl_r = get_cyl_r(self.sample_volume)
        tank_r: float = self.cyl_r  # current outer layer [cm]
        mavric_output: str = self._deck_head(
            title_line=f'DoseEstimatorSquareTank, {self.sample_weight} g, layers {self.layers_thicknesses}')
        mavric_output += f'''
read geometry
global unit 1
    cylinder 1 {self.cyl_r} 2p {self.cyl_r}
    media 1 1 1
'''
        x_planes: list[float] = []  # list of cylinder boundaries for gridgeometry
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_r}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            x_planes.append(tank_r)
        self._set_computed_detectors(tank_r + self.det_standoff_distance)  # Detector is next to the tank
        box_a: float = self.box_a + tank_r  # problem half-width; self.box_a stays the offset
        x_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in x_planes])

        mavric_output += f'''
    cuboid 99999  6p {box_a}
    media 0 1 99999 -{self.outer_body}
boundary 99999
end geometry

read definitions
     location 1
        position {self.det_x} 0 0
    end location
    response 1
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
'''

        mavric_output += self._deck_distributions()
        mavric_output += f'''
    gridGeometry 1
        title="Grid over the problem, location at +x"
        xLinear {self.N_planes_box} {self.cyl_r} {box_a}
'        yLinear {self.N_planes_box} -{box_a} {box_a}
'        zLinear {self.N_planes_box} -{box_a} {box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        xPlanes {x_planes_str} -{box_a} {box_a} {self.det_x + self.planes_xy_around_det} {self.det_x - self.planes_xy_around_det} end
        yPlanes {x_planes_str} -{box_a} {box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
        zPlanes {x_planes_str} -{box_a} {box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
    end gridGeometry
end definitions
'''
        mavric_output += self._deck_tail(
            neutron_source_cylinder=f'{self.cyl_r} {self.cyl_r} -{self.cyl_r}',
            photon_source_cylinder=f'{self.cyl_r} {self.cyl_r} -{self.cyl_r}')
        return mavric_output


class DoseEstimatorStorageTank(DoseEstimatorSquareTank):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a storage tank made of materials
    The storage tank is 4:1  H:R cylinder, which volume is 20% larger than that of the sample.
    Tank plenum is filled with hot helium.
    Geometry center is the middle of the tank, not the sample!
    """

    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        self.plenum_volume_fraction: float = 0.2  # gas plenum above fill is +20% of sample volume
        self.sample_offset_z: (None, float) = None  # tank cylinder z - sample cylinder z [cm]
        self.sample_h2: (None, float) = None  # half-height of the sample
        self.det_z: (None, float) = None  # z location of the detector
        super().__init__(_o)  # Init DoseEstimator

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self._validate_layers()
        tank_inner_volume: float = self.sample_volume + self.sample_volume * self.plenum_volume_fraction
        self.cyl_r = get_cyl_r_4_1(tank_inner_volume)
        tank_r: float = self.cyl_r  # current outer layer radius [cm]
        tank_h2: float = 2.0 * self.cyl_r  # current outer layer half-height [cm]
        self.sample_h2 = get_fill_height_4_1(self.sample_volume, tank_inner_volume) / 2.0
        self.sample_offset_z = 2.0 * (tank_h2 - self.sample_h2)
        sample_z_max: float = tank_h2 - self.sample_offset_z
        sample_z_min: float = -tank_h2

        mavric_output: str = self._deck_head(
            title_line=f'DoseEstimatorStorageTank, {self.sample_weight} g, '
                       f'layers {self.layers_thicknesses}, plenum fraction {self.plenum_volume_fraction}',
            comp_extra='helium 2 end\n')
        mavric_output += f'''
read geometry
global unit 1
    cylinder 1 {tank_r} {sample_z_max} {sample_z_min}
    cylinder 2 {tank_r} {tank_h2} {sample_z_max}
    media 1 1 1
    media 2 1 2
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            tank_h2 += self.layers_thicknesses[k]
            mavric_output += f'    cylinder {k + 3} {tank_r} 2p {tank_h2}\n'
            if k == 0:  # Special case for gas plenum
                mavric_output += f'    media 10 1 -1 -2 3\n'
            else:
                mavric_output += f'    media {k + 10}  1 -{k + 2} {k + 3}\n'
            xy_planes.append(tank_r)
            z_planes.append(tank_h2)
            z_planes.append(-tank_h2)
        # Bodies: sample 1, plenum 2, layer k is body k + 3. Without layers the void is outside sample and plenum.
        outside_tank: str = f'-{len(self.layers_mats) + 2}' if self.layers_mats else '-1 -2'
        self._set_computed_detectors(tank_r + self.det_standoff_distance)  # Detector is next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 {outside_tank}
boundary 99999
end geometry

read definitions
     location 1
        position {self.det_x} 0 {self.det_z}
    end location
    response 1
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
'''

        mavric_output += self._deck_distributions()
        mavric_output += f'''
    gridGeometry 1
        title="Grid over the problem, location at +x"
        xLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a}
'        yLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a} 
'        zLinear {self.N_planes_box} {tank_h2} {tank_h2 + self.box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl}  {sample_z_min} {sample_z_max}
        xPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.det_x - self.planes_xy_around_det} {self.det_x + self.planes_xy_around_det} end
        yPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
        zPlanes {z_planes_str} {tank_h2 + self.box_a} {-tank_h2 - self.box_a} end
    end gridGeometry
end definitions
'''
        mavric_output += self._deck_tail(
            neutron_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}',
            photon_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}')
        return mavric_output


class DoseEstimatorGenericTank(DoseEstimatorSquareTank):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a storage tank made of materials
    The storage tank is a cylinder, filled with the sample.
    """

    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        self.sample_h2: (None, float) = None  # half-height of the sample
        self.det_z: (None, float) = None  # z location of the detector
        self.cyl_r: (None, float) = None  # inner radius of the tank cylinder; h is calculated from V & r
        super().__init__(_o)  # Init DoseEstimator

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self._validate_layers()
        self.sample_h2 = get_cyl_h(self.sample_volume, self.cyl_r) / 2.0
        tank_r: float = self.cyl_r
        tank_h2: float = self.sample_h2
        sample_z_max: float = self.sample_h2
        sample_z_min: float = - self.sample_h2

        mavric_output: str = self._deck_head(
            title_line=f'DoseEstimatorGenericTank, {self.sample_weight} g, layers {self.layers_thicknesses}',
            comp_extra="' helium 2 end\n")
        mavric_output += f'''
read geometry
global unit 1
    cylinder 1 {tank_r} 2p {tank_h2}
    media 1 1 1
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            tank_h2 += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_h2}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            xy_planes.append(tank_r)
            z_planes.append(tank_h2)
            z_planes.append(-tank_h2)
        self._set_computed_detectors(tank_r + self.det_standoff_distance)  # Detector is next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 -{self.outer_body}
boundary 99999
end geometry

read definitions
     location 1
        position {self.det_x} 0 {self.det_z}
    end location
    response 1
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
'''

        mavric_output += self._deck_distributions()
        mavric_output += f'''
    gridGeometry 1
        title="Grid over the problem, location at +x"
        xLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a}
'        yLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a} 
'        zLinear {self.N_planes_box} {tank_h2} {tank_h2 + self.box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl}  {sample_z_min} {sample_z_max}
        xPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.det_x - self.planes_xy_around_det} {self.det_x + self.planes_xy_around_det} end
        yPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
        zPlanes {z_planes_str} {tank_h2 + self.box_a} {-tank_h2 - self.box_a} end
    end gridGeometry
end definitions
'''
        mavric_output += self._deck_tail(
            neutron_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}',
            photon_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}')
        return mavric_output


class HandlingContactDoseEstimatorGenericTank(DoseEstimatorSquareTank):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a storage tank made of materials
    The storage tank is a cylinder, filled with the sample.
    Two detector locations sit at the tank mid-height, measured from the outer surface of the last layer:
    location 1 (contact) at det_standoff_distance (0.1 cm) and location 2 (handling) at
    handling_det_standoff_distance (30 cm).
    Responses by point-detector ID [rem/h]: '1' neutron and '2' photon at location 1, '5' neutron and '6' photon
    at location 2. The beta estimates '3' (contact, from '2') and '7' (handling, from '6') are zero unless the
    sample is bare (no layers).
    contact_dose and handling_dose hold the totals at the two locations. total_dose returns contact_dose,
    since location 1 has the same meaning as the single detector of DoseEstimatorGenericTank.
    """

    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        self.sample_h2: (None, float) = None  # half-height of the sample
        self.det_z: (None, float) = None  # z location of the detector
        self.cyl_r: (None, float) = None  # inner radius of the tank cylinder; h is calculated from V & r
        self.handling_det_x: (None, float) = None
        super().__init__(_o)  # Init DoseEstimator
        self.det_standoff_distance = 0.1  # [cm] contact dose
        self.handling_det_standoff_distance = 30.0  # [cm] handing dose
        self.box_a: float = self.handling_det_standoff_distance + 10.0

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self._validate_layers()
        self.sample_h2 = get_cyl_h(self.sample_volume, self.cyl_r) / 2.0
        tank_r: float = self.cyl_r
        tank_h2: float = self.sample_h2
        sample_z_max: float = self.sample_h2
        sample_z_min: float = - self.sample_h2

        mavric_output: str = self._deck_head(
            title_line=f'DoseEstimatorGenericTank, {self.sample_weight} g, layers {self.layers_thicknesses}',
            comp_extra="' helium 2 end\n")
        mavric_output += f'''
read geometry
global unit 1
    cylinder 1 {tank_r} 2p {tank_h2}
    media 1 1 1
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            tank_h2 += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_h2}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            xy_planes.append(tank_r)
            z_planes.append(tank_h2)
            z_planes.append(-tank_h2)
        self._set_computed_detectors(tank_r + self.det_standoff_distance,
                                     tank_r + self.handling_det_standoff_distance)  # next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 -{self.outer_body}
boundary 99999
end geometry

read definitions
     location 1
        position {self.det_x} 0 {self.det_z}
    end location
    location 2
        position {self.handling_det_x} 0 {self.det_z}
    end location
    response 1
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
'''

        mavric_output += self._deck_distributions()
        mavric_output += f'''
    gridGeometry 1
        title="Grid over the problem, location at +x"
        xLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a}
'        yLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a} 
'        zLinear {self.N_planes_box} {tank_h2} {tank_h2 + self.box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl}  {sample_z_min} {sample_z_max}
        xPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.det_x - self.planes_xy_around_det} {self.det_x + self.planes_xy_around_det} end
        xPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.handling_det_x - self.planes_xy_around_det} {self.handling_det_x + self.planes_xy_around_det} end
        yPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
        zPlanes {z_planes_str} {tank_h2 + self.box_a} {-tank_h2 - self.box_a} end
    end gridGeometry
end definitions
'''
        mavric_output += self._deck_tail(
            neutron_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}',
            photon_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}',
            handling_detector=True)
        return mavric_output

    def get_responses(self):
        """ Reads over the MAVRIC output and returns responses for rem/h doses
            Note, this version keys self.responses by point-detector ID and adds 'particle' and 'pid' fields.
            Raises RuntimeError when a detector is missing from the tally summary. """
        out_path: str = os.path.join(self.cwd, self.case_dir, self.MAVRIC_out_file_name)
        if not os.path.isfile(out_path):
            raise FileNotFoundError("Expected decayed sample MAVRIC output file: \n" + out_path)
        with open(out_path, 'r') as f:
            d: str = f.read()
        detectors: dict = parse_point_detector_responses(d)
        missing: list[str] = [det for det in ('1', '2', '5', '6') if det not in detectors]
        if missing:
            raise RuntimeError(f'Failed to parse point detectors {missing} from {out_path}')
        beta_contact: dict = self._beta_response(detectors['2'])
        beta_handling: dict = self._beta_response(detectors['6'])
        self.responses = {
            '1': detectors['1'], '2': detectors['2'], '3': {'particle': 'beta', 'pid': '', **beta_contact},
            '5': detectors['5'], '6': detectors['6'], '7': {'particle': 'beta', 'pid': '', **beta_handling},
        }
        if self.debug > 3:
            print(self.responses)

    @property
    def contact_dose(self) -> dict:
        """ Contact dose [rem/h] at location 1: neutron '1' + photon '2' + beta '3' (zero unless bare) """
        return combine_dose_responses(self.responses, ['1', '2', '3'], correlated=('2', '3'))

    @property
    def handling_dose(self) -> dict:
        """ Handling dose [rem/h] at location 2: neutron '5' + photon '6' + beta '7' (zero unless bare) """
        return combine_dose_responses(self.responses, ['5', '6', '7'], correlated=('6', '7'))

    @property
    def total_dose(self) -> dict:
        """ Dose at location 1, the same detector meaning as DoseEstimatorGenericTank; see contact_dose """
        return self.contact_dose

    def print_response(self):
        """ Prints dose responses """
        if self.responses:
            r1: dict = self.responses['1']  # Neutron dose, contact
            r2: dict = self.responses['2']  # Photon dose, contact
            r5: dict = self.responses['5']  # Neutron dose, handling (30cm)
            r6: dict = self.responses['6']  # Photon dose, handling (30cm)
            print(self.sample_weight, r1['value'], r1['stdev'], r2['value'], r2['stdev'], r5['value'], r5['stdev'],
                  r6['value'], r6['stdev'])


class MHATank(HandlingContactDoseEstimatorGenericTank):
    """ MHA tank, 1/16 pipe
    https://www.analytical-sales.com/product/ss-tubing-1-16-od-0-04-id-x-25/
    OR = 1.0/16.0 * 2.54 / 2.0 = 0.079375 cm
    IR = 0.04 * 2.54 / 2.0 = 0.05080 cm
    length = 10 ft = 304.8 cm
    """

    def __init__(self, _o: Origen = None):
        super().__init__(_o)
        """ This reads decayed sample information from the Origen object """
        self.sample_h2: float = 304.8 / 2.0  # half-height of the sample
        self.cyl_r: float = 0.05080  # inner radius of the tank cylinder; h is calculated from V & r
        self.det_z: (None, float) = None  # z location of the detector
        object.__setattr__(self, 'handling_det_x', None)  # set by mavric_deck(); not a manual det_x
        self.det_standoff_distance = 0.1  # [cm] contact dose
        self.handling_det_standoff_distance = 30.0  # [cm] handing dose
        self.box_a: float = self.handling_det_standoff_distance + 10.0
        self.source_multiplier: (None, float) = None  # Multiplier of ORIGEN source

    def scale_decayed_pipe_material(self, desired_weight: float = 1.0):
        """ Reads atom density and rho from F71 file """
        self.sample_density = get_rho_from_atom_density(self.decayed_atom_dens)
        if self.debug > 2:
            print(f'Initial sample density {self.sample_density} g/cm3')
        os.chdir(self.ORIGEN_dir)
        try:
            sample_gram_dict: dict = get_burned_nuclide_data(self.DECAYED_SAMPLE_F71_file_name,
                                                             self.DECAYED_SAMPLE_F71_position, 'gram')
        finally:
            os.chdir(self.cwd)
        self.sample_weight = 0.0
        for k,v in sample_gram_dict.items():
            self.sample_weight += v

        # input_volume: float = self.sample_weight / self.sample_density
        self.sample_volume = np.pi * self.cyl_r ** 2 * (2.0 * self.sample_h2)  # V = pi r^2 h
        adens_scaling: float = desired_weight / (self.sample_volume * self.sample_density)
        self.decayed_atom_dens = scale_adens(self.decayed_atom_dens, adens_scaling)
        self.sample_density = get_rho_from_atom_density(self.decayed_atom_dens)
        self.source_multiplier = desired_weight / self.sample_weight
        if self.debug > 2:
            print(f'Atomic density scaling factor {adens_scaling}, source scaling factor {self.source_multiplier}')
            print(f'Sample density {self.sample_density} g/cm3, volume {self.sample_volume} cm3')
            nicely_print_atom_dens(self.decayed_atom_dens)

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self._validate_layers()
        tank_r: float = self.cyl_r
        tank_h2: float = self.sample_h2
        sample_z_max: float = self.sample_h2
        sample_z_min: float = - self.sample_h2

        mavric_output: str = self._deck_head(
            title_line=f'DoseEstimatorGenericTank, {self.sample_weight} g, layers {self.layers_thicknesses}',
            comp_pre="' THIS IS BASICALLY MEANINGLESS, SINCE THIS IS A GAS AND WE ARE USING A SOURCE FROM F71 FILE \n",
            comp_extra="' helium 2 end\n")
        mavric_output += f'''
read geometry
global unit 1
    cylinder 1 {tank_r} 2p {tank_h2}
    media 1 1 1
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            tank_h2 += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_h2}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            xy_planes.append(tank_r)
            z_planes.append(tank_h2)
            z_planes.append(-tank_h2)
        self._set_computed_detectors(tank_r + self.det_standoff_distance,
                                     tank_r + self.handling_det_standoff_distance)  # next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 -{self.outer_body}
boundary 99999
end geometry

read definitions
     location 1
        position {self.det_x} 0 {self.det_z}
    end location
    location 2
        position {self.handling_det_x} 0 {self.det_z}
    end location
    response 1
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
'''

        mavric_output += self._deck_distributions()
        mavric_output += f'''
    gridGeometry 1
        title="Grid over the problem, location at +x"
        xLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a}
'        yLinear {self.N_planes_box} {tank_r} {tank_r + self.box_a} 
'        zLinear {self.N_planes_box} {tank_h2} {tank_h2 + self.box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl}  {sample_z_min} {sample_z_max}
        xPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.det_x - self.planes_xy_around_det} {self.det_x + self.planes_xy_around_det} end
        xPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.handling_det_x - self.planes_xy_around_det} {self.handling_det_x + self.planes_xy_around_det} end
        yPlanes {xy_planes_str} {tank_r + self.box_a} {-tank_r - self.box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
        zPlanes {z_planes_str} {tank_h2 + self.box_a} {-tank_h2 - self.box_a} end
    end gridGeometry
end definitions
'''
        mavric_output += self._deck_tail(
            neutron_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}',
            photon_source_cylinder=f'{self.cyl_r} {sample_z_max} {sample_z_min}',
            multiplier=self.source_multiplier,
            handling_detector=True)
        return mavric_output


if __name__ == "__main__":
    pass  # main()
