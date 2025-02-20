#!/bin/env python3
"""
Handling dose [rem/h] in SCALE
Ondrej Chvala <ochvala@utexas.edu>
"""
import os
import shutil
import re
import subprocess
import math
from datetime import datetime
import numpy as np
from bisect import bisect_left
from sample_decay_dose.read_opus import integrate_opus

NOW: str = datetime.now().replace(microsecond=0).isoformat()

SCALE_bin_path: str = os.getenv('SCALE_BIN', '/opt/scale6.3.1/bin/')
ATOM_DENS_MINIMUM: float = 1e-60
MAVRIC_NG_XSLIB: str = 'v7.1-28n19g'

# Predefined lists of atom densities for convenience. Check and make yours :)
# https://www.sandmeyersteel.com/316H.html
ADENS_SS316H_HOT: dict = {'c': 0.0002384, 'n-14': 0.000254609, 'n-15': 9.30164e-07, 'al-27': 5.30544e-05,
                          'si-28': 1.56475e-05, 'si-29': 7.94907e-07, 'si-30': 5.24622e-07, 'p-31': 6.93324e-05,
                          's-32': 4.23788e-05, 's-33': 3.34604e-07, 's-34': 1.89609e-06, 's-36': 4.46139e-09,
                          'ti-46': 3.28985e-06, 'ti-47': 2.96684e-06, 'ti-48': 2.93973e-05, 'ti-49': 2.15734e-06,
                          'ti-50': 2.06562e-06, 'cr-50': 0.000677912, 'cr-52': 0.0130729, 'cr-53': 0.00148236,
                          'cr-54': 0.00036899, 'mn-55': 1.73977e-05, 'fe-54': 0.00390032, 'fe-56': 0.0612267,
                          'fe-57': 0.00141399, 'fe-58': 0.000188176, 'ni-58': 6.64309e-05, 'ni-60': 2.55891e-05,
                          'ni-61': 1.11234e-06, 'ni-62': 3.54662e-06, 'ni-64': 9.03221e-07, 'mo-92': 0.00018364,
                          'mo-94': 0.00011476, 'mo-95': 0.00019769, 'mo-96': 0.000207388, 'mo-97': 0.000118863,
                          'mo-98': 0.000300762, 'mo-100': 0.00012023}
ADENS_SS316H_COLD: dict = {'c': 0.000311167, 'si-28': 0.00153402, 'si-29': 7.79293e-05, 'si-30': 5.14317e-05,
                           'p-31': 6.78722e-05, 's-32': 4.15187e-05, 's-33': 3.27814e-07, 's-34': 1.85761e-06,
                           's-36': 4.37085e-09, 'cr-50': 0.000663653, 'cr-52': 0.0127979, 'cr-53': 0.00145118,
                           'cr-54': 0.000361229, 'mn-55': 0.00170071, 'fe-54': 0.0031951, 'fe-56': 0.0501563,
                           'fe-57': 0.00115833, 'fe-58': 0.000154152, 'co-59': 7.93501e-05, 'ni-58': 0.00650228,
                           'ni-60': 0.00250467, 'ni-61': 0.000108876, 'ni-62': 0.000347145, 'ni-64': 8.84075e-05,
                           'mo-92': 0.000179806, 'mo-94': 0.000112364, 'mo-95': 0.000193563, 'mo-96': 0.000203058,
                           'mo-97': 0.000116381, 'mo-98': 0.000294483, 'mo-100': 0.00011772}
ADENS_HDPE_COLD: dict = {'c-12': 3.992647e-02, 'c-13': 4.318339e-04, 'h-poly': 8.071660e-02}
ADENS_HELIUM_COLD: dict = {'he-3': 2.680585e-11 , 'he-4': 2.693156e-05}
ADENS_HELIUM_HOT: dict = {'he-3': 8.27507e-12, 'he-4': 8.27506e-06}
ADENS_CONCRETE_COLD: dict = {'h-1': 0.01373939, 'o-16': 0.04606872, 'na-23': 0.001747024, 'al-27': 0.001745235,
                             'si-28': 0.01532717, 'si-29': 0.0007786319, 'si-30': 0.0005138805, 'ca-40': 0.001474023,
                             'ca-42': 9.837866e-06, 'ca-43': 2.052723e-06, 'ca-44': 3.171838e-05, 'ca-46': 6.082143e-08,
                             'ca-48': 2.843401e-06, 'fe-54': 2.029158e-05, 'fe-56': 0.0003185345, 'fe-57': 7.356349e-06,
                             'fe-58': 9.789951e-07, 'co-59': 2.350275e-07, 'eu-151': 4.357685e-08,
                             'eu-153': 4.756903e-08}
ADENS_KAOWOOL_COLD: dict = {'b-10': 4.39982e-07, 'b-11': 1.77098e-06, 'o-16': 0.00300423, 'o-17': 1.14439e-06,
                            'o-18': 6.17368e-06, 'al-27': 0.000850507, 'si-28': 0.000770832, 'si-29': 3.91588e-05,
                            'si-30': 2.5844e-05, 'ca-40': 1.66602e-06, 'ca-42': 1.11193e-08, 'ca-43': 2.32009e-09,
                            'ca-44': 3.58497e-08, 'ca-46': 6.87435e-11, 'ca-48': 3.21376e-09, 'ti-46': 1.69203e-06,
                            'ti-47': 1.5259e-06, 'ti-48': 1.51195e-05, 'ti-49': 1.10956e-06, 'ti-50': 1.06239e-06,
                            'fe-54': 7.05374e-07, 'fe-56': 1.10729e-05, 'fe-57': 2.55721e-07, 'fe-58': 3.40317e-08}
ADENS_LEAD_COLD: dict = {'pb-204': 4.615515e-04, 'pb-206': 7.945277e-03, 'pb-207': 7.285918e-03, 'pb-208': 1.727521e-02}


def nicely_print_atom_dens(adens: dict, n_top_nuc: int = 20, n_per_row: int = 5):
    """ Prints atom density of top nuclides """
    nuclide_list: list = list(adens.items())[:n_top_nuc]
    i_nuc: int = 0
    output: str = ''
    while i_nuc < min(n_top_nuc, len(nuclide_list)):
        ele: str = nuclide_list[i_nuc][0]
        ade: float = nuclide_list[i_nuc][1]
        output += f'{ele:>7s} {ade:.4e} '
        i_nuc += 1
        if i_nuc % n_per_row == 0:
            output += '\n'
    print(output.rstrip())


def get_rho_from_atom_density(adens: dict) -> float:
    """ Calculate mass density [g/cm3] from atom densities [at/b-cm] """
    from sample_decay_dose.isotopes import rel_iso_mass, m_Da
    my_rho = 0.0
    for k, v in adens.items():
        if v > 0.0:
            iso_name = k.lower()
            iso_name = re.sub('m$', '', iso_name)  # strip "m"
            iso_name = re.sub('([a-zA-Z]+)(\d+)', '\g<1>-\g<2>', iso_name)  # add dash
            my_rho += rel_iso_mass[iso_name] * v * m_Da * 1e24
    return my_rho


def get_f71_positions_index(f71file: str) -> dict:
    """ Read info of SCALE's F71 file """
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=info", f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    f71_idx = {}
    skip_data = ['pos', '(-)']  # sip records starting with this
    end_data = 'state definition present'  # stop reading after reaching this
    for line in output:
        data = line.split()
        if data[0].strip() in skip_data:
            continue
        if line.count(end_data) > 0:
            break
        # print(data)
        f71_idx[int(data[0])] = {
            'time': data[1],
            'power': data[2],
            'flux': data[3],
            'fluence': data[4],
            'energy': data[5],
            'initialhm': data[6],
            'libpos': data[7],
            'case': data[8],
            'step': data[9],
            'DCGNAB': data[10]
        }
    return f71_idx


def get_burned_material_atom_dens(f71file: str, position: int) -> dict:
    """ Read atom density of nuclides from SCALE's F71 file """
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'",
         f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    densities = {}  # densities[nuclide] = (density at position of f71 file)
    skip = ["case", "step", "time", "power", "flux", "volume"]
    regexp = re.compile(r"(?P<elem>[a-zA-Z]+)(?P<num>\d+)(?P<meta>m)?")
    for line in output:
        data = line.split(',')
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            dummy = re.search(regexp, data[0].strip())
            elem = dummy.group("elem").lower()  # convert to all lower cases
            num = int(dummy.group("num"))  # to cut off leading zeros
            if dummy.group("meta"):
                nuclide = elem + "-" + str(num) + "m"  # for metastable isotopes
            else:
                nuclide = elem + "-" + str(num)
            if float(data[position]) > ATOM_DENS_MINIMUM:
                # The [x] here is what causes the code to return only the densities at position x of the f71 file
                densities[nuclide] = float(data[position])
    sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
    return sorted_densities


def get_burned_material_total_mass_dens(f71file: str, position: int) -> float:
    """ Read mass density of nuclides from SCALE's F71 file, calculate total \rho [g/cm^3] """
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=gper", "-idform='{:Ee}{:AAA}{:m}'",
         f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    my_rho: float = 0
    skip = ["case", "step", "time", "power", "flux", "volume"]
    for line in output:
        data = line.split(',')
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            my_rho = my_rho + float(data[position])
    return my_rho


def get_F33_num_sets(f33file: str) -> int:
    """ Returns the number of numSets in SCALE F33 file """
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "info", f33file], capture_output=True)
    output = output.stdout.decode().split("\n")
    for line in output:
        data = line.split()
        if data[0] in f33file:
            return int(data[2])
    return -1


def get_cyl_r(cyl_volume: float) -> float:
    """ Radius of a square cylinder from its volume
        V = pi r^2 h = pi r^2 2 r = 2 pi r^3
        r = (V / 2 pi) ^ 1/3
    """
    return (cyl_volume / (2.0 * math.pi)) ** (1.0 / 3.0)


def get_cyl_r_4_1(cyl_volume: float) -> float:
    """ Radius of a 4:1 H:R cylinder from its volume
        V = pi r^2 h = pi r^2 4 r = 4 pi r^3
        r = (V / 4 pi) ^ 1/3
        h_cyl = 4 * r
    """
    return (cyl_volume / (4.0 * math.pi)) ** (1.0 / 3.0)


def get_fill_height_4_1(fill_volume: float, cyl_volume: float) -> float:
    """ Fill height of 4:1 H:R cylinder
        V_cyl = pi r^2 h_cyl = pi r^2 4 r = 4 pi r^3
        r = (V_cyl / 4 pi) ^ 1/3;  h_cyl = 4 * r
        V_fill = pi r^2 h_fill
        h_fill = V_fill / (pi r^2)
    """
    if not fill_volume < cyl_volume:
        raise ValueError("Fill volume is larger than the cylinder")
    r: float = (cyl_volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    return fill_volume / (math.pi * r ** 2)


def get_cyl_h(cyl_volume: float, cyl_r: float) -> float:
    """ Radius of a square cylinder from its volume
        V = pi r^2 h
        h = V / (pi r^2)
    """
    if not cyl_r > 0:
        raise ValueError(f"Cylinder radius has to be positive, {cyl_r}")
    return cyl_volume / math.pi / cyl_r ** 2


def run_scale(deck_file: str):
    """ Run a SCALE deck """
    scale_out = subprocess.run([f"{SCALE_bin_path}/scalerte", "-m", deck_file], capture_output=True)
    scale_out = scale_out.stdout.decode().split("\n")
    if scale_out.count('Error') > 0:
        print('Failed run: ', deck_file)
        return False
    else:
        print('OK run: ', deck_file)
        return True


def atom_dens_for_origen(adens: dict) -> str:
    """ Print atom densities in ORIGEN input format """
    output = ''
    for k, v in adens.items():
        if v > 0.0:
            output += f'{k} = {v} \n'
    return output


def atom_dens_for_mavric(adens: dict, mix_number: int = 1, tempK: float = 873.0) -> str:
    """ Print atom densities in SCALE CSAS/MAVRIC input format """
    output = ''
    for k, v in adens.items():
        if v > 0.0:
            k = re.sub('m$', '', k)
            output += f'{k} {mix_number} 0 {v} {tempK} end\n'
    return output


def extract_flux_values(scale_output_file):
    """ Extract total flux values for mixtures """
    flux_data = {}
    with open(scale_output_file, 'r') as file:
        for line in file:
            match = re.match(r'\s*(\d+)\s+[\d.]+\s+[\d.]+\s+N/A\s+N/A\s+[\d.e+-]+\s+([\d.e+-]+)', line)
            if match:
                mixture = int(match.group(1))
                total_flux = float(match.group(2))
                flux_data[mixture] = total_flux
    return flux_data


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
        self.SAMPLE_F71_position: int = 12  # Sample decay steps
        self.SAMPLE_DECAY_days: float = 30.0  # Sample decay time [days]
        self.sample_weight: float = np.nan  # Mass of the sample [g]
        self.sample_density: float = np.nan  # Mass density of the sample [g/cm3]
        self.sample_volume: float = np.nan  # Sample volume [cm3]

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
        return beta_spectrum_integral / gamma_spectrum_integral

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
        print(f'--> Closest F71 position found at slot {pos_number}, {times[pos_number]} seconds, '
              f'{times[pos_number]/(60.0*60.0*24)} days.')

        self.BURNED_MATERIAL_F71_position = pos_number

    def read_burned_material(self):
        """ Reads atom density and rho from F71 file """
        if self.debug > 0:
            print(f'ORIGEN: reading nuclides from {self.BURNED_MATERIAL_F71_file_name}, '
                  f'position {self.BURNED_MATERIAL_F71_position}')
        self.sample_density = get_burned_material_total_mass_dens(self.BURNED_MATERIAL_F71_file_name,
                                                                  self.BURNED_MATERIAL_F71_position)
        self.sample_volume = self.sample_weight / self.sample_density
        self.burned_atom_dens = get_burned_material_atom_dens(self.BURNED_MATERIAL_F71_file_name,
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
        os.chdir(self.case_dir)

        with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
            f.write(atom_dens_for_origen(self.burned_atom_dens))

        with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.origen_deck())

        if self.debug > 0:
            print(f'ORIGEN: decaying sample for {self.SAMPLE_DECAY_days} days')
            print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

        self.decayed_atom_dens = get_burned_material_atom_dens(self.SAMPLE_F71_file_name,
                                                               self.SAMPLE_F71_position)
        os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def origen_deck(self) -> str:
        """ Sample decay Origen deck """
        time_interp_steps = self.SAMPLE_F71_position - 3
        if time_interp_steps < 1:
            raise ValueError("Too few time steps")

        origen_output = f'''
=shell
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_Origen} .
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
<{self.SAMPLE_ATOM_DENS_file_name_Origen}
        ]
        units=ATOMS-PER-BARN-CM
        volume={self.sample_volume}
    }}
    time {{
        units=DAYS
        t=[{time_interp_steps}L 0.0001 {self.SAMPLE_DECAY_days}]
        start=0
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Neutrons'
typarams=nspectrum
units=intensity
time=days
tmin={self.SAMPLE_DECAY_days}
tmax={self.SAMPLE_DECAY_days}
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Gamma'
typarams=gspectrum
units=intensity
time=days
tmin={self.SAMPLE_DECAY_days}
tmax={self.SAMPLE_DECAY_days}
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Beta'
typarams=bspectrum
units=intensity
time=days
tmin={self.SAMPLE_DECAY_days}
tmax={self.SAMPLE_DECAY_days}
end
'''
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
        os.chdir(self.case_dir)

        self.sample_density = get_rho_from_atom_density(atom_dens)
        with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write at-dens input for Origen irradiation
            f.write(atom_dens_for_origen(atom_dens))
        os.chdir(self.cwd)

    def run_irradiate_decay_sample(self):
        """  Writes Origen input file, runs Origen to irradiate and decay it, and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir + '/' + self.SAMPLE_ATOM_DENS_file_name_Origen):
            raise FileNotFoundError("Write atom density for Origen first")
        os.chdir(self.case_dir)

        with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.origen_deck())

        if self.debug > 0:
            print(f'ORIGEN: burning sample for {self.irradiate_days} days at {self.irradiate_flux} n/cm2/s, '
                  f'then decaying for {self.SAMPLE_DECAY_days} days')
            print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

        self.decayed_atom_dens = get_burned_material_atom_dens(self.SAMPLE_F71_file_name,
                                                               self.SAMPLE_F71_position)
        os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def read_irradiated_material_density(self):
        """ Reads rho from irradiated F71 file """
        self.sample_density = get_burned_material_total_mass_dens(self.irradiate_F71_file_name, 1)
        self.sample_volume = self.sample_weight / self.sample_density

    def origen_deck(self) -> str:
        """ Sample irradiation and decay Origen deck """
        self.sample_volume = self.sample_weight / self.sample_density
        time_interp_steps = self.SAMPLE_F71_position - 3
        if time_interp_steps < 1:
            raise ValueError("Too few time steps")
        # The self.F33_file_name potentially includes path to that file.
        # It gets copied to temp directory, where it is just F33_file
        F33_path, F33_file = os.path.split(self.F33_file_name)
        origen_output = f'''
=shell
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_Origen} .
cp -r {self.cwd}/{self.F33_file_name} .
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
        t=[{self.irradiate_steps - 2}L 0.001 {self.irradiate_days}]
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
        t=[{time_interp_steps}L 0.001 {self.SAMPLE_DECAY_days}]
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Neutrons'
typarams=nspectrum
units=intensity
'time=days
'tmin={self.SAMPLE_DECAY_days}
'tmax={self.SAMPLE_DECAY_days}
npos={self.SAMPLE_F71_position} end
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Gamma'
typarams=gspectrum
units=intensity
'time=days
'tmin={self.SAMPLE_DECAY_days}
'tmax={self.SAMPLE_DECAY_days}
npos={self.SAMPLE_F71_position} end
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Beta'
typarams=bspectrum
units=intensity
'time=days
'tmin={self.SAMPLE_DECAY_days}
'tmax={self.SAMPLE_DECAY_days}
npos={self.SAMPLE_F71_position} end
'''
        return origen_output


def read_cvs_atom_dens(csv_file: str, volume: float = 1.0) -> dict:
    """ Reads CVS file with a nuclide and number of atoms per row:
    <nuc1>, <# of atoms>
    <nuc2>, <# of atoms>
    ...
        Returns atom density
        volume is in cm^3
    """
    import csv
    my_atom_density: dict = {}
    volume_barn_cm: float = volume * 1e24  # Volume [cm^3] -> [barn-cm]
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            my_atoms: float = float(row[1])
            if my_atoms > 0:
                iso_name = row[0].lower()
             #   iso_name = re.sub('([a-zA-Z]+)(\d+)', '\g<1>-\g<2>', iso_name)  # add dash
                my_atom_density[iso_name] = my_atoms / volume_barn_cm
    return my_atom_density


class OrigenDecayBox(Origen):
    """ Origen decay from a simple dict of atom density and volume [cm] """

    def __init__(self, _adens: (None, dict) = None, _vol: float = 0):
        Origen.__init__(self)
        self.SAMPLE_ATOM_DENSITY: (None, dict) = _adens
        self.sample_volume = _vol
        self.case_dir: str = f'_decaybox_{_vol:.5}_g'  # Directory to run the case

    def set_decay_days(self, decay_days: float = 30.0):
        """ Use this to change decay time, as it also updates the case directory """
        self.SAMPLE_DECAY_days = decay_days

    def write_atom_dens(self):
        """ Writes atom density of the sample to decay """
        if self.SAMPLE_ATOM_DENSITY is None or self.sample_volume == 0:
            raise ValueError(f'Adens is None or Volume is zero: {self.SAMPLE_ATOM_DENSITY}, {self.sample_volume}')
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        self.sample_density = get_rho_from_atom_density(self.SAMPLE_ATOM_DENSITY)
        with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write at-dens input for Origen irradiation
            f.write(atom_dens_for_origen(self.SAMPLE_ATOM_DENSITY))
        os.chdir(self.cwd)

    def run_decay_sample(self):
        """  Writes Origen input file, runs Origen to decay it, and Opus to plot spectra.
        Finally, it reads atom density of the decayed sample, used later as a mixture for Mavric.
        """
        if not os.path.exists(self.case_dir + '/' + self.SAMPLE_ATOM_DENS_file_name_Origen):
            raise FileNotFoundError("Write atom density for Origen first")
        os.chdir(self.case_dir)

        with open(self.ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.origen_deck())

        if self.debug > 0:
            print(f'ORIGEN: decaying for {self.SAMPLE_DECAY_days} days')
            print(f"Running case: {self.case_dir}/{self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

        print(self.SAMPLE_F71_file_name,
                                                               self.SAMPLE_F71_position)
        self.decayed_atom_dens = get_burned_material_atom_dens(self.SAMPLE_F71_file_name,
                                                               self.SAMPLE_F71_position)
        os.chdir(self.cwd)
        if self.debug > 2:
            # print(list(self.decayed_atom_dens.items())[:25])
            nicely_print_atom_dens(self.decayed_atom_dens)

    def origen_deck(self) -> str:
        """ Sample irradiation and decay Origen deck """
        time_interp_steps = self.SAMPLE_F71_position - 3
        if time_interp_steps < 1:
            raise ValueError("Too few time steps")
        origen_output = f'''
=shell
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_Origen} .
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
        t=[{time_interp_steps}L 0.001 {self.SAMPLE_DECAY_days}]
    }}
    save {{
        file="{self.SAMPLE_F71_file_name}"
    }}
}}
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Neutrons'
typarams=nspectrum
units=intensity
'time=days
'tmin={self.SAMPLE_DECAY_days}
'tmax={self.SAMPLE_DECAY_days}
npos={self.SAMPLE_F71_position} end
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Gamma'
typarams=gspectrum
units=intensity
'time=days
'tmin={self.SAMPLE_DECAY_days}
'tmax={self.SAMPLE_DECAY_days}
npos={self.SAMPLE_F71_position} end
end

=opus
data='{self.SAMPLE_F71_file_name}'
title='Beta'
typarams=bspectrum
units=intensity
'time=days
'tmin={self.SAMPLE_DECAY_days}
'tmax={self.SAMPLE_DECAY_days}
npos={self.SAMPLE_F71_position} end
'''
        return origen_output


class DoseEstimator:
    """ MAVRIC calculation of rem/h doses from the decayed sample """

    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        self.debug: int = 3  # Debugging flag
        self.MAVRIC_input_file_name: str = 'my_dose.inp'
        self.SAMPLE_ATOM_DENS_file_name_MAVRIC: str = 'my_sample_atom_dens_mavric.inp'
        self.responses: dict = {}  # Dose responses 1: neutron, 2: gamma, 3: beta
        self.det_x: float = 30.0  # Detector distance [cm]
        self.N_planes_box: int = 5  # Planes per box
        self.N_planes_cyl: int = 8  # Planes per cylinder
        self.histories_per_batch: int = 100000  # Monaco hist per batch
        self.batches: int = 10  # Monaco number of batches in total
        self.box_a: float = np.nan
        self.cyl_r: float = np.nan
        self.sample_temperature_K: float = 873.0  # Sample temperature [K]
        self.decayed_atom_dens: dict = {}  # Atom density of the decayed sample
        self.beta_over_gamma: float = 0  # Beta over gamma spectral ratio
        if _o is not None:
            self.sample_weight: float = _o.sample_weight  # Mass of the sample [g]
            self.sample_density: float = _o.sample_density  # Mass density of the sample [g/cm3]
            self.sample_volume: float = _o.sample_volume  # Sample volume [cm3]
            self.DECAYED_SAMPLE_F71_file_name: str = _o.SAMPLE_F71_file_name
            self.DECAYED_SAMPLE_F71_position: int = _o.SAMPLE_F71_position
            self.DECAYED_SAMPLE_days: float = _o.SAMPLE_DECAY_days  # Sample decay time [days]
            self.decayed_atom_dens: dict = _o.decayed_atom_dens  # Atom density of the decayed sample
            self.beta_over_gamma: float = _o.get_beta_to_gamma()  # Beta over gamma spectral ratio
            self.neutron_intensity: float = _o.get_neutron_integral()  # Integral of neutron spectra
            self.ORIGEN_dir: str = _o.case_dir  # Directory to run the case
            self.case_dir: str = self.ORIGEN_dir + '_MAVRIC'
            self.cwd: str = _o.cwd  # Current running directory

    @property
    def MAVRIC_out_file_name(self) -> str:
        return self.MAVRIC_input_file_name.replace('inp', 'out')

    def run_mavric(self):
        """ Writes Mavric inputs and runs the case """
        if not os.path.isfile(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name):
            raise FileNotFoundError("Expected decayed sample F71 file: \n" +
                                    self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name)
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        shutil.copy2(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name,
                     self.cwd + '/' + self.case_dir)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(self.SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
            f.write(atom_dens_for_mavric(self.decayed_atom_dens, 1, self.sample_temperature_K))

        with open(self.MAVRIC_input_file_name, 'w') as f:  # write MAVRIC input deck
            f.write(self.mavric_deck())

        if self.debug > 0:
            print(f"MAVRIC: running case {self.case_dir}/{self.MAVRIC_input_file_name}")
        run_scale(self.MAVRIC_input_file_name)
        os.chdir(self.cwd)

    def get_responses(self):
        """ Reads over the MAVRIC output and returns responses for rem/h doses """
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name):
            raise FileNotFoundError("Expected decayed sample MAVRIC output file: \n" +
                                    self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)

        tally_sep: str = 'Final Tally Results Summary'
        is_in_tally: bool = False

        s: list[str] = ['']
        with open(self.MAVRIC_out_file_name, 'r') as f:
            for line in f.read().splitlines():
                if line.count(tally_sep) > 0:
                    is_in_tally = True
                if is_in_tally:
                    if line.count('response') == 1:
                        s = line.split()
                        if float(s[2]) == 0:
                            self.responses[s[1]] = {'value': float(s[2]), 'stdev': 0.0}
                        else:
                            self.responses[s[1]] = {'value': float(s[2]), 'stdev': float(s[3])}

        if not s:  # This should handle MONACO crashes
            s[1] = -1.0
            s[2] = -1.0
            s[3] = -1.0
        self.responses['3'] = {'value': self.beta_over_gamma * float(s[2]), 'stdev': self.beta_over_gamma * float(s[3])}

        os.chdir(self.cwd)
        if self.debug > 3:
            print(self.responses)

    def print_response(self):
        """ Prints dose responses """
        if self.responses is not {}:
            r1: dict = self.responses['1']
            r2: dict = self.responses['2']
            r3: dict = self.responses['3']
            print(self.sample_weight, r1['value'], r1['stdev'], r2['value'], r2['stdev'], r3['value'], r3['stdev'])

    @property
    def total_dose(self) -> dict:
        """ Return total dose, which is a sum of all responses """
        d: dict = {'value': 0, 'stdev': -1}
        relsig: float = 0
        for k, v in self.responses.items():
            d['value'] += v['value']
            if v['value'] > 0:
                relsig += (v['stdev'] / v['value']) ** 2
        d['stdev'] = d['value'] * math.sqrt(relsig)
        return d

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self.cyl_r = get_cyl_r(self.sample_volume)
        self.box_a: float = self.det_x + 10.0  # Problem box distance [cm]
        adjoint_flux_file = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
end

=mavric parm=(   )
{NOW} Sample dose, {self.sample_weight} g, at x={self.det_x} cm
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=1
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
<{self.SAMPLE_ATOM_DENS_file_name_MAVRIC}
end comp

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

        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    distribution 1
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, neutrons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 1 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution'''

        mavric_output += f'''
    distribution 2
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 5 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution

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

read sources'''
        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cylinder {self.cyl_r} {self.cyl_r} -{self.cyl_r}
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cylinder {self.cyl_r} {self.cyl_r} -{self.cyl_r}
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1
'   adjointFluxes="{adjoint_flux_file}"
   adjointSource 1
        locationID=1
        responseID=1
   end adjointSource
   adjointSource 2
        locationID=1
        responseID=2
   end adjointSource
end importanceMap

read tallies
    pointDetector 1
        title="neutron detector"
        neutron
        locationID=1
        responseID=1
    end pointDetector
    pointDetector 2
        title="photon detector"
        photon
        locationID=1
        responseID=2
    end pointDetector
end tallies

end data
end
'''
        return mavric_output


class DoseEstimatorSquareTank(DoseEstimator):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a square tank made of materials """

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
        if len(self.layers_thicknesses) != len(self.layers_mats):
            raise ValueError("There needs to be the same amount of layers in both lists.")

    def run_mavric(self):
        """ Writes Mavric inputs and runs the case """
        self.case_dir += "".join([f'_{s:.3f}' for s in self.layers_thicknesses])  # unique IDs for parallel run
        if not os.path.isfile(self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name):
            raise FileNotFoundError("Expected decayed sample F71 file: \n" +
                                    self.cwd + '/' + self.ORIGEN_dir + '/' + self.DECAYED_SAMPLE_F71_file_name)
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

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
        run_scale(self.MAVRIC_input_file_name)
        os.chdir(self.cwd)

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self.cyl_r = get_cyl_r(self.sample_volume)
        tank_r: float = self.cyl_r  # current outer layer [cm]
        adjoint_flux_file = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
end

=mavric parm=(   )
{NOW} DoseEstimatorSquareTank, {self.sample_weight} g, layers {self.layers_thicknesses}
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=1
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
<{self.SAMPLE_ATOM_DENS_file_name_MAVRIC}
end comp

read geometry
global unit 1
    cylinder 1 {self.cyl_r} 2p {self.cyl_r}
    media 1 1 1
'''
        x_planes: list[float] = []  # list of cylinder boundaries for gridgeometry
        k: int = 0
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_r}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            x_planes.append(tank_r)
        self.det_x = tank_r + self.det_standoff_distance  # Detector is next to the tank
        self.box_a += tank_r
        x_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in x_planes])

        mavric_output += f'''
    cuboid 99999  6p {self.box_a}
    media 0 1 99999 -{k + 2}
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

        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    distribution 1
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, neutrons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 1 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution'''

        mavric_output += f'''
    distribution 2
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 5 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution

    gridGeometry 1
        title="Grid over the problem, location at +x"
        xLinear {self.N_planes_box} {self.cyl_r} {self.box_a}
'        yLinear {self.N_planes_box} -{self.box_a} {self.box_a}
'        zLinear {self.N_planes_box} -{self.box_a} {self.box_a}
        xLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        yLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        zLinear {self.N_planes_cyl} -{self.cyl_r} {self.cyl_r}
        xPlanes {x_planes_str} -{self.box_a} {self.box_a} {self.det_x + self.planes_xy_around_det} {self.det_x - self.planes_xy_around_det} end
        yPlanes {x_planes_str} -{self.box_a} {self.box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
        zPlanes {x_planes_str} -{self.box_a} {self.box_a} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
    end gridGeometry
end definitions

read sources'''
        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cylinder {self.cyl_r} {self.cyl_r} -{self.cyl_r}
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cylinder {self.cyl_r} {self.cyl_r} -{self.cyl_r}
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1
'   adjointFluxes="{adjoint_flux_file}"
   adjointSource 1
        locationID=1
        responseID=1
   end adjointSource
   adjointSource 2
        locationID=1
        responseID=2
   end adjointSource
end importanceMap

read tallies
    pointDetector 1
        title="neutron detector"
        neutron
        locationID=1
        responseID=1
    end pointDetector
    pointDetector 2
        title="photon detector"
        photon
        locationID=1
        responseID=2
    end pointDetector
end tallies

end data
end
'''
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
        tank_inner_volume: float = self.sample_volume + self.sample_volume * self.plenum_volume_fraction
        self.cyl_r = get_cyl_r_4_1(tank_inner_volume)
        adjoint_flux_file = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        tank_r: float = self.cyl_r  # current outer layer radius [cm]
        tank_h2: float = 2.0 * self.cyl_r  # current outer layer half-height [cm]
        self.sample_h2 = get_fill_height_4_1(self.sample_volume, tank_inner_volume) / 2.0
        self.sample_offset_z = 2.0 * (tank_h2 - self.sample_h2)
        sample_z_max: float = tank_h2 - self.sample_offset_z
        sample_z_min: float = -tank_h2

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
end

=mavric parm=(   )
{NOW} DoseEstimatorStorageTank, {self.sample_weight} g, layers {self.layers_thicknesses}, plenum fraction {self.plenum_volume_fraction}
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=1
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
<{self.SAMPLE_ATOM_DENS_file_name_MAVRIC}
helium 2 end
end comp

read geometry
global unit 1
    cylinder 1 {tank_r} {sample_z_max} {sample_z_min}
    cylinder 2 {tank_r} {tank_h2} {sample_z_max} 
    media 1 1 1
    media 2 1 2
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        k: int = 0
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
        self.det_x = tank_r + self.det_standoff_distance  # Detector is next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 -{k + 3}
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

        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    distribution 1
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, neutrons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 1 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution'''

        mavric_output += f'''
    distribution 2
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 5 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution

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

read sources'''
        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cylinder {self.cyl_r} {sample_z_max} {sample_z_min}
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cylinder {self.cyl_r} {sample_z_max} {sample_z_min}
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1
'   adjointFluxes="{adjoint_flux_file}"
   adjointSource 1
        locationID=1
        responseID=1
   end adjointSource
   adjointSource 2
        locationID=1
        responseID=2
   end adjointSource
end importanceMap

read tallies
    pointDetector 1
        title="neutron detector"
        neutron
        locationID=1
        responseID=1
    end pointDetector
    pointDetector 2
        title="photon detector"
        photon
        locationID=1
        responseID=2
    end pointDetector
end tallies

end data
end
'''
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
        adjoint_flux_file: str = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        self.sample_h2 = get_cyl_h(self.sample_volume, self.cyl_r) / 2.0
        tank_r: float = self.cyl_r
        tank_h2: float = self.sample_h2
        sample_z_max: float = self.sample_h2
        sample_z_min: float = - self.sample_h2

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
end

=mavric parm=(   )
{NOW} DoseEstimatorGenericTank, {self.sample_weight} g, layers {self.layers_thicknesses}
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=1
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
<{self.SAMPLE_ATOM_DENS_file_name_MAVRIC}
' helium 2 end
end comp

read geometry
global unit 1
    cylinder 1 {tank_r} 2p {tank_h2}
    media 1 1 1
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        k: int = 0
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            tank_h2 += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_h2}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            xy_planes.append(tank_r)
            z_planes.append(tank_h2)
            z_planes.append(-tank_h2)
        self.det_x = tank_r + self.det_standoff_distance  # Detector is next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 -{k + 2}
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

        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    distribution 1
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, neutrons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 1 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution'''

        mavric_output += f'''
    distribution 2
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 5 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution

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

read sources'''
        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cylinder {self.cyl_r} {sample_z_max} {sample_z_min}
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cylinder {self.cyl_r} {sample_z_max} {sample_z_min}
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1
'   adjointFluxes="{adjoint_flux_file}"
   adjointSource 1
        locationID=1
        responseID=1
   end adjointSource
   adjointSource 2
        locationID=1
        responseID=2
   end adjointSource
end importanceMap

read tallies
    pointDetector 1
        title="neutron detector"
        neutron
        locationID=1
        responseID=1
    end pointDetector
    pointDetector 2
        title="photon detector"
        photon
        locationID=1
        responseID=2
    end pointDetector
end tallies

end data
end
'''
        return mavric_output


class HandlingContactDoseEstimatorGenericTank(DoseEstimatorSquareTank):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a storage tank made of materials
    The storage tank is a cylinder, filled with the sample.
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
        adjoint_flux_file: str = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        self.sample_h2 = get_cyl_h(self.sample_volume, self.cyl_r) / 2.0
        tank_r: float = self.cyl_r
        tank_h2: float = self.sample_h2
        sample_z_max: float = self.sample_h2
        sample_z_min: float = - self.sample_h2

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
end

=mavric parm=(   )
{NOW} DoseEstimatorGenericTank, {self.sample_weight} g, layers {self.layers_thicknesses}
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=1
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
<{self.SAMPLE_ATOM_DENS_file_name_MAVRIC}
' helium 2 end
end comp

read geometry
global unit 1
    cylinder 1 {tank_r} 2p {tank_h2}
    media 1 1 1
'''
        xy_planes: list[float] = []  # list of XY cylinder boundaries for gridgeometry
        z_planes: list[float] = [sample_z_max, sample_z_min]  # list of Z cylinder boundaries for gridgeometry
        k: int = 0
        for k in range(len(self.layers_mats)):
            tank_r += self.layers_thicknesses[k]
            tank_h2 += self.layers_thicknesses[k]
            mavric_output += f'''
    cylinder {k + 2} {tank_r} 2p {tank_h2}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            xy_planes.append(tank_r)
            z_planes.append(tank_h2)
            z_planes.append(-tank_h2)
        self.det_x = tank_r + self.det_standoff_distance  # Detector is next to the tank, contact dose
        self.handling_det_x = tank_r + self.handling_det_standoff_distance  # Detector is next to the tank, handling dose
        xy_planes_str: str = " ".join([f' {x:.5f} -{x:.5f}' for x in xy_planes])
        self.det_z = (sample_z_max + sample_z_min) / 2.0
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {tank_r + self.box_a} 2p {tank_h2 + self.box_a}  
    media 0 1 99999 -{k + 2}
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

        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    distribution 1
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, neutrons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 1 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution'''

        mavric_output += f'''
    distribution 2
        title="Decayed sample after {self.DECAYED_SAMPLE_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 5 end
        filename="{self.DECAYED_SAMPLE_F71_file_name}"
    end distribution

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

read sources'''
        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cylinder {self.cyl_r} {sample_z_max} {sample_z_min}
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cylinder {self.cyl_r} {sample_z_max} {sample_z_min}
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1
'   adjointFluxes="{adjoint_flux_file}"
   adjointSource 1
        locationID=1
        responseID=1
   end adjointSource
   adjointSource 2
        locationID=1
        responseID=2
   end adjointSource
   adjointSource 5
        locationID=2
        responseID=1
   end adjointSource
   adjointSource 6
        locationID=2
        responseID=2
   end adjointSource
end importanceMap

read tallies
    pointDetector 1
        title="neutron detector"
        neutron
        locationID=1
        responseID=1
    end pointDetector
    pointDetector 2
        title="photon detector"
        photon
        locationID=1
        responseID=2
    end pointDetector
    pointDetector 5
        title="neutron detector"
        neutron
        locationID=2
        responseID=1
    end pointDetector
    pointDetector 6
        title="photon detector"
        photon
        locationID=2
        responseID=2
    end pointDetector

end tallies

end data
end
'''
        return mavric_output

    def get_responses(self):
        """ Reads over the MAVRIC output and returns responses for rem/h doses
            Note, this version uses different format of self.responses """
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name):
            raise FileNotFoundError("Expected decayed sample MAVRIC output file: \n" +
                                    self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(self.MAVRIC_out_file_name, 'r') as f:
            d: str = f.read()
        rs: list = re.findall(r'Point Detector (\d).  (\w+) detector\n.*\n.*\n.*\n.*\n.*\s+response (\d)\s+'
                              r'([+-]?\d+\.\d+[Ee]?[+-]?\d+)\s+([+-]?\d+\.\d+[Ee]?[+-]?\d+)?'
                              r'\s+([+-]?\d+\.\d+[Ee]?[+-]?\d+)?', d)
        self.responses = {t[0]: {'particle': t[1], 'pid': t[2], 'value': float(t[3]),
            'stdev': 0.0 if t[4] == '' else float(t[4])} for t in rs}

        os.chdir(self.cwd)
        if self.debug > 3:
            print(self.responses)

    def print_response(self):
        """ Prints dose responses """
        if self.responses is not {}:
            r1: dict = self.responses['1']  # Neutron dose, contact
            r2: dict = self.responses['2']  # Photon dose, contact
            r5: dict = self.responses['5']  # Neutron dose, handling (30cm)
            r6: dict = self.responses['6']  # Photon dose, handling (30cm)
            print(self.sample_weight, r1['value'], r1['stdev'], r2['value'], r2['stdev'],
                  r5['value'], r5['stdev'], r6['value'], r6['stdev'])


if __name__ == "__main__":
    pass
    # main()
