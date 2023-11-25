#!/bin/env python3
"""
Handling dose [rem/h] in SCALE
Ondrej Chvala <ochvala@utexas.edu>
"""
import os
import re
import subprocess
import math
import numpy as np
from bisect import bisect_left
from sample_decay_dose.read_opus import integrate_opus

SCALE_bin_path: str = os.getenv('SCALE_BIN', '/opt/scale6.3.1/bin/')
ATOM_DENS_MINIMUM: float = 1e-60
MAVRIC_NG_XSLIB: str = 'v7.1-28n19g'

RHO: float = 8.0  # TODO: replace with calculation from atomic density
# https://www.sandmeyersteel.com/316H.html
ATOM_DENS_SS316H_HOT: str = """c = 0.0002384
n-14 = 0.000254609
n-15 = 9.30164e-07
al-27 = 5.30544e-05
si-28 = 1.56475e-05
si-29 = 7.94907e-07
si-30 = 5.24622e-07
p-31 = 6.93324e-05
s-32 = 4.23788e-05
s-33 = 3.34604e-07
s-34 = 1.89609e-06
s-36 = 4.46139e-09
ti-46 = 3.28985e-06
ti-47 = 2.96684e-06
ti-48 = 2.93973e-05
ti-49 = 2.15734e-06
ti-50 = 2.06562e-06
cr-50 = 0.000677912
cr-52 = 0.0130729
cr-53 = 0.00148236
cr-54 = 0.00036899
mn-55 = 1.73977e-05
fe-54 = 0.00390032
fe-56 = 0.0612267
fe-57 = 0.00141399
fe-58 = 0.000188176
ni-58 = 6.64309e-05
ni-60 = 2.55891e-05
ni-61 = 1.11234e-06
ni-62 = 3.54662e-06
ni-64 = 9.03221e-07
mo-92 = 0.00018364
mo-94 = 0.00011476
mo-95 = 0.00019769
mo-96 = 0.000207388
mo-97 = 0.000118863
mo-98 = 0.000300762
mo-100 = 0.00012023"""


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


def atom_dens_for_origen(dens: dict) -> str:
    """ Print atom densities in ORIGEN input format """
    output = ''
    for k, v in dens.items():
        if v > 0.0:
            output += f'{k} = {v} \n'
    return output


def atom_dens_for_mavric(dens: dict) -> str:
    """ Print atom densities in SCALE CSAS/MAVRIC input format """
    output = ''
    for k, v in dens.items():
        if v > 0.0:
            k = re.sub('m$', '', k)
            output += f'{k} 1 0 {v} 873 end\n'
    return output


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
        self.SAMPLE_F71_position: int = 12
        self.SAMPLE_DECAY_days: float = 30.0  # Sample decay time [days]
        self.sample_weight: float = np.NaN  # Mass of the sample [g]
        self.sample_density: float = RHO  # np.NaN  # Mass density of the sample [g/cm3]
        self.sample_volume: float = np.NaN  # Sample volume [cm3]

    def set_decay_days(self, decay_days: float = 30.0):
        """ Use this to change decay time, as it also updates the case directory """
        self.SAMPLE_DECAY_days = decay_days
        self.case_dir: str = f'run_{self.sample_weight:.5}_g-{decay_days:.5}_days'  # Directory to run the case

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
        #    return None
        if t > t_max:
            print(f"Error: Time {t} seconds is longer than {t_max} s, the maximum time in records.")
        #    return None
        self.BURNED_MATERIAL_F71_position = bisect_left(times, t)

    def read_burned_material(self):
        """ Reads atom density and rho from F71 file """
        self.sample_density = get_burned_material_total_mass_dens(self.BURNED_MATERIAL_F71_file_name,
                                                                  self.BURNED_MATERIAL_F71_position)
        self.sample_volume = self.sample_weight / self.sample_density
        self.burned_atom_dens = get_burned_material_atom_dens(self.BURNED_MATERIAL_F71_file_name,
                                                              self.BURNED_MATERIAL_F71_position)
        if self.debug > 2:
            print(list(self.burned_atom_dens.items())[:25])

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

        print(f"\nRUNNING {self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

        self.decayed_atom_dens = get_burned_material_atom_dens(self.SAMPLE_F71_file_name,
                                                               self.SAMPLE_F71_position)
        os.chdir(self.cwd)
        if self.debug > 2:
            print(list(self.decayed_atom_dens.items())[:25])

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

    def write_atom_dens(self, str_atom_dens: str = ATOM_DENS_SS316H_HOT):
        """ Writes atom density of the fresh material to be irradiated """
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        with open(self.SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write at-dens input for Origen irradiation
            f.write(str_atom_dens)
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

        print(f"\nRUNNING {self.ORIGEN_input_file_name}")
        run_scale(self.ORIGEN_input_file_name)

        self.decayed_atom_dens = get_burned_material_atom_dens(self.SAMPLE_F71_file_name,
                                                               self.SAMPLE_F71_position)
        os.chdir(self.cwd)
        if self.debug > 2:
            print(list(self.decayed_atom_dens.items())[:25])

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
cp -r ${{INPDIR}}/{self.F33_file_name} .
end

=origen
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


class DoseEstimator:
    """ MAVRIC calculation of rem/h doses from the decayed sample """
    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        self.debug: int = 3  # Debugging flag
        self.sample_weight: float = _o.sample_weight  # Mass of the sample [g]
        self.sample_density: float = _o.sample_density  # Mass density of the sample [g/cm3]
        self.sample_volume: float = _o.sample_volume  # Sample volume [cm3]
        self.DECAYED_SAMPLE_F71_file_name: str = _o.SAMPLE_F71_file_name
        self.DECAYED_SAMPLE_F71_position: int = _o.SAMPLE_F71_position
        self.DECAYED_SAMPLE_days: float = _o.SAMPLE_DECAY_days  # Sample decay time [days]
        self.decayed_atom_dens: dict = _o.decayed_atom_dens  # Atom density of the decayed sample
        self.beta_over_gamma: float = _o.get_beta_to_gamma()  # Beta over gamma spectral ratio
        self.neutron_intensity: float = _o.get_neutron_integral()  # Integral of neutron spectra
        self.case_dir: str = _o.case_dir  # Directory to run the case
        self.cwd: str = _o.cwd  # Current running fir
        self.MAVRIC_input_file_name: str = 'my_dose.inp'
        self.MAVRIC_out_file_name: str = self.MAVRIC_input_file_name.replace('inp', 'out')
        self.SAMPLE_ATOM_DENS_file_name_MAVRIC: str = 'my_sample_atom_dens_mavric.inp'
        self.responses: dict = {}  # Dose responses 1: neutron, 2: gamma, 3: beta
        self.det_x: float = 30.0  # Detector distance [cm]
        self.N_planes_box: int = 5  # Planes per box
        self.N_planes_cyl: int = 8  # Planes per cylinder
        self.box_a: float = np.NaN
        self.cyl_r: float = np.NaN

    def run_mavric(self):
        """ Writes Mavric inputs and runs the case """
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.DECAYED_SAMPLE_F71_file_name):
            raise FileNotFoundError("Expected decayed sample F71 file: \n" +
                                    self.cwd + '/' + self.case_dir + '/' + self.DECAYED_SAMPLE_F71_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(self.SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
            f.write(atom_dens_for_mavric(self.decayed_atom_dens))

        with open(self.MAVRIC_input_file_name, 'w') as f:  # write MAVRICinput deck
            f.write(self.mavric_deck())

        print(f"\nRUNNING {self.MAVRIC_input_file_name}")
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
Sample dose, {self.sample_weight} g, at x={self.det_x} cm
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=00003ecd7b4e3e8b
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=1
    perBatch=100000 batches=10
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
'   adjointFLuxes="{adjoint_flux_file}"
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


if __name__ == "__main__":
    pass
    # main()
