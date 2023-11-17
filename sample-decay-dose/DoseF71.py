#!/bin/env python3
"""
Handling dose in SCALE
Ondrej Chvala <ochvala@utexas.edu>
"""

# import sys
import os
import re
# import argparse
import subprocess
import math
from bisect import bisect_left

SCALE_bin_path = os.getenv('SCALE_BIN', '/opt/scale6.3.1/bin/')
ORIGEN_input_file_name: str = 'fuel_salt_decay.inp'  # SCALE deck name for sample decay
MAVRIC_input_file_name: str = 'my_sample.inp'  # SCALE deck name for MAVRIC/Monaco transport
SAMPLE_ATOM_DENS_file_name_Origen: str = 'my_sample_atom_dens_origen.inp'
SAMPLE_ATOM_DENS_file_name_MAVRIC: str = 'my_sample_atom_dens_mavric.inp'
ATOM_DENS_MINIMUM: float = 1e-60
MAVRIC_NG_XSLIB: str = 'v7.1-28n19g'


def get_cyl_r(cyl_volume: float) -> float:
    """ Radius of a square cylinder from its volume
        V = pi r^2 h = pi r^2 2 r = 2 pi r^3
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


def get_positions_index(f71file: str) -> dict:
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


def get_burned_salt_atom_dens(f71file: str, position: int) -> dict:
    """ Read atom density of nuclides from SCALE's F71 file
    """
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


def atom_dens_for_origen(dens: dict) -> str:
    """ Print atom densities in ORIGEN-input format
    """
    output = ''
    for k, v in dens.items():
        if v > 0.0:
            output += f'{k} = {v} \n'
    return output


def atom_dens_for_mavric(dens: dict) -> str:
    """ Print atom densities in SCALE-input format
    """
    output = ''
    for k, v in dens.items():
        if v > 0.0:
            k = re.sub('m$', '', k)
            output += f'{k} 1 0 {v} 873 end\n'
    return output


class DoseEstimator:
    """ """
    def __init__(self, _f71: str = './SCALE_FILE.f71', _mass: float = 0.1, _density: float = 2.406):
        self.debug: int = 3  # Debugging flag
        self.BURNED_SALT_F71_file_name: str = _f71  # Burned core F71 file from TRITON
        self.SAMPLE_MASS: float = _mass  # Mass of the sample [g]
        self.SAMPLE_DENSITY: float = _density  # Mass density of the sample [g/cm3]. TODO: Calculate from obiwan data
        self.BURNED_SALT_F71_index: dict = get_positions_index(self.BURNED_SALT_F71_file_name)
        self.BURNED_SALT_F71_position: int = 16
        self.decks = None  # ScaleInput instance
        self.case_dir: str = f'run_{_mass:.5}_g'  # Directory to run the case
        self.cwd: str = os.getcwd()  # Current running fir
        self.burned_atom_dens: dict = {}  # Atom density of the burned material from F71 file
        self.DECAYED_SALT_F71_file_name: str = MAVRIC_input_file_name.replace('inp', 'f71')
        self.DECAYED_SALT_out_file_name: str = MAVRIC_input_file_name.replace('inp', 'out')
        self.DECAYED_SALT_F71_position: int = 12
        self.decayed_atom_dens: dict = {}  # Atom density of the decayed sample
        self.responses: dict = {}  # Dose responses

    def set_f71_pos(self, t: float = 5184000.0, case: str = '1'):
        """ Returns closest position in the F71 file for a case """
        pos_times = [(k, float(v['time'])) for k, v in self.BURNED_SALT_F71_index.items() if v['case'] == case]
        times = [x[1] for x in pos_times]
        t_min = min(times)
        t_max = max(times)
        if t < t_min:
            print(f"Error: Time {t} seconds is less than {t_min} s, the minimum time in records.")
        #    return None
        if t > t_max:
            print(f"Error: Time {t} seconds is longer than {t_max} s, the maximum time in records.")
        #    return None
        self.BURNED_SALT_F71_position = bisect_left(times, t)

    def define_decks(self):
        """ Initializes Origen and Mavric deck writers """
        self.decks = ScaleInput(self.SAMPLE_MASS, self.SAMPLE_DENSITY)
        self.decks.DECAYED_SALT_F71_file_name = self.DECAYED_SALT_F71_file_name
        self.decks.DECAYED_SALT_F71_position = self.DECAYED_SALT_F71_position

    def read_burned_material(self):
        """ Reads atom density from F71 file """
        self.burned_atom_dens = get_burned_salt_atom_dens(self.BURNED_SALT_F71_file_name, self.BURNED_SALT_F71_position)
        if self.debug > 2:
            print(list(self.burned_atom_dens.items())[:25])

    def run_decay_sample(self):
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(SAMPLE_ATOM_DENS_file_name_Origen, 'w') as f:  # write Origen at-dens sample input
            f.write(atom_dens_for_origen(self.burned_atom_dens))

        with open(ORIGEN_input_file_name, 'w') as f:  # write ORIGEN input deck
            f.write(self.decks.origen_deck())

        print(f"\nRUNNING {ORIGEN_input_file_name}")
        run_scale(ORIGEN_input_file_name)

        self.decayed_atom_dens = get_burned_salt_atom_dens(self.DECAYED_SALT_F71_file_name,
                                                           self.DECAYED_SALT_F71_position)
        if self.debug > 2:
            print(list(self.decayed_atom_dens.items())[:25])

    def run_mavric(self):
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.DECAYED_SALT_F71_file_name):
            raise FileNotFoundError("Expected decayed sample F71 file: \n" +
                                    self.cwd + '/' + self.case_dir + '/' + self.DECAYED_SALT_F71_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
            f.write(atom_dens_for_mavric(self.decayed_atom_dens))

        with open(MAVRIC_input_file_name, 'w') as f:  # write MAVRICinput deck
            f.write(self.decks.mavric_deck())

        print(f"\nRUNNING {MAVRIC_input_file_name}")
        run_scale(MAVRIC_input_file_name)

    def get_responses(self):
        """ Reads over the MAVRIC output and returns responses for rem/h doses """
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.DECAYED_SALT_out_file_name):
            raise FileNotFoundError("Expected decayed sample MAVRIC output file: \n" +
                                    self.cwd + '/' + self.case_dir + '/' + self.DECAYED_SALT_out_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)

        tally_sep: str = 'Final Tally Results Summary'
        is_in_tally: bool = False

        with open(self.DECAYED_SALT_out_file_name, 'r') as f:
            for line in f.read().splitlines():
                if line.count(tally_sep) > 0:
                    is_in_tally = True
                if is_in_tally:
                    if line.count('response') == 1:
                        s = line.split()
                        self.responses[s[1]] = {'value': s[2], 'stdev': s[3]}

        if self.debug > 3:
            print(self.responses)

    def print_response(self):
        if self.responses is not {}:
            r1 = self.responses['1']
            r2 = self.responses['2']
            print(self.SAMPLE_MASS, r1['value'], r1['stdev'], r2['value'], r2['stdev'])


# -------------------------------------------------------------------------------------------------------------
class ScaleInput:
    """ This class writes SCALE decks for
    1. Origen - to decay sample
    2. MAVRIC - to calculate dose at det_x distance
    """

    def __init__(self, _sample_weight: float = 0.1, _sample_density: float = 2.406):
        self.SAMPLE_ATOM_DENS_file_name_Origen = None
        self.DECAYED_SALT_F71_file_name: str = 'my_sample.f71'
        self.DECAYED_SALT_F71_position: int = 12
        self.sample_weight: float = _sample_weight  # [g]
        self.sample_density: float = _sample_density  # [g/cm^3]
        self.sample_volume = self.sample_weight / self.sample_density
        self.cyl_r = get_cyl_r(self.sample_volume)
        self.det_x: float = 30.0  # Detector distance [cm]
        self.box_a: float = self.det_x + 10.0  # Problem box distance [cm]
        self.N_planes_box: int = 10  # Planes per box
        self.N_planes_cyl: int = 8  # Planes per cylinder

    def origen_deck(self) -> str:
        """ Salt decay Origen deck
        """
        origen_output = f'''
=shell
cp -r ${{INPDIR}}/{SAMPLE_ATOM_DENS_file_name_Origen} .
end

=origen
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
}}
case {{
    gamma=yes
    neutron=yes
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
        t=[0.0001 0.001 0.01 0.1 0.5 1 2 5 10 15 30]
        start=0
    }}
    save {{
        file="{self.DECAYED_SALT_F71_file_name}"
    }}
}}
end
'''
        return origen_output

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation
        """
        adjoint_flux_file = MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SALT_F71_file_name} .
cp -r ${{INPDIR}}/{SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
end

=mavric parm=(   )
Sample dose, {self.sample_weight} g, at x={self.det_x} cm
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=00003ecd7b4e3e8b
    ceLibrary="ce_v7.1_endf.xml"
    neutrons  photons
    fissionMult=1  secondaryMult=0
    perBatch=1000000 batches=10
end parameters

read comp
<{SAMPLE_ATOM_DENS_file_name_MAVRIC}
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
        title="ANSI standard (1991) flux-to-dose-rate factors, neutrons"
        doseData=9031
    end response
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors, photons"
        doseData=9505
    end response

    distribution 1
        title="Fuel salt after 1 month, neutrons"
        special="origensBinaryConcentrationFile"
        parameters 1 1 end
        filename="{self.DECAYED_SALT_F71_file_name}"
    end distribution
    distribution 2
        title="Fuel salt after 1 month, photons"
        special="origensBinaryConcentrationFile"
        parameters 1 5 end
        filename="{self.DECAYED_SALT_F71_file_name}"
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

read sources
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cylinder {self.cyl_r} {self.cyl_r} -{self.cyl_r}
        eDistributionID=1
    end src
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


# class OutputReaderMonaco:
#     def __init__(self):
#         self.file_in_name = 'my_sample.out'
#         self.responses: dict = None
#
#     def get_responses(self) -> dict:
#         tally_sep: str = 'Final Tally Results Summary'
#         is_in_tally: bool = False
#
#         with open(self.file_in_name, 'r') as f:
#             for l in f.read().splitlines():
#                 if l.count(tally_sep) > 0:
#                     is_in_tally = True
#                 if is_in_tally:
#                     if l.count('response') == 1:
#                         s = l.split()
#                         self.responses[s[1]] = {'value': s[2], 'stdev': s[3]}
#
#         # print(responses)
#         # return responses
#
#     def read_response(self):
#         cwd = os.getcwd()
#         with os.scandir(cwd) as it:
#             for entry in it:
#                 if entry.name.startswith('g_') and entry.is_dir():
#                     # print(entry.name)
#                     mass_g = entry.name.replace('g_', '')
#                     resps = self.get_responses(entry.name + '/' + self.file_in_name)
#                     r1 = resps['1']
#                     r2 = resps['2']
#                     print(mass_g, r1['value'], r1['stdev'], r2['value'], r2['stdev'])


if __name__ == "__main__":
    pass
    # main()
