import io
import re
import subprocess
import numpy as np
import pandas as pd
from sample_decay_dose.SampleDose import SCALE_bin_path, ATOM_DENS_MINIMUM


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


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
            iso_name = re.sub(r'm$', '', iso_name)  # strip "m"
            iso_name = re.sub(r'([a-zA-Z]+)(\d+)', '\g<1>-\g<2>', iso_name)  # add dash
            my_rho += rel_iso_mass[iso_name] * v * m_Da * 1e24
    return my_rho


def scale_adens(adens: dict, scalef: float = 1.0) -> dict:
    """" Scale atom densities by a factor scalef """
    assert not np.isnan(scalef)
    new_adens: dict = {}
    for k, v in adens.items():
        new_adens[k] = v * scalef
    return new_adens


def get_f71_positions_index(f71file: str) -> dict:
    """ Read info of SCALE's F71 file """
    output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=info", f71file], capture_output=True)
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
        try:
            f71_idx[int(data[0])] = {'time': data[1], 'power': data[2], 'flux': data[3], 'fluence': data[4],
                'energy': data[5], 'initialhm': data[6], 'libpos': data[7], 'case': data[8], 'step': data[9],
                'DCGNAB': data[10]}
        except (IndexError, ValueError):
            raise IndexError
    return f71_idx


def get_last_position_for_case(f71file: str, case: int = 1) -> int:
    iii = get_f71_positions_index(f71file)
    if not iii:
        raise RuntimeError("Got no data from get_f71_positions_index()")
    maxi = max([i for i, rec in iii.items() if rec['case'] == str(case)])
    return maxi


def get_f71_nuclide_case(f71file: str, f71units: str = 'atom', my_cases=None) -> pd.DataFrame:
    """ Read atom nuclide data from SCALE's F71 file
    f71units:   abso|fiss|capt|airm|apel|atom|becq|curi|gamw|gamm|gato|gper|gram|h2om|
                kilo|wpel|watt|mevs|part|inte|ener """
    if my_cases is None:
        my_cases = [1]
    cases_str: str = f'-cases={my_cases}'
    output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=20", f"-units={f71units}",
        "-idform='{:Ee}{:AAA}{:m}'", cases_str, f71file], capture_output=True)
    return pd.read_csv(io.StringIO(output.stdout.decode()), skipinitialspace=True, index_col=0)


def get_f71_elements_case(f71file: str, f71units: str = 'atom', my_cases=None) -> pd.DataFrame:
    """ Read atom nuclide data from SCALE's F71 file
    f71units:   abso|fiss|capt|airm|apel|atom|becq|curi|gamw|gamm|gato|gper|gram|h2om|
                kilo|wpel|watt|mevs|part|inte|ener """
    if my_cases is None:
        my_cases = [1]
    cases_str: str = f'-cases={my_cases}'
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=20", f"-units={f71units}", "-idform='{:Ee}'",
            cases_str, f71file], capture_output=True)
    return pd.read_csv(io.StringIO(output.stdout.decode()), skipinitialspace=True, index_col=0)


def get_burned_nuclide_atom_dens(f71file: str, position: int, my_cases: list[int] = None) -> dict:
    """ Read atom density of nuclides from SCALE's F71 file """
    runlist: list = [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom",
        "-idform='{:Ee}{:AAA}{:m}'", f71file]
    if my_cases:
        cases_str: str = f'-cases={my_cases}'
        runlist = [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom",
            "-idform='{:Ee}{:AAA}{:m}'", cases_str, f71file]

    output = subprocess.run(runlist, capture_output=True)
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


def get_burned_nuclide_data(f71file: str, position: int, f71units: str = 'atom', my_cases: list[int] = None) -> dict:
    """ Read atom nuclide data from SCALE's F71 file
    f71units:   abso|fiss|capt|airm|apel|atom|becq|curi|gamw|gamm|gato|gper|gram|h2om|
                kilo|wpel|watt|mevs|part|inte|ener """
    runlist: list = [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", f"-units={f71units}",
        "-idform='{:Ee}{:AAA}{:m}'", f71file]
    if my_cases:
        cases_str: str = f'-cases={my_cases}'
        runlist = [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", f"-units={f71units}",
            "-idform='{:Ee}{:AAA}{:m}'", cases_str, f71file]
    output = subprocess.run(runlist, capture_output=True)
    output = output.stdout.decode().split("\n")
    f71unit_data = {}
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
                f71unit_data[nuclide] = float(data[position])
    sorted_f71unit_data = {k: v for k, v in sorted(f71unit_data.items(), key=lambda item: -item[1])}
    return sorted_f71unit_data


def get_single_nuclide_case(f71file: str, position: int, my_nuclide: str) -> float:
    my_nuclide = my_nuclide.lower()
    output = subprocess.run(
        [f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=becq", "-idform='{:Ee}{:AAA}{:m}'",
            f71file], capture_output=True)
    output = output.stdout.decode().split("\n")
    # bqs = {}  # densities[nuclide] = (density at position of f71 file)
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
                # bqs[nuclide] = float(data[position])
                if nuclide == my_nuclide:
                    # print(f"{my_nuclide}: {data[position]}")
                    return float(data[position])
    # sorted_densities = {k: v for k, v in sorted(bqs.items(), key=lambda item: -item[1])}
    return -1.0


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
    output = subprocess.run([f"{SCALE_bin_path}/obiwan", "info", f33file], capture_output=True)
    output = output.stdout.decode().split("\n")
    for line in output:
        data = line.split()
        if data[0] in f33file:
            return int(data[2])
    return -1


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


def get_cyl_r(cyl_volume: float) -> float:
    """ Radius of a square cylinder from its volume
        V = pi r^2 h = pi r^2 2 r = 2 pi r^3
        r = (V / 2 pi) ^ 1/3
    """
    return (cyl_volume / (2.0 * np.pi)) ** (1.0 / 3.0)


def get_cyl_r_4_1(cyl_volume: float) -> float:
    """ Radius of a 4:1 H:R cylinder from its volume
        V = pi r^2 h = pi r^2 4 r = 4 pi r^3
        r = (V / 4 pi) ^ 1/3
        h_cyl = 4 * r
    """
    return (cyl_volume / (4.0 * np.pi)) ** (1.0 / 3.0)


def get_fill_height_4_1(fill_volume: float, cyl_volume: float) -> float:
    """ Fill height of 4:1 H:R cylinder
        V_cyl = pi r^2 h_cyl = pi r^2 4 r = 4 pi r^3
        r = (V_cyl / 4 pi) ^ 1/3;  h_cyl = 4 * r
        V_fill = pi r^2 h_fill
        h_fill = V_fill / (pi r^2)
    """
    if fill_volume > cyl_volume:
        raise ValueError("Fill volume is larger than the cylinder")
    r: float = (cyl_volume / (4.0 * np.pi)) ** (1.0 / 3.0)
    return fill_volume / (np.pi * r ** 2)


def get_cyl_h(cyl_volume: float, cyl_r: float) -> float:
    """ Radius of a square cylinder from its volume
        V = pi r^2 h
        h = V / (pi r^2)
    """
    if not cyl_r > 0:
        raise ValueError(f"Cylinder radius has to be positive, {cyl_r}")
    return cyl_volume / np.pi / cyl_r ** 2


def run_scale(deck_file: str, nmpi: int = 1):
    """ Run a SCALE deck """
    scale_out = subprocess.run([f"{SCALE_bin_path}/scalerte", "-N", str(nmpi), "-m", deck_file], capture_output=True)
    scale_out = scale_out.stdout.decode().split("\n")
    if any("error" in line.lower() for line in scale_out):
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


def read_scale_keff(scale_out_file: str) -> tuple[float, float]:
    """ Returns a k_eff from SCALE output file """
    with open(scale_out_file, 'r') as file:
        content = file.read()
    subtext1: str = re.findall(r'best estimate system k-eff \s+(.*)', content)[0]
    k_list: list = re.findall(r'(\d+\.\d+)', subtext1)
    assert len(k_list) == 2
    return float(k_list[0]), float(k_list[1])


def rho(k: float) -> float:
    """ Returns pcm reactivity """
    return 1e5 * (k - 1.0) / k


def rho_err(k: tuple[float, float]) -> tuple[float, float]:
    """ Returns pcm with error: d rho / d k = 1/k^2, sig_rho = sqrt [ (d rho/d k)^2 * sig_k^2 ]  """
    my_rho: float = 1e5 * (k[0] - 1.0) / k[0]
    my_rho_err: float = 1e5 * k[1] / (k[0] ** 2)
    return my_rho, my_rho_err
