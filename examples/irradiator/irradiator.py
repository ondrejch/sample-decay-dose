#!/bin/env python3
"""
Wastewater irradiator using SampleDose - irradiation of SS316
Ondrej Chvala <ochvala@utexas.edu>
"""

from sample_decay_dose import SampleDose
import os
import shutil
import numpy as np
import json5

# SS-316H with 0.1 w% cobalt, 600C
SS316_Co1ppm_ADENS_HOT: dict = {'c': 0.000311167, 'si-28': 0.00153402, 'si-29': 7.79293e-05, 'si-30': 5.14317e-05,
    'p-31': 6.78722e-05, 's-32': 4.15187e-05, 's-33': 3.27814e-07, 's-34': 1.85761e-06, 's-36': 4.37085e-09,
    'cr-50': 0.000663653, 'cr-52': 0.0127979, 'cr-53': 0.00145118, 'cr-54': 0.000361229, 'mn-55': 0.00170071,
    'fe-54': 0.0031951, 'fe-56': 0.0501563, 'fe-57': 0.00115833, 'fe-58': 0.000154152, 'co-59': 7.93501e-05,
    'ni-58': 0.00650228, 'ni-60': 0.00250467, 'ni-61': 0.000108876, 'ni-62': 0.000347145, 'ni-64': 8.84075e-05,
    'mo-92': 0.000179806, 'mo-94': 0.000112364, 'mo-95': 0.000193563, 'mo-96': 0.000203058, 'mo-97': 0.000116381,
    'mo-98': 0.000294483, 'mo-100': 0.00011772}
# SS-316H with 0.1 w% cobalt, 20C
SS316_Co1ppm_ADENS_COLD: dict = {'c': 0.000320573, 'si-28': 0.00158039, 'si-29': 8.02849e-05, 'si-30': 5.29863e-05,
    'p-31': 6.99238e-05, 's-32': 4.27737e-05, 's-33': 3.37723e-07, 's-34': 1.91376e-06, 's-36': 4.50297e-09,
    'cr-50': 0.000683713, 'cr-52': 0.0131847, 'cr-53': 0.00149504, 'cr-54': 0.000372148, 'mn-55': 0.00175212,
    'fe-54': 0.00329168, 'fe-56': 0.0516724, 'fe-57': 0.00119334, 'fe-58': 0.000158812, 'co-59': 8.17487e-05,
    'ni-58': 0.00669882, 'ni-60': 0.00258038, 'ni-61': 0.000112167, 'ni-62': 0.000357638, 'ni-64': 9.10798e-05,
    'mo-92': 0.000185241, 'mo-94': 0.00011576, 'mo-95': 0.000199414, 'mo-96': 0.000209196, 'mo-97': 0.000119899,
    'mo-98': 0.000303385, 'mo-100': 0.000121279}

# How much cobalt is in steel
w_co = 0.1 * 1e-2  # convert to percent

#
# def cobalt_steel(wf_co: float = 0.1e-2) -> dict:
#     """ Calculate atom density of SS-316 with cobalt """
#     from MSRRpy.mat.material_types import Solid
#     ss_composition = [['c', 0.000800], ['mn', 0.020000], ['p', 0.000450], ['s', 0.000300], ['si', 0.010000],
#         ['cr', 0.170000], ['ni', 0.120000], ['mo', 0.025000], ['fe', 0.653450]]
#     my_steel = Solid(name='StainlessSteel', temperature=20.0 + 273.0, composition=ss_composition,
#                      composition_mode='weight', alpha=17.2e-6, ref_temperature=20.0 + 273.0, ref_density=8.0)
#     my_steel.add_impurity(composition=[['co', 1.0]], composition_mode='weight', fraction=wf_co, fraction_mode='weight')
#     return dict(my_steel.mixing_table)


def square_plate_volume(a: float = 100.0, h: float = 2.54) -> float:
    return a * a * h


class IrradiatorDose(SampleDose.Origen):
    """ Irradiator Dose using MAVRIC """

    def __init__(self, _o: SampleDose.Origen = None):
        self.debug: int = 3  # Debugging flag
        self.MAVRIC_input_file_name: str = 'my_dose.inp'
        self.IRRPLATE_ATOM_DENS_file_name_MAVRIC: str = 'my_sample_atom_dens_mavric.inp'
        self.histories_per_batch: int = 100000  # Monaco hist per batch
        self.batches: int = 10  # Monaco number of batches in total
        self.irr_plate_temperature_K: float = 293.0  # Sample temperature [K]
        self.irr_plate_atom_dens: dict = {}  # Atom density of the irradiating plate
        self.plate_a: float = 0.0  # side of irradiation plate [cm]
        self.plate_h: float = 0.0  # height of irradiation plate [cm]
        self.water_h: float = 0.0  # height of irradiated water above and below the irradiated plate [cm]
        self.tally_planes_xy: int = 1
        self.tally_planes_z: int = 1
        self.denovo_planes_xy: int = 1  # Planes inside the irradiation plate
        self.denovo_planes_z: int = 1
        if _o is not None:
            self.irr_plate_weight: float = _o.sample_weight  # Mass of the sample [g]
            self.irr_plate_density: float = _o.sample_density  # Mass density of the sample [g/cm3]
            self.irr_plate_volume: float = _o.sample_volume  # Sample volume [cm3]
            self.irr_plate_F71_file_name: str = _o.SAMPLE_F71_file_name
            self.irr_plate_F71_position: int = _o.SAMPLE_F71_position
            self.irr_plate_decay_days: float = _o.SAMPLE_DECAY_days  # Sample decay time [days]
            self.irr_plate_atom_dens: dict = _o.decayed_atom_dens  # Atom density of the decayed sample
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
        if not os.path.isfile(self.cwd + '/' + self.ORIGEN_dir + '/' + self.irr_plate_F71_file_name):
            raise FileNotFoundError(
                "Expected decayed sample F71 file: \n" + self.cwd + '/' + self.ORIGEN_dir + '/' + self.irr_plate_F71_file_name)
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)

        shutil.copy2(self.cwd + '/' + self.ORIGEN_dir + '/' + self.irr_plate_F71_file_name,
                     self.cwd + '/' + self.case_dir)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(self.IRRPLATE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
            f.write(SampleDose.atom_dens_for_mavric(self.irr_plate_atom_dens, 1, self.irr_plate_temperature_K))

        with open(self.MAVRIC_input_file_name, 'w') as f:  # write MAVRIC input deck
            f.write(self.mavric_deck())

        if self.debug > 0:
            print(f"MAVRIC: running case {self.case_dir}/{self.MAVRIC_input_file_name}")
        SampleDose.run_scale(self.MAVRIC_input_file_name)
        os.chdir(self.cwd)

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        adjoint_flux_file = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        plate_a_half: float = self.plate_a / 2.0
        plate_h_half: float = self.plate_h / 2.0
        box_z_half: float = plate_h_half + self.water_h  # half-height of the box in Z
        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.irr_plate_F71_file_name} .
cp -r ${{INPDIR}}/{self.IRRPLATE_ATOM_DENS_file_name_MAVRIC} .
end

=mavric parm=(   )
{SampleDose.NOW} Irradiator dose, {self.irr_plate_weight} g,
{SampleDose.MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    photons
    secondaryMult=1
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
<{self.IRRPLATE_ATOM_DENS_file_name_MAVRIC}
water 2  1.0  end
end comp

read geometry
global unit 1
' irradiated plate
cuboid 1 4p {plate_a_half} 2p {plate_h_half}
media 1 1  1
' water above
cuboid 2 4p {plate_a_half} {box_z_half} 0  
media 2 1  2 -1
' water below
cuboid 3 4p {plate_a_half} {0} {-box_z_half}
media 2 1  3 -1
' boundary  
cuboid 99  4p {plate_a_half} 2p {box_z_half} 
media 0 1 99 -1 -2 -3
boundary 99
end geometry

read definitions
    response 102
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response 
    response 112 
        title="ICRU-57 Table A.21 (air) Kerma rad/h, photons"
        photon
        doseData=9507
    end response 
    response 122 
        title="Henderson conversion factors rad/h, photons"
        photon
        doseData=9502
    end response 
'''

        mavric_output += f'''
    distribution 2
        title="Decayed irradiator after {self.irr_plate_decay_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.irr_plate_F71_position} 5 end
        filename="{self.irr_plate_F71_file_name}"
    end distribution

    gridGeometry 1
        title="Grid over water"
        xLinear {self.tally_planes_xy} {-plate_a_half} {plate_a_half}
        yLinear {self.tally_planes_xy} {-plate_a_half} {plate_a_half}
        zLinear {self.tally_planes_z} {plate_h_half} {box_z_half}
        zLinear {self.tally_planes_z} {-box_z_half} {-plate_h_half} 
    end gridGeometry

    gridGeometry 9
        title="Denovo grid"
        xLinear {self.tally_planes_xy} {-plate_a_half} {plate_a_half}
        yLinear {self.tally_planes_xy} {-plate_a_half} {plate_a_half}
        zLinear {self.tally_planes_z} {plate_h_half} {box_z_half}
        zLinear {self.tally_planes_z} {-box_z_half} {-plate_h_half} 
        zLinear {self.denovo_planes_z} {-plate_h_half} {plate_h_half} 
    end gridGeometry
end definitions

read sources
    src 2
        title="Sample photons"
        photon
        useNormConst
        cuboid {-plate_a_half} {plate_a_half} {-plate_a_half} {plate_a_half} {-plate_a_half} {plate_a_half} 
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1
'   adjointFluxes="{adjoint_flux_file}"
   adjointSource 2
        boundingBox {-plate_a_half} {plate_a_half} {-plate_a_half} {plate_a_half} {-box_z_half} {box_z_half} 
        responseID=102
        mixture=2
   end adjointSource
end importanceMap

read tallies
    meshTally 1 
        title="Photon ANSI dose [rem/h]]" 
        photon 
        gridGeometryID=1 
        responseID=102 
        noGroupFluxes
    end meshTally 
    meshTally 2 
        title="Photon ICRU-57 (air) Kerma [rad/h]" 
        photon 
        gridGeometryID=1 
        responseID=112 
        noGroupFluxes
    end meshTally 
    meshTally 3 
        title="Henderson conversion factors [rad/h]]" 
        photon 
        gridGeometryID=1 
        responseID=122 
        noGroupFluxes
    end meshTally  
end tallies
end data
end
'''
        return mavric_output


def run_analysis():
    # Irradiator plate
    plate_density: float = SampleDose.get_rho_from_atom_density(SS316_Co1ppm_ADENS_COLD)  # [g/cm3]
    plate_size: float = 100.0  # [cm]
    plate_height: float = 2.54  # [cm]
    plate_mass: float = square_plate_volume(plate_size, plate_height) * plate_density  # [g]

    irr = SampleDose.OrigenIrradiation('EIRENE.mix0002.f33', plate_mass)
    irr.irradiate_days = 7.0 * 365.24  # 7 years
    irr.irradiate_flux = 7.37e12  # n/s/cm2
    irr.set_decay_days(0.5 * 365.24)  # 6 months
    irr.set_decay_days(1. * 365.24)
    irr.SAMPLE_F71_position = 50  # How many decay steps
    irr.case_dir = 'run_irradiator_1m'
    irr.write_atom_dens(SS316_Co1ppm_ADENS_COLD)
    irr.run_irradiate_decay_sample()

    dose_irr = IrradiatorDose(irr)
    dose_irr.MAVRIC_input_file_name = 'irradiator.inp'
    dose_irr.batches = 40
    dose_irr.histories_per_batch = 500000
    # dose_irr.batches = 10
    # dose_irr.histories_per_batch = 50000
    dose_irr.plate_h = plate_height  # Geometry
    dose_irr.plate_a = plate_size
    dose_irr.water_h = 200.0
    dose_irr.tally_planes_z = 40  # Meshing
    dose_irr.tally_planes_xy = 20
    dose_irr.denovo_planes_xy = 20
    dose_irr.denovo_planes_z = 4
    dose_irr.run_mavric()


if __name__ == "__main__":
    run_analysis()
