import os
import numpy as np

from sample_decay_dose.SampleDose import NOW, MAVRIC_NG_XSLIB, Origen, run_scale, atom_dens_for_mavric, \
    get_f71_positions_index, OrigenFromTriton
from sample_decay_dose.HotCell import HotCellDoses
from typing import TypedDict


class RadiatorGeometry(TypedDict):
    length_I: float     # Length of inner fluid
    n_X: int            # Number of pins along X
    n_Z: int            # Number of pins along Z
    pin_ID: float       # Inside diameter of a tube pin
    pin_OD: float       # Outside diameter of a tube pin
    pin_pitch: float    # Tube pin pitch


class RadiatorBox(HotCellDoses):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a cubical radiator
        1st material in "shielding" layer, mix=10, is the tube material
        2nd material in "shielding" layer, mix=11, is in-between tubes (air)
    """
    def __init__(self, _o: OrigenFromTriton = None):
        """ This reads decayed sample information from the Origen object """
        super().__init__(_o)  # Init DoseEstimator
        self.origen_from_triton: OrigenFromTriton = _o
        self.ORIGEN_dir = 'run'
        self.decay_days = 0.0               # no decay, just raw F71, copy the data from F71 file
        self.DECAYED_SAMPLE_F71_file_name = _o.BURNED_MATERIAL_F71_file_name
        self.DECAYED_SAMPLE_F71_index = get_f71_positions_index(self.DECAYED_SAMPLE_F71_file_name)
        self.handling_det_x: (None, float) = None
        self.radiator_geometry: RadiatorGeometry = {  # 10ft MSRE-like model
            "length_I": 10 * 30.48,
            "n_X": 10,
            "n_Z": 13,
            "pin_ID":  0.75 * 2.54,
            "pin_OD":  (0.75 + 2.0 * 0.072) * 2.54,
            "pin_pitch": 1.5 * 2.54
            # "length_I": 50,
            # "n_X": 10,
            # "n_Z": 10,
            # "pin_ID": 2,
            # "pin_OD":  3,
            # "pin_pitch": 5
        }
        # Add steel layer around the problem to simulate back-scatter, TBD
        self.radiator_surrounded_by_steel_box_width: (None, float) = None

    @property
    def pin_IR(self) -> float:
        return self.radiator_geometry['pin_ID'] / 2.0

    @property
    def pin_OR(self) -> float:
        return self.radiator_geometry['pin_OD'] / 2.0

    @property
    def half_pitch(self) -> float:
        return self.radiator_geometry['pin_pitch'] / 2.0

    @property
    def half_pin_length(self) -> float:
        return self.radiator_geometry['length_I'] / 2.0
    @property
    def wall_thickness(self) -> float:
        return self.pin_OR - self.pin_IR

    @property
    def pin_volume(self) -> float:
        return self.radiator_geometry["length_I"] * np.pi * self.pin_IR ** 2

    @property
    def all_pins_volume(self) -> float:
        return self.radiator_geometry['n_X'] * self.radiator_geometry['n_Z'] * self.pin_volume

    def run_mavric(self, nmpi: int = 1):
        """ Writes Mavric inputs and runs the case """
        self.case_dir: str = f'run_MAVRIC_{NOW}_{self.all_pins_volume:.5}_cm-{self.decay_days:.5}_days'  # Directory to run the case
        self.DECAYED_SAMPLE_F71_position = self.origen_from_triton.BURNED_MATERIAL_F71_position

        if not os.path.isfile(os.path.join(self.cwd, self.DECAYED_SAMPLE_F71_file_name)):
            raise FileNotFoundError(
                "Expected decayed sample F71 file: \n" + os.path.join(self.cwd, self.DECAYED_SAMPLE_F71_file_name))
        if not os.path.exists(self.case_dir):
            os.mkdir(self.case_dir)
        os.chdir(self.case_dir)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(self.SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
            f.write(atom_dens_for_mavric(self.origen_from_triton.burned_atom_dens, 1, self.sample_temperature_K))
            for k in range(len(self.layers_mats)):
                f.write(atom_dens_for_mavric(self.layers_mats[k], k + 10, self.layers_temperature_K[k]))

        with open(self.MAVRIC_input_file_name, 'w') as f:  # write MAVRICinput deck
            f.write(self.mavric_deck())

        if self.debug > 0:
            print(f"MAVRIC: running case {self.case_dir}/{self.MAVRIC_input_file_name}")

        run_scale(self.MAVRIC_input_file_name, nmpi)
        os.chdir(self.cwd)

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        if not self.det_z:
            self.det_z = 0.0
        adjoint_flux_file: str = os.path.join(self.cwd, self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff'))
        # self.cyl_r = get_cyl_r(self.sample_volume)
        # self.sample_h2 = self.cyl_r             # sample is a square cylinder
        # sample_r: float = self.cyl_r            # sample outer layer [cm]
        # sample_h2: float = self.sample_h2
        box_y2: float = self.half_pin_length + self.wall_thickness          # X box around
        box_x2: float = self.radiator_geometry['n_X'] * self.half_pitch     # Y box around
        box_z2: float = self.radiator_geometry['n_Z'] * self.half_pitch     # Z box around

        place_x: float = -box_x2 + self.half_pitch  # array placement vector for array unit 1 1 1
        place_y: float = 0.0
        place_z: float = -box_z2 + self.half_pitch

        nX = self.radiator_geometry['n_X']
        nZ = self.radiator_geometry['n_Z']

        self.det_x = box_x2 + self.det_standoff_distance  # Detector is next to the tank, contact dose
        self.handling_det_x = box_x2 + self.handling_det_standoff_distance  # Detector is next to the tank, handling dose

        x_planes: list[float] = list(np.linspace(-box_x2, box_x2, nX))      # boundaries for gridgeometry
        y_planes: list[float] = list(np.linspace(-box_y2, box_y2, self.N_planes_box))
        z_planes: list[float] = list(np.linspace(-box_z2, box_z2, nX))

        x_planes_str: str = " ".join([f' {x:.5f}' for x in x_planes])
        y_planes_str: str = " ".join([f' {x:.5f}' for x in y_planes])
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/../{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
end

=mavric parm=(   )
{NOW} HotCellDoses, {self.sample_weight} g, layers {self.layers_thicknesses}
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
unit 10 
    ycylinder 10 {self.pin_IR} 2p{self.half_pin_length}
    ycylinder 20 {self.pin_IR + self.wall_thickness} 2p{self.half_pin_length + self.wall_thickness}
    cuboid 30 2p{self.half_pitch} 2p{self.half_pin_length + self.wall_thickness} 2p{self.half_pitch}

' pin with activated material
    media  1 1  10
' tube around the pin     
    media  10 1  20 -10
' dry air around the tubes    
    media 11 1  30 -20
boundary 30
    
global unit 1
    cuboid 1 2p{box_x2} 2p{box_y2} 2p{box_z2}
    array 7 1  place 1 1 1 {place_x} {place_y} {place_z}
    cuboid 99999  2p{box_x2 + self.box_a} 2p{box_y2 + self.box_a} 2p{box_z2 + self.box_a}  
    media 0 1 99999 -1
boundary 99999
end geometry

read array
    ara=7 nux={nX} nuy=1 nuz={nZ}  
        fill f10 end fill
end array

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
        xPlanes {x_planes_str} 
            {self.det_x - self.planes_xy_around_det} {self.det_x + self.planes_xy_around_det} 
        end
        yPlanes {y_planes_str} 
            {self.planes_xy_around_det} {-self.planes_xy_around_det} 
        end
        zPlanes {z_planes_str} end
    end gridGeometry
end definitions

read sources'''
        if self.neutron_intensity > 0.0:
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        cuboid  -{box_x2} {box_x2} -{box_y2} {box_y2} -{box_z2} {box_z2}
        mixture=1
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cuboid  -{box_x2} {box_x2} -{box_y2} {box_y2} -{box_z2} {box_z2}
        mixture=1
        eDistributionID=2
    end src
end sources

read importanceMap
   gridGeometryID=1'''
        if self.reuse_adjoint_flux:
            mavric_output += f'''
   adjointFluxes="{adjoint_flux_file}"'''
        else:
            mavric_output += f'''
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
   beckerMethod=2
   respWeighting'''
        mavric_output += f'''
'   xblocks=4
'   yblocks=4
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
