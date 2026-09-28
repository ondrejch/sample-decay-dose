import os
import shutil
import numpy as np

from sample_decay_dose.SampleDose import NOW, MAVRIC_NG_XSLIB, DAY_IN_SECONDS, OrigenFromTriton, run_scale_or_raise
from sample_decay_dose.utils import get_f71_volume, atom_dens_for_mavric
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
        layers_thicknesses is not used; the tube geometry comes from radiator_geometry.
    The source is the salt state at F71 position origen_from_triton.BURNED_MATERIAL_F71_position, read directly
    from the burned-material F71 without an ORIGEN decay step.
    Monaco's origensBinaryConcentrationFile distribution reads the emission spectra stored on the F71.
    The selected position must therefore carry gamma spectra (and neutron spectra for the neutron source),
    i.e., flags G and N in obiwan's DCGNAB column. ORIGEN writes them for cases run with gamma=yes and
    neutron=yes. TRITON depletion positions carry concentrations only (DC----), and mavric_deck() raises for them.
    Source strength assumption: the F71 spectra are totals over the F71 material volume (obiwan "volume" row).
    The radiator holds all_pins_volume of that salt, so the sources are scaled by
    all_pins_volume / V_F71 unless source_multiplier is set.
    """
    def __init__(self, _o: OrigenFromTriton = None):
        """ This reads decayed sample information from the Origen object """
        super().__init__(_o)  # Init DoseEstimator
        self.origen_from_triton: OrigenFromTriton = _o
        self.ORIGEN_dir = 'run'
        self.decay_days = 0.0               # no decay, just raw F71, copy the data from F71 file
        self.DECAYED_SAMPLE_days = self.decay_days
        self.DECAYED_SAMPLE_F71_file_name = _o.BURNED_MATERIAL_F71_file_name
        self.DECAYED_SAMPLE_F71_index: dict = _o.BURNED_MATERIAL_F71_index  # already read by OrigenFromTriton
        self.DECAYED_SAMPLE_F71_position = _o.BURNED_MATERIAL_F71_position
        # The OPUS integrals of an ORIGEN decay case do not describe the burned F71 used here. The neutron source
        # is included unless neutron_intensity is set to 0.0, which omits neutron transport explicitly.
        self.neutron_intensity = None
        self.source_multiplier: float | None = None  # None: all_pins_volume / F71 volume at the source position
        self.check_f71_spectra: bool = True  # raise if the F71 position lacks the needed emission spectra
        object.__setattr__(self, 'handling_det_x', None)  # set by mavric_deck(); not a manual det_x
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

    @property
    def beta_applies(self) -> bool:
        """ The salt sits inside the tube walls, so the bare-sample beta estimate does not apply """
        return False

    @property
    def f71_basename(self) -> str:
        """ The F71 file is copied into case_dir and referenced by its basename """
        return os.path.basename(self.DECAYED_SAMPLE_F71_file_name)

    def _f71_path(self) -> str:
        return os.path.join(self.cwd, self.DECAYED_SAMPLE_F71_file_name)

    def _validate_layers(self):
        """ The radiator uses exactly two materials: tube (mix 10) and the space between tubes (mix 11) """
        if len(self.layers_mats) != 2 or len(self.layers_temperature_K) != 2:
            raise ValueError(f"RadiatorBox needs two layers_mats and two layers_temperature_K (tube, between tubes), "
                             f"got {len(self.layers_mats)} and {len(self.layers_temperature_K)}")

    def _set_source_position(self):
        """ Source position and its time from the OrigenFromTriton object, which may have changed after __init__ """
        self.DECAYED_SAMPLE_F71_position = self.origen_from_triton.BURNED_MATERIAL_F71_position
        record: dict | None = self.DECAYED_SAMPLE_F71_index.get(self.DECAYED_SAMPLE_F71_position)
        if record is None:
            raise ValueError(f"Position {self.DECAYED_SAMPLE_F71_position} is not in the index of "
                             f"{self.DECAYED_SAMPLE_F71_file_name}")
        return record

    def _check_spectra(self, record: dict):
        """ Monaco reads the stored emission spectra, see the class docstring """
        if not self.check_f71_spectra:
            return
        flags: str = record.get('DCGNAB', '')
        needed: list[tuple[str, str]] = [('G', 'gamma')]
        if self._include_neutron_source():
            needed.append(('N', 'neutron'))
        missing: list[str] = [name for flag, name in needed if flag not in flags]
        if missing:
            raise ValueError(
                f"F71 position {self.DECAYED_SAMPLE_F71_position} of {self.DECAYED_SAMPLE_F71_file_name} has no "
                f"{' or '.join(missing)} emission spectra (DCGNAB flags '{flags}'), which the MAVRIC source needs. "
                f"Select a position from an ORIGEN case run with gamma=yes neutron=yes, or decay the material with "
                f"OrigenFromTriton.run_decay_sample(). Set neutron_intensity = 0.0 to omit the neutron source.")

    def _source_multiplier(self) -> float:
        """ Scales the F71 spectra (totals over the F71 material volume) to the radiator salt volume """
        if self.source_multiplier is not None:
            return self.source_multiplier
        f71_volume: float = get_f71_volume(self._f71_path(), self.DECAYED_SAMPLE_F71_position)
        if not f71_volume > 0.0:
            raise ValueError(f"Non-positive material volume {f71_volume} at F71 position "
                             f"{self.DECAYED_SAMPLE_F71_position}; set source_multiplier explicitly")
        return self.all_pins_volume / f71_volume

    def run_mavric(self, nmpi: int = 1):
        """ Writes Mavric inputs and runs the case """
        self._validate_layers()
        self.case_dir: str = f'run_MAVRIC_{NOW}_{self.all_pins_volume:.5}_cm-{self.decay_days:.5}_days'  # Directory to run the case
        self._set_source_position()

        if not os.path.isfile(self._f71_path()):
            raise FileNotFoundError("Expected decayed sample F71 file: \n" + self._f71_path())
        if not self.origen_from_triton.burned_atom_dens:
            raise ValueError("No burned-material composition, call read_burned_material() first")
        os.makedirs(os.path.join(self.cwd, self.case_dir), exist_ok=True)
        shutil.copy2(self._f71_path(), os.path.join(self.cwd, self.case_dir, self.f71_basename))
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            with open(self.SAMPLE_ATOM_DENS_file_name_MAVRIC, 'w') as f:  # write MAVRIC at-dens sample input
                f.write(atom_dens_for_mavric(self.origen_from_triton.burned_atom_dens, 1, self.sample_temperature_K))
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
        if not self.det_z:
            self.det_z = 0.0
        record: dict = self._set_source_position()
        self._check_spectra(record)
        source_time_days: float = float(record['time']) / DAY_IN_SECONDS
        multiplier: float = self._source_multiplier()
        adjoint_flux_file: str = self._reused_adjoint_flux_path(self.case_dir) if self.reuse_adjoint_flux else ''
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

        self._set_computed_detectors(box_x2 + self.det_standoff_distance,
                                     box_x2 + self.handling_det_standoff_distance)  # next to the tank

        # Grid planes over the radiator, plus the problem boundary and both detectors: Monaco stops when a particle
        # leaves the importance map inside the geometry
        x_planes: list[float] = list(np.linspace(-box_x2, box_x2, nX))      # boundaries for gridgeometry
        x_planes += [-box_x2 - self.box_a, box_x2 + self.box_a,
                     self.handling_det_x - self.planes_xy_around_det, self.handling_det_x + self.planes_xy_around_det]
        y_planes: list[float] = list(np.linspace(-box_y2, box_y2, self.N_planes_box))
        y_planes += [-box_y2 - self.box_a, box_y2 + self.box_a]
        z_planes: list[float] = list(np.linspace(-box_z2, box_z2, nZ))
        z_planes += [-box_z2 - self.box_a, box_z2 + self.box_a]

        x_planes_str: str = " ".join([f' {x:.5f}' for x in x_planes])
        y_planes_str: str = " ".join([f' {x:.5f}' for x in y_planes])
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.f71_basename} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
end

=mavric parm=(   )
{NOW} RadiatorBox, F71 position {self.DECAYED_SAMPLE_F71_position} (t = {source_time_days:.6g} d), salt {self.all_pins_volume:.6g} cm3
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

        if self._include_neutron_source(note=True):
            mavric_output += f'''
    distribution 1
        title="Salt at F71 position {self.DECAYED_SAMPLE_F71_position}, t = {source_time_days:.6g} days, neutrons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 1 end
        filename="{self.f71_basename}"
    end distribution'''

        mavric_output += f'''
    distribution 2
        title="Salt at F71 position {self.DECAYED_SAMPLE_F71_position}, t = {source_time_days:.6g} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.DECAYED_SAMPLE_F71_position} 5 end
        filename="{self.f71_basename}"
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
        if self._include_neutron_source():
            mavric_output += f'''
    src 1
        title="Sample neutrons"
        neutron
        useNormConst
        multiplier={multiplier}
        cuboid  -{box_x2} {box_x2} -{box_y2} {box_y2} -{box_z2} {box_z2}
        mixture=1
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        multiplier={multiplier}
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
   macromaterial
       mmsubcell=8
   end macromaterial 
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
