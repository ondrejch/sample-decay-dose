import os

from sample_decay_dose.SampleDose import NOW, MAVRIC_NG_XSLIB, HandlingContactDoseEstimatorGenericTank, Origen
from sample_decay_dose.utils import get_cyl_r


class HotCellDoses(HandlingContactDoseEstimatorGenericTank):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a cubical hotcell.
    The sample is a square cylinder at the origin. The layers are nested cuboids centred on the sample.
    layers_thicknesses[0] is the half-width of the innermost cuboid, i.e., the cell interior measured from the
    sample centre. It must be at least the sample radius. Each later entry is the thickness of the next cuboid shell.
    With no layers the sample is bare, and the detectors are measured from the sample surface.
    Detectors and responses are those of HandlingContactDoseEstimatorGenericTank.
    reuse_adjoint_flux=True skips the Denovo adjoint calculation and reads the adjoint flux from adjoint_flux_file,
    by default <case_dir>/my_dose.adjoint.dff written by an earlier run of the same case directory.
    The file must come from the same geometry, layers, detector positions, and grid.
    """
    def __init__(self, _o: Origen = None):
        """ This reads decayed sample information from the Origen object """
        self.sample_h2: (None, float) = None  # half-height of the sample
        self.det_z: float = 0  # z location of the detector
        self.handling_det_x: (None, float) = None
        super().__init__(_o)  # Init DoseEstimator
        self.det_standoff_distance = 0.1  # [cm] contact dose
        self.handling_det_standoff_distance = 30.0  # [cm] handing dose
        self.box_a: float = self.handling_det_standoff_distance + 10.0  # outer bounding box [cm]
        self.reuse_adjoint_flux: bool = False
        self.adjoint_flux_file: str | None = None  # adjoint flux to reuse; None: the MAVRIC output in case_dir

    def _reused_adjoint_flux_path(self, case_dir: str) -> str:
        """ Absolute path of the adjoint flux file to reuse. MAVRIC writes <input name>.adjoint.dff next to
        its output, in case_dir. Raises FileNotFoundError when the file does not exist. """
        default_file: str = self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff')
        path: str = os.path.abspath(self.adjoint_flux_file or os.path.join(self.cwd, case_dir, default_file))
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"reuse_adjoint_flux is set, but the adjoint flux file {path} does not exist. Run the case once with "
                f"reuse_adjoint_flux=False, or set adjoint_flux_file to the {default_file} of a run with the same "
                f"geometry, layers, and detector positions.")
        return path

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        self._validate_layers()
        if not self.det_z:
            self.det_z = 0.0
        adjoint_flux_file: str = self._reused_adjoint_flux_path(self._layer_case_dir()) \
            if self.reuse_adjoint_flux else ''
        self.cyl_r = get_cyl_r(self.sample_volume)
        self.sample_h2 = self.cyl_r             # sample is a square cylinder
        sample_r: float = self.cyl_r            # sample outer layer [cm]
        sample_h2: float = self.sample_h2
        if self.layers_thicknesses and self.layers_thicknesses[0] < max(sample_r, sample_h2):
            raise ValueError(f"HotCellDoses: layers_thicknesses[0] = {self.layers_thicknesses[0]} cm is the half-width "
                             f"of the innermost cuboid and must be at least the sample radius {sample_r:.5g} cm")
        box_xy2: float = 0.0                    # XY box around
        box_z2: float = 0.0                     # Z box around

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
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
' helium 2 end
end comp

read geometry
global unit 1
    cylinder 1 {sample_r} 2p {sample_h2}
    media 1 1 1
'''
        xy_planes: list[float] = [-sample_r, sample_r]      # list of XY boundaries for gridgeometry
        z_planes: list[float] = [-sample_h2, sample_h2]     # list of Z boundaries for gridgeometry
        for k in range(len(self.layers_mats)):
            box_xy2 += self.layers_thicknesses[k]
            box_z2 += self.layers_thicknesses[k]
            mavric_output += f'''
    cuboid {k + 2} 4p {box_xy2} 2p {box_z2}   
    media {k + 10}  1 -{k + 1} {k + 2}'''
            xy_planes.append(box_xy2)
            xy_planes.append(-box_xy2)
            z_planes.append(box_z2)
            z_planes.append(-box_z2)
        if not self.layers_mats:  # bare sample: the detectors and the boundary are measured from the sample surface
            box_xy2 = sample_r
            box_z2 = sample_h2
        self._set_computed_detectors(box_xy2 + self.det_standoff_distance,
                                     box_xy2 + self.handling_det_standoff_distance)  # next to the tank
        xy_planes_str: str = " ".join([f' {x:.5f}' for x in xy_planes])
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {box_xy2 + self.box_a} 2p {box_z2 + self.box_a}  
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

        if self._include_neutron_source(note=True):
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
        xLinear {self.N_planes_box}  {sample_r}  {box_xy2 + self.box_a}
        xLinear {self.N_planes_box} -{sample_r} -{box_xy2 + self.box_a}
        yLinear {self.N_planes_box}  {sample_r}  {box_xy2 + self.box_a} 
        yLinear {self.N_planes_box} -{sample_r} -{box_xy2 + self.box_a} 
        zLinear {self.N_planes_box}  {sample_h2}  {box_z2 + self.box_a}
        zLinear {self.N_planes_box} -{sample_h2} -{box_z2 + self.box_a}
        xLinear {self.N_planes_cyl} -{sample_r} {sample_r}          
        yLinear {self.N_planes_cyl} -{sample_r} {sample_r}
        zLinear {self.N_planes_cyl} -{sample_h2} {sample_h2}
        xPlanes {xy_planes_str} {self.det_x - self.planes_xy_around_det} {self.det_x + self.planes_xy_around_det} end
        yPlanes {xy_planes_str} {self.planes_xy_around_det} {-self.planes_xy_around_det} end
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
        cylinder {self.cyl_r} -{sample_h2} {sample_h2}
        eDistributionID=1
    end src'''

        mavric_output += f'''
    src 2
        title="Sample photons"
        photon
        useNormConst
        cylinder {self.cyl_r} -{sample_h2} {sample_h2}
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
