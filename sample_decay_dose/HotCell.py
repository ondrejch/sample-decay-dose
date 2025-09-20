import os
import re

from sample_decay_dose.SampleDose import NOW, MAVRIC_NG_XSLIB, HandlingContactDoseEstimatorGenericTank, Origen, get_cyl_r


class HotCellDoses(HandlingContactDoseEstimatorGenericTank):
    """ MAVRIC calculation of rem/h doses from the decayed sample in a cubical hotcell """
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

    def mavric_deck(self) -> str:
        """ MAVRIC dose calculation input file """
        if not self.det_z:
            self.det_z = 0.0
        adjoint_flux_file: str = os.path.join(self.cwd,
                                              self.MAVRIC_input_file_name.replace('.inp', '.adjoint.dff'))
        self.cyl_r = get_cyl_r(self.sample_volume)
        self.sample_h2 = self.cyl_r / 2.0       # sample is a square cylinder
        sample_r: float = self.cyl_r            # sample outer layer [cm]
        sample_h2: float = self.sample_h2
        box_xy2: float = 0.0                    # XY box around
        box_z2: float = 0.0                     # Z box around

        mavric_output = f'''
=shell
cp -r ${{INPDIR}}/{self.DECAYED_SAMPLE_F71_file_name} .
cp -r ${{INPDIR}}/{self.SAMPLE_ATOM_DENS_file_name_MAVRIC} .
'cp -r ${{INPDIR}}/{adjoint_flux_file} .
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
        k: int = 0
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
        self.det_x = box_xy2 + self.det_standoff_distance  # Detector is next to the tank, contact dose
        self.handling_det_x = box_xy2 + self.handling_det_standoff_distance  # Detector is next to the tank, handling dose
        xy_planes_str: str = " ".join([f' {x:.5f}' for x in xy_planes])
        z_planes.append(self.det_z + self.planes_xy_around_det)
        z_planes.append(self.det_z - self.planes_xy_around_det)
        z_planes_str: str = " ".join([f' {x:.5f}' for x in z_planes])
        mavric_output += f'''
    cuboid 99999  4p {box_xy2 + self.box_a} 2p {box_z2 + self.box_a}  
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
        if self.neutron_intensity > 0.0:
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

    def get_responses(self):
        """ Reads over the MAVRIC output and returns responses for rem/h doses
            Note, this version uses different format of self.responses """
        if not os.path.isfile(self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name):
            raise FileNotFoundError(
                "Expected decayed sample MAVRIC output file: \n" + self.cwd + '/' + self.case_dir + '/' + self.MAVRIC_out_file_name)
        os.chdir(self.cwd + '/' + self.case_dir)

        with open(self.MAVRIC_out_file_name, 'r') as f:
            d: str = f.read()
        rs: list = re.findall(r'Point Detector (\d).  (\w+) detector\n.*\n.*\n.*\n.*\n.*\s+response (\d)\s+'
                              r'([+-]?\d+\.\d+[Ee]?[+-]?\d+)\s+([+-]?\d+\.\d+[Ee]?[+-]?\d+)?'
                              r'\s+([+-]?\d+\.\d+[Ee]?[+-]?\d+)?', d)
        self.responses = {
            t[0]: {'particle': t[1], 'pid': t[2], 'value': float(t[3]), 'stdev': 0.0 if t[4] == '' else float(t[4])} for
            t in rs}

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
            print(self.sample_weight, r1['value'], r1['stdev'], r2['value'], r2['stdev'], r5['value'], r5['stdev'],
                  r6['value'], r6['stdev'])
