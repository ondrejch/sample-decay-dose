#!/bin/env python3
"""
Prototype: activation and contact dose of a reactor-cavity concrete lid.

The class chains the two-step workflow from concrete_irrad/README.md:

1. ORIGEN irradiation (cavity neutron spectrum via an alpha-library F33) and
   decay of the two modeled lid regions, the bottom concrete slab and the
   homogenized concrete-rebar layer.
2. MAVRIC fixed-source photon transport with point-detector dose at the top
   surface of the lid, using ORIGEN spectral distributions read from the F71
   files.

Deck structure follows the conventions of sample_decay_dose.SampleDose
(OrigenIrradiation chaining, DoseEstimator useNormConst sources and response
parsing). The ORIGEN deck runs under SCALE 6.3.3, checked with a stand-in F33.
The MAVRIC deck has not been run yet, so check its first output before relying
on the dose.
"""
import os
import re
import shutil

from sample_decay_dose.SampleDose import NOW, MAVRIC_NG_XSLIB, run_scale_or_raise
from sample_decay_dose.utils import atom_dens_for_mavric, atom_dens_for_origen, \
    get_burned_nuclide_atom_dens

from concrete_irrad.concrete_rebar_mixer import ConcreteRebarMixer, expand_to_nuclides, \
    wt_fractions_to_atom_densities


class ConcreteLidContactDose:
    """
    Activation and contact gamma dose for a three-slice reactor-cavity lid.

    The modeled regions are the bottom concrete slab and the homogenized
    concrete-rebar mixed layer; the upper concrete is excluded per the README.
    Each region is activated in its own ORIGEN case pair. Within one ORIGEN
    input the cases chain: each irradiation case carries an explicit fresh
    material, and each following decay case continues from the preceding case.
    """

    def __init__(self,
                 name: str = 'lid',
                 length_cm: float = 100.0,
                 width_cm: float = 100.0,
                 bottom_thickness_cm: float = 5.0,
                 mixer: ConcreteRebarMixer | None = None,
                 irradiation_lib_f33: str = '',
                 irradiation_lib_pos: int | None = None,
                 irradiation_flux: float = 1.0e10,
                 irradiation_days: float = 365.24,
                 decay_days: float = 30.0,
                 f71_position: int = 12,
                 contact_standoff_cm: float = 0.1,
                 sample_temperature_K: float = 300.0,
                 histories_per_batch: int = 50000,
                 batches: int = 10):
        if length_cm <= 0 or width_cm <= 0 or bottom_thickness_cm <= 0:
            raise ValueError("Lid dimensions must be positive")
        if irradiation_days <= 0:
            raise ValueError(f"irradiation_days must be positive, got {irradiation_days}")
        if decay_days < 0:
            raise ValueError(f"decay_days must be non-negative, got {decay_days}")
        self.name: str = name
        self.length_cm: float = float(length_cm)
        self.width_cm: float = float(width_cm)
        self.mixer: ConcreteRebarMixer = mixer \
            if mixer is not None else ConcreteRebarMixer(rebar_diameter_cm=1.59,
                                                         spacing_x_cm=30.48, spacing_y_cm=15.24)
        self.bottom_thickness_cm: float = float(bottom_thickness_cm)
        self.mixed_thickness_cm: float = self.mixer.rebar_diameter_cm
        self.irradiation_lib_f33: str = irradiation_lib_f33
        self.irradiation_lib_pos: int | None = irradiation_lib_pos
        self.irradiation_flux: float = float(irradiation_flux)
        self.irradiation_days: float = float(irradiation_days)
        self.irradiation_steps: int = max(3, f71_position - 2)
        self.decay_days: float = float(decay_days)
        self.f71_position: int = int(f71_position)
        self.contact_standoff_cm: float = float(contact_standoff_cm)
        self.sample_temperature_K: float = float(sample_temperature_K)
        self.histories_per_batch: int = int(histories_per_batch)
        self.batches: int = int(batches)

        self.debug: int = 3
        self.cwd: str = os.getcwd()
        self.case_dir: str = f"run_{name}_{irradiation_days:.5}d-{decay_days:.5}d"
        self.case_dir_mavric: str = self.case_dir + '_MAVRIC'
        self.ORIGEN_input_file_name: str = f'{name}_activation.inp'
        self.MAVRIC_input_file_name: str = f'{name}_dose.inp'

        self.regions: dict = {
            'bottom_slab': {
                'thickness': self.bottom_thickness_cm,
                'volume': self.length_cm * self.width_cm * self.bottom_thickness_cm,
                'adens_file': f'{name}_bottom_adens_origen.inp',
                'f71': f'{name}_bottom.f71',
                'mavric_comp_file': f'{name}_bottom_adens_mavric.inp',
                'mix_id': 2},
            'mixed_layer': {
                'thickness': self.mixed_thickness_cm,
                'volume': self.length_cm * self.width_cm * self.mixed_thickness_cm,
                'adens_file': f'{name}_mixed_adens_origen.inp',
                'f71': f'{name}_mixed.f71',
                'mavric_comp_file': f'{name}_mixed_adens_mavric.inp',
                'mix_id': 1}}
        self.fresh_adens: dict = {
            'bottom_slab': expand_to_nuclides(
                wt_fractions_to_atom_densities(self.mixer.concrete_wt_fractions,
                                               self.mixer.concrete_density_g_cc)),
            'mixed_layer': self.mixer.get_atom_densities()}
        self.decayed_adens: dict = {}
        self.responses: dict = {}

    @property
    def MAVRIC_out_file_name(self) -> str:
        return self.MAVRIC_input_file_name.replace('inp', 'out')

    @property
    def decayed_f71_position(self) -> int:
        """ F71 position of the decayed state on each decay-case file

        Position 1 is the decay-case start (step 0), which ORIGEN stores with
        its gamma spectrum. It holds the end-of-irradiation inventory, so zero
        cooling reads it directly. Positive cooling reads the last position,
        f71_position, at t = decay_days.
        """
        return 1 if self.decay_days == 0 else self.f71_position

    def _decay_time_grid(self) -> str:
        """ Strictly increasing ORIGEN decay-case time list [d] with f71_position - 3 interpolated steps """
        interp_steps: int = self.f71_position - 3
        if self.decay_days > 0:
            return f'{interp_steps}L {min(0.0001, self.decay_days / 1e3)} {self.decay_days}'
        # ORIGEN requires times after start=0. Zero cooling therefore keeps a nominal
        # 1-day grid and reads the decay-case start instead (decayed_f71_position).
        return f'{interp_steps}L 0.0001 1.0'

    def origen_deck(self) -> str:
        """ Chained ORIGEN irradiation-deck cases for both lid regions

        Each irradiation case converts the ATOMS-PER-BARN-CM input to moles with
        `volume`, the full region volume (length x width x thickness). The decay
        case continues that whole-region inventory, so the F71 photon spectrum
        is the emission rate of the full region [photons/s].
        """
        irr_min: float = min(0.0001, self.irradiation_days / 1e3)
        if self.f71_position - 3 < 1:
            raise ValueError("Too few time steps for the requested f71 position")
        lib_pos_line: str = '' if self.irradiation_lib_pos is None else f'pos={self.irradiation_lib_pos} '
        if not self.irradiation_lib_f33:
            raise ValueError("irradiation_lib_f33 with the cavity neutron spectrum is required")
        lib_name: str = os.path.basename(self.irradiation_lib_f33)

        deck: str = f'''=shell
cp -r ${{INPDIR}}/*.inp .
cp -r ${{INPDIR}}/{lib_name} .
end

=origen
' {self.__class__.__name__} {NOW}
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
'''
        for region in ('bottom_slab', 'mixed_layer'):
            r: dict = self.regions[region]
            deck += f'''case(irr_{region}) {{
    lib {{
        file="{lib_name}"
        {lib_pos_line.rstrip()}
    }}
    mat{{
        iso=[
<{r['adens_file']}
]
        units=ATOMS-PER-BARN-CM
        volume={r['volume']:.6e}
    }}
    time {{
        units=DAYS
        start=0
        t=[{self.irradiation_steps - 2}L {irr_min} {self.irradiation_days}]
    }}
    flux=[{self.irradiation_steps}R {self.irradiation_flux}]
    save {{
        file="{r['f71'].replace('.f71', '_irr.f71')}"
    }}
}}
case(dec_{region}) {{
    gamma=yes
    neutron=yes
    beta=yes
    lib {{
        file="end7dec"
    }}
    time {{
        units=DAYS
        start=0
        t=[{self._decay_time_grid()}]
    }}
    save {{
        file="{r['f71']}"
    }}
}}
'''
        deck += '''end
'''
        return deck

    def _stage_library(self):
        """ Copies the spectrum F33 into the ORIGEN case directory, where the deck names it by basename """
        src: str = os.path.join(self.cwd, self.irradiation_lib_f33)  # absolute paths pass through unchanged
        if not os.path.isfile(src):
            raise FileNotFoundError(f"Spectrum library F33 not found: {src}")
        dst: str = os.path.join(self.cwd, self.case_dir, os.path.basename(src))
        if not (os.path.exists(dst) and os.path.samefile(src, dst)):
            shutil.copy2(src, dst)

    def write_inputs(self):
        """ Writes fresh atom densities, the ORIGEN deck, and a copy of the F33 into the case directory """
        deck: str = self.origen_deck()
        os.makedirs(os.path.join(self.cwd, self.case_dir), exist_ok=True)
        self._stage_library()
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            for region, adens in self.fresh_adens.items():
                with open(self.regions[region]['adens_file'], 'w') as f:
                    f.write(atom_dens_for_origen(adens))
            with open(self.ORIGEN_input_file_name, 'w') as f:
                f.write(deck)
        finally:
            os.chdir(self.cwd)

    def run_activation(self, nmpi: int = 1):
        """ Runs the ORIGEN activation-decay chain and reads back decayed atom densities """
        if not os.path.isfile(os.path.join(self.cwd, self.case_dir, self.ORIGEN_input_file_name)):
            self.write_inputs()
        os.chdir(os.path.join(self.cwd, self.case_dir))
        try:
            if self.debug > 0:
                print(f"ORIGEN: irradiating {self.irradiation_days} d at "
                      f"{self.irradiation_flux:.3e} n/cm2/s, decaying {self.decay_days} d")
            run_scale_or_raise(self.ORIGEN_input_file_name, nmpi=nmpi, context_dir=self.case_dir)
            for region in self.regions:
                self.decayed_adens[region] = get_burned_nuclide_atom_dens(
                    self.regions[region]['f71'], self.decayed_f71_position)
        finally:
            os.chdir(self.cwd)

    def mavric_deck(self) -> str:
        """ MAVRIC photon contact-dose deck over the two activated lid slabs

        Each src carries useNormConst, so its strength equals the normalization
        constant of its F71 photon distribution. That constant is the total
        photon emission rate [photons/s] of the region inventory at the decayed
        position. The F71 inventory covers the full region volume (see
        origen_deck), and each src cuboid spans that same region, so no
        multiplier is applied. A SCALE 6.3.3 ORIGEN check confirmed that doubling
        the case `volume` doubles the F71 photon total.

        Geometry bodies and source cuboids list bounds max-first,
        +X -X +Y -Y +Z -Z, per the KENO-VI body definition and the Monaco
        source-shape table.
        """
        hx: float = self.length_cm / 2.0
        hy: float = self.width_cm / 2.0
        z_bot_lo, z_bot_hi = 0.0, self.bottom_thickness_cm
        z_mix_hi: float = z_bot_hi + self.mixed_thickness_cm
        box_pad: float = 10.0 + self.contact_standoff_cm
        bx, by, bz = hx + box_pad, hy + box_pad, z_mix_hi + box_pad  # boundary half-extents
        det_z: float = z_mix_hi + self.contact_standoff_cm

        comp_includes: str = ''
        distributions: str = ''
        sources: str = ''
        src_id: int = 1
        for region in ('mixed_layer', 'bottom_slab'):
            r: dict = self.regions[region]
            comp_includes += f'<{r["mavric_comp_file"]}\n'
            dist_id: int = src_id + 10
            distributions += f'''
    distribution {dist_id}
        title="Decayed {region} after {self.decay_days} days, photons"
        special="origensBinaryConcentrationFile"
        parameters {self.decayed_f71_position} 5 end
        filename="{r['f71']}"
    end distribution'''
            z_lo, z_hi = (z_bot_lo, z_bot_hi) if region == 'bottom_slab' else (z_bot_hi, z_mix_hi)
            sources += f'''
    src {src_id}
        title="Decay photons, {region}"
        photon
        useNormConst
        cuboid {hx} {-hx} {hy} {-hy} {z_hi} {z_lo}
        eDistributionID={dist_id}
    end src'''
            src_id += 1

        return f'''=shell
cp -r ${{INPDIR}}/*.f71 .
cp -r ${{INPDIR}}/*_adens_mavric.inp .
end

=mavric parm=(   )
{NOW} {self.name}, concrete lid contact dose after {self.decay_days} d
{MAVRIC_NG_XSLIB}

read parameters
    randomSeed=0000000100000001
    ceLibrary="ce_v7.1_endf.xml"
    photons
    perBatch={self.histories_per_batch} batches={self.batches}
end parameters

read comp
{comp_includes}end comp

read geometry
global unit 1
    cuboid 11 {hx} {-hx} {hy} {-hy} {z_bot_hi} {z_bot_lo}
    media 2 1 11
    cuboid 12 {hx} {-hx} {hy} {-hy} {z_mix_hi} {z_bot_hi}
    media 1 1 12
    cuboid 99999 {bx} {-bx} {by} {-by} {bz} {-bz}
    media 0 1 99999 -11 -12
boundary 99999
end geometry

read definitions
    location 1
        position 0.0 0.0 {det_z}
    end location
    response 2
        title="ANSI standard (1991) flux-to-dose-rate factors [rem/h], photons"
        doseData=9505
    end response
{distributions}

    gridGeometry 1
        title="Grid over the lid, detector at top center"
        xLinear 5 {-hx} {hx}
        yLinear 5 {-hy} {hy}
        zLinear 8 {z_bot_lo} {z_mix_hi + box_pad}
        xPlanes {-hx} {hx} end
        yPlanes {-hy} {hy} end
        zPlanes {z_bot_lo} {z_bot_hi} {z_mix_hi} {det_z} end
    end gridGeometry
end definitions

read sources{sources}
end sources

read importanceMap
    gridGeometryID=1
    adjointSource 1
        locationID=1
        responseID=2
    end adjointSource
    beckerMethod=2
end importanceMap

read tallies
    pointDetector 2
        title="photon contact dose at lid top center"
        photon
        locationID=1
        responseID=2
    end pointDetector
end tallies

end data
end
'''

    def run_mavric(self, nmpi: int = 1):
        """ Writes MAVRIC inputs from the decayed inventories and runs the transport case """
        if not self.decayed_adens:
            raise RuntimeError("Run run_activation() before run_mavric()")
        mavric_dir: str = os.path.join(self.cwd, self.case_dir_mavric)
        os.makedirs(mavric_dir, exist_ok=True)
        os.chdir(mavric_dir)
        try:
            for region, r in self.regions.items():
                src: str = os.path.join(self.cwd, self.case_dir, r['f71'])
                if not os.path.isfile(src):
                    raise FileNotFoundError(f"Expected decayed F71 file: {src}")
                shutil.copy2(src, os.path.join(mavric_dir, r['f71']))  # the deck copies ${INPDIR}/*.f71
                with open(r['mavric_comp_file'], 'w') as f:
                    f.write(atom_dens_for_mavric(self.decayed_adens[region], r['mix_id'],
                                                 self.sample_temperature_K))
            with open(self.MAVRIC_input_file_name, 'w') as f:
                f.write(self.mavric_deck())
            if self.debug > 0:
                print(f"MAVRIC: running case {self.case_dir_mavric}/{self.MAVRIC_input_file_name}")
            run_scale_or_raise(self.MAVRIC_input_file_name, nmpi=nmpi, context_dir=self.case_dir_mavric)
        finally:
            os.chdir(self.cwd)

    def get_responses(self):
        """ Parses the photon point-detector response and its standard deviation [rem/h] from the MAVRIC output

        A MAVRIC summary line reads "response 2  <average value>  <standard
        deviation>  <relative uncertainty> ...". The first two numbers are
        stored as 'value' and 'stdev', both in rem/h.
        """
        out_path: str = os.path.join(self.cwd, self.case_dir_mavric, self.MAVRIC_out_file_name)
        if not os.path.isfile(out_path):
            raise FileNotFoundError(f"Expected MAVRIC output file: {out_path}")
        tally_sep: str = 'Final Tally Results Summary'
        is_in_tally: bool = False
        response_re = re.compile(
            r'\bresponse\s+(\d+)\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)?')
        self.responses = {}
        with open(out_path, 'r') as f:
            for line in f.read().splitlines():
                if tally_sep in line:
                    is_in_tally = True
                    continue
                if not is_in_tally:
                    continue
                m = response_re.search(line)
                if m:
                    self.responses[m.group(1)] = {'value': float(m.group(2)),
                                                  'stdev': 0.0 if m.group(3) in (None, '') else float(m.group(3))}
        if '2' not in self.responses:
            raise RuntimeError(f"Failed to parse photon response from {out_path}")

    @property
    def contact_dose(self) -> dict:
        """ Contact photon dose at the lid top center [rem/h] with its absolute 1-sigma 'stdev' [rem/h] """
        if '2' not in self.responses:
            raise RuntimeError("No responses parsed; call get_responses() first")
        r: dict = self.responses['2']
        return {'value': r['value'], 'stdev': r['stdev']}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Concrete-lid activation and contact dose prototype")
    parser.add_argument('--run', action='store_true', help="Execute ORIGEN and MAVRIC via scalerte")
    parser.add_argument('--lib', default='',
                        help='Alpha-library F33 with the cavity neutron spectrum (required with --run)')
    args = parser.parse_args()
    if args.run and not args.lib:
        parser.error("--run requires --lib <f33>")

    lib: str = args.lib if args.lib else 'cavity_spectrum.f33'
    calc = ConcreteLidContactDose(name='lid', irradiation_lib_f33=lib)
    print(f"Case directory : {calc.case_dir}")
    print(f"Spectrum F33   : {lib}{'' if args.lib else ' (placeholder; pass --lib for real runs)'}")
    if args.run:
        calc.write_inputs()
        calc.run_activation()
        print("Decayed nuclides (top 10 per region):")
        for region, adens in calc.decayed_adens.items():
            print(f"  {region}: {list(adens.items())[:10]}")
        calc.run_mavric()
        calc.get_responses()
        dose: dict = calc.contact_dose
        print(f"\nContact photon dose after {calc.decay_days} d: "
              f"{dose['value']:.4e} +/- {dose['stdev']:.2e} rem/h")
    else:
        print(calc.origen_deck())
        print(calc.mavric_deck())
