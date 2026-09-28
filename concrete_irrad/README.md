# How-to: Dose from an irradiated reactor-cavity concrete lid

This directory collects a workflow that computes the dose rate from
the concrete lid of a reactor cavity after irradiation. The lid is activated by
neutrons from the cavity during operation, and the dose is produced by decay
gamma rays from the activated concrete and rebar.

## Workflow

1. **Activation (ORIGEN).** Irradiate the lid material under the cavity neutron
   spectrum and follow the radioactive decay of the activation products to the
   time of interest. The result is a set of nuclide atom densities for each
   modeled region.
2. **Dose estimate (MAVRIC).** Place the decayed source terms in a fixed-source
   Monaco/MAVRIC model of the lid and compute the gamma dose rate at the point
   of interest, for example at the lid surface or in a hot-cell above it.

## Start here

1. Export the SCALE installation: `export SCALE_BIN=/opt/scale6.3.2-mpi/bin` (adjust to your site).
2. Edit the `EDIT ME` block in `run_contact_dose.py`: lid span, bottom slab
   thickness, rebar diameter and grid spacings, cavity flux, and irradiation
   duration.
3. Obtain the cavity neutron spectrum as an ORIGEN alpha-library F33 and pass
   it with `--lib`.
4. Run one cooling time:
   `PYTHONPATH=. python concrete_irrad/run_contact_dose.py --lib cavity_spectrum.f33 --decay-days 30`
   Add `--nmpi N` to pass N MPI tasks to each SCALE run.
5. Read the contact photon dose [rem/h] and its 1-sigma standard deviation
   [rem/h] from stdout, or scan several cooling times with
   `--scan 1,30,90,365`, which writes `lid_contact_doses.csv`. A cooling time
   of 0 gives the dose at shutdown.

Case directories `run_*/` hold the ORIGEN decks, outputs, F71 files, and the
MAVRIC deck and output for inspection.

## Getting the spectrum library (F33)

ORIGEN irradiation needs a transition-matrix library weighted by the neutron
spectrum at the lid. Two standard routes exist. COUPLE builds an alpha library
from a group spectrum extracted from the cavity shielding model; see the SCALE
manual chapter on material specification and cross-section processing for the
COUPLE input format. A TRITON or Polaris F33 is reusable only when its
spectrum represents the cavity environment. Check any candidate file with
`obiwan info <file>.f33` before use.

## Geometry

The lid is a rectangular brick. The model splits the brick into three
horizontal slices stacked along the vertical axis:

1. **Bottom concrete slab.** Pure concrete between the cavity-facing surface
   and the first layer of reinforcement.
2. **Concrete-rebar mixed layer.** A single slice containing the first mat of
   reinforcement bars running in two orthogonal directions. The slice thickness
   equals the bar diameter, so the mat fits inside one slice. The bars and the
   surrounding concrete are homogenized into a single material whose atom
   densities come from `concrete_rebar_mixer.py`.
3. **Upper concrete slab.** The remaining concrete above the mixed layer. This
   region contributes negligibly to the dose rate and is excluded from the
   quantitative model.

## Material model

The mixed-layer composition follows from simple volume averaging. In a unit
cell of area `spacing_x * spacing_y` and height equal to the bar diameter `d`,
the cell contains one x-running bar segment of length `spacing_x` and one
y-running bar segment of length `spacing_y`. The two segments cross once per
cell. Their shared volume is the Steinmetz bicylinder of two perpendicular
cylinders of diameter `d`, which is `2 d^3 / 3`. The steel volume counts it
once:

```
V_steel = pi * d^2 * (spacing_x + spacing_y) / 4 - 2 * d^3 / 3
f_steel = V_steel / (spacing_x * spacing_y * d)
```

Parallel bars overlap when `d` exceeds a spacing, and the mixer raises an
error for that geometry. At the defaults (`d = 1.59` cm on a 30.48 x 15.24 cm
grid) the steel volume fraction is 0.1193 and the mixture density is 2.962
g/cm^3. Counting the crossing twice would give 0.1229 and 2.982 g/cm^3.

Each constituent keeps its own density and elemental composition. The mixture
mass density is `f_steel * rho_steel + (1 - f_steel) * rho_concrete`, and the
element number densities are the volume-weighted sums converted through natural
isotopic abundances.

The concrete composition is PNNL-15870 Rev. 1 material #96, "Concrete,
Ordinary (NBS 04)". Its tabulated weight fractions sum to 0.999993, and the
code normalizes them to unity. The concrete density is 2.30 g/cm^3, the
ANSI/ANS-6.4.3 value for ordinary concrete. PNNL-15870 lists 2.35 g/cm^3 for
the same composition.

## Activation-relevant impurities

Long-lived gamma emitters dominate the decay dose from irradiated concrete, and
they originate from trace elements rather than from the major constituents.
Kinno et al. (2002) showed that Co-60, Eu-152, Eu-154, and Cs-134 account for
99--100% of the total activity of ordinary concrete shields. The mixer
therefore carries explicit trace impurities in both materials.

The rebar steel carries a Co-59 impurity that produces the Co-60 source term.
It is set through `steel_co59_wt_fraction` with a default of 0.25 wt%, the
middle of the 0.1--0.4 wt% design band. Measured carbon steels contain less
cobalt: 93--151 ppm across the carbon-steel samples surveyed by Radulescu and
Banerjee (2020), and about 200 micro-g/g in Hiroshima structural steels (Kerr
et al., DS02). Values in the design band are therefore conservative for
ordinary rebar. The mixer accepts any value from 0 to 1 wt%, so measured
cobalt levels can be entered directly. It prints a note for values outside the
design band.

In the Japanese shielding-concrete survey that Alhajali et al. (2016) compiled
from Suzuki et al. (2001), ordinary concrete contains 0.16--21.9 ppm cobalt and
0.049--1.08 ppm europium. Additional measurements give Co 4.5 ppm with Eu 0.098
ppm (Kakinuma et al., 2007) and Co 3.5 ppm with Eu 0.85 ppm for barite
concrete (Gaudry and Delmas, 2007). A United States decommissioning study used
about 1 ppm Co-59 and 0.01 ppm total Eu as natural levels (NSTec, 2007).
Yoshida et al. (2020) measured facility-to-facility variations at the same
order of magnitude by neutron activation analysis. The defaults are Co 10 ppm
and Eu 1.0 ppm. The Co default lies inside the surveyed range and above the
single measurements listed here. The Eu default lies near the top of the
surveyed range, which puts the Eu-152 and Eu-154 source above most of the
listed concretes. Both defaults are overridable through
`concrete_impurities_wt`. Natural europium expands to Eu-151 and Eu-153, the
precursors of Eu-152 and Eu-154.

## Source term and file staging

- `write_inputs()` copies the F33 named by `--lib` into the ORIGEN case
  directory. The deck's shell block copies it into the SCALE working
  directory, and the deck names it by basename.
- Each ORIGEN irradiation case enters the atom densities of one region
  together with the full region volume (length x width x thickness). ORIGEN
  converts them to moles for the whole region, and the decay case continues
  that inventory. The F71 photon spectrum is therefore the emission rate of the
  whole region in photons/s. A SCALE 6.3.3 ORIGEN check confirmed that doubling
  the case volume doubles the F71 photon total.
- `run_mavric()` copies the decayed F71 files into the MAVRIC case directory.
  Each MAVRIC source uses `useNormConst`, which sets its strength to the
  normalization constant of its F71 photon distribution. Each source cuboid
  spans its full region, so each source carries the absolute emission rate of
  its region.
- ORIGEN requires decay times after the case start. For zero cooling the decay
  case therefore keeps a nominal 1-day grid, and the model reads F71 position
  1. That position is the decay-case start, which holds the end-of-irradiation
  inventory and its photon spectrum.
- Geometry bodies and source cuboids list their bounds max-first,
  `+X -X +Y -Y +Z -Z`. This follows the KENO-VI body definition and the Monaco
  source-shape table.
- The reported `stdev` is the absolute 1-sigma standard deviation in rem/h,
  read from the MAVRIC tally summary.

## Files

- `concrete_rebar_mixer.py`: computes the homogenized atom densities,
  weight fractions, and mass density of the concrete-rebar layer from the bar
  diameter and the bar spacings in x and y.
- `concrete_lid_contact_dose.py`: prototype `ConcreteLidContactDose` class
  that chains the whole workflow. It builds a four-case ORIGEN deck
  (irradiation and decay for each of the two modeled regions), reads back the
  decayed inventories from the F71 files, and builds a MAVRIC photon deck with
  point-detector contact dose at the lid top center. The ORIGEN deck runs
  under SCALE 6.3.3, checked with a stand-in F33. The MAVRIC deck has not been
  run yet, so check its first output before relying on the dose.
- `run_contact_dose.py`: driver for daily use. Edit the parameter block, then
  run single doses or cooling-time scans; see "Start here" above.

## Usage

```bash
PYTHONPATH=. python concrete_irrad/concrete_rebar_mixer.py
```

```python
from concrete_irrad.concrete_rebar_mixer import ConcreteRebarMixer

mixer = ConcreteRebarMixer(rebar_diameter_cm=1.59, spacing_x_cm=30.48, spacing_y_cm=15.24)
atom_dens = mixer.get_atom_densities()      # atoms/barn-cm, isotope level
rho_mix = mixer.mixture_density()           # g/cm3

low_co = ConcreteRebarMixer(1.59, 30.48, 15.24, steel_co59_wt_fraction=0.001)
```

## Assumptions

- Bars run along both x and y in a single plane, on rectangular grids.
- The mixed-layer thickness equals the bar diameter; concrete displaced by
  steel inside the slice belongs to the mixture volume only.
- Default compositions are ASTM A615 Grade 60 reinforcing steel and ordinary
  concrete per PNNL-15870 Rev. 1 material #96 at the ANSI/ANS-6.4.3 density of
  2.30 g/cm^3. Both are overridable constructor arguments.
- Impurity weight fractions are specified relative to their own constituent;
  the base composition is scaled down to conserve the unit mass.

## References

- Alhajali, S., Yousef, S., Naoum, B. (2016). Appropriate concrete for nuclear
  reactor shielding. *Applied Radiation and Isotopes* 107, 29--32.
  doi:10.1016/j.apradiso.2015.09.001.
- Gaudry, A., Delmas, M.-C. (2007). As cited in Alhajali et al. (2016).
- Kakinuma, S. et al. (2007). As cited in Alhajali et al. (2016).
- Kinno, M. et al. (2002). Raw materials for low-activation concrete shields:
  studies of Mn and rare earth elements. *Journal of Nuclear Science and
  Technology* 39, 1275--1280. As cited in Alhajali et al. (2016).
- Kerr, G.D. et al. (2005). Activation measurements for thermal neutrons,
  Part A: Cobalt-60 activation. In *DS02: Reassessment of the Atomic Bomb
  Radiation Dosimetry for Hiroshima and Nagasaki*, Vol. 1, RERF, Chapter 8.
- McConn, R.J. Jr., Gesh, C.J., Pagh, R.T., Rucker, R.A., Williams, R.G. III
  (2011). Compendium of Material Composition Data for Radiation Transport
  Modeling. PNNL-15870 Rev. 1, Pacific Northwest National Laboratory.
- NSTec (2007). Cost Effective Decommissioning of Shield Wall Structures.
  OSTI 908403.
- Radulescu, G., Banerjee, K. (2020). Best Practices for Shielding Analyses of
  Activated Metals and Spent Resins from Reactor Operation.
  ORNL/SPR-2020/1586. doi:10.2172/1669765.
- Suzuki, T. et al. (2001). As cited in Alhajali et al. (2016).
- Yoshida, G., Nishikawa, K., Nakamura, H., Yashima, H., Sekimoto, S.,
  Miura, T., Masumoto, K., Toyoda, A., Matsumura, H. (2020). Investigation of
  variations in cobalt and europium concentrations in concrete to prepare for
  accelerator decommissioning. *Journal of Radioanalytical and Nuclear
  Chemistry* 325, 801--806. doi:10.1007/s10967-020-07212-7.
