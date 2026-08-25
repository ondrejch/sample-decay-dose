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
y-running bar segment of length `spacing_y`. The steel volume fraction is

```
f_steel = pi * d * (spacing_x + spacing_y) / (4 * spacing_x * spacing_y)
```

Each constituent keeps its own density and elemental composition. The mixture
mass density is `f_steel * rho_steel + (1 - f_steel) * rho_concrete`, and the
element number densities are the volume-weighted sums converted through natural
isotopic abundances.

## Activation-relevant impurities

Long-lived gamma emitters dominate the decay dose from irradiated concrete, and
they originate from trace elements rather than from the major constituents.
Kinno et al. (2002) showed that Co-60, Eu-152, Eu-154, and Cs-134 account for
99--100% of the total activity of ordinary concrete shields. The mixer
therefore carries explicit trace impurities in both materials.

The rebar steel carries a Co-59 impurity that produces the Co-60 source term.
The design band is 0.1--0.4 wt%, set through `steel_co59_wt_fraction` with a
default of 0.25 wt%. Measured carbon steels contain less cobalt: 93--151 ppm
across the carbon-steel samples surveyed by Radulescu and Banerjee (2020), and
about 200 micro-g/g in Hiroshima structural steels (Kerr et al., DS02). Values
in the design band are therefore conservative for ordinary rebar.

Ordinary concrete contains cobalt at the 0.16--21.9 ppm level (average 21.9
ppm) and europium at the 0.049--1.08 ppm level (average 1.08 ppm) in the
Japanese shielding-concrete survey compiled by Alhajali et al. (2016) from
Suzuki et al. (2001). Additional measurements give Co 4.5 ppm with Eu 0.098
ppm (Kakinuma et al., 2007) and Co 3.5 ppm with Eu 0.85 ppm for barite
concrete (Gaudry and Delmas, 2007). A United States decommissioning study used
about 1 ppm Co-59 and 0.01 ppm total Eu as natural levels (NSTec, 2007).
Yoshida et al. (2020) measured facility-to-facility variations at the same
order of magnitude by neutron activation analysis. The defaults are Co 10 ppm
and Eu 1.0 ppm, mid-range across these sources; both are overridable through
`concrete_impurities_wt`. Natural europium expands to Eu-151 and Eu-153, the
precursors of Eu-152 and Eu-154.

## Files

- `concrete_rebar_mixer.py`: computes the homogenized atom densities,
  weight fractions, and mass density of the concrete-rebar layer from the bar
  diameter and the bar spacings in x and y.
- `concrete_lid_contact_dose.py`: prototype `ConcreteLidContactDose` class
  that chains the whole workflow. It builds a four-case ORIGEN deck
  (irradiation and decay for each of the two modeled regions), reads back the
  decayed inventories from the F71 files, and builds a MAVRIC photon deck with
  point-detector contact dose at the lid top center. Preview both decks with
  `PYTHONPATH=. python concrete_irrad/concrete_lid_contact_dose.py`; execute
  the chain with `--run --lib <alpha_library.f33>`. The irradiation requires an
  alpha-library F33 carrying the cavity neutron spectrum. The deck structure
  follows the SampleDose conventions but has not yet been exercised against a
  real SCALE installation; verify on first run.

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
  concrete per ANSI/ANS-6.4.3; both are overridable constructor arguments.
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
