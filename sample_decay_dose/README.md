# sample_decay_dose

Core library for sample decay and dose calculations.

**Key modules**
- `SampleDose.py` main workflow and high-level API.
- `HotCell.py` sample in a cubical hot cell; `Radiator.py` salt in a tube-bundle radiator.
- `read_opus.py` OPUS reader.
- `utils.py`, `constants.py`, `isotopes.py` shared utilities.
- `data/` packaged reference data.

## ORIGEN decay grids

Every ORIGEN deck decays the sample with `start=0`, so F71 position 1 holds the state at the start of the decay case.
A positive decay time writes `SAMPLE_F71_position` positions, and the last one holds the decayed state.
The first decay point is capped at 1/1000 of the decay time, which keeps the time grid strictly increasing for short decays.
A zero decay time writes one short step (`ZERO_DAY_STEP_DAYS`) and reads the zero-day state at position 1.
`Origen.decayed_F71_position` returns the position that holds the decayed state.
The OPUS spectra are selected by that position (`npos`).
For `OrigenFromTritonMHA`, ORIGEN saves position 1 after the `retained` processing, so it holds the retained elements only.

## Dose responses

Responses are point-detector doses in rem/h (ANSI-1991 flux-to-dose factors) with absolute 1-sigma uncertainties.
- `DoseEstimator` models a bare square cylinder of the sample. Its detector sits at `det_x` from the sample centre (30 cm by default).
  The distance from the sample surface is therefore `det_x` minus the sample radius. `mavric_deck()` raises `ValueError` if the detector lies inside the sample.
- The tank estimators place the detector `det_standoff_distance` from the outer surface of the last layer.
  The three layer lists (`layers_thicknesses`, `layers_mats`, `layers_temperature_K`) must have equal lengths, and every deck checks them.
  Empty lists give a bare sample.
- `HandlingContactDoseEstimatorGenericTank`, `MHATank`, and `HotCellDoses` have a contact detector (location 1, 0.1 cm) and a handling detector (location 2, 30 cm).
  Their responses are keyed by point-detector ID: `'1'`/`'2'` neutron/photon at contact and `'5'`/`'6'` at handling.
  `contact_dose` and `handling_dose` hold the two totals. `total_dose` returns `contact_dose`.

**Beta estimate.** Monaco does not transport electrons. The beta response is `beta_over_gamma` times the photon response, where the ratio comes from the ORIGEN beta and gamma spectral integrals.
This estimate applies to a bare sample, and `beta_applies` is true only then.
Behind any layer, and in the radiator, the beta response is reported as zero.
The contact/handling estimators add beta as `'3'` (contact) and `'7'` (handling).
If the ratio is unknown for a bare sample, `get_responses()` raises `RuntimeError`.

**Totals.** `total_dose`, `contact_dose`, and `handling_dose` sum the responses.
The neutron and photon tallies are independent, so their sigmas add in quadrature.
The beta response is proportional to the photon response, so the photon and beta sigmas add linearly before the quadrature sum.

**ORIGEN spectra not in memory.** A `DoseEstimator` built from an `Origen` object whose case was not run in the current process reads the spectral integrals from the OPUS files in the ORIGEN case directory.
If those files are missing, `beta_over_gamma` and `neutron_intensity` are `None`.
An unknown neutron intensity includes the neutron source. A known zero intensity omits it.

## Radiator

`RadiatorBox` reads the salt state directly from the burned-material F71 position, without an ORIGEN decay.
Monaco's `origensBinaryConcentrationFile` distribution uses the emission spectra stored on the F71.
The position must carry gamma spectra, and neutron spectra when the neutron source is included (flags `G` and `N` in the `DCGNAB` column of `obiwan view -format=info`).
ORIGEN writes these for cases run with `gamma=yes neutron=yes`. TRITON depletion positions carry concentrations only, and `mavric_deck()` raises for them.
Setting `neutron_intensity = 0.0` omits the neutron source explicitly.

Source strength assumption: the F71 spectra are totals over the F71 material volume (the `volume` row of `obiwan view`).
The radiator holds `all_pins_volume` of salt, so both sources are scaled by `all_pins_volume / V_F71` unless `source_multiplier` is set.

## Hot cell

In `HotCellDoses`, `layers_thicknesses[0]` is the half-width of the innermost cuboid (the cell interior, measured from the sample centre).
It must be at least the sample radius. Each later entry is the thickness of the next cuboid shell.
With `reuse_adjoint_flux = True`, MAVRIC skips the adjoint calculation and reads `adjoint_flux_file`.
The default file is `<case_dir>/my_dose.adjoint.dff`, which an earlier MAVRIC run of the same case directory writes.
The deck generator raises `FileNotFoundError` if the file does not exist.

## Utilities

- `integrate_opus()` integrates an OPUS histogram written as (E_low, y), (E_high, y) point pairs.
  Adjacent bins with equal y count as separate bins.
- `run_scale()` reports a failure for a non-zero exit code or for SCALE's failure markers in the messages or the output file.
  The markers are `***Error:` lines, `terminated due to errors`, and `<module> failed.`.
  scalerte can exit with code 0 after a module aborts, which is why the output is checked.
- `utils` reads `constants.ATOM_DENS_MINIMUM` at call time, so setting it on `sample_decay_dose.constants` takes effect.
- `get_rho_from_atom_density()` maps SCALE bound-scatterer IDs such as `h-poly` to the nuclide mass, and gives metastable states the ground-state mass.
