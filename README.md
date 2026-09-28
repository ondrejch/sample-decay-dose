# sample-decay-dose

This is a Python framework that calculates a handling dose of a sample from either
(A) a nuclide composition (SCALE F71 file), such as fuel salt sample; or
(B) an atom density and a transition matrix (SCALE F33 file), such as an irradiated coupon.

## Workflow:

1. SCALE/ORIGEN to (irradiate and) decay a sample, generate sources spectra and intensities. 
This is implemented in **Origen** class, which has a child class **OrigenFromTriton** for reading the F71 file
for the use case (A), and a child class **OrigenIrradiation** for the use case (B).
2. MAVRIC/Monaco to calculate ANSI-1991 neutron and gamma dose responses [rem/h] at a point detector.
This is implemented in the **DoseEstimator** class, which places the detector at `det_x` = 30 cm from the sample centre
and raises an error if that point lies inside the sample.
3. Since Monaco does not transport electrons, the ratio of beta/gamma spectral integrals is used for the beta dose estimate. 
This estimate applies to a bare sample. The estimators with shielding layers report the beta response `'3'` as zero.
To model a sample container or other shielding, modify the MAVRIC input deck, **mavric_deck()** method of the
**DoseEstimator** class, or use one of the derived classes below.
4. An example of how to modify the **DoseEstimator** class is its child **DoseEstimatorSquareTank** class, where 
the sample is surrounded by layers of shielding. Thicknesses and atom densities of the shielding layers can be easily specified. 
Three layers of shielding are included as a default and an example.  
The tank estimators place the detector `det_standoff_distance` (1 cm by default) from the outer surface of the last
shielding layer. They recompute `det_x` from it in **mavric_deck()**, so set `det_standoff_distance` to move the detector.
5. **DoseEstimatorStorageTank** is another derived class that models a storage tank with a gas plenum above the sample.
**DoseEstimatorGenericTank** models a cylindrical tank filled with the sample.
6. **HandlingContactDoseEstimatorGenericTank** and **MHATank** compute a contact dose 0.1 cm and a handling dose 30 cm
from the outer surface of the last layer. Their `contact_dose` and `handling_dose` properties hold the two totals.
For a bare sample (no layers) the beta estimates are responses `'3'` (contact) and `'7'` (handling).

## Repository structure   

Module scripts are in the **sample\_decay\_dose** directory
* **SampleDose.py** - The main module
* **read\_opus.py** - OPUS file reader for spectral integrals
* **isotopes.py** - Relative isotopic masses for conversion from atom to mass density 

Scripts showing how to use the framework are in the **examples** directory
* **calc\*py** - scripts showing example calculations 
* **plot\_doses.py** - results plotter
* **storage_tank_doses_scan_parallel.py**  - an example how to parallelize several MAVRIC calculations for the same decayed sample. 
Run from `examples/`, it decays the sample once and runs a 5 x 5 steel/concrete shielding scan in parallel.

Leaky-box ORIGEN utilities are in **leaky\_box\_origen**. Scratch scripts live in **play**.
The reactor-cavity concrete lid activation and contact-dose workflow, with the concrete/rebar layer homogenizer,
is in **concrete\_irrad**.
Unit tests are in **test**.

See the per-directory READMEs for details:
* concrete_irrad/README.md
* examples/README.md
* leaky_box_origen/README.md
* play/README.md
* sample_decay_dose/README.md
* test/README.md

## Installation

1. Clone this repository
2. Python 3.10+ is required.
3. Install Python dependencies:
   - `pip install -r requirements.txt`
   - or `pip install ./`
   - Some examples need extra packages (`scipy`, `MSRRpy`); see examples/README.md.
4. Set SCALE binaries location, for example:
   - `export SCALE_BIN=/opt/scale6.3.2-mpi/bin`

## Runtime Requirements

- SCALE/ORIGEN/MAVRIC executables available under `SCALE_BIN`.
- The dose conversion factor tables ship in `leaky_box_origen/data/`. Regenerating them needs the source PDFs
  under `PDF/` (not tracked) and the extras `pip install '.[dcf-extract]'`:
  - `leaky_box_origen/extract_fgr11_dcf.py` reads FGR-11 (`PDF/6294233.pdf`, with `PDF/EPA 1988_FGR11_0.pdf`
    as a second scan). It needs `pdftoppm`, `pdftotext`, `pdfinfo` (poppler) and `tesseract`.
  - `leaky_box_origen/extract_icrp72_dcf.py` reads ICRP-72 (`PDF/ANIB_26_1.pdf`). It needs `pdftoppm` and
    `tesseract`.
  - See leaky_box_origen/README.md for what each table contains.

## Quick Validation

- Run tests: `pytest -q`
- Typical generated outputs:
  - Leaky-box runs: `leaky_box_origen/run_YYYY-MM-DD/` (contains `_box_*`, `box*.json5`, `leaky_boxes*.xlsx/.json/.csv/.png`)
  - Dose post-processing: `leaky_box_origen/run_YYYY-MM-DD/leaky_boxes_dose*.csv/.png`

## SSH-native scalerte job API

Gateway executables and docs were moved to:
`sample_decay_dose/gw_exe/` (see `sample_decay_dose/gw_exe/README.md`).

Ondrej Chvala <ochvala@utexas.edu>
