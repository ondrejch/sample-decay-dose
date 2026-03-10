# sample-decay-dose

This is a Python framework that calculates a handling dose of a sample from either
(A) a nuclide composition (SCALE F71 file), such as fuel salt sample; or
(B) an atom density and a transition matrix (SCALE F33 file), such as an irradiated coupon.

## Workflow:

1. SCALE/ORIGEN to (irradiate and) decay a sample, generate sources spectra and intensities. 
This is implemented in **Origen** class, which has a child class **OrigenFromTriton** for reading the F71 file
for the use case (A), and a child class **OrigenIrradiation** for the use case (B).
2. MAVRIC/Monaco to calculate ANSI-1991 neutron and gamma dose responses [rem/h] at 30 cm distance away from the sample.
This is implemented in the **DoseEstimator** class.
3. Since Monaco does not transport electrons, the ratio of beta/gamma spectral integrals is used for the beta dose estimate. 
This assumes an unshielded sample. If a sample container or other shielding is present, the user shall modify  
the MAVRIC input deck, **mavric_deck()** method of the **DoseEstimator** class.
4. An example of how to modify the **DoseEstimator** class is its child **DoseEstimatorSquareTank** class, where 
the sample is surrounded by layers of shielding. Thicknesses and atom densities of the shielding layers can be easily specified. 
Three layers of shielding are included as a default and an example.  
5. **DoseEstimatorStorageTank** is another derived class that models a storage tank with a gas plenum above teh sample.

## Repository structure   

Module scripts are in the **sample\_decay\_dose** directory
* **SampleDose.py** - The main module
* **read\_opus.py** - OPUS file reader for spectral integrals
* **isotopes.py** - Relative isotopic masses for conversion from atom to mass density 

Scripts showing how to use the framework are in the **examples** directory
* **calc\*py** - scripts showing example calculations 
* **plot\_doses.py** - results plotter
* **storage_tank_doses_scan_parallel.py**  - an example how to parallelize several MAVRIC calculations for the same decayed sample. 

Leaky-box ORIGEN utilities are in **leaky\_box\_origen**. Scratch scripts live in **play**.
Unit tests are in **test**.

See the per-directory READMEs for details:
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
4. Set SCALE binaries location, for example:
   - `export SCALE_BIN=/opt/scale6.3.2-mpi/bin`

## Runtime Requirements

- SCALE/ORIGEN/MAVRIC executables available under `SCALE_BIN`.
- For FGR-11 OCR extraction (`leaky_box_origen/extract_fgr11_dcf.py`):
  - source PDF default: `PDF/EPA 1988_FGR11_0.pdf`
  - `pdftoppm`
  - `tesseract`
  - generated DCF CSVs are written to `leaky_box_origen/data/`

## Quick Validation

- Run tests: `pytest -q`
- Typical generated outputs:
  - Leaky-box runs: `leaky_box_origen/run_YYYY-MM-DD/` (contains `_box_*`, `box*.json5`, `leaky_boxes*.xlsx/.json/.csv/.png`)
  - Dose post-processing: `leaky_box_origen/run_YYYY-MM-DD/leaky_boxes_dose*.csv/.png`

## SSH-native scalerte job API

Gateway executables and docs were moved to:
`sample_decay_dose/gw_exe/` (see `sample_decay_dose/gw_exe/README.md`).

Ondrej Chvala <ochvala@utexas.edu>
