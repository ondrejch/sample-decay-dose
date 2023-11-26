# sample-decay-dose

This is a Python framework that calculates a handling dose of a sample from either
(A) a nuclide composition (SCALE F71 file), such as fuel salt sample; or
(B) an atom density and a transition matrix (SCALE F33 file), such as an irradiated coupon.

## Workflow:

1. SCALE/ORIGEN to (irradiate and) decay a sample, generate sources spectra and intensities. 
This is implemented in Origen class, which has one child class for reading the F71 file OrigenFromTriton
for use case (A), and OrigenIrradiation for user case (B).
2. Use MAVRIC/Monaco to calculate dose response for neutrons and gamma at 30 cm distance away from the sample.
3. Since Monaco does not transport electrons, use the ratio of beta/gamma spectral integrals for beta dose estimate. 
This assumes unshielded sample. If a sample container or other shielding is present, the user should modify  
the MAVRIC input deck in the DoseEstimator class.

## Repository structure   

Module scripts are in sample\_decay\_dose directory
* SampleDose.py - The main module
* read\_opus - OPUS file reader for spectral integrals
* isotopes.py - Relative isotopic masses for conversion from atom to mass density 

Scripts showing how to use the framework are in examples directory
* calc\*py - scripts showing example calculations 
* plot\_doses.py - results plotter

## Installation

1. Clone the repository
2. $ cd sample-decay-dose && pip install ./

Ondrej Chvala <ochvala@utexas.edu>
