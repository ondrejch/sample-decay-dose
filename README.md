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

## Repository structure   

Module scripts are in the **sample\_decay\_dose** directory
* **SampleDose.py** - The main module
* **read\_opus.py** - OPUS file reader for spectral integrals
* **isotopes.py** - Relative isotopic masses for conversion from atom to mass density 

Scripts showing how to use the framework are in the **examples** directory
* **calc\*py** - scripts showing example calculations 
* **plot\_doses.py** - results plotter
* ** ** 

## Installation

1. Clone this repository
2. $ cd sample-decay-dose && pip install ./

Ondrej Chvala <ochvala@utexas.edu>
