# sample-decay-dose

This is a Python framework that calculates a handling dose of a sample from SCALE F71 file.

Workflow:
1. Extract nuclide composition from F71 file
2. Decay in ORIGEN and get gamma and beta spectra
3. Use MAVRIC/Monaco to calculate dose response for neutrons and gamma at 30 cm distance away from the sample
4. Since Monaco does not transport electrons, use the ratio of beta/gamma spectral integrals for beta dose estimate. 

All scripts are in sample\_decay\_dose directory
* DoseD71.py - main script
* read\_opus - OPUS file reader for spectral integrals 
* calc\*py - scripts showing example calculations 
* plot\_doses.py - results plotter

Ondrej Chvala <ochvala@utexas.edu>
