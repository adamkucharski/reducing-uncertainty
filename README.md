# reducing-uncertainty

### Summary

This repository contains code used to generate the figures in:

Kucharski AJ, Riley S. Reducing uncertainty about flavivirus infections. Lancet Infect Dis. 2016

Assay sensitivity (defined as the probability of avoiding false-negative, given previous infection with a homotypic virus) is assumed to be 95% by default. Specificity (ie. probability of avoiding false-positive, given heterotypic infection) is assumed to be 70%.


### Guide to files

`simulation_sera_run.R` Main model code - calls following source file:

> `simulation_sera_model.R` Generate simulated data sets and infer parameters
