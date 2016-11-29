# Reducing uncertainty about flavivirus infections using serological assays and novel analytics


This repository contains code used to generate the figures in:

Kucharski AJ, Riley S. Reducing uncertainty about flavivirus infections. Lancet Infect Dis. 2016

### Guide to files

`simulation_sera_run.R` Main model code - calls following source file:

`simulation_sera_model.R` Functions to generate simulated data sets, infer parameters and plot outputs

### Simulation model

We considered a single population that had experienced recent introduction of both dengue and zika virus. We simulated annual dengue attack rates for four DENV serotypes, using a lognormal distribution with mean 5% and standard deviation 1.5. In 2016, we assumed there was a ZIKV epidemic with 50% attack rate. For each age group, we calculated the probability that an individual would have previously been infected with each serotype, based on the simulated attack rates. 

We assumed the serological assay had imperfect sensitivity and specificity in identifying individuals that were seropositive. Assay sensitivity (defined as the probability of avoiding false-negative, given previous infection with a homotypic virus) was assumed to be 95% by default. Specificity (ie. probability of avoiding false-positive, given heterotypic infection) was assumed to be 70%.

### Simulation study to compare inference approaches

There is potential for cross-reaction between viruses, which could inflate the level of seropositivity. To investigate this bias, we simulated 100 randomly generated DENV and ZIKV time series and converted the annual attack rates into age-stratified seroprevalence data using a binomial distribution with 50 samples for each age group. 

We then used four different models to infer the annual force of infection through Markov chain Monte Carlo. We then used the estimates of annual force of infection to estimate the probability an individual aged 20 had previously been infected with ZIKV. The four models were as follows:

* **Model 1:** estimate the ZIKV annual force of infection using only ZIKV seroprevalence data
* **Model 2:** estimate force of infection for all five serotypes, plus sensitivity and specificity of assay. 
* **Model 3:** as in model 2, but only allow non-zero ZIKV force of infection in 2016
* **Model 4:** as in model 2, but only allow non-zero force of infection for DENV serotype or ZIKV in years where attack rate for that virus was >1%
