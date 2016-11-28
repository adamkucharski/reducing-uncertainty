# Code to accompany:
# Reducing uncertainty about flavivirus infections. Lancet Infect Dis. 2016
# (c) Kucharski AJ, Riley S

# - - - - - - - - - -
# Main code run file

library(foreach)
library(doMC)
library(mvtnorm)
library(MASS)
registerDoMC(4)  #change the 2 to your number of CPU cores
getDoParWorkers()

rm(list=ls(all=TRUE))

setwd("~/Documents/reducing-uncertainty/model_code/")
source("simulation_sera_model.R")

dir.create(file.path(getwd(), "posterior_runs"))
dir.create(file.path(getwd(), "plot_simulations"))

# - - - - - - - - - - - - - - - - - 
# Simulate serology for ZIKV/DENV

# inf.years.sera = range of years for outbreaks
# zikv.attack = attack rate in final year for strain 5
# p.inf.in = mean attack for each strain (drawn from lognormal with sd.val.in=1 by default)
# dmatrix.in = impose cross-reaction structure

# By default cross reaction structure assumes type 1 error (given infection) and type 2 (from cross-reaction)
# error = probability negative given infection (1-sensitivity)
# sigma = probability positive from cross-reaction (1-specificity)

theta.serology=c(error=0.05,sigma=0.3) # Model parameters
per_sample0=50 # Number of paticipants in each year group
seedRuns=100 # Number of repeated simulations

# Generate data
for(seedK in 1:seedRuns){

  simulate_sera_data(strains=5,inf.years.sera=c(1985:2016),time.series.in=NULL,theta=theta.serology,
                     p.inf.in=0.05*c(1,1,1,1,1),sd.val.in=1.5,seedi=seedK,roundv=F,dmatrix.in=NULL,zikv.attack=0.5,per_sample=per_sample0)

}


# Plot results
plot_simulated_sera_data(strains=5,seedi=1)


# - - - - - - - - - - - - - - - - - 
# Run MCMC inference for each simulated data set for 4 scenarios
#Model 1: estimated force of infection for ZIKV Zika virus only. 
#Model 2: estimated force of infection for all five serotypes, plus sensitivity and specificity of assay. 
#Model 3: as in model 2, but only including non-zero ZIKV force of infection in 2016
#Model 4: as in model 2, but only including non-zero force of infection for a serotype in years where attack rate for that serotype was >1%

for(scenario in c(1,4)){
  foreach(seedK=c(1:seedRuns)) %dopar% {
    inference_model(seedK,strains=5,runsMCMC=1e5,scenario,per_sample=per_sample0,switch0=10)
  }
}

# - - - - - - - - - - - - - - - - - 
# Plot output of MCMC

plot.posteriors(per_sample=per_sample0,strains=5,scenario=1,seedK=4)
plot.performance(per_sample=per_sample0,age_out=20,strains=5,scenarioN=4,runs=seedRuns)

