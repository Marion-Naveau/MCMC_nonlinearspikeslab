rm(list = ls())

source("Simul_donnees_Marion.R")
source("nimble_fit_functions.R")
library(nimble)

nimbleMCMC_object = OneVariableParamSpikeAndSlabDelay_marginal_delta_setup(data_input)

## Time the execution for MCMC sampling (1 chain)
temps=TimeOneVariableParamSpikeAndSlabDelay_marginal_delta(Cmcmc = nimbleMCMC_object, nu_0 = 0.04, niter = 3000)

save(temps,file="temps3_p5000.Rdata")

## Actually get a sample from posterior.
#FitOneVariableParamSpikeAndSlabDelay_marginal_delta(Cmcmc = nimbleMCMC_object, nu_0 = 0.01, niter = 5000)

load("temps.Rdata") #p=5000, nu0=0.04
load("temps2_p5000.Rdata") #p=5000, nu0=0.04
load("temps3_p5000.Rdata") #p=5000, nu0=0.04

load("temps1.Rdata") #nu0=0.01
temps1=temps
load("temps2.Rdata") #nu0=0.01
temps2=temps
load("temps3.Rdata") #nu0=0.01
temps3=temps
load("temps4.Rdata") #nu0=0.01
temps4=temps
load("temps5.Rdata") #nu0=0.01
temps5=temps

load("temps6.Rdata") #nu0=0.04
temps6=temps
load("temps7.Rdata") #nu0=0.04
temps7=temps
load("temps8.Rdata") #nu0=0.04
temps8=temps
load("temps9.Rdata") #nu0=0.04
temps9=temps
load("temps10.Rdata") #nu0=0.04
temps10=temps

load("temps11.Rdata") #nu0=2
temps11=temps
load("temps12.Rdata") #nu0=2
temps12=temps
load("temps13.Rdata") #nu0=2
temps13=temps
load("temps14.Rdata") #nu0=2
temps14=temps
load("temps15.Rdata") #nu0=2
temps15=temps
