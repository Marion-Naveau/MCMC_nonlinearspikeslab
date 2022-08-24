rm(list = ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(nimble)
library(microbenchmark)

source("Simul_donnees_Marion.R")
source("nimble_fit_functions.R")

#nimbleMCMC_object = OneVariableParamSpikeAndSlabDelay_marginal_delta_setup(data_input)

temps=microbenchmark(FitOneVariableParamSpikeAndSlabDelay_marginal_delta_slow(data = data_input, nu_0 = 0.04, niter = 3000),times=5,unit="s")

save(temps,file="Temps_MCMC_004.Rdata")

load("Temps_MCMC.Rdata") #nu0=0.01
temps

load("Temps_MCMC_004.Rdata") #nu0=0.04
temps
