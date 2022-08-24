rm(list = ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(nimble)
library(microbenchmark)
library(nonlinearspikeslab)
library(tidyverse)

source("Simul_donnees_Marion2.R")
source("nimble_fit_functions.R")

PlotSimulatedData(data_input, maxspecies = 30)

nu_0 = 0.04
set.seed(1)
mcmc_out=FitOneVariableParamSpikeAndSlabDelay_marginal_delta_slow(data = data_input, nu_0 = nu_0, niter = 6000)

save(mcmc_out,file="Test_MCMC_p5000.Rdata")

#load("Test_MCMC_p2000.Rdata")
load("Test_cv_MCMC_p1000_2.Rdata") #avec bons hyperparam et Init Gamma2=500

load("Test_MCMC_p5000.Rdata")
traceplot(mcmc.out = mcmc_out, parnames = c("Gamma2", "sigma2")) + scale_y_log10()
PriorPosteriorPlotSmallVectorParams(mcmc_out = mcmc_out, parnames = c("sigma2", "Gamma2")) + scale_x_log10()

traceplot(mcmc.out = mcmc_out, parnames = c("psi"))
PriorPosteriorPlotSmallVectorParams(mcmc_out = mcmc_out, parnames = c("psi"))

traceplot(mcmc.out = mcmc_out, parnames = c("alpha"))
PriorPosteriorPlotSmallVectorParams(mcmc_out = mcmc_out, parnames = c("alpha"))

traceplot(mcmc.out = mcmc_out, parnames = c("mu"))
PriorPosteriorPlotSmallVectorParams(mcmc_out = mcmc_out, parnames = c("mu"))

traceplot(mcmc.out = mcmc_out, parnames = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]"))
#PriorPosteriorPlotSmallVectorParams(mcmc_out = mcmc_out, parnames = c("beta"))

traceplot(mcmc.out = mcmc_out,parnames = paste("phi[", 1:9, "]", sep = ""))
posterior_predictive_estimate(data = data_input, mcmc_out = mcmc_out, maxspecies = 20)
