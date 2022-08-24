rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)
library(microbenchmark)
library(tidyverse)
library(nonlinearspikeslab)
library(tidyverse)

source("nimble_fit_functions.R")

S=5                         #nb de simulations
n <- 200 # nb d'individus
J <- 10 # nb de repetitions
p <- 2500 # nb de covariables
sigma2 <- 30 # variance residuelle des y
Gamma2 <- 200 # variance inter-individuelle des phi_i
mu <- 1200 # intercept
beta <- c(100, 50, 20, rep(0, p - 3)) # vecteur d'effets fixes
betatild <- c(mu, beta)
c0 <- 200 # parametre de la fct g
c2 <- 300 # parametre de la fct g, doit etre different de 0
psi <- c(c0, c2)

g <- function(phi_i, psi, t_ij) {
  return(psi[1] / (1 + exp(-(t_ij - phi_i) / psi[2])))
}

ncore = 10
cl = makeCluster(ncore)
registerDoParallel(cl)

nb_iter_finish=0
save(nb_iter_finish,file="Nb_iter_Test_cvMCMC_p2500.Rdata")

Res_Tot<-foreach(s = 1:S, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","microbenchmark","tidyverse","nimble")) %dopar% {

  print(s)
  # Simulation du modèle
  set.seed(s)
  phi <- rep(0, n)
  Vtild <- matrix(NA, nrow = n, ncol = p + 1)
  Vtild[, 1] <- rep(1, n)
  for (i in 1:n) {
    # simulation de V_i
    Vtild[i, 2:(p + 1)] <- rmvn(1, mu = rep(0, p), sigma = diag(1, p))
  }
  Vtild[, c(2:(p + 1))] <- scale(Vtild[, c(2:(p + 1))])
  # simulation de phi_i
  for (i in 1:n) {
    phi[i] <- sum(betatild * Vtild[i, ]) + rnorm(1, mean = 0, sd = sqrt(Gamma2))
  }
  V <- Vtild[, -1]
  # Simulation des données observées y

  Y <- rep(0, n * J) # vecteur des y_ij avec y_ij=Y[(i-1)*J+j]
  t <- rep(seq(150, 3000, length = J), n) # vecteur des t_ij avec t_ij=t[(i-1)*J+j]
  for (i in 1:n) {
    for (j in 1:J) {
      Y[(i - 1) * J + j] <- g(phi_i = phi[i], psi = psi, t_ij = t[(i - 1) * J + j]) + rnorm(1, mean = 0, sd = sqrt(sigma2))
    }
  }

  Id <- rep(c(1:n), each = J)
  data <- data.frame(Id, Y, t)
  id <- as.matrix(Id)
  data2 <- data[id <= 30, ]
  ggplot(data2, aes(x = t, y = Y, group = Id)) +
    geom_point() +
    geom_line() +
    theme_bw()

  data_input = data %>%
    as_tibble() %>%
    left_join(V %>%
                as_tibble() %>%
                rowid_to_column("Id")) %>%
    rename(times = t, i = Id)

  nu0=0.04

  mcmc_out=FitOneVariableParamSpikeAndSlabDelay_marginal_delta_slow(data = data_input, nu_0 = nu0, niter = 3000)

  load("Nb_iter_Test_cvMCMC_p2500.Rdata")
  nb_iter_finish=nb_iter_finish+1
  save(nb_iter_finish,file="Nb_iter_Test_cvMCMC_p2500.Rdata")

  mcmc_out
}
stopCluster(cl)

save(Res_Tot,file="Res_Test_cvMCMC_p2500_2.Rdata")

#load("Res_Test_cvMCMC_p2500.Rdata") #niter=10000

load("Res_Test_cvMCMC_p2500_2.Rdata") #niter=3000
mcmc_out=Res_Tot[[5]]
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

