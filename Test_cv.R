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

source("nimble_fit_functions.R")

###### Simulation des donnees ######

n <- 200 # nb d'individus
J <- 10 # nb de repetitions
p <- 500 # nb de covariables
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

# Simulation des parametres individuels phi_i
s <- 1
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


# Valeurs initiales des paramètres à estimer avec le MCMC-SAEM
betatild_init <- c(1400, rep(100, 10), rep(1, p - 10))
sigma2_init <- 100
Gamma2_init <- 5000

eta_init <- rep(400, 2)
Omega_init <- rep(20, 2) # Omega=(w_1^2,w_2^2)

# Hyperparamètres
M <- 20
Nu <- 10^(seq(-2, 2, length.out = M)) # grille de nu0 tq nu0/nu1 dans [10^(-6),10^(-3)]
nu1 <- 12000 # parametre slab
nu_Gamma <- 1 # parametre prior Gamma2
lb_Gamma <- 1 # parametre prior Gamma2
nu_sigma <- 1 # parametre prior sigma2
lb_sigma <- 1 # parametre prior sigma2
a <- 1 # parametre prior theta
b <- p # parametre prior theta
sigma2_mu <- 3000^2 # parametre prior mu
rho2 <- rep(1200, 2) # parametre prior eta
nu_omega <- 1 # parametre prior omega_m
lb_omega <- 1 # parametre prior omega_m


data_input = data %>%
  as_tibble() %>%
  left_join(V %>%
              as_tibble() %>%
              rowid_to_column("Id")) %>%
  rename(times = t, i = Id)


nu_0 = 0.01
set.seed(1)
mcmc_out=FitOneVariableParamSpikeAndSlabDelay_marginal_delta_slow(data = data_input, nu_0 = nu_0, niter = 4000)

save(mcmc_out,file="Test_cv_MCMC_p500_2.Rdata")

# load("Test_cv_MCMC_p500.Rdata")
load("Test_cv_MCMC_p500_2.Rdata") #avec bons hyperparam et init

#load("Test_cv_MCMC_p1000.Rdata")

#load("Test_cv_MCMC_p1500.Rdata")

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

