# This file tests the Gibbs sampling code to sample from the joint
# posterior f(beta, gamma, sigma2, rho | y).
# Random effects assume a CAR model.
# Rho (spatial percentage of total variation parameter) prior is uniform on [0, 1].
# Alpha (CAR spatial smoothing parameter) is also assumed fixed.

library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
rm(list = ls())
set.seed(1130)

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))

Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_flatbeta_MCMC.cpp"))

# Priors
a_sigma <- 1
b_sigma <- 1
lambda_rho <- 0.2
BYM2Sampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled, a_sigma, b_sigma, lambda_rho)
BYM2Sampler$initOLS()

BYM2Sampler$BurnMCMCSample(1000)
BYM2Sampler$sigma2
samps <- BYM2Sampler$MCMCSample(100)

n_burnin <- 2000
n_sim <- 3000
system.time({
  BYM2Sampler$BurnMCMCSample(n_sim)
  samps <- BYM2Sampler$MCMCSample(n_sim)
})

beta_sim <- samps$beta
gamma_sim <- samps$gamma
sigma2_sim <- samps$sigma2
rho_sim <- samps$rho
YFit_sim <- samps$YFit

hist(rho_sim)
hist(sigma2_sim)
