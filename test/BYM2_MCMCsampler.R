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
a_sigma <- 2.1
b_sigma <- 0.1
lambda_rho <- 1.5
BYM2Sampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled, a_sigma, b_sigma, lambda_rho)
BYM2Sampler$initOLS()
# cheating
BYM2Sampler$rho <- 0.9901583
BYM2Sampler$burnMCMCSample(1000)

BYM2Sampler$updateBeta()
BYM2Sampler$beta
BYM2Sampler$updateGamma()
BYM2Sampler$gamma
BYM2Sampler$vec_N
BYM2Sampler$burnMCMCSample(1000)
BYM2Sampler$sigma2
BYM2Sampler$rho
samps <- BYM2Sampler$MCMCSample(100)
BYM2Sampler$logPostd(TRUE)

rho_chain_analysis <- BYM2Sampler$rhoChainAnalysis(50, 2000, 100)
plot(rho_chain_analysis$log_postd, type = "l")
plot(rho_chain_analysis$rho, type = "l")
plot(rho_chain_analysis$rtr, type = "l")
# examine posterior density for multi-modality
e <- y - X %*% BYM2Sampler$beta
ete <- t(e) %*% e
r <- y - X %*% BYM2Sampler$beta - BYM2Sampler$gamma
rtr <- t(r) %*% r
log_d <- -N / 2 * log((1 - BYM2Sampler$rho) * 2 * pi * BYM2Sampler$sigma2) -
  rtr / (2 * BYM2Sampler$sigma2) / (1 - BYM2Sampler$rho)
# seems like it is very easy for sampler to get stuck in local minima...
log_d2 <- -N / 2 * log((1 - rho) * 2 * pi * BYM2Sampler$sigma2) -
  rtr / (2 * BYM2Sampler$sigma2) / (1 - rho)

#check true mean
rho <- 0.5
XtX <- t(X) %*% X
XtX_inv <- chol2inv(chol(XtX))
beta <- XtX_inv %*% t(X) %*% y
H <- X %*% XtX_inv %*% t(X)
B <- Q_scaled * (1 / rho - 1) + (diag(N) - H)
Binv <- chol2inv(chol(B))
U <- BYM2Sampler$U
Binv2 <- U %*% diag(1 / ((1/ rho - 1) + BYM2Sampler$D[,1])) %*% t(U)
Binv[1:5, 1:5]
Binv2[1:5, 1:5]
v <- XtX %*% beta + t(X) %*% Binv %*% y
C <- XtX + t(X) %*% Binv %*% X
Cinv <- chol2inv(chol(C))
v <- Cinv %*% v
u <- y - X %*% v
gamma_mean <- Binv %*% u

chol(C)
BYM2Sampler$matrix_pp
v <- XtX %*% beta + t(X) %*% Binv %*% y
v <- solve(chol(C), solve(t(chol(C)), v))

v <- XtX %*% beta + t(X) %*% Binv %*% y
v <- solve(BYM2Sampler$matrix_pp, solve(t(BYM2Sampler$matrix_pp), v))
BYM2Sampler$vec_p
v

v <- XtX %*% beta + t(X) %*% Binv %*% y
solve(t(BYM2Sampler$matrix_pp), v)
BYM2Sampler$vec2_p

# check gamma variance
temp <- chol(diag(N) / (1 - rho) + Q_scaled / rho)
chol2inv(temp)[1:5, 1:5]
temp2 <- 1 / sqrt(BYM2Sampler$Lambda / rho + 1 / (1 - rho))
temp2 <- BYM2Sampler$P %*% diag(temp2[,1])

Q_recovered <- BYM2Sampler$P %*% diag(BYM2Sampler$Lambda[,1]) %*% t(BYM2Sampler$P)
Q_recovered[1:5, 1:5]
Q_scaled[1:5, 1:5]
z <- rnorm(N)
error <- sqrt(sigma2) * solve(temp, z)
summary(error)
error2 <- sqrt(sigma2) * temp2 %*% z
summary(error2[,1])
summary(BYM2Sampler$vec_N[,1])

summary(gamma_mean + error)
summary(BYM2Sampler$gamma - BYM2Sampler$vec_N)



# DRAW POSTERIOR SAMPLES
n_sim <- 50000
system.time({
  samps <- BYM2Sampler$MCMCSample(n_sim)
})

beta_sim <- samps$beta
gamma_sim <- samps$gamma
sigma2_sim <- samps$sigma2
rho_sim <- samps$rho
YFit_sim <- samps$YFit

hist(rho_sim)
hist(sigma2_sim)

plot(rho_sim, type = "l")
plot(sigma2_sim, type = "l")
