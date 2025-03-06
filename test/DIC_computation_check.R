library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
rm(list = ls())
set.seed(122)

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Exact sampling and Gibbs sampling
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_flatbeta_MCMC.cpp"))
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "NGRegression.cpp"))

# Set priors
a_0 <- 2.1
b_0 <- 0.1
M_0inv <- diag(rep(1e-4, p))
M_0inv <- diag(rep(1e-12), p)
m_0 <- rep(0, p)
n_sim <- 1000
target_rho <- .99999
BYM2Sampler <- new(BYM2ExactSampler, X, y, Q_scaled, target_rho)
BYM2Sampler$SetPriors(M_0inv, m_0, a_0, b_0)
exact_samps <- BYM2Sampler$ExactSample(n_sim)
m_n <- BYM2Sampler$m_nstar
M_n_invchol <- BYM2Sampler$M_nstar_inv
v <- solve(t(M_n_invchol), m_n)
M1m1 <- solve(M_n_invchol, v)
t(v) %*% v
b_0 + (t(y) %*% y / (1 - target_rho) - t(v) %*% v) / 2

beta_ols <- chol2inv(chol(t(X) %*% X)) %*% t(X) %*% y
Q_eigen <- eigen(Q_scaled)
P <- Q_eigen$vectors
Lambda <- Q_eigen$values
Q_neghalf <- P %*% diag(Lambda^(-0.5)) %*% t(P)
I_H <- diag(N) - X %*% chol2inv(chol(t(X) %*% X)) %*% t(X)
temp_eigen <- eigen(Q_neghalf %*% I_H %*% Q_neghalf)
O <- temp_eigen$vectors
D <- temp_eigen$values
U <- Q_neghalf %*% O
e <- I_H %*% y
v <- t(U) %*% e
target_rho <- .99999999999
SSE <- t(e) %*% e
B_inv <- U %*% diag(1 / (D + (1 - target_rho) / target_rho)) %*% t(U)
# another equation for b_n
b_0 + (t(e) %*% e - sum(v^2 / (D + (1 - target_rho) / target_rho))) / (2 * (1 - target_rho))
# ||y - WM_1m_1||^2
W <- cbind(X, diag(N))
z <- y - W %*% M1m1
quantity <- t(z) %*% z
z - I_H %*% (y - B_inv %*% e)
