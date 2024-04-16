library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
rm(list = ls())

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Exact sampling
sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
set.seed(122)
model_rho <- .99999999999
BYM2Sampler <- new(BYM2ExactSampler, X, y, Q_scaled, model_rho)
# Set priors
a_sigma <- 0.1
b_sigma <- 0.1
M_0inv <- diag(rep(1e-10, p))
m_0 <- rep(0, p)
BYM2Sampler$SetPriors(M_0inv, m_0, a_sigma, b_sigma)

n_sim <- 1000
system.time({
  exact_samps <- BYM2Sampler$ExactSample(n_sim)
})

beta_sim <- exact_samps$beta
gamma_sim <- exact_samps$gamma
sigma2_sim <- exact_samps$sigma2
YFit_sim <- exact_samps$YFit
denom <- sqrt(sigma2_sim * model_rho)
phi_sim <- apply(gamma_sim, MARGIN = 2, function(x) x / denom)
rm(denom)

summary(beta_sim)

N <- nrow(X)
XtX <- t(X) %*% X
XtX_inv <- chol2inv(chol(XtX))
XtX_invXt <- XtX_inv %*% t(X)
Q_eigen <- eigen(Q_scaled)
if (any(Q_eigen$values <= 0)) stop("Q not positive definite")
Q_neghalf <- Q_eigen$vectors %*% diag(Q_eigen$values^(-0.5)) %*% t(Q_eigen$vectors)
H <- chol(XtX)
H <- solve(t(H), t(X))
H <- t(H) %*% H
I_H <- diag(N) - H
B <- Q_neghalf %*% I_H %*% Q_neghalf
B_eigen <- eigen(B)
D <- B_eigen$values
O <- B_eigen$vectors
U <- Q_neghalf %*% O
U_inv <- solve(U)
D_star <- ifelse(D > 0, 1, 0)
Eb <- U_inv %*% I_H %*% y
Eb <- D_star * Eb
Eb <- y - U %*% Eb
Eb <- XtX_invXt %*% Eb
UDsUt <- U %*% diag(D_star) %*% t(U)
Varb <- XtX_invXt %*% UDsUt %*% t(XtX_invXt)

beta_wls <- chol2inv(chol(t(X) %*% Q_scaled %*% X))
beta_wls <- beta_wls %*% t(X) %*% Q_scaled %*% y
