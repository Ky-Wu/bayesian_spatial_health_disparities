library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
rm(list = ls())

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "CA_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Exact sampling
sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
set.seed(122)
model_rho <- .999999
BYM2Sampler <- new(BYM2ExactSampler, X, Y, Q_scaled, model_rho)
# Set priors
a_sigma <- 0.1
b_sigma <- 0.1
M_0inv <- diag(rep(1e-10, p))
m_0 <- rep(0, p)
BYM2Sampler$SetPriors(M_0inv, m_0, a_sigma, b_sigma)

n_sim <- 100000
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
D_star <- c(rep(1, N - p), rep(0, p))

B <- Q_scaled + diag(N) - X %*% XtX_invXt
B_eigen <- eigen(B)
summary(B_eigen$values)

# U^{-1}(I - H)y
Eb_limit <- U_inv %*% I_H %*% Y
# D^{\star}U^{-1}(I - H)y
Eb_limit <- D_star * Eb_limit
# (y - UD^{\star}U^{-1}e)
Eb_limit <- Y - U %*% Eb_limit
# (X'X)^{-1}X'(y - UD^{\star}U^{-1}e)
Eb_limit <- XtX_invXt %*% Eb_limit
UDsUt <- U %*% diag(D_star) %*% t(U)
Varb <- XtX_invXt %*% UDsUt %*% t(XtX_invXt)

beta_wls <- chol2inv(chol(t(X) %*% Q_scaled %*% X))
beta_wls <- beta_wls %*% t(X) %*% Q_scaled %*% Y
beta_ols <- XtX_invXt %*% Y

beta_hat <- apply(beta_sim, MARGIN = 2, mean)
beta_hat - beta_ols
Eb_limit
beta_hat
apply(beta_sim, MARGIN = 2, var)

summary(lm(Y ~ X - 1))

n <- length(Y)
A <- diag(N) - O %*% diag(D_star) %*% t(O)
sum(diag(t(A) %*% A))

max_change <- -XtX_invXt %*% U %*% diag(D_star) %*% U_inv %*% I_H %*% Y
max_change

B <- diag(N) - U %*% diag(D_star) %*% U_inv
sum(diag(t(B) %*% B))

e <- I_H %*% Y
z <- U %*% diag(D_star) %*% U_inv %*% e
z2 <- U %*% diag(1 - D_star) %*% U_inv %*% e
summary(lm(z ~ X - 1))
summary(lm(z2 ~ X - 1))
summary(lm(e ~ X - 1))

n <- nrow(X)
p <- ncol(X)
XtX <- t(X) %*% X
XtX_inv <- chol2inv(chol(XtX))
loss <- function(a) {
  a <- matrix(a, nrow = p, ncol = p)
  M0_inv <- t(a) %*% a
  Ainv <- chol2inv(chol(M0_inv + XtX))
  B <- Q_scaled + rho / (1 - rho) * (diag(n) - X %*% Ainv %*% t(X))
  B_inv <- chol2inv(chol(B))
  XtBinv <- t(X) %*% B_inv
  LHS <- Ainv %*% XtBinv %*% Q_scaled
  RHS <- XtX_inv %*% t(X)
  C <- LHS - RHS
  sqrt(sum(C^2))
}

init <- rep(1, p^2)
loss_optim <- optim(init, loss, control = list(abstol = .Machine$double.eps))
loss_optim
A <- loss_optim$par


