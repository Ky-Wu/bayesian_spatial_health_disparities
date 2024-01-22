library(Rcpp)
library(RcppArmadillo)
set.seed(1130)
rm(list = ls())

Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_MCMCsampler.cpp"))

# Generate Data
N <- 4000
p <- 6
L <- matrix(rnorm(N^2), N, N)
Sigma <- L %*% t(L)
Q <- chol2inv(chol(Sigma))
scaling_factor <- exp(mean(log(diag(Sigma))))
Q_scaled <- Q * scaling_factor
Sigma_scaled <- Sigma / scaling_factor
Sigma_eigen <- eigen(Sigma_scaled)
P <- Sigma_eigen$vectors
Lambda <- 1 / Sigma_eigen$values
X <- matrix(rnorm(N * p), N, p)
truth_beta <- runif(p)
sigma2 <- 4
rho <- 0.78
sigma2Sp <- rho * sigma2
sigma2NSp <- (1 - rho) * sigma2
Q_scaled_cholR <- chol(Q_scaled)
# Generate random effects
z <- rnorm(N, 0, 1)
phi <- solve(Q_scaled_cholR, z)
X <- cbind(1, matrix(rnorm(N * (p - 1), 0, 1), ncol = p - 1))
mu <- X %*% truth_beta
y <- rnorm(N, mean = mu + sqrt(sigma2Sp) * phi, sd = sqrt(sigma2NSp))

BYM2Sampler <- new(BYM2ModelMCMC, X, y, Q_scaled)
BYM2Sampler$initOLS()
BYM2Sampler$beta
BYM2Sampler$phi
BYM2Sampler$updateBeta()
BYM2Sampler$beta
BYM2Sampler$updatePhi()
BYM2Sampler$phi
