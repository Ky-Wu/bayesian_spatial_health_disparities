library(data.table)
library(ggplot2)
library(sf)
library(spdep)
library(stringr)
library(maps)
library(maptools)
library(magrittr)
library(gridExtra)
library(rstan)
library(pROC)
rm(list = ls())
set.seed(1130)

# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# FITTING BYM2 MODEL with CAR SPATIAL COMPONENT VIA NIMBLE
source(file.path(getwd(), "src", "R", "simulation", "nimble_codes.R"))
# Hyperparameters
a0_sigma <- 2
b0_sigma <- .01
N <- 3000
W <- matrix(0, N, N)
n_nb_samp <- 2
ij_list <- lapply(seq_len(N - 1), function(n) {
  data.frame(i = rep(n, n_nb_samp),
             j = sample(seq(n + 1, N), n_nb_samp, replace = T))
}) %>%
  do.call(rbind, .)
ij_list <- unique(ij_list)
alpha <- 0.99
for (r in seq_len(nrow(ij_list))) {
  target_i <- ij_list[r,]$i
  target_j <- ij_list[r,]$j
  W[target_i, target_j] <- 1
  W[target_j, target_i] <- 1
}
D <- diag(rowSums(W))
alpha <- 0.99
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
scaling_factor <- exp(mean(log(diag(Sigma))))
Q_scaled <- Q * scaling_factor
Q_scaled_cholR <- Q_cholR * sqrt(scaling_factor)
Sigma_scaled <- Sigma / scaling_factor
# p = proportion of overall variance attributed to spatial variance
# sigma2 = overall variance
sigma2 <- 4
rho <- 0.93
sigma2Sp <- rho * sigma2
sigma2NSp <- (1 - rho) * sigma2
# Generate random effects
z <- rnorm(N, 0, 1)
phi <- solve(Q_scaled_cholR, z)
beta <- c(2, 5)
p <- length(beta)
X <- cbind(1, matrix(rnorm(N * (p - 1), 0, 1), ncol = p - 1))
mu <- X %*% beta
# K replicates in each region
K <- 1
Y <- replicate(K, rnorm(N, mean = mu + sqrt(sigma2Sp) * phi, sd = sqrt(sigma2NSp)))


# SAMPLING
a0_tau2 <- 2
b0_tau2 <- .01
nimble_constants <- list(X = X,
                         N = N,
                         p = p,
                         Q_scaled = Q_scaled,
                         mu_phi = rep(0, N),
                         a0_tau2 = a0_tau2,
                         b0_tau2 = b0_tau2)
nimble_data <- list(Y = Y)
# sampling parameters
n_samp <- 10 ## Number of posterior samples to be used for inference
n_chains <- 5 ## Number of different chains
n_burn <- 10 ## Number of initial runs ("burn" period) for MCMC to converge.

# Fit Nimble Model
base_init <- list(tau2 = 0.5,
                  rho = 0.5)
nimble_inits <- list(
  c(base_init, list(beta = rnorm(p, 1)), list(phi = rnorm(N, 1))),
  c(base_init, list(beta = rnorm(p, 1)), list(phi = rnorm(N, 1))),
  c(base_init, list(beta = rnorm(p, 1)), list(phi = rnorm(N, 1)))
)

model_parameters = c("beta",
                     "sigma2Sp",
                     "sigma2NSp",
                     "Yfit",
                     "sigma",
                     "rho",
                     "phi")

# sampling parameters
n_samp <- 150 ## Number of posterior samples to be used for inference
n_chains <- 3 ## Number of different chains
n_burn <- 50 ## Number of initial runs ("burn" period) for MCMC to converge.

rModel <- nimbleModel(code = BYM2_code,
                      constants = nimble_constants,
                      data = nimble_data)
mcmc.out <- nimbleMCMC(model = rModel,
                       nchains = n_chains,
                       inits = nimble_inits,
                       niter = n_burn + n_samp,
                       nburnin = n_burn,
                       monitors = model_parameters)
