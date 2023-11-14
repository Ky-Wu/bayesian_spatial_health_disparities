library(data.table)
library(ggplot2)
library(sf)
library(spdep)
library(stringr)
library(maps)
library(maptools)
library(magrittr)
library(nimble)
rm(list = ls())
set.seed(1130)

source(file.path(getwd(), "R", "simulation", "helper_functions.R"))

# Import CA counties
county_poly <- maps::map("county","california", fill=TRUE, plot=FALSE)
county_state <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[1]]))
county_names <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[2]]))
sf_use_s2(FALSE)
county_sp <- maptools::map2SpatialPolygons(county_poly, IDs = county_names)
county_nbs <- poly2nb(county_sp)
W <- nb2mat(county_nbs, style="B")
D <- diag(rowSums(W))
alpha <- 0.8
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
N <- nrow(W)
# noise to spatial variance ratio
delta2 <- 1 / 0.3
# spatial variance parameter
sigma2 <- 2
# Generate random effects
z <- rnorm(N, 0, 1)
phi <- solve(Q_cholR, z) * sqrt(sigma2)

beta <- c(2, 0.2, 0.3, 0.4, 0.5)
p <- length(beta)
X <- cbind(1, matrix(runif(N * (p - 1), -1, 1), ncol = p - 1))
mu <- X %*% beta
Y <- rnorm(N, mean = mu + phi, sd = sqrt(sigma2 * delta2))

plot(cbind(Y, st_as_sf(county_sp)))

# FITTING CAR MODEL VIA RJAGS


# Hyperparameters
C <- CAR_calcC(unlist(county_nbs), diag(D))
M <- 1 / diag(D)
L <- length(C)
tau2_beta <- .00001
a0_tau2 <- .01
b0_tau2 <- .01
a0_delta2 <- 1
b0_delta2 <- 1
nimble_constants <- list(X = X,
                         N = N,
                         p = p,
                         L = L,
                         C = C,
                         adj = unlist(county_nbs),
                         num = diag(D),
                         M = M,
                         mu_phi = rep(0, N),
                         mu_beta = rep(0, p),
                         tau2_beta = tau2_beta,
                         a0_tau2 = a0_tau2,
                         b0_tau2 = b0_tau2,
                         a0_delta2 = a0_delta2,
                         b0_delta2 = b0_delta2)
nimble_data <- list(Y = Y)


# Fit Nimble Model
base_init <- list(tau2 = 1.0, delta2 = 0.5, phi = rep(0, times = N), yFit = rep(0, times = N),
                  alpha = 0.8)
nimble_inits <- list(
  c(base_init, list(beta = rep(0, times = p))),
  c(base_init, list(beta = rep(-5, times = p))),
  c(base_init, list(beta = rep(5, times = p)))
)

# sampling parameters
n_samp <- 10000 ## Number of posterior samples to be used for inference
n_chains <- 3 ## Number of different chains
n_burn <- 1000 ## Number of initial runs ("burn" period) for MCMC to converge.

model_parameters = c("beta","sigma2", "sigma2Sp", "alpha", "delta2", "phi", "yFit")

spCode <- nimbleCode({
  for (i in 1:N) {
    Y[i] ~ dnorm(mu[i], tau2Sp * delta2)
    mu[i] <- inprod(X[i, 1:p], beta[1:p]) + phi[i]
    yFit[i] ~ dnorm(mu[i], tau2Sp * delta2) ## Posterior predictive model fit
  }
  phi[1:N] ~ dcar_proper(mu_phi[1:N], C[1:L], adj[1:L], num[1:N], M[1:N], tau2Sp, gamma = alpha)
  for (i in 1:p) {
    beta[i] ~ dnorm(mu_beta[i], tau2_beta)
  }
  sigma2 <- 1 / tau2
  sigma2Sp <- 1 / tau2Sp
  tau2Sp <- tau2 / delta2
  tau2 ~ dgamma(a0_tau2, b0_tau2)
  delta2 ~ dbeta(a0_delta2, b0_delta2)
  alpha ~ dunif(0.5, 1)
})
rModel <- nimbleModel(code = spCode, constants = nimble_constants, data = nimble_data)
mcmc.out <- nimbleMCMC(model = rModel,
                       nchains = n_chains,
                       inits = nimble_inits,
                       niter = n_burn + n_samp,
                       nburnin = n_burn,
                       monitor = model_parameters)
hist(mcmc.out$chain1[, "beta[1]"])
summary(mcmc.out$chain1[, "sigma2"])
hist(mcmc.out$chain1[, "sigma2Sp"])
summary(mcmc.out$chain1[, "sigma2Sp"])

samps <- rbind(mcmc.out$chain1, mcmc.out$chain2, mcmc.out$chain3)
dim(samps)

credible_intervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))


# COMPUTE SIMULATED STD DIFFERENCES OF NEIGHBORS
ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL

sim_phi <- samps[, grepl("phi", colnames(samps))]
sim_alpha <- samps[, "alpha"]
sim_delta2 <- samps[, "delta2"]
M0_inv <- diag(tau2_beta, ncol = p, nrow = p)
phi_std_diffs <- CARStdDiff(sim_phi, sim_alpha, sim_delta2, D, W, X, M0_inv,
                            ij_list, mc.cores = 8)
