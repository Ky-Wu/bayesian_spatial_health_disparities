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
alpha <- 0.99
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
Q_scaled <-
N <- nrow(W)
# spatial variance to noise ratio
alpha <- 0.5
delta2 <- alpha / (1 - alpha)
# spatial variance parameter
sigma2 <- 4
sigma2Sp <- delta2 * sigma2
# Generate random effects
z <- rnorm(N, 0, 1)
phi <- solve(Q_cholR, z) * sqrt(sigma2Sp)

beta <- c(2, 0, .25, .5, .75)
p <- length(beta)
X <- cbind(1, matrix(runif(N * (p - 1), -1, 5), ncol = p - 1))
mu <- X %*% beta
Y <- rnorm(N, mean = mu + phi, sd = sqrt(sigma2))

plot(cbind(Y, st_as_sf(county_sp)))
plot(cbind(phi, st_as_sf(county_sp)))
# FITTING CAR MODEL VIA NIMBLE

# Hyperparameters
C <- CAR_calcC(unlist(county_nbs), diag(D))
M <- 1 / diag(D)
L <- length(C)
tau2_beta <- .0000001
a0_tau2Sp <- .1
b0_tau2Sp <- .1
a0_tau2 <- .001
b0_tau2 <- .001
a0_alpha <- 1
b0_alpha <- 1
nimble_constants <- list(X = X,
                         N = N,
                         p = p,
                         Q = Q,
                         mu_phi = rep(0, N),
                         a0_alpha = a0_alpha,
                         b0_alpha = b0_alpha,
                         a0_tau2 = a0_tau2,
                         b0_tau2 = b0_tau2)
nimble_data <- list(Y = Y)


# Fit Nimble Model
base_init <- list(tau2 = 0.5,
                  alpha = 0.5,
                  phi = rep(0, times = N))
nimble_inits <- list(
  c(base_init, list(beta = rep(0, times = p))),
  c(base_init, list(beta = rep(-5, times = p))),
  c(base_init, list(beta = rep(5, times = p)))
)

# sampling parameters
n_samp <- 10000 ## Number of posterior samples to be used for inference
n_chains <- 3 ## Number of different chains
n_burn <- 1000 ## Number of initial runs ("burn" period) for MCMC to converge.

model_parameters = c("beta","sigma", "sigmaSp", "phi", "delta2", "alpha")

rModel <- nimbleModel(code = BYM_Code,
                      constants = nimble_constants,
                      data = nimble_data)
mcmc.out <- nimbleMCMC(model = rModel,
                       nchains = n_chains,
                       inits = nimble_inits,
                       niter = n_burn + n_samp,
                       nburnin = n_burn,
                       monitors = model_parameters)

samps <- rbind(mcmc.out$chain1, mcmc.out$chain2, mcmc.out$chain3)
dim(samps)
hist(samps[, "sigma"])
summary(samps[, "sigma"])
hist(samps[, "sigmaSp"])
summary(samps[, "sigmaSp"])
hist(samps[, "delta2"])
summary(samps[, "delta2"])
hist(samps[, "alpha"])
summary(samps[, "alpha"])
credible_intervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))
phi_credible_intervals <- credible_intervals[grepl("phi", rownames(credible_intervals)),]


# COMPUTE SIMULATED STD DIFFERENCES OF NEIGHBORS
ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL

sim_phi <- samps[, grepl("phi", colnames(samps))]
sim_delta2 <- samps[, "delta2"]
M0_inv <- diag(0, ncol = p, nrow = p)
phi_diffs <- sapply(seq_len(nrow(ij_list)), function(r) {
  target_i <- ij_list[r, ]$i
  target_j <- ij_list[r, ]$j
  sim_phi[, target_i] - sim_phi[, target_j]
})
phi_std_diffs <- apply(phi_diffs, MARGIN = 2, function(x) x / sd(x))
# phi_std_diffs <- CARStdDiff(sim_phi, sim_delta2, D, W, X, M0_inv,
#                             ij_list, mc.cores = 4)
# true_phi_std_diffs <- ComputeSimSTDDifferences(matrix(phi, ncol = 1), Sigma,
#                                                ij_list = ij_list)
true_phi_std_diffs <- sapply(seq_len(nrow(ij_list)), function(r) {
  target_i <- ij_list[r, ]$i
  target_j <- ij_list[r, ]$j
  phi[target_i] - phi[target_j]
})
true_phi_std_diffs <- true_phi_std_diffs / sd(true_phi_std_diffs)

plot(density(abs(phi_std_diffs[,5])), xlim = c(0, 6))

## loss function as heuristic to choose epsilon threshold
loss_function <- function(V, epsilon) GeneralLoss(V, epsilon, I, log)
eps_optim <- optim(median(abs(phi_std_diffs)), function(e) {
  e_vij <- ComputeSimVij(phi_std_diffs, ij_list, epsilon = e)
  loss_function(e_vij, epsilon = e)
}, method = "Brent", lower = 0, upper = max(abs(phi_std_diffs)))
optim_e <- eps_optim$par
mean(abs(true_phi_std_diffs) > optim_e)


optim_e_vij <- ComputeSimVij(phi_std_diffs, ij_list,
                             epsilon = optim_e)
# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
eta <- .05
t_seq_length <- 1000
t_seq <- seq(0, max(optim_e_vij) - .001, length.out = t_seq_length)
t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t), numeric(1))
t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t), numeric(1))
optim_t <- min(c(t_seq[t_FDR <= eta], 1))


e2 <- round(optim_e / 3, digits = 3)
e3 <- round(optim_e / 1.5, digits = 3)
e4 <- round(optim_e * 1.5, digits = 3)
e5 <- round(optim_e * 3, digits = 3)

e2_vij <- ComputeSimVij(phi_std_diffs, ij_list,
                        epsilon = e2)
e3_vij <- ComputeSimVij(phi_std_diffs, ij_list,
                        epsilon = e3)
e4_vij <- ComputeSimVij(phi_std_diffs, ij_list,
                        epsilon = e4)
e5_vij <- ComputeSimVij(phi_std_diffs, ij_list,
                        epsilon = e5)

optim_e_vij_order <- order(optim_e_vij, decreasing = F)
e2_vij_order <- order(e2_vij[optim_e_vij_order], decreasing = F)
e3_vij_order <- order(e3_vij[optim_e_vij_order], decreasing = F)
e4_vij_order <- order(e4_vij[optim_e_vij_order], decreasing = F)
e5_vij_order <- order(e5_vij[optim_e_vij_order], decreasing = F)
rejection_path <- data.table(
  optim_e_vij = seq_along(optim_e_vij),
  e2_vij_order = e2_vij_order,
  e3_vij_order = e3_vij_order,
  e4_vij_order = e4_vij_order,
  e5_vij_order = e5_vij_order,
  true_diff = (abs(true_phi_std_diffs) > optim_e)[optim_e_vij_order]
)
rejection_path <- melt(rejection_path,
                      id.vars = c("optim_e_vij", "true_diff"),
                      variable.name = "order_type",
                      value.name = "order")
rejection_path[, order_type := fcase(
  order_type == "e2_vij_order", paste0("eps = ", e2),
  order_type == "e3_vij_order", paste0("eps = ", e3),
  order_type == "e4_vij_order", paste0("eps = ", e4),
  order_type == "e5_vij_order", paste0("eps = ", e5)
)]
sim_vij_pval_order_graph <- ggplot() +
  geom_point(data = rejection_path,
             aes(x = optim_e_vij, y = order, color = true_diff),
             alpha = 0.3, size = 1) +
  facet_grid(~order_type) +
  labs(x = "Rank of Optimal E V_ij", y = "Other Eps order") +
  theme_minimal() +
  scale_color_manual(name = "",
                     labels = c("No Difference", "True Difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue"))
sim_vij_pval_order_graph

true_diff <- (abs(true_phi_std_diffs) > optim_e)
decisions <- logical(L / 2)
decisions[optim_e_vij_order[seq_len(sum(true_diff))]] <- TRUE
ComputeClassificationMetrics(decisions, true_diff)

