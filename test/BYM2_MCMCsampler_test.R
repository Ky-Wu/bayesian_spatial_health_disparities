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
a_rho <- .1
b_rho <- .1
# limits on possible values of rho
lower_rho <- 0.00
upper_rho <- 1.0
BYM2Sampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled)
BYM2Sampler$setPriors(a_sigma, b_sigma, a_rho, b_rho, lower_rho, upper_rho)
BYM2Sampler$initRandom()

for(j in seq_len(100)) {
  print(paste0(j, "/100"))
  print(paste0("current rho: ", BYM2Sampler$rho))
  print(paste0("current sigma2: ", BYM2Sampler$sigma2))
  BYM2Sampler$burnMCMCSample(100)
}
# DRAW POSTERIOR SAMPLES
n_sim <- 10000
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
hist(beta_sim[,1])
hist(beta_sim[,2])

phi_sim <- apply(gamma_sim, MARGIN = 2, function(x) {
  x / sqrt(sigma2_sim * rho_sim)
})
hist(phi_sim[,130], xlim = c(-5, 5), breaks = 50)

phi_diffs <- BYM2_StdDiff(phi_sim, rho_sim, Q_scaled, X, ij_list)
# phi_diffs <- vapply(seq_len(k), function(pair_indx) {
#   i <- ij_list[pair_indx,]$i
#   j <- ij_list[pair_indx,]$j
#   sim_phi[,i] - sim_phi[,j]
# }, numeric(n_sim))
dev.off()
par(mfrow = c(2,2))
plot(density(abs(phi_diffs[,5])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,15])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,55])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,130])), xlim = c(0, 6))

phi_truediffs <- BYM2_StdDiff(matrix(phi, nrow = 1),
                              rho, Q_scaled, X, ij_list)

# Maximize entropy of posterior distribution with respect to epsilon
loss_function <- function(V, epsilon) -Entropy(V)
eps_optim <- optim(median(abs(phi_diffs)), function(e) {
  e_vij <- ComputeSimVij(phi_diffs, ij_list, epsilon = e)
  loss_function(e_vij, epsilon = e)
}, method = "Brent", lower = 0.0001, upper = 5.0)
optim_e <- eps_optim$par
eps_optim

true_diff <- abs(phi_truediffs) > optim_e
#true_diff <- (abs(true_phi_diffs) > optim_e)
mean(true_diff)
# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
optim_e_vij <- ComputeSimVij(phi_diffs, ij_list,
                             epsilon = optim_e)
optim_e_vij_order <- order(optim_e_vij, decreasing = F)
eta <- .15
t_seq_length <- 10000
t_seq <- seq(0, max(optim_e_vij) - .001, length.out = t_seq_length)
t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t, e = 0), numeric(1))
t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t, e = 0), numeric(1))
optim_t <- min(c(t_seq[t_FDR <= eta], 1))
optim_t

sum(optim_e_vij >= optim_t)

decisions <- logical(nrow(ij_list))
decisions[optim_e_vij >= optim_t] <- TRUE
print(paste0("Optimal epsilon: ", optim_e))
print(paste0("Optimal t: ", optim_t))
print(ComputeClassificationMetrics(decisions, true_diff))
mean(decisions)

rej_indx <- optim_e_vij_order[optim_e_vij_order %in%
                                which(optim_e_vij >= optim_t)]
detected_borders <- ij_list[rej_indx,]
node1_all <- county_sf[detected_borders$i,]
node2_all <- county_sf[detected_borders$j,]
intersections <- lapply(seq_len(sum(decisions)), function(i) {
  node1 <- node1_all[i,]
  node2 <- node2_all[i,]
  suppressMessages(st_intersection(st_buffer(node1, 0.0003),
                                   st_buffer(node2, 0.0003)))

}) %>%
  do.call(rbind, .)
intersections$rank <- rev(seq_len(sum(decisions)))

x1 <- cbind("y" = y, county_sf)
yfit_pmeans <- colMeans(YFit_sim)
yfit_pmeans_df <- cbind(y_pmeans = yfit_pmeans, county_sf)
p1 <- ggplot() +
  geom_sf(data = x1, aes(fill = y)) +
  scale_fill_viridis_c() +
  geom_sf(data = intersections, col = "red", alpha = 0) +
  theme_minimal() +
  theme(legend.position = "bottom")
p3 <- ggplot() +
  geom_sf(data = yfit_pmeans_df, aes(fill = y_pmeans)) +
  scale_fill_viridis_c() +
  geom_sf(data = intersections, col = "red", alpha = 0) +
  theme_minimal() +
  theme(legend.position = "bottom")
p3


## rejection order graph
e2 <- round(optim_e / 3, digits = 3)
e3 <- round(optim_e / 1.5, digits = 3)
e4 <- round(optim_e * 1.5, digits = 3)
e5 <- round(optim_e * 3, digits = 3)

e2_vij <- ComputeSimVij(phi_diffs, ij_list,
                        epsilon = e2)
e3_vij <- ComputeSimVij(phi_diffs, ij_list,
                        epsilon = e3)
e4_vij <- ComputeSimVij(phi_diffs, ij_list,
                        epsilon = e4)
e5_vij <- ComputeSimVij(phi_diffs, ij_list,
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
  true_diff = true_diff[optim_e_vij_order]
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
sim_vij_order_graph <- ggplot() +
  geom_point(data = rejection_path,
             aes(x = optim_e_vij, y = order, color = true_diff),
             alpha = 0.3, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of v_ij(", round(optim_e, digits = 3), ")"),
       y = "Rank of v_ij(eps)") +
  theme_minimal() +
  scale_color_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No Difference", "True Difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue"))
sim_vij_order_graph
