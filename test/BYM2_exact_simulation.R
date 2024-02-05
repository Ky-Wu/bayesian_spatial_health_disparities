# This file tests the exact sampling theory provided by the augmented linear
# system applied to the linear mixed effects model to sample from the joint
# posterior f(beta, gamma, sigma2 | y).
# Random effects assume a CAR model.
# Rho (spatial percentage of total variation parameter) is assumed fixed.
# Alpha (CAR spatial smoothing parameter) is also assumed fixed.

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
# Exact sampling
sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))

BYM2Sampler <- new(BYM2ExactSampler, X, y, Q_scaled, rho)
# Set priors
a_0 <- 2.1
b_0 <- 0.1
M_0inv <- diag(rep(1e-4, p))
m_0 <- rep(0, p)
BYM2Sampler$SetPriors(M_0inv, m_0, a_0, b_0)

n_sim <- 10000
system.time({
  samps <- BYM2Sampler$ExactSample(n_sim)
})

beta_sim <- samps$beta
gamma_sim <- samps$gamma
sigma2_sim <- samps$sigma2
YFit_sim <- samps$YFit
denom <- sqrt(sigma2_sim * rho)
phi_sim <- apply(gamma_sim, MARGIN = 2, function(x) x / denom)
rm(denom)
plot(density(phi_sim[, 105]))
hist(phi_sim[, 20])

phi_diffs <- BYM2_StdDiff(phi_sim, rep(rho, n_sim), Q_scaled, X, ij_list)
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

true_phi_diffs <- ComputeSimSTDDifferences(matrix(phi, ncol = 1),
                                           Sigma_scaled,
                                           ij_list = ij_list)
mean(abs(true_phi_diffs) > optim_e)
# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
optim_e_vij <- ComputeSimVij(phi_diffs, ij_list,
                             epsilon = optim_e)
eta <- .05
t_seq_length <- 10000
t_seq <- seq(0, max(optim_e_vij) - .001, length.out = t_seq_length)
t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t, e = 0), numeric(1))
t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t, e = 0), numeric(1))
optim_t <- min(c(t_seq[t_FDR <= eta], 1))
optim_t

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
true_diff <- (abs(true_phi_diffs) > optim_e)
# true_diff <- vapply(seq_len(nrow(ij_list)), function(r) {
#   phi[ij_list[r,]$i] != phi[ij_list[r,]$j]
# }, logical(1))
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

decisions <- logical(nrow(ij_list))
decisions[optim_e_vij >= optim_t] <- TRUE
print(paste0("Optimal epsilon: ", optim_e))
print(paste0("Optimal t: ", optim_t))
print(ComputeClassificationMetrics(decisions, true_diff))
mean(decisions)

x1 <- cbind("y" = y, st_as_sf(county_sp))
x2 <- cbind("phi" = phi, st_as_sf(county_sp))
yfit_pmeans <- colMeans(YFit_sim)
yfit_pmeans_df <- cbind(y_pmeans = yfit_pmeans, st_as_sf(county_sp))
p1 <- ggplot() +
  geom_sf(data = x1, aes(fill = y)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom")
p2 <- ggplot() +
  geom_sf(data = x2, aes(fill = phi)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom")
p3 <- ggplot() +
  geom_sf(data = yfit_pmeans_df, aes(fill = y_pmeans)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom")
p1
p2
p3

ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_data.png"),
       width = 6, height = 4.5, units = "in", p1)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_phi.png"),
       width = 6, height = 4.5, units = "in", p2)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_yfitpmean.png"),
       width = 6, height = 4.5, units = "in", p3)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "vij_order.png"),
       width = 6, height = 4.5, units = "in", sim_vij_order_graph)

