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

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Exact sampling
sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
set.seed(122)
model_rho <- rho
BYM2Sampler <- new(BYM2ExactSampler, X, y, Q_scaled, model_rho)
# Set priors
a_sigma <- 0.1
b_sigma <- 0.1
M_0inv <- diag(rep(1e-4, p))
m_0 <- rep(0, p)
BYM2Sampler$SetPriors(M_0inv, m_0, a_sigma, b_sigma)

n_sim <- 10000
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


phi_diffs <- BYM2_StdDiff(phi_sim, rep(model_rho, n_sim), Q_scaled, X, ij_list)
phi_truediffs <- BYM2_StdDiff(matrix(phi, nrow = 1),
                              rho, Q_scaled, X, ij_list)
# Maximize entropy of posterior distribution with respect to epsilon
loss_function <- function(V, epsilon) -ConditionalEntropy(V)
eps_optim <- optim(median(abs(phi_diffs)), function(e) {
  e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
  loss_function(e_vij, epsilon = e)
}, method = "Brent", lower = 0.0001, upper = 5.0)
optim_e <- eps_optim$par
eps_optim

# number of true disparities at e_CE level
sum(abs(phi_truediffs) > optim_e)

# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
optim_e_vij <- ComputeSimVij(phi_diffs, epsilon = optim_e)

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

e2_vij <- ComputeSimVij(phi_diffs, epsilon = e2)
e3_vij <- ComputeSimVij(phi_diffs, epsilon = e3)
e4_vij <- ComputeSimVij(phi_diffs, epsilon = e4)
e5_vij <- ComputeSimVij(phi_diffs, epsilon = e5)

optim_e_vij_order <- order(optim_e_vij, decreasing = F)
e2_vij_order <- order(e2_vij[optim_e_vij_order], decreasing = F)
e3_vij_order <- order(e3_vij[optim_e_vij_order], decreasing = F)
e4_vij_order <- order(e4_vij[optim_e_vij_order], decreasing = F)
e5_vij_order <- order(e5_vij[optim_e_vij_order], decreasing = F)
true_diff <- (abs(phi_truediffs) > optim_e)
rejection_path <- data.table(
  optim_e_vij = as.numeric(optim_e_vij),
  e2_vij_order = as.numeric(e2_vij),
  e3_vij_order = as.numeric(e3_vij),
  e4_vij_order = as.numeric(e4_vij),
  e5_vij_order = as.numeric(e5_vij),
  true_diff = as.vector(true_diff)
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
             aes(x = optim_e_vij, y = order, color = true_diff,
                 shape = true_diff),
             alpha = 0.5, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of h_ij(", round(optim_e, digits = 3), ")"),
       y = "Rank of h_ij(eps)") +
  theme_bw() +
  scale_color_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = 2, "TRUE" = 7))
sim_vij_order_graph

decisions <- logical(nrow(ij_list))
decisions[optim_e_vij >= optim_t] <- TRUE
print(paste0("Optimal epsilon: ", optim_e))
print(paste0("Optimal t: ", optim_t))
print(ComputeClassificationMetrics(decisions, true_diff))
mean(decisions)

county_sf$y <- y
county_sf$phi <- phi
county_sf$y_pmeans <- colMeans(YFit_sim)

y_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = y)) +
  scale_fill_viridis_c(name = "Simulated Y") +
  theme_bw() +
  coord_sf(crs = st_crs(5070)) +
  theme(legend.position = "bottom")
# alternative: quantile map
# cut_pts <- quantile(y, seq(0, 1, length = 6))
# y_sf <- cbind(y_cut = cut(y, cut_pts, right = FALSE,
#                       include.lowest = TRUE),
#               x1)
# st_crs(y_sf) <- st_crs(county_sf)
# y_map <- ggplot() +
#   geom_sf(data = county_sf) +
#   geom_sf(data = y_sf, aes(fill = y_cut)) +
#   scale_fill_viridis_d(name = "Y quantile", drop = FALSE) +
#   theme_minimal() +
#   theme(legend.position = "bottom", legend.title=element_text(size=10))
cut_pts <- quantile(phi, seq(0, 1, length = 6))
phi_cut_sf <- cbind(phi_cut = cut(phi, cut_pts, right = FALSE,
                                  include.lowest = TRUE),
                        county_sf)
phi_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = phi)) +
  scale_fill_viridis_c(name = "Simulated phi") +
  theme_bw() +
  coord_sf(crs = st_crs(5070)) +
  theme(legend.position = "bottom")
phi_cut_map <- ggplot() +
  geom_sf(data = phi_cut_sf, aes(fill = phi_cut)) +
  scale_fill_viridis_d(name = "Simulated Phi") +
  theme_bw() +
  coord_sf(crs = st_crs(5070)) +
  theme(legend.position = "bottom")
pmeans_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = y_pmeans)) +
  scale_fill_viridis_c(name = "Y posterior mean") +
  theme_bw() +
  coord_sf(crs = st_crs(5070)) +
  theme(legend.position = "bottom")

#eps_seq <- seq(0.001, optim_e * 2.5, length.out = 100)
eps_seq <- seq(0.001, max(abs(phi_diffs)) / 2, length.out = 100)
sim_vij <- ComputeSimVij(phi_diffs, epsilon = eps_seq)
sim_loss <- loss_function(sim_vij, epsilon = eps_seq)
eps_loss_graph <- ggplot() +
  geom_line(data = data.frame(sim_loss = sim_loss, epsilon = eps_seq),
            aes(x = epsilon, y = sim_loss), color = "dodgerblue") +
  labs(x = "Epsilon", y = "Loss") +
  geom_vline(xintercept = optim_e, lwd = 0.8, linetype = "dotted",
             color = "red") +
  theme_bw()

eps_seq <- c(e2, e3, optim_e, e4, e5)
e_vijs <- cbind(e2_vij, e3_vij, optim_e_vij, e4_vij, e5_vij)
FDR_FNR_curves <- lapply(seq_len(ncol(e_vijs)), function(i) {
  v <- e_vijs[, i]
  FDR <- vapply(t_seq, function(t) FDR_estimate(v, t, e = 0), numeric(1))
  FNR <- vapply(t_seq, function(t) FNR_estimate(v, t, e = 0.1), numeric(1))
  sensitivity <-
    data.table(epsilon = round(eps_seq[i], digits = 3),
               t = t_seq,
               FDR = FDR,
               FNR = FNR)
}) %>% do.call(rbind, .)
FDR_FNR_curves[, epsilon := factor(epsilon, levels = sort(unique(epsilon)),
                                   ordered = T)]
FDRt_graph <- ggplot() +
  geom_line(data = FDR_FNR_curves, aes(x = t, y = FDR, group = epsilon,
                                        color = epsilon, linetype = epsilon)) +
  labs(x = "t", y = "Bayesian FDR") +
  scale_color_hue(l = 30) +
  theme_bw()
FNRt_graph <- ggplot() +
  geom_line(data = FDR_FNR_curves, aes(x = t, y = FNR, group = epsilon,
                                        color = epsilon, linetype = epsilon)) +
  labs(x = "t", y = "Bayesian FNR") +
  scale_color_hue(l = 30) +
  theme_bw()
FDRFNR_eps_graph <- ggplot() +
  geom_line(data = FDR_FNR_curves, aes(x = FNR, y = FDR, group = epsilon,
                                        color = epsilon, linetype = epsilon)) +
  labs(x = "Bayesian FNR", y = "Bayesian FDR") +
  scale_color_hue(l = 30) +
  theme_bw()

optim_e_vij_hist <- ggplot() +
  geom_histogram(aes(x = optim_e_vij), fill = "dodgerblue", color = "black",
                 breaks = seq(0, 1, by = .05)) +
  lims(x = c(0, 1)) +
  labs(x = paste0("h_ij(", round(optim_e, digits = 3), ")")) +
  theme_bw()

data_maps <- ggarrange(y_map, pmeans_map, common.legend = FALSE, nrow = 1, legend = "bottom")
epsloss_rhohist <- ggarrange(eps_loss_graph, optim_e_vij_hist, nrow = 1)

ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_data.png"),
       width = 6, height = 4.5, units = "in", y_map, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_phi.png"),
       width = 6, height = 4.5, units = "in", phi_map, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_phi_cut.png"),
       width = 7, height = 4.5, units = "in", phi_cut_map, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "US_sim_yfitpmean.png"),
       width = 6, height = 4.5, units = "in", pmeans_map, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "vij_order.png"),
       width = 6, height = 5.5, units = "in", sim_vij_order_graph, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "eps_loss_graph.png"),
       width = 6, height = 4.5, units = "in", eps_loss_graph, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "optim_e_vij_hist.png"),
       width = 6, height = 4.5, units = "in", optim_e_vij_hist, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "FDRFNR_eps_graph.png"),
       width = 6, height = 4.5, units = "in", FDRFNR_eps_graph, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "y_postpredmeans.png"),
       width = 12, height = 4.5, units = "in", data_maps, dpi = 900)
ggsave(file.path(getwd(), "output", "US_exact_sample_sim", "epsloss_rhohist.png"),
       width = 8, height = 4.5, units = "in", epsloss_rhohist, dpi = 900)

saveRDS(phi_diffs, file.path(getwd(), "output", "US_exact_sample_sim",
                             "phi_diffs.Rds"))
