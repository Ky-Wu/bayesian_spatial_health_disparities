# This file tests the Gibbs sampling code to sample from the joint
# posterior f(beta, gamma, sigma2, rho | y).
# Random effects assume a CAR model.
# Rho (spatial percentage of total variation parameter) prior is a PC prior on [0, 1].
# Alpha (CAR spatial smoothing parameter) is also assumed fixed.

library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(ggpubr)
rm(list = ls())

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Gibbs sampler
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_flatbeta_MCMC.cpp"))
# Fixed Rho Exact sampling
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
set.seed(51624)
# Priors
a_sigma <- 0.1
b_sigma <- 0.1
# PC prior, ignore IG prior parameters
a_rho <- 0.0
b_rho <- 0.0
# limits on possible values of rho
lower_rho <- 0
upper_rho <- 1.0
# PC prior: pi(rho) = lambda * exp(-lambda * d(rho)), d(rho) = sqrt(2 * KLD(rho))
lambda_rho <- .0335
BYM2Sampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled)
BYM2Sampler$setPriors(a_sigma, b_sigma, rho_prior_type = "pc",
                      a_rho, b_rho,
                      lower_rho, upper_rho, lambda_rho)
BYM2Sampler$initOLS()

warm_up_cycles <- 100
for (j in seq_len(warm_up_cycles)) {
  print(paste0(j, "/", warm_up_cycles))
  print(paste0("current rho: ", BYM2Sampler$rho))
  print(paste0("current sigma2: ", BYM2Sampler$sigma2))
  BYM2Sampler$burnMCMCSample(100)
}
# DRAW POSTERIOR SAMPLES
n_sim <- 30000
if (!file.exists(file.path(getwd(), "output", "US_gibbs_sample_sim", "samps.Rds"))) {
  system.time({
    samps <- BYM2Sampler$MCMCSample(n_sim)
  })
  saveRDS(samps, file.path(getwd(), "output", "US_gibbs_sample_sim", "samps.Rds"))
} else {
  samps <- readRDS(file.path(getwd(), "output", "US_gibbs_sample_sim", "samps.Rds"))
}

beta_sim <- samps$beta
gamma_sim <- samps$gamma
sigma2_sim <- samps$sigma2
rho_sim <- samps$rho
YFit_sim <- samps$YFit

hist(rho_sim)
hist(sigma2_sim)
hist(beta_sim[,1])
hist(beta_sim[,2])
# trace plots for rho and sigma2
plot(rho_sim, type = "l")
plot(sigma2_sim, type = "l")

phi_sim <- apply(gamma_sim, MARGIN = 2, function(x) {
  x / sqrt(sigma2_sim * rho_sim)
})

phi_diffs <- BYM2_StdDiff(phi_sim, rho_sim, Q_scaled, X, ij_list)

# look at some differences
dev.off()
par(mfrow = c(2,2))
plot(density(abs(phi_diffs[,5])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,15])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,55])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,130])), xlim = c(0, 6))
# compute the true differences
phi_truediffs <- BYM2_StdDiff(matrix(phi, nrow = 1),
                              rho, Q_scaled, X, ij_list)
# Maximize entropy of posterior distribution with respect to epsilon
loss_function <- function(v, epsilon) -ConditionalEntropy(v)
system.time({
  eps_optim <- optim(1, function(e) {
    e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0.0001, upper = 2.0)
})

optim_e <- eps_optim$par

true_diff <- abs(phi_truediffs) > optim_e
#true_diff <- (abs(true_phi_diffs) > optim_e)
mean(true_diff)
# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
optim_e_vij <- ComputeSimVij(phi_diffs, epsilon = optim_e)
optim_e_vij_order <- order(optim_e_vij, decreasing = F)
t_seq_length <- 10000
t_seq <- seq(0, max(optim_e_vij) - .001, length.out = t_seq_length)
t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t, e = 0), numeric(1))
t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t, e = 0), numeric(1))


eta <- .05
optim_t <- min(c(t_seq[t_FDR <= eta], 1))
optim_t
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

### DATA MAPS
x1 <- cbind("y" = y, county_sf)
yfit_pmeans <- colMeans(YFit_sim)
yfit_pmeans_df <- cbind(y_pmeans = yfit_pmeans, county_sf)
y_map <- ggplot() +
  geom_sf(data = x1, aes(fill = y)) +
  scale_fill_viridis_c() +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom")
y_pmeans_map <- ggplot() +
  geom_sf(data = yfit_pmeans_df, aes(fill = y_pmeans)) +
  scale_fill_viridis_c(name = "Y posterior mean") +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom")
y_pmeans_map
y_map

## rejection order graph
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

rejection_path <- data.table(
  optim_e_vij = seq_along(optim_e_vij),
  e2_vij_order = e2_vij_order,
  e3_vij_order = e3_vij_order,
  e4_vij_order = e4_vij_order,
  e5_vij_order = e5_vij_order,
  true_diff = true_diff[optim_e_vij_order]
)
indx <- seq(8820, 9119)
print(paste0("Correlation between top 300 rankings for eps = ", optim_e, " and eps = ", e5, ": ",
             rejection_path[indx, cor(optim_e_vij, e5_vij_order)]))

rejection_path <- melt(rejection_path[indx,],
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
             alpha = 0.3, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of v_ij(", round(optim_e, digits = 3), ")"),
       y = "Rank of v_ij(eps)") +
  theme_bw() +
  scale_color_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = 2, "TRUE" = 7)) +
  theme(legend.position = "bottom")
sim_vij_order_graph

# Bayesian FDR vs. Sensitivity graphs: comparison to exact model
phi_diffs2 <- readRDS(file.path(getwd(), "output", "US_exact_sample_sim",
                                "phi_diffs.Rds"))
optim2_e_vij <- ComputeSimVij(phi_diffs2, epsilon = optim_e)

sim2_FDR_sensitivity <- data.table(
  t = t_seq,
  FDR = vapply(t_seq, function(target_t) FDR_estimate(optim2_e_vij, target_t, e = 0), numeric(1)),
  true_sensitivity = vapply(t_seq, function(target_t) {
    decisions <- logical(nrow(ij_list))
    decisions[optim2_e_vij >= target_t] <- TRUE
    ComputeClassificationMetrics(decisions, true_diff)["sensitivity"]
  }, numeric(1))
)
sim_FDR_sensitivity <- data.table(
  t = t_seq,
  FDR = vapply(t_seq, function(target_t) FDR_estimate(optim_e_vij, target_t, e = 0), numeric(1)),
  true_sensitivity = vapply(t_seq, function(target_t) {
    decisions <- logical(nrow(ij_list))
    decisions[optim_e_vij >= target_t] <- TRUE
    ComputeClassificationMetrics(decisions, true_diff)["sensitivity"]
  }, numeric(1))
)
FDR_sensitivity_df <- rbind(
  cbind(sim_FDR_sensitivity, model = paste0("rho ~ PC(", lambda_rho, ")")),
  cbind(sim2_FDR_sensitivity, model = paste0("rho = ", rho))
)
FDR_sensitivity_df$model <- factor(FDR_sensitivity_df$model,
                                   levels = c(paste0("rho ~ PC(", lambda_rho, ")"),
                                              paste0("rho = ", rho)))
FDR_sensitivity_plot <- ggplot() +
  geom_line(data = FDR_sensitivity_df,
            aes(x = true_sensitivity, y = FDR, col = model, linetype = model)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Sensitivity", y = "Bayesian FDR")
FDR_sensitivity_plot

# Histogram of difference probabilities
optim_e_vij_hist <- ggplot() +
  geom_histogram(aes(x = optim_e_vij), fill = "dodgerblue", color = "black",
                 breaks = seq(0, 1, by = .05)) +
  lims(x = c(0, 1)) +
  labs(x = paste0("v_ij(", round(optim_e, digits = 3), ")")) +
  theme_bw()
# Histogram of posterior samples of rho
rho_hist <- ggplot() +
  geom_histogram(aes(x = rho_sim), fill = "dodgerblue", color = "black") +
  labs(x = paste0("rho")) +
  theme_bw()

roc_list <- list(
  pROC::roc(as.vector(true_diff), as.vector(optim_e_vij)),
  pROC::roc(as.vector(true_diff), as.vector(optim2_e_vij))
)
names(roc_list) <- c(
  paste0("rho ~ PC(", lambda_rho, ")"),
  paste0("rho = ", rho)
)
auc <- lapply(roc_list, function(r) round(pROC::auc(r), 3))
roc_plot <- pROC::ggroc(roc_list, aes = c("colour", "linetype"), linewidth = 0.8) +
  geom_abline(intercept = 1, slope = 1, color = "darkgrey", linetype = "dotted") +
  scale_color_discrete(name = "Model") +
  scale_linetype_discrete(name = "Model") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Specificity", y = "Sensitivity")

ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "vij_order_graph.png"),
       width = 8, height = 5.5, units = "in", sim_vij_order_graph)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "FDR_sensitivity_plot.png"),
       width = 4.5, height = 6, units = "in", FDR_sensitivity_plot)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "optim_e_vij_hist.png"),
       width = 6, height = 4.5, units = "in", optim_e_vij_hist)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "rho_hist.png"),
       width = 6, height = 4.5, units = "in", rho_hist)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "auc.png"),
       width = 4.5, height = 6, units = "in", roc_plot)

auc_FDR <- ggpubr::ggarrange(roc_plot, FDR_sensitivity_plot, common.legend = TRUE,
                             nrow = 1, legend = "bottom")
auc_FDR_vijorder <- ggpubr::ggarrange(auc_FDR, sim_vij_order_graph)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "auc_FDR.png"),
       width = 9, height = 6, units = "in", auc_FDR)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "auc_FDR_vijorder.png"),
       width = 12, height = 6, units = "in", auc_FDR_vijorder)


# look at top 300 rankings
rejection_path <- data.table(
  optim_e_vij = seq_along(optim_e_vij),
  e2_vij_order = e2_vij_order,
  e3_vij_order = e3_vij_order,
  e4_vij_order = e4_vij_order,
  e5_vij_order = e5_vij_order,
  true_diff = true_diff[optim_e_vij_order]
)
indx <- seq(8820, 9119)
rejection_path <- melt(rejection_path[indx,],
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
  geom_point(data = rejection_path, color = "dodgerblue",
             aes(x = optim_e_vij, y = order, color = true_diff),
             alpha = 1, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of v_ij(", round(optim_e, digits = 3), ")"),
       y = "Rank of v_ij(eps)") +
  theme_bw() +
  theme(legend.position = "bottom")
sim_vij_order_graph

ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "top300_vij_order_graph.png"),
       width = 8, height = 5.5, units = "in", sim_vij_order_graph)

# what if we dont divide by conditional posterior variance?
phi_diffs2 <- vapply(seq_len(nrow(ij_list)), function(pair_indx) {
  i <- ij_list[pair_indx,]$i
  j <- ij_list[pair_indx,]$j
  phi_sim[, i] - phi_sim[,j]
}, numeric(nrow(phi_sim)))
phi_truediffs <- vapply(seq_len(nrow(ij_list)), function(pair_indx) {
  i <- ij_list[pair_indx,]$i
  j <- ij_list[pair_indx,]$j
  phi[i] - phi[j]
}, numeric(1))

loss_function <- function(v, epsilon) -ConditionalEntropy(v)
system.time({
  eps_optim2 <- optim(1, function(e) {
    e_vij <- ComputeSimVij(phi_diffs2, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0.0001, upper = 3.0)
})
e <- eps_optim2$par
e_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e)
e2 <- signif(e / 3, 3)
e3 <- signif(e / 1.5, 3)
e4 <- signif(e * 1.5, 3)
e5 <- signif(e * 3, 3)

e2_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e2)
e3_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e3)
e4_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e4)
e5_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e5)
true_diff <- abs(phi_truediffs) > e

e_vij_order2 <- order(e_vij2, decreasing = F)
e2_vij_order2 <- order(e2_vij2[e_vij_order2], decreasing = F)
e3_vij_order2 <- order(e3_vij2[e_vij_order2], decreasing = F)
e4_vij_order2 <- order(e4_vij2[e_vij_order2], decreasing = F)
e5_vij_order2 <- order(e5_vij2[e_vij_order2], decreasing = F)
DescTools::KendallW(cbind(seq_along(e_vij2),
                          e3_vij_order2), test = TRUE)

rejection_path2 <- data.table(
  e_vij = seq_along(e_vij2),
  e2_vij_order = e2_vij_order2,
  e3_vij_order = e3_vij_order2,
  e4_vij_order = e4_vij_order2,
  e5_vij_order = e5_vij_order2,
  true_diff = true_diff[e_vij_order2]
)
print(paste0("Correlation between top 300 rankings for eps = ", e, " and eps = ", e5, ": ",
             rejection_path2[indx, cor(e_vij, e5_vij_order)]))
rejection_path2 <- melt(rejection_path2[indx,],
                        id.vars = c("e_vij", "true_diff"),
                        variable.name = "order_type",
                        value.name = "order")
rejection_path2[, order_type := fcase(
  order_type == "e2_vij_order", paste0("eps = ", e2),
  order_type == "e3_vij_order", paste0("eps = ", e3),
  order_type == "e4_vij_order", paste0("eps = ", e4),
  order_type == "e5_vij_order", paste0("eps = ", e5)
)]
sim_vij_order_graph2 <- ggplot() +
  geom_point(data = rejection_path2, color = "dodgerblue",
             aes(x = e_vij, y = order, color = true_diff),
             alpha = 1, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of v_ij(", round(e, digits = 3), ")"),
       y = "Rank of v_ij(eps)") +
  theme_bw() +
  theme(legend.position = "bottom")
sim_vij_order_graph2

roc_list <- list(
  pROC::roc(as.vector(true_diff), as.vector(e_vij2))
)
names(roc_list) <- c(
  paste0("rho ~ PC(", lambda_rho, ")")
)
auc <- lapply(roc_list, function(r) round(pROC::auc(r), 3))
roc_plot <- pROC::ggroc(roc_list, aes = c("colour", "linetype"), linewidth = 0.8) +
  geom_abline(intercept = 1, slope = 1, color = "darkgrey", linetype = "dotted") +
  scale_color_discrete(name = "Model") +
  scale_linetype_discrete(name = "Model") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Specificity", y = "Sensitivity")

ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "top300_unstd_vij_order_graph.png"),
       width = 8, height = 5.5, units = "in", sim_vij_order_graph2)
