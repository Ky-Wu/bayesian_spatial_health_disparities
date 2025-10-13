# This file tests the Gibbs sampling code to sample from the joint
# posterior f(beta, gamma, sigma2, rho | y).
# Random effects assume a CAR model.
# Rho (spatial percentage of total variation parameter) prior is a PC prior on [0, 1].
# Alpha (CAR spatial smoothing parameter) is also assumed fixed.
# Data is generated over a map of US counties.

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

delta <- .05
optim_t <- computeFDRCutoff(optim_e_vij, delta = delta)
optim_t <- optim_t$cutoff
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
e2 <- round(optim_e * 1/3, digits = 3)
e3 <- round(optim_e * 2/3, digits = 3)
e4 <- round(optim_e * 4/3, digits = 3)
e5 <- round(optim_e * 5/3, digits = 3)

e2_vij <- ComputeSimVij(phi_diffs, epsilon = e2)
e3_vij <- ComputeSimVij(phi_diffs, epsilon = e3)
e4_vij <- ComputeSimVij(phi_diffs, epsilon = e4)
e5_vij <- ComputeSimVij(phi_diffs, epsilon = e5)



rej_dt <- mapply(function(vij, e) {
  t_star <- computeFDRCutoff(vij, delta = delta)
  decisions <- vij >= t_star$cutoff
  n_boundaries <- sum(decisions)
  true_diff <- abs(phi_truediffs) > e
  n_true_diff <- sum(true_diff)
  res <- ComputeClassificationMetrics(as.vector(decisions), as.vector(true_diff))
  data.frame(epsilon = e,
             t_star = t_star$cutoff,
             n_boundaries = n_boundaries,
             n_true_diff = n_true_diff,
             sensitivity = res["sensitivity"],
             specificity = res["specificity"],
             accuracy = res["accuracy"],
             fdr = res["fdr"],
             F1 = res["F1"])
}, list(e2_vij, e3_vij, e4_vij, e5_vij, optim_e_vij),
c(e2, e3, e4, e5, optim_e), SIMPLIFY = FALSE) %>%
  do.call(rbind, .)
rocs <- mapply(function(vij, e) {
  true_diff <- abs(phi_truediffs) > e
  pROC::roc(as.vector(true_diff), as.vector(vij))
}, list(e2_vij, e3_vij, optim_e_vij, e4_vij, e5_vij),
c(e2, e3, optim_e, e4, e5), SIMPLIFY = FALSE)
names(rocs) <- paste0("epsilon = ", signif(c(e2, e3, optim_e, e4, e5), 4))
auc <- lapply(rocs, function(r) round(pROC::auc(r), 3))
roc_plot <- pROC::ggroc(rocs, aes = c("colour", "linetype"), linewidth = 0.8) +
  geom_abline(intercept = 1, slope = 1, color = "darkgrey", linetype = "dotted") +
  scale_color_discrete(name = "Difference Threshold") +
  scale_linetype_discrete(name = "Difference Threshold") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Specificity", y = "Sensitivity")
roc_plot

rownames(rej_dt) <- NULL
rej_dt <- rej_dt[order(rej_dt$epsilon),]
rej_dt

colnames(rej_dt) <- c("$\\epsilon$", "$t_{\\star}$", "Detected Boundaries", "True Boundaries",
                      "Sensitivity", "Specificity", "Accuracy", "True FDR", "F1")
rej_dt <- rej_dt[, c(1, 2, 3, 4, 5, 6, 7)]
print(xtable::xtable(rej_dt, caption = "Cutoff probabilities, number of detected/true difference boundaries, and
                     classification performance with Bayesian FDR tolerance $\\delta = 0.05$
                     across different values of difference threshold $\\epsilon$.",
                     label = "tab:epsilon_sensitivity",
                     align = rep("c", ncol(rej_dt) + 1), digits = 3),
      type = "latex", include.rownames = FALSE, sanitize.text.function = function(x) {x},
      file.path(getwd(), "output", "US_gibbs_sample_sim", "epsilon_sensitivity_table.tex"))

optim_e_vij_order <- order(optim_e_vij, decreasing = F)
e2_vij_order <- order(e2_vij[optim_e_vij_order], decreasing = F)
e3_vij_order <- order(e3_vij[optim_e_vij_order], decreasing = F)
e4_vij_order <- order(e4_vij[optim_e_vij_order], decreasing = F)
e5_vij_order <- order(e5_vij[optim_e_vij_order], decreasing = F)

rejection_path <- data.table(
  optim_e_vij = as.numeric(optim_e_vij),
  e2_vij_order = as.numeric(e2_vij),
  e3_vij_order = as.numeric(e3_vij),
  e4_vij_order = as.numeric(e4_vij),
  e5_vij_order = as.numeric(e5_vij),
  true_diff = as.vector(true_diff)#[optim_e_vij_order]
)
indx <- seq(9020, 9119)
print(paste0("Correlation between top 300 rankings for eps = ", optim_e, " and eps = ", e5, ": ",
             rejection_path[indx, cor(optim_e_vij, e5_vij_order)]))

rejection_path <- melt(rejection_path,#[indx,],
                       id.vars = c("optim_e_vij", "true_diff"),
                       variable.name = "order_type",
                       value.name = "order")
rejection_path[, order_type := fcase(
  order_type == "e2_vij_order", paste0("eps = ", e2),
  order_type == "e3_vij_order", paste0("eps = ", e3),
  order_type == "e4_vij_order", paste0("eps = ", e4),
  order_type == "e5_vij_order", paste0("eps = ", e5)
), ]
rejection_path[, hline := fcase(
  order_type == "e2_vij_order", paste0("eps = ", e2),
  order_type == "e3_vij_order", paste0("eps = ", e3),
  order_type == "e4_vij_order", paste0("eps = ", e4),
  order_type == "e5_vij_order", paste0("eps = ", e5)
), ]
sim_vij_order_graph <- ggplot() +
  geom_point(data = rejection_path,
             aes(x = optim_e_vij, y = order, color = true_diff,
                 shape = true_diff),
             alpha = 0.3, size = 1) +
  #geom_hline(yintercept = optim_t, color = "darkgreen", linetype = 2) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("v_ij(", round(optim_e, digits = 3), ")"),
       y = "v_ij(eps)") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
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
system.time({
  eps_optim2 <- optim(1, function(e) {
    e_vij <- ComputeSimVij(phi_diffs2, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0.0001, upper = 3.0)
})
optim2_e <- eps_optim2$par
optim2_e_vij <- ComputeSimVij(phi_diffs2, epsilon = optim2_e)
optim_t <- computeFDRCutoff(optim2_e_vij, delta = 0.05)$cutoff
decisions <- logical(nrow(optim2_e_vij))
decisions[as.vector(optim2_e_vij) >= optim_t] <- TRUE
sum(decisions)
true_diff <- abs(phi_truediffs) > optim_e
true_diff2 <- abs(phi_truediffs) > optim2_e
ComputeClassificationMetrics(decisions, true_diff2)

comparison_df <- data.frame(
  prior_vij = optim_e_vij,
  exact_vij = optim2_e_vij,
  true_diff = as.vector(true_diff)
)
prob_comparison <- ggplot(data = comparison_df) +
  geom_point(aes(x = prior_vij, y = exact_vij, color = true_diff, shape = true_diff), alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  labs(x = "Difference Probability (PC Prior Model)", y = "Difference Probability (Conditioned Model)") +
  scale_color_discrete(name = "County Pair", labels = c("No Disparity", "Disparity")) +
  scale_shape_discrete(name = "County Pair", labels = c("No Disparity", "Disparity"))


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
  pROC::roc(as.vector(true_diff2), as.vector(optim2_e_vij))
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
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "prob_comparison.png"),
       width = 8, height = 5.5, units = "in", prob_comparison)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "optim_e_vij_hist.png"),
       width = 6, height = 4.5, units = "in", optim_e_vij_hist)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "rho_hist.png"),
       width = 6, height = 4.5, units = "in", rho_hist)
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "auc.png"),
       width = 4.5, height = 6, units = "in", roc_plot)

roc_probs <- ggpubr::ggarrange(roc_plot, prob_comparison, legend = "bottom")
ggsave(file.path(getwd(), "output", "US_gibbs_sample_sim", "roc_probs.png"),
       width = 9, height = 6, units = "in", roc_probs)
