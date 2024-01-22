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
# Simualate data Y with CAR model
source(file.path(getwd(), "src", "R", "simulation", "CAR_data_generation.R"))
# Hyperparameters
a0_sigma <- 2
b0_sigma <- .01
nbs_df <- data.frame(
  node1 = rep(seq_len(N), times = lapply(county_nbs, length)),
  node2 = unlist(county_nbs)
)
nbs_df <- nbs_df[nbs_df$node1 < nbs_df$node2,]
data <- list(Y = Y[,1],
             X = X,
             N = N,
             Sigma_scaled = Sigma_scaled,
             mu_phi = rep(0, N),
             N_edges = nrow(nbs_df),
             p = p,
             a0_sigma = a0_sigma,
             b0_sigma = b0_sigma,
             node1 = nbs_df$node1,
             node2 = nbs_df$node2)
inits <- list(list(rho = 0.5),
              list(rho = 0.5),
              list(rho = 0.5),
              list(rho = 0.5),
              list(rho = 0.5))

# sampling parameters
n_samp <- 10000 ## Number of posterior samples to be used for inference
n_chains <- 5 ## Number of different chains
n_burn <- 5000 ## Number of initial runs ("burn" period) for MCMC to converge.
rstan_fit <- stan(file.path(getwd(), "src", "stan", "bym2.stan"),
                  data = data, warmup = n_burn, iter = n_samp + n_burn,
                  chains = n_chains, init = inits,
                  pars = c("beta" , "sigma2", "rho", "phi", "YFit"))
samps <- as.matrix(rstan_fit, pars = "lp__", include = FALSE)
apply(samps, MARGIN = 2, summary)
apply(samps, MARGIN = 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))

par(mfrow = c(1,1))
plot(samps[, "rho"], type = "l")
hist(samps[, "rho"])

phi_sim <- samps[, grepl("phi", colnames(samps))]
rho_sim <- samps[, "rho"]

phi_diffs <- BYM2_StdDiff(phi_sim, rho_sim, Q_scaled, X, ij_list)
diff_credible_intervals <- t(apply(phi_diffs, 2, function(x) {
  quantile(x, c(0.50, 0.025, 0.975))}))
diff_credible_intervals

dev.off()
par(mfrow = c(2,2))
plot(density(abs(phi_diffs[,5])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,15])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,55])), xlim = c(0, 6))
plot(density(abs(phi_diffs[,130])), xlim = c(0, 6))

loss_function <- function(V, epsilon) GeneralLoss(V, epsilon,
                                                  sqrt,
                                                  log)
eps_optim <- optim(median(abs(phi_diffs)), function(e) {
  e_vij <- ComputeSimVij(phi_diffs, ij_list, epsilon = e)
  loss_function(e_vij, epsilon = e)
}, method = "Brent", lower = 0, upper = max(abs(phi_diffs)))
optim_e <- eps_optim$par
true_phi_diffs <- ComputeSimSTDDifferences(matrix(phi, ncol = 1),
                                           Sigma_scaled,
                                           ij_list = ij_list)
mean(abs(true_phi_diffs) > optim_e)
# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
optim_e_vij <- ComputeSimVij(phi_diffs, ij_list,
                             epsilon = optim_e)
eta <- .15
t_seq_length <- 10000
t_seq <- seq(0, max(optim_e_vij) - .001, length.out = t_seq_length)
t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t), numeric(1))
t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t), numeric(1))
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

x1 <- cbind("y" = Y[,1], st_as_sf(county_sp))
x2 <- cbind("phi" = phi, st_as_sf(county_sp))
yfit_sim <- samps[, grepl("YFit", colnames(samps))]
yfit_pmeans <- colMeans(yfit_sim)
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
y_phi_pfitmeans_plot <- grid.arrange(p1, p2, p3, nrow = 1)
y_phi_pfitmeans_plot

# SAVE GGPLOT2 GRAPHS
ggsave(file.path(getwd(), "output", "rstan_BYM2_sim", "BYM2_maps.png"),
       width = 6, height = 4.5, units = "in",
       y_phi_pfitmeans_plot)
ggsave(file.path(getwd(), "output", "rstan_BYM2_sim", "vij_order.png"),
       width = 6, height = 4.5, units = "in",
       sim_vij_order_graph)

# SAVE ROC GRAPH

dev.off()
par(mfrow = c(1, 1))
x <- roc(true_diff, optim_e_vij)
vijs <- list(e2_vij, e3_vij, e4_vij, e5_vij)
png(file.path(getwd(), "output", "rstan_BYM2_sim", "roc.png"))
cols <- c("#000004CC", "#56106ECC", "#BB3754CC", "#F98C0ACC", "#FCFFA4CC")
plot(x, legacy.axes = TRUE, print.auc = TRUE,
     col = cols[3], add = FALSE, lwd = 2, xbreaks = 5)
for (i in c(1:2, 4:5)) {
  y <- roc(true_diff, vijs[[i]])
  col <- cols[i]
  plot(y, add = T, col = col, lwd = 2)
}
legend("right", col = cols, lwd = 2,
       lty = 1, legend = round(c(e2, e3, optim_e, e4, e5), digits = 3),
       title = "epsilon")
dev.off()

