library(rstan)
library(parallel)
library(data.table)
library(sf)
library(spdep)
library(maps)
library(maptools)
library(magrittr)
library(stringr)
library(ggplot2)
library(fields)
rm(list = ls())
set.seed(113001)

source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))

# Import US counties
county_poly <- maps::map("county", "california", fill = TRUE, plot = FALSE)
county_state <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[1]]))
county_names <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[2]]))
sf_use_s2(TRUE)
county_sp <- maptools::map2SpatialPolygons(county_poly, IDs = county_poly$names)
county_nbs <- poly2nb(county_sp)
no_neighbors <- vapply(county_nbs, function(x) identical(x, 0L), logical(1))
# restrict to connected county map
county_sp <- county_sp[!no_neighbors,]
county_state <- county_state[!no_neighbors]
county_names <- county_names[!no_neighbors]
county_nbs <- poly2nb(county_sp)
county_sf <- st_as_sf(county_sp)
rownames(county_sf) <- NULL
st_crs(county_sf) <- st_crs(st_as_sf(county_poly))

# data generation spatial variance: Matern covariance
county_cent <- st_centroid(st_as_sf(county_sp))
st_crs(county_cent) <- st_crs(county_sf)
dist_matrix <- matrix(st_distance(county_cent), nrow = nrow(county_sf),
                      ncol = nrow(county_sf)) / 1000
#dist_matrix <- st_distance(county_cent[1,], county_cent[2,])
Sigma <- Matern(dist_matrix, range = 0.5 * 100, phi = 1, smoothness = 0.5, nu = 0.5)
Sigma_chol <- chol(Sigma)
Q <- chol2inv(Sigma_chol)
N <- nrow(Q)

adj_df <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
adj_df <- adj_df[adj_df$i < adj_df$j, ]
rownames(adj_df) <- NULL

beta <- c(-5, 0.5)
cent_coords <- st_coordinates(county_cent)
mean_lat <- mean(cent_coords[,2])
#x <- numeric(N)
#high_risk <- cent_coords[,2] > mean_lat
#x[high_risk] <- rnorm(sum(high_risk), mean = 2, sd = 1)
#x[!high_risk] <- rnorm(sum(!high_risk), mean = -2, sd = 1)
x <- rnorm(N, mean = 2, sd = 1)
county_sf$x <- x
X <- cbind(1, x)
E <- ceiling(runif(N, 10000, 5e5))
#E[high_risk] <- ceiling(runif(sum(high_risk), 10000, 50000))
#E[!high_risk] <- ceiling(runif(sum(!high_risk), 50000, 5e5))
sigma2 <- 2
rho <- 0.93
#phi <- solve(Q_scaled_cholR, rnorm(N))
phi <- t(Sigma_chol) %*% rnorm(N)
eps <- rnorm(N)
total_err <- sqrt(sigma2) * (sqrt(rho) * phi + sqrt(1 - rho) * eps)
Y <- rpois(N, exp(log(E) + X %*% beta + total_err))
county_sf$y <- Y
county_sf$E <- E
# analysis parameters
W <- nb2mat(county_nbs, style="B")
D <- diag(rowSums(W))
alpha <- 0.99
Q_analysis <- D - alpha * W
Sigma_analysis <- chol2inv(chol(Q_analysis))
scaling_factor <- exp(mean(log(diag(Sigma_analysis))))
Sigma_analysis <- Sigma_analysis / scaling_factor
Sigma_analysis_chol <- chol(Sigma_analysis)

a0_sigma <- 0.1
b0_sigma <- 0.1
data <- list(
  N = N,
  Sigma_chol = t(Sigma_analysis_chol),
  mu_phi = rep(0, N),
  Y = Y,
  E = E,
  p = ncol(X),
  X = X,
  a0_sigma = a0_sigma,
  b0_sigma = b0_sigma,
  Lambda = eigen(Q_analysis * scaling_factor)$values,
  lambda_rho = 0.2
)
plot(county_sf)

fit1 <- stan(
  file = file.path(getwd(), "src", "stan", "bym2_poisson.stan"),
  pars = c("beta", "phi", "sigma2", "rho", "alpha"),
  data = data,
  chains = 4,
  warmup = 40000,
  iter = 60000,
  cores = 4
)

print(fit1)

samps <- as.matrix(fit1)
#HDInterval::hdi(samps[,"rho"])
phi_samps <- samps[, paste0("phi[", seq_len(N), "]")]
sigma2_samps <- samps[, "sigma2"]
summary(sigma2_samps)
rho_samps <- samps[, "rho"]
summary(rho_samps)

HDInterval::hdi(rho_samps)
V_est <- cov(phi_samps)
n_s <- nrow(phi_samps)
k <- nrow(adj_df)
phi_pmeans <- colMeans(phi_samps)
phi_diffs <- vapply(seq_len(k), function(pair_indx) {
  i <- adj_df[pair_indx,]$i
  j <- adj_df[pair_indx,]$j
  var <- V_est[i, i] + V_est[j, j] - 2 * V_est[i, j]
  (phi_samps[,i] - phi_samps[,j]) / sqrt(var)
}, numeric(n_s))
g_phi_diffs <-  vapply(seq_len(k), function(pair_indx) {
  i <- adj_df[pair_indx,]$i
  j <- adj_df[pair_indx,]$j
  d <- exp(abs(phi_samps[,i] - phi_samps[,j]))
  d
}, numeric(n_s))

phi_truediff <- vapply(seq_len(k), function(pair_indx) {
  i <- adj_df[pair_indx,]$i
  j <- adj_df[pair_indx,]$j
  #var <- V_est[i, i] + V_est[j, j] - 2 * V_est[i, j]
  var <- Sigma[i, i] + Sigma[j, j] - 2 * Sigma[i, j]
  (phi[i] - phi[j]) / sqrt(var)
}, numeric(1))
g_phi_truediff <- vapply(seq_len(k), function(pair_indx) {
  i <- adj_df[pair_indx,]$i
  j <- adj_df[pair_indx,]$j
  #var <- V_est[i, i] + V_est[j, j] - 2 * V_est[i, j]
  log_var <- Sigma[i, i] + Sigma[j, j] - 2 * Sigma[i, j]
  d <- exp(abs(phi[i] - phi[j]))
  d
}, numeric(1))

loss_function <- function(v, epsilon) -ConditionalEntropy(v)
system.time({
  eps_optim <- optim(1, function(e) {
    e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0.0001, upper = 2.0)
})
optim_e <- eps_optim$par
optim_e_vij <- ComputeSimVij(phi_diffs, epsilon = optim_e)
optim_e_vij_order <- order(optim_e_vij, decreasing = F)

true_diff <- abs(phi_truediff) > optim_e
#true_diff <- (abs(true_phi_diffs) > optim_e)
mean(true_diff)
# number of true difference boundaries
sum(true_diff)

# indx <- abs(phi_truediff)  > median(abs(phi_truediff))
indx <- optim_e_vij >= sort(optim_e_vij, decreasing = TRUE)[20]
detected_borders <- adj_df[indx,]
county_sf2 <- county_sf
county_sf2$x <- NULL
node1_all <- county_sf2[detected_borders$i,]
node2_all <- county_sf2[detected_borders$j,]
sf_use_s2(FALSE)
intersections <- lapply(seq_len(sum(indx)), function(i) {
  #print(i)
  node1 <- node1_all[i,]
  node2 <- node2_all[i,]
  suppressMessages(st_intersection(st_buffer(node1, 0.001),
                                   st_buffer(node2, 0.001)))

}) %>%
  do.call(rbind, .)
rates <- Y / E
rates_boundaries_df <- data.frame(node1_rate = rates[adj_df[indx,]$i],
                                  node2_rate = rates[adj_df[indx,]$j])
sum(apply(rates_boundaries_df, 1, function(x) all(x < 0.05)))
county_sf$phi <- phi
rate_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = y / E), color = "black") +
  geom_sf(data = intersections, color = "red", linewidth = 1) +
  scale_fill_viridis_c(name = "Simulated Rate") +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10))
lograte_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = log(y / E)), color = "black") +
  geom_sf(data = intersections, color = "red", linewidth = 1) +
  scale_fill_viridis_c(name = "Simulated Log(Rate)") +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10))
#county_sf$mu <- exp(X %*% beta)
county_sf$phi <- phi
county_sf$gamma <- phi * sqrt(sigma2 * rho)
county_sf$exp_phi <- exp(phi)
county_sf$exp_gamma <- exp(phi * sqrt(sigma2 * rho))
county_sf$phi_pmean <- phi_pmeans
gamma_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = gamma), color = "black") +
  geom_sf(data = intersections, color = "red", linewidth = 1) +
  scale_fill_viridis_c(name = "gamma") +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  labs(title = "Top 20 Log-Rate Disparities") +
  theme(legend.position = "bottom", legend.title=element_text(size=10))
gamma_map


system.time({
  eps_optim2 <- optim(1, function(e) {
    e_vij <- ComputeSimVij(g_phi_diffs, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0.001, upper = 5.0)
})
optim_e_exp <- eps_optim2$par
true_diff_exp <- abs(g_phi_truediff) > optim_e_exp
optim_e_vij_exp <- ComputeSimVij(g_phi_diffs, epsilon = optim_e_exp)

# indx <- abs(phi_truediff)  > median(abs(phi_truediff))
indx <- optim_e_vij_exp >= sort(optim_e_vij_exp, decreasing = TRUE)[20]
detected_borders <- adj_df[indx,]
county_sf2 <- county_sf
county_sf2$x <- NULL
node1_all <- county_sf2[detected_borders$i,]
node2_all <- county_sf2[detected_borders$j,]
sf_use_s2(FALSE)
intersections2 <- lapply(seq_len(sum(indx)), function(i) {
  #print(i)
  node1 <- node1_all[i,]
  node2 <- node2_all[i,]
  suppressMessages(st_intersection(st_buffer(node1, 0.001),
                                   st_buffer(node2, 0.001)))

}) %>%
  do.call(rbind, .)
gamma_map2 <- ggplot() +
  geom_sf(data = county_sf, aes(fill = gamma), color = "black") +
  geom_sf(data = intersections2, color = "red", linewidth = 1) +
  scale_fill_viridis_c(name = "gamma") +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10)) +
  labs(title = "Top 20 IRR Disparities")
gamma_map2


## rejection order graph
e2 <- round(optim_e / 3, digits = 3)
e3 <- round(optim_e * 1.5, digits = 3)
e4 <- round(optim_e * 2, digits = 3)
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
             alpha = 1, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of tau_ij(", round(optim_e, digits = 3), ")"),
       y = "Rank of tau_ij(eps)") +
  theme_bw() +
  scale_color_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = 2, "TRUE" = 7)) +
  theme(legend.position = "bottom")
sim_vij_order_graph

# compute unstandardized difference probabilities
phi_diffs2 <- vapply(seq_len(k), function(pair_indx) {
  i <- adj_df[pair_indx,]$i
  j <- adj_df[pair_indx,]$j
  (phi_samps[,i] - phi_samps[,j])
}, numeric(n_s))

phi_truediff2 <- vapply(seq_len(k), function(pair_indx) {
  i <- adj_df[pair_indx,]$i
  j <- adj_df[pair_indx,]$j
  (phi[i] - phi[j])
}, numeric(1))
system.time({
  eps_optim <- optim(1, function(e) {
    e_vij <- ComputeSimVij(phi_diffs2, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0.0001, upper = 4.0)
})
e <- eps_optim$par
optim_e_vij2 <- ComputeSimVij(phi_diffs2, epsilon = optim_e)
e2 <- round(e / 3, digits = 3)
e3 <- round(e * 1.5, digits = 3)
e4 <- round(e * 2, digits = 3)
e5 <- round(e * 3, digits = 3)
true_diff2 <- abs(phi_truediff2) > e

e2_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e2)
e3_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e3)
e4_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e4)
e5_vij2 <- ComputeSimVij(phi_diffs2, epsilon = e5)
optim_e_vij_order <- order(optim_e_vij2, decreasing = F)
e2_vij2_order <- order(e2_vij2[optim_e_vij_order], decreasing = F)
e3_vij2_order <- order(e3_vij2[optim_e_vij_order], decreasing = F)
e4_vij2_order <- order(e4_vij2[optim_e_vij_order], decreasing = F)
e5_vij2_order <- order(e5_vij2[optim_e_vij_order], decreasing = F)
rejection_path <- data.table(
  optim_e_vij = seq_along(optim_e_vij),
  e2_vij_order = e2_vij2_order,
  e3_vij_order = e3_vij2_order,
  e4_vij_order = e4_vij2_order,
  e5_vij_order = e5_vij2_order,
  true_diff = true_diff2[optim_e_vij_order]
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
sim_vij2_order_graph <- ggplot() +
  geom_point(data = rejection_path,
             aes(x = optim_e_vij, y = order, color = true_diff,
                 shape = true_diff),
             alpha = 1, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of tau_ij(", round(e, digits = 3), ")"),
       y = "Rank of tau_ij(eps)") +
  theme_bw() +
  scale_color_manual(name = paste0("eps = ", round(e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = paste0("eps = ", round(e, digits = 3)),
                     labels = c("No difference", "True difference"),
                     values = c("FALSE" = 2, "TRUE" = 7)) +
  theme(legend.position = "bottom")
sim_vij2_order_graph

# examine top 40 rankings
indx1 <- optim_e_vij >= sort(optim_e_vij)[100]
sum(indx1)
rank_stability_score <- cor(rank(e2_vij[indx1]), rank(e5_vij[indx1]))
indx2 <- optim_e_vij2 >= sort(optim_e_vij2)[100]
rank_stability_score2 <- cor(rank(e2_vij2[indx2]), rank(e5_vij2[indx2]))

# indx3 <- optim_e_vij_exp >= sort(optim_e_vij_exp)[100]
# rank_stability_score3 <- cor(
#   rank(ComputeSimVij(g_phi_diffs, epsilon = optim_e_exp / 3)[indx3]),
#   rank(ComputeSimVij(g_phi_diffs, epsilon = optim_e_exp * 3)[indx3])
# )

Wt2 <- DescTools::KendallW(cbind(rank(e2_vij2[indx2]), rank(e5_vij2[indx2])),
                           correct = TRUE, test = TRUE)
rank_stability_df <- data.table(
  "Difference Type" = c("Standardized Difference", "Unstandardized Difference"),
  "Rank Stability Score" = c(rank_stability_score, rank_stability_score2)
)
rank_stability_df

roc_list <- list(
  "Standardized Difference" = pROC::roc(as.vector(true_diff), as.vector(optim_e_vij)),
  "Unstandardized Difference" = pROC::roc(as.vector(true_diff2), as.vector(optim_e_vij2)),
  "IRR Probability" = pROC::roc(as.vector(true_diff_exp), as.vector(optim_e_vij_exp))
)
auc_values <- vapply(roc_list, function(x) x$auc, numeric(1))
auc_values

tau_df <- data.table(
  diff_prob = rep(c("Standardized", "Unstandardized"), each = nrow(adj_df)),
  e1 = rep(c(round(optim_e / 3, digits = 3),
             round(e / 3, digits = 3)), each = nrow(adj_df)),
  e2 = rep(c(round(optim_e * 3, digits = 3),
             round(e * 3, digits = 3))),
  tau1_rank = c(rank(e2_vij), rank(e2_vij2)),
  tau2_rank = c(rank(e5_vij), rank(e5_vij2)),
  true_diff = c(true_diff, true_diff2)
)

stability_plot <- ggplot(data = tau_df) +
  geom_point(aes(x = tau1_rank, y = tau2_rank, color = true_diff, shape = true_diff)) +
  facet_grid(~diff_prob) +
  scale_color_manual(name = "e_CE Difference",
                     labels = c("No disparity", "True disparity"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = "e_CE Difference",
                     labels = c("No disparity", "True disparity"),
                     values = c("FALSE" = 2, "TRUE" = 7)) +
  labs(x = "Ranks of tau_ij(e_CE / 3)",
       y = "Ranks of tau_ij(e_CE * 3)") +
  theme_bw() +
  theme(legend.position = "bottom")
stability_plot

stability_plot2 <- ggplot(data = tau_df[diff_prob == "Standardized",]) +
  geom_point(aes(x = tau1_rank, y = tau2_rank, color = true_diff, shape = true_diff)) +
  scale_color_manual(name = "e_CE Difference",
                     labels = c("No disparity", "True disparity"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = "e_CE Difference",
                     labels = c("No disparity", "True disparity"),
                     values = c("FALSE" = 2, "TRUE" = 7)) +
  labs(x = "Ranks of tau_ij(e_CE / 3)",
       y = "Ranks of tau_ij(e_CE * 3)") +
  theme_bw() +
  theme(legend.position = "bottom")

roc_plot <- pROC::ggroc(roc_list, aes = c("colour", "linetype"), linewidth = 0.8) +
  geom_abline(intercept = 1, slope = 1, color = "darkgrey", linetype = "dotted") +
  scale_color_discrete(name = "Model") +
  scale_linetype_discrete(name = "Model") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Specificity", y = "Sensitivity")
roc_plot

roc_plot_s <- pROC::ggroc(roc_list[1], linewidth = 0.8) +
  geom_abline(intercept = 1, slope = 1, color = "darkgrey", linetype = "dotted") +
  #scale_color_discrete(name = "Model") +
  #scale_linetype_discrete(name = "Model") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Specificity", y = "Sensitivity")

fig <- ggpubr::ggarrange(gamma_map, gamma_map2, roc_plot,
                         nrow = 1)
ggsave(file.path(getwd(), "output", "poisson_sim", "std_diff_results.png"),
       width = 9, height = 4, units = "in", fig, dpi = 400, scale = 1.7)
