library(data.table)
library(stringr)
library(magrittr)
library(spdep)
library(maps)
library(maptools)
library(ggplot2)
library(rstan)
library(Rcpp)
library(RcppArmadillo)
rm(list = ls())
set.seed(1130)

# BYM2 Model Sampling functions
source(file.path(getwd(), "src", "R", "bym2_sampling.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "vij_computation.R"))
sourceCpp(file.path(getwd(), "src", "rcpp", "epsilonBinarySearch.cpp"))
# Helper functions for epsilon-loss bayesian FDR control procedure
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
# Read in and setup lung data
source(file.path(getwd(), "src", "R", "RDA", "CA_data_setup.R"))

plot(county_sf[, c("young", "old", "highschool", "poverty", "unemployed",
                   "black_percent", "male_percent", "uninsured_percent",
                   "smoking", "standard_ratio")])
plot(county_sf[, "standard_ratio"])
Y <- county_sf[, "standard_ratio", drop = T]
N <- length(Y)
X <- county_sf[, c("young", "old", "highschool", "poverty", "unemployed",
                   "black_percent", "male_percent", "uninsured_percent",
                   "smoking")]
X <- county_sf[, c("smoking")]
# X <- county_sf[, c("smoking", "unemployed")]
X <- cbind(1, X)
colnames(X)[1] <- "Intercept"
X <- as.matrix(st_drop_geometry(X))
p <- ncol(X)

# Test for spatial autocorrelation
ols_fit <- lm(Y ~ X)
e <- residuals(ols_fit)
lw <- nb2listw(county_nbs, style = "W", zero.policy = TRUE)
e_lag <- lag.listw(lw, e)
plot(e_lag ~ e, pch = 16, asp = 1)
M1 <- lm(e_lag ~ e)
abline(M1, col = "blue")
coef(M1)[2]
county_sf$ols_e <- e
moran.test(e, lw, alternative = "greater")
moran.test(Y, lw, alternative = "greater")

# Set hyperparameters
a0_sigma <- 2
b0_sigma <- .01
# Collect data and hyperparameters
# sampling parameters
n_samp <- 10000 ## Number of posterior samples to be used for inference
n_chains <- 5 ## Number of different chains
n_burn <- 5000 ## Number of initial runs ("burn" period) for MCMC to converge.
samps <- DrawBYM2Samples(X, Y, Sigma_scaled, a0_sigma, b0_sigma,
                         node1 = county_nbs$node1,
                         node2 = county_nobs$node2,
                         n_samp, n_chains, n_burn)
par(mfrow = c(1,1))
hist(samps$rho_sim)

phi_sim_pmeans <- apply(samps$phi_sim, 2, mean)
yfit_pmeans <- apply(samps$yfit_sim, 2, mean)
county_sf$yfit_pmean <- yfit_pmeans
county_sf$phi_pmean <- phi_sim_pmeans
rho_sim <- samps$rho_sim
phi_diffs <- BYM2_StdDiff(samps$phi_sim, samps$rho_sim, Q_scaled, X, ij_list)
params <- with(samps, cbind(beta_sim, sigma2_sim, rho_sim))
params_confint <- apply(params, 2, function(x) {
  pmean <- round(mean(x), digits = 3)
  low <- round(quantile(x, .025), digits = 3)
  high <- round(quantile(x, .975), digits = 3)
  paste0(pmean, " (", low, ", ", high, ")")
})
names(params_confint) <- NULL
params_confint <- data.table(term = colnames(X) %>%
                             str_replace_all("_", " ") %>%
                             str_to_title() %>%
                             c("sigma2", "rho"),
                           confint = params_confint)

# Pick maximum epsilon such that Bayesian FDR can be controlled
# under threshold eta with T positive results
eta <- .10
min_pos <- 50
optim_e <- epsilonBinarySearch(abs(phi_diffs),
                               eta, n_positive = 50,
                               lower = 0, upper = max(phi_diffs), toler = 1e-10)
optim_e_vij <- ComputeSimVij(phi_diffs, ij_list,
                             epsilon = optim_e)
optim_e_vij_order <- order(optim_e_vij, decreasing = F)
optim_t <- sort(optim_e_vij, decreasing = T)[min_pos]

# Relative rankings at surrounding values of epsilon
vij_order_graph <- VijOrderGraph(phi_diffs, ij_list, optim_e, T_line = min_pos)
vij_order_graph

decisions <- logical(nrow(ij_list))
decisions[optim_e_vij >= optim_t] <- TRUE
sum(decisions)
rej_indx <- optim_e_vij_order[optim_e_vij_order %in%
                                which(optim_e_vij >= optim_t)]
print(paste0("Optimal epsilon: ", optim_e))
print(paste0("Optimal t: ", optim_t))
decisions
detected_borders <- ij_list[rej_indx,]
node1_all <- county_sf[detected_borders$i,]
node2_all <- county_sf[detected_borders$j,]
intersections <- lapply(seq_len(sum(decisions)), function(i) {
  node1 <- node1_all[i,]
  node2 <- node2_all[i,]
  st_intersection(st_buffer(node1, 0.01), st_buffer(node2, .01))
}) %>%
  do.call(rbind, .)
intersections$rank <- rev(seq_len(sum(decisions)))
lower <- with(county_sf, min(c(standard_ratio, yfit_pmeans)))
upper <- with(county_sf, max(c(standard_ratio, yfit_pmeans)))
beta_sim <- samps$beta_sim
Xb <- X %*% t(beta_sim)
mus <- apply(Xb, 1, mean)
Xb_plus_phi <- Xb + t(samps$phi_sim * sqrt(samps$sigma2_sim) * sqrt(samps$rho_sim))
county_sf$mu_pmean <- mus
county_sf$Xb_plus_phi <- apply(Xb_plus_phi, 1, mean)

SIR_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = standard_ratio)) +
  scale_fill_viridis_c(name = "SIR", limits = c(lower, upper)) +
  scale_linewidth_identity() +
  theme_minimal()
SIR_map_boundaries <- SIR_map +
  geom_sf(data = intersections, col = "red", alpha = 0,
          aes(linewidth = rank / 80 + .2))
y_pmeans_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = yfit_pmeans)) +
  geom_sf(data = intersections, col = "red", alpha = 0,
          aes(linewidth = rank / 80 + .2)) +
  scale_fill_viridis_c(name = "Y Posterior mean") +
  scale_linewidth_identity() +
  theme_minimal()
Xb_pmeans_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = mu_pmean)) +
  geom_sf(data = intersections, col = "red", alpha = 0,
          aes(linewidth = rank / 80 + .2)) +
  scale_fill_viridis_c(name = "Xb posterior mean") +
  scale_linewidth_identity() +
  theme_minimal()
phi_postmeans_map <- ggplot() +
  geom_sf(data = county_sf, aes(fill = phi_pmean), col = "black") +
  geom_sf(data = intersections, col = "red", alpha = 0,
          aes(linewidth = rank / 80 + .2)) +
  scale_fill_viridis_c(name = "Phi Posterior Mean") +
  scale_linewidth_identity() +
  theme_minimal()
phi_postmeans_map

# ggsave(file.path(getwd(), "output", "RDA", "SIR_map.png"),
#        width = 6, height = 4.5, units = "in",
#        SIR_map)
# ggsave(file.path(getwd(), "output", "RDA", "SIR_map_boundaries.png"),
#        width = 6, height = 4.5, units = "in",
#        SIR_map_boundaries)
# ggsave(file.path(getwd(), "output", "RDA", "y_pmeans_map.png"),
#        width = 6, height = 4.5, units = "in",
#        y_pmeans_map)
# ggsave(file.path(getwd(), "output", "RDA", "phi_postmeans.png"),
#        width = 6, height = 4.5, units = "in",
#        phi_postmeans_map)
# ggsave(file.path(getwd(), "output", "RDA", "vij_order.png"),
#        width = 6, height = 4.5, units = "in",
#        vij_order_graph)
# fwrite(beta_confint, file.path(getwd(), "output", "RDA", "beta_confint.csv"))
