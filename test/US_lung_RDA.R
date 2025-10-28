library(data.table)
library(stringr)
library(magrittr)
library(spdep)
library(maps)
library(maptools)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(ggpubr)
library(xtable)
rm(list = ls())
set.seed(1969)

# BYM2 Model Sampling functions
source(file.path(getwd(), "src", "R", "bym2_sampling.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Helper functions for epsilon-loss bayesian FDR control procedure
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
# Read in and setup lung + smoking data
source(file.path(getwd(), "src", "R", "RDA", "US_data_setup.R"))
# Exact sampling and Gibbs sampling
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_flatbeta_MCMC.cpp"))

X <- with(cancer_smoking_sf, {
  cbind(1, total_mean_smoking,
        unemployed_2014, SVI_2014, physical_inactivity_2014, uninsured_2012_2016)
})
y <- cancer_smoking_sf$mortality2014
# Test for spatial autocorrelation
ols_fit <- lm(y ~ X)
e <- residuals(ols_fit)
lw <- nb2listw(county_nbs, style = "W", zero.policy = TRUE)
e_lag <- lag.listw(lw, e)
plot(e_lag ~ e, pch = 16, asp = 1)
M1 <- lm(e_lag ~ e)
abline(M1, col = "blue")
coef(M1)[2]
moran.mc(e, lw, nsim = 10000, alternative = "greater")
geary.mc(e, lw, nsim = 10000, alternative = "greater")

# Priors
a_sigma <- 0.1
b_sigma <- 0.1
a_rho <- 0
b_rho <- 0
# limits on possible values of rho
lower_rho <- 0.00
upper_rho <- 1.0
lambda_rho <- 0.0335
gibbsSampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled)
gibbsSampler$setPriors(a_sigma, b_sigma, rho_prior_type = "pc",
                       a_rho, b_rho,
                       lower_rho, upper_rho, lambda_rho)
gibbsSampler$initOLS()

print("10000 Burn in:")
for(j in seq_len(100)) {
  print(paste0(j, "/100"))
  print(paste0("current rho: ", gibbsSampler$rho))
  print(paste0("current sigma2: ", gibbsSampler$sigma2))
  gibbsSampler$burnMCMCSample(100)
}

# DRAW POSTERIOR SAMPLES
n_sim <- 30000
system.time({
  samps <- gibbsSampler$MCMCSample(n_sim)
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

# Maximize conditional entropy with respect to epsilon
loss_function <- function(V, epsilon) -ConditionalEntropy(V)
eps_optim <- optim(median(abs(phi_diffs)), function(e) {
  e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
  loss_function(e_vij, epsilon = e)
}, method = "Brent", lower = 0.0001, upper = 5.0)
optim_e <- eps_optim$par

optim_e_vij <- ComputeSimVij(phi_diffs,
                             epsilon = optim_e)
# select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
delta <- .05
FDR_res <- computeFDRCutoff(optim_e_vij, delta = delta)
optim_t <- FDR_res$cutoff

sum(optim_e_vij >= optim_t)

decisions <- logical(nrow(ij_list))
decisions[optim_e_vij >= optim_t] <- TRUE
print(paste0("Optimal epsilon: ", optim_e))
print(paste0("Optimal t: ", optim_t))
mean(decisions)

## rejection order graph
optim_e_vij_order <- order(optim_e_vij, decreasing = F)
rej_indx <- optim_e_vij_order[optim_e_vij_order %in%
                                which(optim_e_vij >= optim_t)]
detected_borders <- ij_list[rej_indx,]
node1_all <- cancer_smoking_sf[detected_borders$i,]
node2_all <- cancer_smoking_sf[detected_borders$j,]
intersections <- lapply(seq_len(sum(decisions)), function(i) {
  node1 <- node1_all[i,]
  node2 <- node2_all[i,]
  suppressMessages(st_intersection(st_buffer(node1, 0.0003),
                                   st_buffer(node2, 0.0003)))

}) %>%
  do.call(rbind, .)
intersections$rank <- rev(seq_len(sum(decisions)))

e2 <- round(optim_e / 3, digits = 3)
e3 <- round(optim_e / 1.5, digits = 3)
e4 <- round(optim_e * 1.5, digits = 3)
e5 <- round(optim_e * 3, digits = 3)

e2_vij <- ComputeSimVij(phi_diffs, epsilon = e2)
e3_vij <- ComputeSimVij(phi_diffs, epsilon = e3)
e4_vij <- ComputeSimVij(phi_diffs, epsilon = e4)
e5_vij <- ComputeSimVij(phi_diffs, epsilon = e5)

# Note: the R order function returns permutation index to sort vector, NOT vector of ranks
optim_e_vij_order <- order(optim_e_vij, decreasing = F)
e2_vij_order <- order(e2_vij[optim_e_vij_order], decreasing = F)
e3_vij_order <- order(e3_vij[optim_e_vij_order], decreasing = F)
e4_vij_order <- order(e4_vij[optim_e_vij_order], decreasing = F)
e5_vij_order <- order(e5_vij[optim_e_vij_order], decreasing = F)

rejection_path <- data.table(
  optim_e_order = seq_along(optim_e_vij),
  e2_vij_order = e2_vij_order,
  e3_vij_order = e3_vij_order,
  e4_vij_order = e4_vij_order,
  e5_vij_order = e5_vij_order,
  classified_diff = decisions[optim_e_vij_order])

rejection_path <- melt(rejection_path,
                       id.vars = c("optim_e_order", "classified_diff"),
                       variable.name = "order_type",
                       value.name = "order")
rejection_path[, order_type := fcase(
  order_type == "e2_vij_order", paste0("eps = ", e2),
  order_type == "e3_vij_order", paste0("eps = ", e3),
  order_type == "e4_vij_order", paste0("eps = ", e4),
  order_type == "e5_vij_order", paste0("eps = ", e5)
)]
vij_order_graph <- ggplot() +
  geom_point(data = rejection_path,
             aes(x = optim_e_order, y = order, color = classified_diff,
                 shape = classified_diff),
             alpha = 0.3, size = 1) +
  #geom_vline(xintercept = nrow(ij_list) - sum(optim_e_vij == 1)) +
  facet_grid(~order_type) +
  labs(x = paste0("Rank of v_ij(", round(optim_e, digits = 3), ")"),
       y = "Rank of v_ij(eps)") +
  theme_bw() +
  scale_color_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("Not classified difference", "Classified difference"),
                     values = c("FALSE" = "darkgrey", "TRUE" = "dodgerblue")) +
  scale_shape_manual(name = paste0("eps = ", round(optim_e, digits = 3)),
                     labels = c("Not classified difference", "Classified difference"),
                     values = c(2, 7)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 2, alpha = 1)))
vij_order_graph

optim_e_vij_hist <- ggplot() +
  geom_histogram(aes(x = optim_e_vij), fill = "dodgerblue", color = "black",
                 breaks = seq(0, 1, by = .05)) +
  lims(x = c(0, 1)) +
  labs(x = paste0("v_ij(", round(optim_e, digits = 3), ")")) +
  theme_bw()

yfit_pmeans <- colMeans(YFit_sim)
cut_pts <- quantile(y, seq(0, 1, length = 6))
yfit_pmeans_sf <- cbind(y_pmeans = cut(yfit_pmeans, cut_pts, right = FALSE,
                                       include.lowest = TRUE),
                        cancer_smoking_sf)

lung_map <- ggplot() +
  geom_sf(data = county_sf) +
  geom_sf(data = cancer_smoking_sf[!is.na(cancer_smoking_sf$mortality2014),],
          aes(fill = cut(mortality2014, cut_pts, right = FALSE,
                         include.lowest = TRUE)), col = "gray") +
  scale_fill_viridis_d(name = "Mortality",
                       drop = FALSE) +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10))
smoking_map <- ggplot() +
  geom_sf(data = county_sf) +
  geom_sf(data = cancer_smoking_sf[!is.na(cancer_smoking_sf$mortality2014),],
          aes(fill = cut(total_mean_smoking,
                         quantile(total_mean_smoking, seq(0, 1, length = 6)),
                         right = FALSE, include.lowest = TRUE)),
          col = "gray") +
  scale_fill_viridis_d(name = "Smoking",
                       drop = FALSE) +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10))
ols_e_map <- ggplot() +
  geom_sf(data = county_sf) +
  geom_sf(data = cancer_smoking_sf, col = "gray",
          aes(fill = cut(y - e, quantile(y - e, seq(0, 1, 0.2))))) +
  geom_sf(data = intersections, col = "red", alpha = 0, lwd = 0.4) +
  scale_fill_viridis_d() +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10))

yfit_pmeans_map <- ggplot() +
  geom_sf(data = county_sf, col = "gray") +
  geom_sf(data = yfit_pmeans_sf, aes(fill = y_pmeans), col = "gray") +
  scale_fill_viridis_d(name = "Lung Cancer Mortality Posterior Mean", drop = FALSE) +
  geom_sf(data = intersections, col = "red", alpha = 0, lwd = 0.4) +
  coord_sf(crs = st_crs(5070)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title=element_text(size=10))
yfit_pmeans_map
rho_hist <- ggplot() +
  geom_histogram(aes(x = rho_sim), fill = "dodgerblue", color = "black",
                 breaks = seq(0, 1, by = .05)) +
  lims(x = c(0, 1)) +
  labs(x = paste0("rho")) +
  theme_bw()
data_maps <- ggarrange(lung_map, smoking_map, common.legend = FALSE, nrow = 1, legend = "bottom")
ggsave(file.path(getwd(), "output", "RDA", "US_data", "vij_order_graph.png"),
       width = 8, height = 5.5, units = "in", vij_order_graph, dpi = 900)
ggsave(file.path(getwd(), "output", "RDA", "US_data", "lung_map.png"),
       width = 6, height = 4.5, units = "in", lung_map, dpi = 900)
ggsave(file.path(getwd(), "output", "RDA", "US_data", "smoking_map.png"),
       width = 6, height = 4.5, units = "in", smoking_map, dpi = 900)
ggsave(file.path(getwd(), "output", "RDA", "US_data", "yfit_pmeans_map.png"),
       width = 8, height = 4.5, units = "in", yfit_pmeans_map, dpi = 900)
ggsave(file.path(getwd(), "output", "RDA", "US_data", "optim_e_vij_hist.png"),
       width = 6, height = 4.5, units = "in", optim_e_vij_hist)
ggsave(file.path(getwd(), "output", "RDA", "US_data", "rho_hist.png"),
       width = 4.5, height = 6, units = "in", rho_hist, dpi = 900)
ggsave(file.path(getwd(), "output", "RDA", "US_data", "data_maps.png"),
       width = 12, height = 4.5, units = "in", data_maps, dpi = 900)
params_sim <- cbind(samps$beta, samps$sigma2, samps$rho)
alpha <- 0.05
params_summary <- apply(params_sim, MARGIN = 2, function(x) {
  data.table(mean = mean(x), lower = quantile(x, alpha / 2),
             upper = quantile(x, 1 - alpha / 2))
}) %>%
  do.call(rbind, .)
round_digits <- 3
params_summary[, `:=`(
  Parameter = c("$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$",
                "$\\beta_4$", "$\\beta_5$", "$\\sigma^2$", "$\\rho$"),
  Description = c("Intercept",
                  "Smoking prevalence",
                  "Unemployment Percentile",
                  "SVI Percentile",
                  "Physically Inactive (\\%)",
                  "Uninsured (\\%)",
                  "Total variance",
                  "Spatial proportion of variance"),
  `Mean` = round(mean, digits = round_digits),
  `95\\% Credible Interval` = paste0("(", round(lower, digits = round_digits), ", ",
                    round(upper, digits = round_digits), ")")
)]
params_summary <- params_summary[, .(Parameter, Description, `Mean`, `95\\% Credible Interval`)]
rownames(params_summary) <- NULL
print(xtable(params_summary, type = "latex",
             caption = "Posterior summaries of non-spatial effect parameters.",
             label = "tab:RDA_post_summaries", align = c("c", "c", "c", "c", "c")),
      type = "latex", sanitize.text.function = function(x) {x},
      include.rownames = FALSE,
      file = file.path(getwd(), "output", "RDA", "US_data", "params_summary.tex"))

rankings <- order(optim_e_vij[decisions], decreasing = TRUE)
node1_indx <- ij_list[decisions,]$i[rankings]
node2_indx <- ij_list[decisions,]$j[rankings]
disparity_df <- data.table(
  `County 1` = cancer_smoking_sf$location[node1_indx],
  `County 2` = cancer_smoking_sf$location[node2_indx],
  `$v_k(\\epsilon_{CE})$` = optim_e_vij[decisions][rankings]
)
print(xtable(disparity_df, type = "latex", digits = 3,
             caption = "Reported spatial disparities in lung cancer mortality between neighboring US counties.",
             label = "tab:RDA_disparities_table", align = c("c", "c", "c", "c")),
      type = "latex", sanitize.text.function = function(x) {x},
      include.rownames = TRUE, tabular.environment = "longtable",
      file = file.path(getwd(), "output", "RDA", "US_data", "disparity_tab.tex"))

