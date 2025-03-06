library(data.table)
library(ggplot2)
library(pROC)
set.seed(1130)

rm(list = ls())
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))

############# EXAMPLE 1 ###################
# Simulation parameters
## maximum allowable Bayesian FDR
eta <- .05
## number of observations to generate
n_sim = 10000
## loss function as heuristic to choose epsilon threshold
loss <- function(V, epsilon) GeneralLoss(V, epsilon, sqrt, log)
## Data generation parameters
n <- 1000
sigma2 <- 16
possible_beta <- c(0, 0.5, 1)
p <- 100
beta <- sample(possible_beta, p, replace = TRUE)
print(beta)
p <- length(beta)
# X <- matrix(rep(diag(p), n), ncol = p, byrow = T)
X <- matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p)
# X <- matrix(runif(n * p, -1, 1), n, p)
Y <- rnorm(n, mean = X %*% beta, sd = sqrt(sigma2))
# get pairs to compare
ij_list <- expand.grid(i = seq_len(p), j = seq_len(p))
ij_list <- ij_list[ij_list$i < ij_list$j,]
ij_list <- within(ij_list, {
  difference <- beta[i] != beta[j]
})
# randomly subset ij list
ij_list_no_diff <- ij_list[!ij_list$difference,]
ij_list_no_diff_s <- ij_list_no_diff[sample.int(nrow(ij_list_no_diff),
                                                nrow(ij_list_no_diff) / 2),]
ij_list_diff <- ij_list[ij_list$difference,]
ij_list_diff_s <- ij_list_diff[sample.int(nrow(ij_list_diff),
                                          nrow(ij_list_diff) / 20),]
ij_list_s <- rbind(ij_list_no_diff_s, ij_list_diff_s)
#ij_list_s <- ij_list
mean(ij_list_s$difference)
sum(ij_list_s$difference)
nrow(ij_list_s)
# simulate decisions under epsilon-t procedure
ij_list_s$difference <- NULL

sim_decisions <- SimDecisions(X, beta, sigma2, loss, eta = .05, n_sim,
                              t_seq_length = 300, Y = Y, ij_list = ij_list_s)
eps_loss_graph <- GGPlotEpsLoss(sim_decisions, loss)
FDRFNR_graph <- GGPlotFDRFNR(sim_decisions)
optim_e <- sim_decisions$optim_params$optim_e

# COMPARE TO FREQUENTIST DECISIONS
alpha <- .05
freq_p <- ComputeFrequentistPValues(X, Y, ij_list_s)
freq_p_adj <- p.adjust(freq_p, method = "BY")
freq_decisions <- freq_p_adj <= alpha

# Prior hyperparameters (relatively non-informative priors)
M0 <- diag(p) * 100
M0_inv <- chol2inv(chol(M0))
m0 <- rep(0, p)
a0 <- .01
b0 <- .01
# Compute parameters of posterior
a1 <- a0 + n/2
M1 <- chol2inv(chol(M0_inv + t(X) %*% X))
m1 <- m0 + t(X) %*% Y
b1 <- b0 + 0.5 * (t(Y) %*% Y + t(m0) %*% M0 %*% m0 - t(m1) %*% M1 %*% m1)

results <- sim_decisions$ij_decisions
results$freq_decision <- freq_decisions
eps_differences <- vapply(seq_len(nrow(results)), function(r) {
  i <- results[r,]$i
  j <- results[r,]$j
  sd <- sqrt(M1[i, i] + M1[j, j] - 2 * M1[i, j])
  (abs(beta[i] - beta[j]) / sd) > optim_e
}, logical(1))

print("frequentist decision metrics:")
ComputeClassificationMetrics(results$freq_decision, results$difference_truth)
print("bayesian epsilon method decision metrics (any difference):")
ComputeClassificationMetrics(results$decision, results$difference_truth)
print("bayesian epsilon method decision metrics (epsilon difference):")
ComputeClassificationMetrics(results$decision, eps_differences)
# compare rejection paths
# frequentist: low p-value => reject H_0 => beta_i != beta_j
# bayesian procedure: high v_ij => "reject" |beta_i - beta_j| <= epsilon => declare |beta_i - beta_j| > epsilon
freq_p
bayes_vij <- ComputeSimVij(sim_decisions$sim$std_differences,
                           epsilon = optim_e)
#bayes_vij <- as.vector(sim_decisions$optim_params$optim_e_vij)
bayes_vij_order <- order(bayes_vij, decreasing = F)
freq_pval_order <- order(freq_p[bayes_vij_order], decreasing = F)
rejection_path <- data.frame(
  bayes_vij_order = seq_along(bayes_vij),
  freq_pval_order = freq_pval_order,
  true_diff = results$difference_truth[bayes_vij_order],
  epsilon = paste0("epsilon = ", round(optim_e, digits = 3))
)

e_2 <- optim_e * 4
bayes_vij2 <- ComputeSimVij(sim_decisions$sim$std_differences,
                           epsilon = e_2)
#bayes_vij <- as.vector(sim_decisions$optim_params$optim_e_vij)
bayes_vij2_order <- order(bayes_vij2, decreasing = F)
freq_pval2_order <- order(freq_p[bayes_vij2_order], decreasing = F)
rejection_path2 <- data.frame(
  bayes_vij_order = seq_along(bayes_vij2),
  freq_pval_order = freq_pval2_order,
  true_diff = results$difference_truth[bayes_vij2_order],
  epsilon = paste0("epsilon = ", round(e_2, digits = 3))
)
all_rejection_paths <- rbind(rejection_path, rejection_path2)


sim_vij_pval_order_graph <- ggplot() +
  geom_point(data = all_rejection_paths,
             aes(x = bayes_vij_order, y = freq_pval_order,
                 color = true_diff),
             alpha = 0.3, size = 1) +
  facet_grid(~epsilon) +
  labs(x = "Rank of Bayesian V_ij", y = "Rank of Frequentist P-value") +
  theme_minimal() +
  scale_color_manual(name = "",
                     labels = c("No Difference", "True Difference"),
                     values = c("FALSE" = "red", "TRUE" = "dodgerblue"))
sim_vij_pval_order_graph

# COMPARE FDR AND FNR TO TRUTH

mean(sim_decisions$ij_decisions$difference_truth)
sum(sim_decisions$ij_decisions$difference_truth)
with(sim_decisions, {
  optim_params$optim_e_vij[ij_decisions$decision == T & ij_decisions$difference_truth == F]
})

beta_sim <- sim_decisions$sim$beta_sim
beta_diff <- ComputeSimSTDDifferences(beta_sim, M1, ij_list = ij_list_s)
eps_seq <- c(seq(0, max(beta_diff) / 2, length.out = 200), optim_e)
V <- ComputeSimVij(beta_diff, epsilon = eps_seq)
t_seq <- seq(0, 1, length = 1000)
FDRFNR_curves <- ComputeFDRFNRCurves(V, t_seq, eps_seq)
FDRFNR_curves_s <- FDRFNR_curves[with(FDRFNR_curves,{epsilon %in%
    c(eps_seq[c(3, 15, 50, 80, 100)], optim_e)}), ]
FDRFNR_curves_s <- within(FDRFNR_curves_s, { epsilon <- factor(round(epsilon, digits = 3))})
FDRt_graph <- ggplot() +
  geom_line(data = FDRFNR_curves_s, aes(x = t, y = FDR, group = epsilon,
                                        color = epsilon)) +
  labs(x = "t", y = "Bayesian FDR") +
  theme_minimal()
FNRt_graph <- ggplot() +
  geom_line(data = FDRFNR_curves_s, aes(x = t, y = FNR, group = epsilon,
                                        color = epsilon)) +
  labs(x = "t", y = "Bayesian FNR") +
  theme_minimal()
FDRFNR_eps_graph <- ggplot() +
  geom_line(data = FDRFNR_curves_s, aes(x = FNR, y = FDR, group = epsilon,
                                        color = epsilon)) +
  labs(x = "Bayesian FNR", y = "Bayesian FDR") +
  theme_minimal()
FDRt_graph
FNRt_graph
FDRFNR_eps_graph

# SAVE GRAPHS
ggsave(file.path(getwd(), "output", "FE_sim", "eps_loss.png"), eps_loss_graph,
       width = 7, height = 5.2, units = "in")
ggsave(file.path(getwd(), "output", "FE_sim", "FDRFNR_graph.png"), FDRFNR_graph,
       width = 7, height = 5.2, units = "in")
ggsave(file.path(getwd(), "output", "FE_sim", "FDRt_graph.png"), FDRt_graph,
       width = 7, height = 5.2, units = "in")
ggsave(file.path(getwd(), "output", "FE_sim", "FNRt_graph.png"), FNRt_graph,
       width = 7, height = 5.2, units = "in")
ggsave(file.path(getwd(), "output", "FE_sim", "FDRFNR_eps_graph.png"), FDRFNR_eps_graph,
       width = 7, height = 5.2, units = "in")
ggsave(file.path(getwd(), "output", "FE_sim", "vij_pval_order.png"),
       width = 6, height = 4.5, units = "in",
       sim_vij_pval_order_graph)

dev.off()
par(mfrow = c(1, 1))
x <- pROC::roc(sim_decisions$ij_decisions$difference_truth, bayes_vij)
png(filename = file.path(getwd(), "output", "FE_sim", "roc.png"))
plot(x, legacy.axes = T, print.auc = TRUE, col = "dodgerblue", add = FALSE)
dev.off()
