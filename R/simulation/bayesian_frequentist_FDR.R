library(data.table)
library(ggplot2)
set.seed(1130)

rm(list = ls())
source(file.path(getwd(), "helper_functions.R"))

#### DOES CONTROLLING BAYESIAN FDR CONTROL TRUE FDR?
n_sim <- 5000
eta <- .05
set.seed(100)
beta <- rnorm(500, -1, sd = 2)
beta_truth <- beta >= 0
sigma2 <- 1
mean(beta_truth)
n <- 200
p <- length(beta)
X <- matrix(runif(n * p, -2, 2), n, p)
Y <- rnorm(n, X %*% beta, sd = sqrt(sigma2))
M0 <- diag(p) * 100
m0 <- rep(0, p)
a0 <- .01
b0 <- .01
# Compute parameters of posterior
a1 <- a0 + n/2
M1 <- chol2inv(chol(chol2inv(chol(M0)) + t(X) %*% X))
m1 <- m0 + t(X) %*% Y
b1 <- b0 + 0.5 * (t(Y) %*% Y + t(m0) %*% M0 %*% m0 - t(m1) %*% M1 %*% m1)
sigma2_sim <- 1 / rgamma(n_sim, a1, rate = b1)
L <- t(chol(M1))
beta_sim <- vapply(seq_len(n_sim), function(i) {
  z <- rnorm(p)
  as.vector(L %*% (sigma2_sim[i] * z + t(L) %*% m1))
}, numeric(p))
beta_vij <- apply(beta_sim, MARGIN = 1, function(x) mean(x >= 0))
t_seq <- seq(.0001, 1, length.out = 500)
t_FDR <- vapply(t_seq, function(target_t) FDR_estimate(beta_vij, target_t), numeric(1))
optimal_t <- min(c(t_seq[t_FDR <= eta], 1))
bayes_decision <- beta_vij >=  optimal_t
# CONTROLLING BAYESIAN FDR DOES NOT CONTROL FDR | beta IN GENERAL
ComputeClassificationMetrics(bayes_decision, beta_truth)
bayes_FDR <- FDR_estimate(beta_vij, t = optimal_t, e = 0)
bayes_FDR
bayes_FNR <- FNR_estimate(beta_vij, t = optimal_t, e = 0)
bayes_FNR
