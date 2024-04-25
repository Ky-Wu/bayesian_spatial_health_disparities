SimulateBetaDifferences <- function(X, beta, sigma2, n_sim,
                                    Y = NULL, ij_list = NULL) {
  n <- nrow(X)
  p <- length(beta)
  if (p != ncol(X)) stop("X and beta non-conformable")
  # Simulate Data, V_y = I_n
  if (is.null(Y)) {
    Y <- rnorm(n, X %*% beta, sd = sqrt(sigma2))
  }
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
  # simulate from posterior
  sigma2_sim <- 1 / rgamma(n_sim, a1, rate = b1)
  L <- t(chol(M1))
  beta_sim <- vapply(seq_len(n_sim), function(i) {
    z <- rnorm(p)
    as.vector(L %*% (sqrt(sigma2_sim[i]) * z + t(L) %*% m1))
  }, numeric(p))
  if (is.null(ij_list)) {
    ij_list <- expand.grid(i = seq_len(p), j = seq_len(p))
    ij_list <- ij_list[ij_list$i < ij_list$j,]
  }
  sim_std_differences <- ComputeSimSTDDifferences(beta_sim, M1, ij_list, sigma2_sim)
  out <- list(beta_sim = beta_sim,
              Y_sim = Y,
              sigma2_sim = sigma2_sim,
              std_differences = sim_std_differences,
              ij_list = ij_list,
              beta_truth = beta,
              sigma2_truth = sigma2,
              X = X)
  out
}

SimDecisions <- function(X, beta, sigma2, loss_function, eta, n_sim,
                         t_seq_length = 100, Y = NULL, ij_list = NULL) {
  # Inputs:
  #   X: numeric nxp
  #   beta: numeric px1
  #   sigma2: numeric(1), simulate data Y ~ N(XB, sigma2 * I_n)
  #   loss_function: function(V, epsilon)
  #   eta: numeric(1), compute optimal t cutoff by controlling FDR <= eta
  #   n_sim: numeric(1)
  n <- nrow(X)
  p <- length(beta)
  # Simulate data Y, compute beta posterior and differences
  sim <- SimulateBetaDifferences(X, beta, sigma2, n_sim = n_sim,
                                 Y = Y, ij_list = ij_list)
  # Compute v_ij = P(|beta_i - beta_j| > epsilon | y)
  eps_optim <- optim(median(sim$std_differences), function(e) {
    e_vij <- ComputeSimVij(sim$std_differences, epsilon = e)
    loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0, upper = max(abs(sim$std_differences)))
  optim_e <- eps_optim$par
  optim_e_vij <- ComputeSimVij(sim$std_differences, epsilon = eps_optim$par)
  # select cutoff t in d(i, j) = I(v_ij > t) to control FDR and minimize FNR
  t_seq <- seq(0, max(optim_e_vij) - .001, length.out = t_seq_length)
  t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t), numeric(1))
  t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t), numeric(1))
  optim_t <- min(c(t_seq[t_FDR <= eta], 1))
  ij_decisions <- sim$ij_list
  ij_decisions <- within(ij_decisions, {
    beta_i <- beta[i]
    beta_j <- beta[j]
    difference_truth = beta[i] != beta[j]
    decision = as.vector(optim_e_vij > optim_t)
  })
  out <- list(
    ij_decisions = ij_decisions,
    optim_params = list(
      optim_e = optim_e,
      optim_e_vij = optim_e_vij,
      optim_t = optim_t,
      optim_t_FDR = t_FDR[t_seq == optim_t],
      optim_t_FNR = t_FNR[t_seq == optim_t]
    ),
    sim = sim,
    eps_optim = eps_optim
  )
  out
}

ComputeFrequentistPValues <- function(X, Y, ij_list) {
  N <- nrow(X)
  p <- ncol(X)
  XtX <- t(X) %*% X
  XtX_chol <- chol(XtX)
  XtX_inv <- chol2inv(XtX_chol)
  beta_hat <- XtX_inv %*% (t(X) %*% Y)
  r <- Y - X %*% beta_hat
  sigma2_hat <- t(r) %*% r / (N - p)
  p_values <- vapply(seq_len(nrow(ij_list)), function(r) {
    c <- numeric(p)
    target_i <- ij_list[r,]$i
    target_j <- ij_list[r,]$j
    c[target_i] <- 1
    c[target_j] <- -1
    ctXtXc_chol <- t(chol(t(c) %*% XtX_inv %*% c))
    V <- solve(ctXtXc_chol, t(c) %*% beta_hat)
    F_statistic <- t(V) %*% V / sigma2_hat
    pf(F_statistic, df1 = 1, df2 = N - p, lower.tail = FALSE)
  }, numeric(1))
  p_values
}

ComputeClassificationMetrics <- function(predict, truth, e = 0L) {
  # predict: k logical vector
  # truth: k logical vector
  # e: numeric(1), small number to avoid zero denominator
  TP <- sum(predict & truth)
  FP <- sum(predict & !truth)
  TN <- sum(!predict & !truth)
  FN <- sum(!predict & truth)
  c("fdr" = FP / (TP + FP + e),
    "fnr" = FN / (FN + TN + e),
    "sensitivity" = TP / (TP + FN + e),
    "specificity" = TN / (TN + FP + e),
    "F1" = TP / (TP + 0.5 * (FP + FN)),
    "accuracy" = (TP + TN) / length(predict))
}
