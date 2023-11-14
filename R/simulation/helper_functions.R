library(parallel)

ComputeSimSTDDifferences <- function(sim_values, C, ij_list) {
  k <- nrow(ij_list)
  n_sim <- ncol(sim_values)
  std_differences <- vapply(seq_len(k), function(pair_indx) {
    i <- ij_list[pair_indx,]$i
    j <- ij_list[pair_indx,]$j
    sd <- sqrt(C[i, i] + C[j, j] - 2 * C[i, j])
    (sim_values[i,] - sim_values[j,]) / sd
  }, numeric(n_sim))
  std_differences
}

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
  sim_std_differences <- ComputeSimSTDDifferences(beta_sim, M1, ij_list)
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

ComputeSimVij <- function(d, ij_list, epsilon) {
  # simulations indexed by row
  k <- nrow(ij_list)
  V <- vapply(epsilon, function(target_epsilon) {
    v <- apply(d, MARGIN = 2, function(col_d) {
      mean(abs(col_d) > target_epsilon)
    })
    v
  }, numeric(k))
  V
}
GetConstantFunction <- function(c, ...) {
  function(...) c
}

GeneralLoss <- function(V, epsilon, g, h) {
  if (length(epsilon) != ncol(V))  {
    stop("length of epsilon sequence does not match number of columns of V")
  }
  # number of pairs
  k <- nrow(V)
  loss <- vapply(seq_along(epsilon), function(i) {
    v <- V[, i]
    target_epsilon <- epsilon[i]
    # both reward and penalty scales with epsilon, want to balance risk and reward
    - sum(v) * g(target_epsilon) + h(target_epsilon) * sum(1 - v)
  }, numeric(1))
  loss
}

FDR_estimate <- function(v, t, e = 0.5) {
  v_s <- v[v > t]
  sum(1 - v_s) / (length(v_s) + e)
}

FNR_estimate <- function(v, t, e = 0.5) {
  v_s <- v[v <= t]
  sum(v_s) / (length(v_s) + e)
}

FD_estimate <- function(v, t) {
  v_s <- v[v > t]
  sum(1 - v_s)
}
FN_estimate <- function(v, t) {
  v_s <- v[v <= t]
  sum(v_s)
}

ComputeFDRFNRCurves <- function(V, t_seq, eps_seq, denom_e = 0.1) {
  # V: numeric m x r matrix, columns match v_ij for epsilon sequence eps_seq
  # t_seq: numeric(k) vector
  # eps_seq: numeric(r) vector
  out <- lapply(seq_len(ncol(V)), function(i) {
    v <- V[, i]
    FDR <- vapply(t_seq, function(t) FDR_estimate(v, t, e = denom_e), numeric(1))
    FNR <- vapply(t_seq, function(t) FNR_estimate(v, t, e = denom_e), numeric(1))
    sensitivity <-
    data.frame(epsilon = eps_seq[i],
               t = t_seq,
               FDR = FDR,
               FNR = FNR)
  })
  out <- do.call(rbind, out)
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
     e_vij <- ComputeSimVij(sim$std_differences, sim$ij_list, epsilon = e)
     loss_function(e_vij, epsilon = e)
  }, method = "Brent", lower = 0, upper = max(sim$std_differences))
  optim_e <- eps_optim$par
  optim_e_vij <- ComputeSimVij(sim$std_differences, sim$ij_list,
                               epsilon = eps_optim$par)
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

GGPlotEpsLoss <- function(sim_results, loss_function, eps_seq_length = 50) {
  differences <- sim_results$sim$std_differences
  eps_seq <- seq(0.001, max(differences) / 2, length.out = eps_seq_length)
  sim_vij <- ComputeSimVij(differences, sim_results$sim$ij_list, epsilon = eps_seq)
  sim_loss <- loss_function(sim_vij, epsilon = eps_seq)
  optim_e <- sim_results$optim_params$optim_e
  ggplot() +
    geom_line(data = data.frame(sim_loss = sim_loss, epsilon = eps_seq),
              aes(x = epsilon, y = sim_loss), color = "dodgerblue") +
    labs(x = "Epsilon", y = "loss", title = "",
         subtitle = paste0("Optimal Epsilon = ", round(optim_e, digits = 4))) +
    geom_vline(xintercept = optim_e, lwd = 0.8, linetype = "dotted",
               color = "red") +
    theme_minimal()
}

GGPlotFDRFNR <- function(sim_results, t_seq_length = 100) {
  optim_e_vij <- sim_results$optim_params$optim_e_vij
  t_seq <- seq(0, max(optim_e_vij) - .001,
               length.out = t_seq_length)
  t_FDR <- vapply(t_seq, function(t) FDR_estimate(optim_e_vij, t), numeric(1))
  t_FNR <- vapply(t_seq, function(t) FNR_estimate(optim_e_vij, t), numeric(1))
  data <- data.frame(t = rep(t_seq, 2),
                     statistic = rep(c("Bayesian FDR", "Bayesian FNR"), each = t_seq_length),
                     value = c(t_FDR, t_FNR))
  ggplot() +
    geom_line(data = data, aes(x = t, y = value), color = "dodgerblue") +
    geom_vline(xintercept = sim_results$optim_params$optim_t,
               color = "red", lwd = 0.8, linetype = "dotted") +
    facet_wrap(~statistic)
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


CARStdDiff <- function(sim_phi, alpha, delta2, D, W, X, M0_inv, ij_list,
                       mc.cores = 1) {
  XtX <- t(X) %*% X
  D_half <- sqrt(D)
  D_neghalf <- diag(diag(D^(-1/2)))
  Wstar <- D_neghalf %*% W %*% D_neghalf
  Wstar_eigen <- eigen(Wstar)
  P <- Wstar_eigen$vectors
  lambda <- Wstar_eigen$values
  N <- nrow(X)
  M0_inv_zero <- all(M0_inv == 0)
  if (M0_inv_zero) {
    H <- chol(XtX)
    H <- solve(t(H), t(X))
    H <- t(H) %*% H
  }
  std_diffs <- mcmapply(function(target_alpha, target_delta2, target_sim_phi_rindx) {
    target_sim_phi <- matrix(sim_phi[target_sim_phi_rindx,], ncol = 1)
     if (!M0_inv_zero) {
       H <- XtX + target_delta2 * M0_inv
       H <- chol(H)
       H <- solve(t(H), t(X))
       H <- t(H) %*% H
     }
     B <- diag(N) - H
     B_star <- D_half %*% t(P) %*% B %*% P %*% D_half / target_delta2
     A_inv <- diag(1 / (1 - target_alpha * lambda))
     phi_cov_core <- A_inv + B_star
     phi_cov_core <- chol(phi_cov_core)
     phi_cov <- solve(t(phi_cov_core), t(P) %*% D_half)
     phi_cov <- t(phi_cov) %*% phi_cov
     std_diffs <- ComputeSimSTDDifferences(target_sim_phi, phi_cov, ij_list)
  }, alpha, delta2, seq_len(nrow(sim_phi)), mc.cores = mc.cores)
  std_diffs
}
