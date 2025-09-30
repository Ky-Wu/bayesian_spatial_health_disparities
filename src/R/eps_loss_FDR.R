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

ConditionalEntropy <- function(V) {
  vapply(seq_len(ncol(V)), function(i) {
    v <- V[, i]
    entropy <- ifelse(v == 0 | v == 1, 0.0, -v * log(v) - (1 - v) * log(1 - v))
    sum(entropy)
  }, numeric(1))
}

FDR_estimate <- function(v, t, e = 0.5) {
  v_s <- v[v >= t]
  sum(1 - v_s) / (length(v_s) + e)
}

FNR_estimate <- function(v, t, e = 0.5) {
  v_s <- v[v < t]
  sum(v_s) / (length(v_s) + e)
}

FD_estimate <- function(v, t) {
  v_s <- v[v >= t]
  sum(1 - v_s)
}
FN_estimate <- function(v, t) {
  v_s <- v[v < t]
  sum(v_s)
}

computeFDRCutoff <- function(diff_prob, delta) {
  t_seq <- sort(unique(diff_prob), decreasing = FALSE)
  t_FDR <- vapply(t_seq, function(t) FDR_estimate(diff_prob, t, e = 0), numeric(1))
  optim_t <- min(c(t_seq[t_FDR <= delta], 1), na.rm = TRUE)
  list(cutoff = optim_t, FDR_estimate = FDR_estimate(diff_prob, t = optim_t, e = 0))
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

GGPlotEpsLoss <- function(sim_results, loss_function, eps_seq_length = 50) {
  differences <- sim_results$sim$std_differences
  eps_seq <- seq(0.001, max(differences) / 2, length.out = eps_seq_length)
  sim_vij <- ComputeSimVij(differences, epsilon = eps_seq)
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
