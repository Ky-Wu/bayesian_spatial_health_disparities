ComputeSimSTDDifferences <- function(sim_values, C, ij_list, sigma2_sim = 1) {
  # sim indexed by column
  k <- nrow(ij_list)
  n_sim <- ncol(sim_values)
  sigma_sim <- sqrt(sigma2_sim)
  std_differences <- vapply(seq_len(k), function(pair_indx) {
    i <- ij_list[pair_indx,]$i
    j <- ij_list[pair_indx,]$j
    sd <- sqrt(C[i, i] + C[j, j] - 2 * C[i, j]) * sigma_sim
    (sim_values[i,] - sim_values[j,]) / sd
  }, numeric(n_sim))
  std_differences
}

ComputeSimScaledDiffs <- function(sim_values, ij_list) {
  k <- nrow(ij_list)
  n_sim <- nrow(sim_values)
  diffs <- vapply(seq_len(k), function(pair_indx) {
    i <- ij_list[pair_indx,]$i
    j <- ij_list[pair_indx,]$j
    sim_values[, i] - sim_values[, j]
  }, numeric(n_sim))
  apply(diffs, MARGIN = 2, function(x) abs(x) / sd(x))
}

ComputeSimVij <- function(d, epsilon) {
  # simulations indexed by row, v_ij indexed by column
  k <- ncol(d)
  V <- vapply(epsilon, function(target_epsilon) {
    over_threshold <- abs(d) > target_epsilon
    colMeans(over_threshold)
  }, numeric(k))
  V
}

BYM2_StdDiff <- function(sim_phi, sim_rho, Q, X, ij_list) {
  # sim indexed by row
  # assume Q positive definite, apply simultaneous diagonalization
  N <- nrow(X)
  k <- nrow(ij_list)

  XtX <- t(X) %*% X
  Q_eigen <- eigen(Q)
  if (any(Q_eigen$values <= 0)) stop("Q not positive definite")
  Q_neghalf <- Q_eigen$vectors %*% diag(Q_eigen$values^(-0.5)) %*% t(Q_eigen$vectors)
  H <- chol(XtX)
  H <- solve(t(H), t(X))
  H <- t(H) %*% H
  I_H <- diag(N) - H
  B <- Q_neghalf %*% I_H %*% Q_neghalf
  B_eigen <- eigen(B)
  D <- B_eigen$values
  O <- B_eigen$vectors
  U <- Q_neghalf %*% O
  n_sim <- nrow(sim_phi)

  U2_contrasts <- vapply(seq_len(k), function(pair_indx) {
    i <- ij_list[pair_indx,]$i
    j <- ij_list[pair_indx,]$j
    (U[i,] - U[j,])^2
  }, numeric(N))
  var_core <- vapply(sim_rho, function(target_rho) {
    1 / (1 + target_rho / (1 - target_rho) * D)
  }, numeric(N))
  sds <- sqrt(t(var_core) %*% U2_contrasts)
  diffs <- vapply(seq_len(k), function(pair_indx) {
    i <- ij_list[pair_indx,]$i
    j <- ij_list[pair_indx,]$j
    sim_phi[,i] - sim_phi[,j]
  }, numeric(n_sim))
  diffs / sds
}

VijOrderGraph <- function(phi_diffs, ij_list, optim_e,
                          other_es = round(optim_e * c(1 / 3, 2 / 3, 3 / 2, 3),
                                           digits = 3),
                          T_line = NULL) {
  optim_e_vij <- ComputeSimVij(phi_diffs, epsilon = optim_e)
  optim_e_vij_order <- order(optim_e_vij, decreasing = F)
  rejection_paths <- lapply(other_es, function(e) {
    e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
    e_vij_order <- order(e_vij[optim_e_vij_order], decreasing = F)
    rejection_path <- data.table(optim_e_vij = seq_along(optim_e_vij),
                                 e_vij_order = e_vij_order)
    rejection_path <- melt(rejection_path, id.vars = "optim_e_vij",
                           variable.name = "eps_value", value.name = "order")
    rejection_path[eps_value == "e_vij_order",
                   eps_value := paste0("eps = ", e)]
    rejection_path
  })
  rejection_paths <- do.call(rbind, rejection_paths)
  vij_order_graph <- ggplot() +
    geom_point(data = rejection_paths,
               aes(x = optim_e_vij, y = order),
               alpha = 0.3, size = 1) +
    facet_grid(~eps_value) +
    labs(x = paste0("Rank of v_ij(", round(optim_e, digits = 3), ")"),
         y = "Rank of v_ij(eps)") +
    theme_minimal()
  if (is.numeric(T_line)) {
    vij_order_graph <- vij_order_graph +
      geom_vline(xintercept = nrow(ij_list) - T_line)
  }
  vij_order_graph
}

ComputeFDREps <- function(phi_diffs, delta, T_edge, tol = 1e-7) {
  lower_eps <- 0.01
  upper_eps <- 5.0
  while((upper_eps - lower_eps) > tol) {
    eps <- (upper_eps + lower_eps) / 2
    e_vij <- ComputeSimVij(phi_diffs, epsilon = eps)
    FDR <- FDR_estimate(e_vij, sort(e_vij, decreasing = TRUE)[T_edge + 1])
    if (FDR > delta) {
      upper_eps <- eps
    } else {
      lower_eps <- eps
    }
  }
  list(eps = lower_eps,
       FDR = FDR_estimate(e_vij, lower_eps),
       e_vij = e_vij)
}

