library(Matrix)

scaled_CARQ <- function(D_W, W, O, G, alpha) {
  # input:
  # D_W: vector of rowsums of W
  # W: adjacency matrix
  # O: simultaneous reduction matrix, O^{\T}WO = G, O^{\T}D_W O = I_n
  # G: nx1 vector, eigenvalues of D_W^{-1/2}WD_W^{-1/2}
  # alpha: CAR spatial smoothing parameter
  Sigma <- O %*% diag(1 / (1 - alpha * G)) %*% t(O)
  scaling_factor <- exp(mean(log(diag(Sigma))))
  Q_scaled <- scaling_factor * alpha * -W
  diag(Q_scaled) <- diag(Q_scaled) + D_W * scaling_factor
  Q_scaled
}

scaled_SARQ <- function(P, Pinv, L, alpha) {
  # input:
  # P: eigenvectors of W_tilde
  # L: eigenvalues of W_tilde
  # alpha: SAR spatial autocorrelation parameter: B = alpha * \tilde{W}
  Q_scaled <- P %*% diag(1 - alpha * L) %*% Pinv
  Q_scaled <- Q_scaled %*% t(Q_scaled)
  Sigma <- t(Pinv) %*% diag(1 / (1 - alpha * L)) %*% t(P)
  Sigma <- Sigma %*% t(Sigma)
  scaling_factor <- exp(mean(log(diag(Sigma))))
  Q_scaled <- Q_scaled * scaling_factor
  Q_scaled
}

DAGAR_Q_ordered <- function(W, alpha) {
  # Inputs:
  # D_W: degree matrix matrix (ordered with W)
  # W: adjacency matrix
  # alpha: spatial autocorrelation parameter
  B <- as(W, "dgCMatrix")
  B[upper.tri(B)] <- 0
  n_i <- rowSums(B)
  nonzero_entries <- which(B != 0, arr.ind = TRUE)
  B[B != 0] <- alpha / (1 + (n_i[nonzero_entries[,1]] - 1) * alpha^2)
  L <- (1 + (n_i - 1) * alpha^2) / (1 - alpha^2)
  I_B <- -B
  diag(I_B) <- diag(I_B) + 1
  Q <- t(I_B) %*% diag(L) %*% I_B
  Q <- as(Q, "dgCMatrix")
  list(Q = Q, I_B = I_B, L = L)
}

scaled_orderedDAGARQ <- function(W, alpha) {
  # Inputs:
  # D_W: degree matrix matrix (ordered with W)
  # W: adjacency matrix
  # alpha: spatial autocorrelation parameter
  DAGAR_Q <- DAGAR_Q_ordered(W, alpha)
  I_Binv <- as(forwardsolve(DAGAR_Q$I_B, diag(N)), "dgCMatrix")
  Sigma <- I_Binv %*% diag(1 / DAGAR_Q$L) %*% t(I_Binv)
  scaling_factor <- exp(mean(log(diag(Sigma))))
  Q_scaled <- scaling_factor * DAGAR_Q$Q
  Q_scaled
}

DAGAR_f <- function(alpha, n_j) {
  r <- seq_len(n_j)
  sum(r / (1 + (r - 1) * alpha^2))
}
DAGAR_g <- function(alpha, n_k) {
  1 / (2 * (n_k - 1)) - 1 / ((n_k - 1) * n_k * (n_k + 1)) * DAGAR_f(alpha, n_k)
}
DAGAR_Q_order_free <- function(W, adj2_nbmat, D_W = rowSums(W), alpha) {
  # Inputs:
  # W: adjacency matrix with zero diagonal
  # adj2_nbmat: 2nd order neighbors matrix
  # D_W: nx1 vector with elements equal to number of neighbors of each region
  # alpha: spatial autocorrelation parameter
  N <- length(D_W)
  Q <- matrix(0, N, N)
  DAGAR_fg_table <- data.frame(n_j = seq(0, max(D_W)))
  DAGAR_fg_table$f <- vapply(DAGAR_fg_table$n_j, function(n_i) DAGAR_f(alpha, n_i),
                             numeric(1))
  DAGAR_fg_table$g <- with(DAGAR_fg_table, {
    1 / (2 * (n_j - 1)) - 1 / ((n_j - 1) * n_j * (n_j + 1)) * f
  })
  diag(Q) <- vapply(seq_len(N), function(i) {
    n_i <- D_W[i]
    nb_indx <- W[i,] != 0
    n_js <- D_W[nb_indx]
    Q_i <- vapply(n_js, function(n_k) {
      with(DAGAR_fg_table, f[n_j == n_k]) / (n_k * (n_k + 1))
    }, numeric(1))
    Q_i <- (alpha^2) / (1 - alpha^2) * sum(Q_i)
    Q_i <- 1 + n_i * alpha^2 / (2 * (1 - alpha^2)) + Q_i
    Q_i
  }, numeric(1))
  for (i in seq_len(N - 1)) {
    for(j in seq(i + 1, N)) {
      Q_ij <- ifelse(W[i, j] != 0, -alpha / (1 - alpha^2), 0)
      if (adj2_nbmat[i, j] != 0) {
        nb_intersection <- (W[i,] != 0) & (W[j,] != 0)
        k_indx <- rowSums(as.matrix(W[,nb_intersection])) != 0
        n_ks <- D_W[k_indx]
        nk_sum <- vapply(n_ks, function(n_k) {
          if (n_k == 1) {
            out <- 0
          } else {
            out <- with(DAGAR_fg_table, g[n_j == n_k])
          }
          out
        }, numeric(1))
        Q_ij <- Q_ij + sum(nk_sum) / (1 - alpha^2)
      }
      Q[i, j] <- Q_ij
      Q[j, i] <- Q_ij
    }
  }
  Q
}

