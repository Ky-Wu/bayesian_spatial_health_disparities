library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(ggpubr)
library(rjags)
library(R2jags)
library(Matrix)
library(rstan)
library(pROC)
rm(list = ls())

outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim")
# Data generation

# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
source(file.path(getwd(), "src", "R", "kernel_matrix_construction.R"))


source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_setup.R"))
N <- nrow(W)
# Priors
a_sigma <- 0.1
b_sigma <- 0.1

T_edge = seq(50, 100, by = 5)
n_sim <- 50
all_vij_df <- data.table()
indx <- seq_len(N)
for (sim_i in seq_len(n_sim)) {
  print(paste0("Running simulation ", sim_i, "/", n_sim, "..."))
  set.seed(sim_i)
  source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_data.R"))
  # epsilon-difference probabilities
  data <- list(
    N = N,
    Sigma_chol = t(chol(Sigma_scaled)),
    #Sigma_scaled = Sigma_scaled,
    Y = y,
    p = ncol(X),
    X = X,
    a0_sigma = a_sigma,
    b0_sigma = b_sigma
  )

  # fit BYM2 model (rstan only for small data)
  fit1 <- stan(
    file = file.path(getwd(), "src", "stan", "bym2.stan"),
    pars = c("beta", "phi", "sigma2", "rho"),
    data = data,
    chains = 4,
    warmup = 40000,
    iter = 60000,
    cores = 4
  )
  samps <- as.matrix(fit1)
  phi_sim <- samps[, paste0("phi[", seq_len(N), "]")]
  rho_sim <- samps[, "rho"]
  sigma2_sim <- samps[, "sigma2"]
  phi_diffs <- BYM2_StdDiff(phi_sim, rho_sim, Q_scaled, X, ij_list)
  # Maximize entropy of posterior distribution with respect to epsilon
  loss_function <- function(v, epsilon) -ConditionalEntropy(v)
  system.time({
    eps_optim <- optim(1, function(e) {
      e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
      loss_function(e_vij, epsilon = e)
    }, method = "Brent", lower = 0.0001, upper = 10.0)
  })
  optim_e <- eps_optim$par
  optim_e_vij <- ComputeSimVij(phi_diffs, epsilon = optim_e)

  # compute ARDP-DAGAR difference probabilities
  model.data1 <- list(k = n_county,  I = diag(ncol(X)),
                      X = X, Y = y,
                      Ik = diag(nrow(X)), index1 = index1,
                      alpha = 1, ncolumn = ncol(X), H = 15, ns = dni,
                      udnei = udnei, cn = c(0, cni))
  model.inits <- rep(list(list(rho1 = 0.1, tau1 = 1, taue1 = 5,
                               taus = 1 / 5, beta = rep(0, ncol(X)))), 2)
  model.param <- c("beta", "rho1", "tau1", "vare1", "taus", "phi", "u")
  run.DAGAR1 <- jags(model.data1, model.inits, model.param,
                     file.path(outputdir, "DAGAR_ARDP.txt"),
                     n.chains = 2,
                     n.iter = 60000,
                     n.burnin = 40000,
                     n.thin = 1)
  gammas <-  run.DAGAR1$BUGSoutput$sims.matrix[,4:61]
  county_sf$ARDP_gamma_pmeans <- apply(gammas, 2, mean)
  plot(county_sf)
  # estimate difference boundaries
  vij_samples <- vapply(seq_len(nrow(ij_list)), function(pair_indx) {
    i <- ij_list[pair_indx,]$i
    j <- ij_list[pair_indx,]$j
    gammas[,i] != gammas[,j]
  }, numeric(nrow(gammas)))
  ARDP_vij <- apply(vij_samples, 2, mean)
  all_vij_df <- rbind(all_vij_df, data.table(
    sim_i = sim_i,
    pair_indx = seq_len(nrow(ij_list)),
    true_diff = true_diff,
    e_vij = optim_e_vij,
    ARDP_vij = ARDP_vij
  ))
  colnames(all_vij_df) <- c("sim_i", "pair_indx", "true_diff", "e_vij", "ARDP_vij")
  fwrite(all_vij_df, file.path(outputdir, "all_vij_df.csv"))
}

