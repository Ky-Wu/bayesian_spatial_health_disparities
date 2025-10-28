library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(ggpubr)
library(rjags)
library(R2jags)
library(Matrix)
#library(rstan)
library(pROC)
rm(list = ls())

outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim")
# Data generation

# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
source(file.path(getwd(), "src", "R", "kernel_matrix_construction.R"))
# Gibbs sampler
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_flatbeta_MCMC.cpp"))

source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_setup.R"))
N <- nrow(W)
# Priors
a_sigma <- 0.1
b_sigma <- 0.1
# PC prior, ignore IG prior parameters
a_rho <- 0.0
b_rho <- 0.0
# limits on possible values of rho
lower_rho <- 0
upper_rho <- 0.99
# PC prior: pi(rho) = lambda * exp(-lambda * d(rho)), d(rho) = sqrt(2 * KLD(rho))
lambda_rho <- .2

T_edge = seq(50, 100, by = 5)
n_sim <- 50
all_vij_df <- data.table()
indx <- seq_len(N)
current_i <- 1
if(file.exists(file.path(outputdir, "all_vij_df.csv"))) {
  all_vij_df <- fread(file.path(outputdir, "all_vij_df.csv"))
  current_i <- max(all_vij_df$sim_i) + 1
}
for (sim_i in seq(current_i, n_sim)) {
  print(paste0("Running simulation ", sim_i, "/", n_sim, "..."))
  set.seed(sim_i)
  source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_data.R"))
  # epsilon-difference probabilities
  our_time <- system.time({
    BYM2Sampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled)
    BYM2Sampler$setPriors(a_sigma, b_sigma, rho_prior_type = "pc",
                          a_rho, b_rho,
                          lower_rho, upper_rho, lambda_rho)
    BYM2Sampler$initOLS()
    BYM2Sampler$burnMCMCSample(40000)
    # fit BYM2 model
    samps <- BYM2Sampler$MCMCSample(20000)
    gamma_sim <- samps$gamma
    rho_sim <- samps$rho
    sigma2_sim <- samps$sigma2
    phi_sim <- apply(gamma_sim, MARGIN = 2, function(x) {
      x / sqrt(sigma2_sim * rho_sim)
    })
    phi_diffs <- BYM2_StdDiff(phi_sim, rho_sim, Q_scaled, X, ij_list)
    # Maximize entropy of posterior distribution with respect to epsilon
    loss_function <- function(v, epsilon) -ConditionalEntropy(v)

    eps_optim <- optim(1, function(e) {
      e_vij <- ComputeSimVij(phi_diffs, epsilon = e)
      loss_function(e_vij, epsilon = e)
    }, method = "Brent", lower = 0.0001, upper = 10.0)
    optim_e <- eps_optim$par
    optim_e_vij <- ComputeSimVij(phi_diffs, epsilon = optim_e)
  })

  # compute ARDP-DAGAR difference probabilities
  ARDP_time <- system.time({
    model.data1 <- list(k = n_county,  I = diag(ncol(X)),
                        X = X, Y = y,
                        Ik = diag(nrow(X)), index1 = index1,
                        alpha = 1, ncolumn = ncol(X), H = 15, ns = dni,
                        udnei = udnei, cn = c(0, cni))
    model.inits <- rep(list(list(rho1 = 0.1, tau1 = 1, taue1 = 5,
                                 taus = 1 / 5, beta = rep(0, ncol(X)))), 4)
    model.param <- c("beta", "rho1", "tau1", "vare1", "taus", "phi", "u")
    run.DAGAR1 <- jags(model.data1, model.inits, model.param,
                       file.path(outputdir, "DAGAR_ARDP.txt"),
                       n.chains = 2,
                       n.iter = 60000,
                       n.burnin = 20000,
                       n.thin = 1)
    name <- paste0("sim", sim_i, "DAGAR_ARDP_samps.RData")
    save(run.DAGAR1, file = file.path(outputdir, "ARDP_DAGAR_samps", name))
    gammas <-  run.DAGAR1$BUGSoutput$sims.matrix[,4:61]
    #county_sf$ARDP_gamma_pmeans <- apply(gammas, 2, mean)
    #plot(county_sf)
    # estimate difference boundaries
    vij_samples <- vapply(seq_len(nrow(ij_list)), function(pair_indx) {
      i <- ij_list[pair_indx,]$i
      j <- ij_list[pair_indx,]$j
      gammas[,i] != gammas[,j]
    }, numeric(nrow(gammas)))
    ARDP_vij <- apply(vij_samples, 2, mean)
  })

  all_vij_df <- rbind(all_vij_df, data.table(
    sim_i = sim_i,
    pair_indx = seq_len(nrow(ij_list)),
    true_diff = true_diff,
    e_vij = as.vector(optim_e_vij),
    ARDP_vij = ARDP_vij,
    analysis_time = our_time["elapsed"],
    ARDP_time = ARDP_time["elapsed"]
  ))
  colnames(all_vij_df) <- c("sim_i", "pair_indx", "true_diff", "e_vij", "ARDP_vij",
                            "analysis_time", "ARDP_time")
  fwrite(all_vij_df, file.path(outputdir, "all_vij_df.csv"))
}

