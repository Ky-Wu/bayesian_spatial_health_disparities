DrawBYM2Samples <- function(X, Y, Sigma_scaled, a0_sigma, b0_sigma, node1, node2,
                            n_samp = 10000, n_chains = 5, n_burn = n_samp / 2) {
  N <- nrow(X)
  p <- ncol(X)
  stan_data <- list(Y = Y,
                    X = X,
                    N = N,
                    Sigma_scaled = Sigma_scaled,
                    mu_phi = rep(0, N),
                    N_edges = nrow(nbs_df),
                    p = p,
                    a0_sigma = a0_sigma,
                    b0_sigma = b0_sigma,
                    node1 = nbs_df$node1,
                    node2 = nbs_df$node2)
  inits <- rep(list(list(rho = 0.5)), n_chains)
  rstan_fit <- stan(file.path(getwd(), "src", "stan", "bym2.stan"),
                    data = stan_data, warmup = n_burn, iter = n_samp + n_burn,
                    chains = n_chains, init = inits,
                    pars = c("beta" , "sigma2", "rho", "phi", "YFit"))
  samps <- as.matrix(rstan_fit, pars = "lp__", include = FALSE)
  list(rho_sim = samps[, "rho"],
       phi_sim = samps[, grepl("phi", colnames(samps))],
       yfit_sim = samps[, grepl("YFit", colnames(samps))],
       beta_sim = samps[, grepl('beta', colnames(samps))],
       sigma2_sim = samps[, "sigma2"])
}




