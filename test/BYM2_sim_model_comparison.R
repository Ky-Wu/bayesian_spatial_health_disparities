library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
rm(list = ls())
set.seed(122)

# Data generation
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))
# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "eps_loss_FDR.R"))
source(file.path(getwd(), "src", "R", "vij_computation.R"))
# Exact sampling and Gibbs sampling
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2ExactSampling.cpp"))
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "BYM2_flatbeta_MCMC.cpp"))
Rcpp::sourceCpp(file.path(getwd(), "src", "rcpp", "NGRegression.cpp"))

# Set priors
a_0 <- 0.1
b_0 <- 0.1
M_0inv <- diag(rep(1e-4, p))
M_0inv <- diag(rep(1e-12), p)
m_0 <- rep(0, p)

candidate_rhos <- c(seq(0.1, 0.9, by = 0.1), rho, .99)
n_sim <- 1000
ICs <- lapply(candidate_rhos, function(target_rho) {
  BYM2Sampler <- new(BYM2ExactSampler, X, y, Q_scaled, target_rho)
  BYM2Sampler$SetPriors(M_0inv, m_0, a_0, b_0)
  exact_samps <- BYM2Sampler$ExactSample(n_sim)
  data.table(rho = target_rho, DIC = exact_samps$DIC, WAIC = exact_samps$WAIC,
             nsp_var_mean = (1 - target_rho) * mean(exact_samps$sigma2),
             p_DIC = exact_samps$p_DIC, p_WAIC2 = exact_samps$p_WAIC2,
             pred_D = exact_samps$pred_D, pred_G = exact_samps$pred_G,
             pred_P = exact_samps$pred_G)
})
ICs <- do.call(rbind, ICs)
setkey(ICs, rho)
ICs

# simpler fits
ols_fit <- new(NGRegression, X, y, diag(N))
ols_fit$SetPriors(M_0inv, m_0, a_0, b_0)
exact_samps <- ols_fit$ExactSample(n_sim)
ICs <- rbind(ICs, data.table(rho = 0, DIC = exact_samps$DIC, WAIC = exact_samps$WAIC,
                             nsp_var_mean = mean(exact_samps$sigma2),
                             p_DIC = exact_samps$p_DIC, p_WAIC2 = exact_samps$p_WAIC2,
                             pred_D = exact_samps$pred_D, pred_G = exact_samps$pred_G,
                             pred_P = exact_samps$pred_G))
wls_fit <- new(NGRegression, X, y, Q_scaled)
wls_fit$SetPriors(M_0inv, m_0, a_0, b_0)
exact_samps <- wls_fit$ExactSample(n_sim)
ICs <- rbind(ICs, data.table(rho = 1, DIC = exact_samps$DIC, WAIC = exact_samps$WAIC,
                             nsp_var_mean = 0,
                             p_DIC = exact_samps$p_DIC, p_WAIC2 = exact_samps$p_WAIC2,
                             pred_D = exact_samps$pred_D, pred_G = exact_samps$pred_G,
                             pred_P = exact_samps$pred_G))
setkey(ICs, rho)
ICs

rm(wls_fit, ols_fit)

# Gibbs sampling: unknown rho
a_rho <- .1
b_rho <- .1
lower_rho <- 0.00
upper_rho <- 1.0
BYM2Sampler <- new(BYM2FlatBetaMCMC, X, y, Q_scaled)
BYM2Sampler$setPriors(a_0, b_0, a_rho, b_rho, lower_rho, upper_rho)
BYM2Sampler$initRandom()

for (j in seq_len(100)) {
  print(paste0(j, "/100"))
  print(paste0("current rho: ", BYM2Sampler$rho))
  print(paste0("current sigma2: ", BYM2Sampler$sigma2))
  BYM2Sampler$burnMCMCSample(100)
}
# DRAW POSTERIOR SAMPLES
n_sim <- 20000
system.time({
  samps <- BYM2Sampler$MCMCSample(n_sim)
})

par(mfrow = c(1, 2))
cap <- paste0("Simulation DIC (true rho = ", rho, ")")
plot(DIC ~ rho, data = ICs[ICs$rho != "gibbs",], pch = 16,
     main = cap)
abline(h = ICs[ICs$rho == "gibbs",]$DIC)
legend("bottomleft", pch = c(16, NA), lty = c(NA, 1),
       legend = c("rho fixed", "rho ~ Unif(0, 1)"))
cap <- paste0("Simulation WAIC (true rho = ", rho, ")")
plot(WAIC ~ rho, data = ICs[ICs$rho != "gibbs",], pch = 16,
     main = cap)
abline(h = ICs[ICs$rho == "gibbs",]$WAIC)
legend("bottomleft", pch = c(16, NA), lty = c(NA, 1),
       legend = c("rho fixed", "rho ~ Unif(0, 1)"))

ICs <- rbind(ICs, data.table(rho = "gibbs", DIC = samps$DIC, WAIC = samps$WAIC,
                             nsp_var_mean = mean((1 - samps$rho) * samps$sigma2),
                             p_DIC = samps$p_DIC, p_WAIC2 = samps$p_WAIC2,
                             pred_D = samps$pred_D, pred_G = samps$pred_G,
                             pred_P = samps$pred_G))



