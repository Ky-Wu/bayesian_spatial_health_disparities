library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
rm(list = ls())

KLD <- function(rho, Lambda) {
  d <- -N / 2
  d <- d + sum(rho / Lambda + (1 - rho)) / 2
  d <- d - sum(log(rho / Lambda + (1 - rho))) / 2
}
rhoPCPrior <- function(rho, Lambda, l) {
  d <- sqrt(2 * KLD(rho, Lambda))
  l * exp(-l * d)
}
# simulation selection of lambda_rho
source(file.path(getwd(), "src", "R", "simulation", "US_data_generation.R"))

Q_eigen <- eigen(Q_scaled)
Lambda <- Q_eigen$values
rho_seq <- seq(0.01, 1, length = 100)
rho_prior_values <- vapply(rho_seq, function(target_rho) {
  rhoPCPrior(target_rho, Lambda, l = .0325)},
  numeric(1))
plot(rho_prior_values ~ rho_seq, type = "l", ylim = c(0, max(rho_prior_values)),
     main = "PC Rho Prior (lamb = 0.325)")
samplePCPrior <- function(n, Lambda, l) {
  M <- rhoPCPrior(0, Lambda, l)
  out <- vapply(seq_len(n), function(i) {
    y <- runif(1, 0, 1)
    u <- runif(1)
    while(u >= rhoPCPrior(y, Lambda, l) / M) {
      y <- runif(1, 0, 1)
      u <- runif(1)
    }
    y
  }, numeric(1))
  out
}

set.seed(1130)
alpha <- 2/3
U <- 0.5
l_seq <- seq(0.02, 0.05, by = .0015)
# alpha = P(rho <= U)
alphas <- vapply(l_seq, function(target_l) {
  mean(samplePCPrior(10000, Lambda, target_l) <= U)
}, numeric(1))
# select lambda such that
min(l_seq[alphas >= alpha])

# US lung mortality RDA selection of lambda_rho
# Read in and setup lung + smoking data
source(file.path(getwd(), "src", "R", "RDA", "US_data_setup.R"))
set.seed(1130)
Q_eigen <- eigen(Q_scaled)
Lambda <- Q_eigen$values
alpha <- 2/3
U <- 0.5
l_seq <- seq(0.02, 0.05, by = .0015)
# alpha = P(rho <= U)
alphas <- vapply(l_seq, function(target_l) {
  mean(samplePCPrior(10000, Lambda, target_l) <= U)
}, numeric(1))
# select lambda such that p(rho <= U) ~= alpha
chosen_lambda <- min(l_seq[alphas >= alpha])
chosen_lambda

samps <- samplePCPrior(100000, Lambda, chosen_lambda)
quantile(samps, c(0.025, 0.975))

