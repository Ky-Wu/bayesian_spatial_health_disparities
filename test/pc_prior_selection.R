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
rho_KLD <- vapply(rho_seq, function(target_rho) {
  KLD(target_rho, Lambda)
}, numeric(1))
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
alphas <- vapply(l_seq, function(target_l) {
  mean(samplePCPrior(10000, Lambda, target_l) <= U)
}, numeric(1))
min(l_seq[alphas >= alpha])

# US lung mortality RDA selection of lambda_rho
# Read in and setup lung + smoking data
source(file.path(getwd(), "src", "R", "RDA", "US_data_setup.R"))
Q_eigen <- eigen(Q_scaled)
Lambda <- Q_eigen$values
l_seq <- seq(0.02, 0.05, by = .0015)
alpha <- 2/3
U <- 0.5
alphas <- vapply(l_seq, function(target_l) {
  mean(samplePCPrior(10000, Lambda, target_l) <= U)
}, numeric(1))
min(l_seq[alphas >= alpha])
