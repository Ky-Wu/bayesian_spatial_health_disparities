KLD <- function(rho, Lambda) {
  # rho: numeric vector, spatial proportion of total variance
  # Lambda: numeric vector, eigenvalues of spatial covariance matrix
  # OUTPUT: Kullback-Leiber divergence between BYM2 model and OLS model
  d <- -N / 2
  d <- d + sum(rho / Lambda + (1 - rho)) / 2
  d <- d - sum(log(rho / Lambda + (1 - rho))) / 2
}
rhoPCPrior <- function(rho, Lambda, l) {
  # rho: numeric vector, spatial proportion of total variance
  # Lambda: numeric vector, eigenvalues of spatial covariance matrix
  # # l: positive numeric, hyperparameter
  d <- sqrt(2 * KLD(rho, Lambda))
  l * exp(-l * d)
}
samplePCPrior <- function(n, Lambda, l) {
  # n: integer, number of samples
  # Lambda: numeric vector, eigenvalues of spatial covariance matrix
  # l: positive numeric, hyperparameter
  # OUTPUT: n samples from PC(l)
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
