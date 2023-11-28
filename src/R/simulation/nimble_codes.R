library(nimble)

BYM_code <- nimbleCode({
  for (i in 1:N) {
    Y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- inprod(X[i, 1:p], beta[1:p]) + phi[i]
    yFit[i] ~ dnorm(mu[i], tau2) ## Posterior predictive model fit
  }
  phi_prec[1:N, 1:N] <- tau2Sp * Q[1:N, 1:N]
  phi[1:N] ~ dmnorm(mu_phi[1:N], phi_prec[1:N, 1:N])
  for (i in 1:p) {
    beta[i] ~ dflat()
  }
  tau2 ~ dgamma(a0_tau2, b0_tau2)
  sigma <- 1 / sqrt(tau2)
  #sigmaSp <- 1 / sqrt(tau2Sp)
  tau2Sp <- 1 / sigmaSp^2
  sigmaSp <- sigma * delta2
  delta2 <- alpha / (1 - alpha)
  alpha ~ dbeta(a0_alpha, b0_alpha)
})

BYM2_code <- nimbleCode({
  for (i in 1:N) {
    for (j in 1:K) {
      Y[i, j] ~ dnorm(mu[i], tau2NSp)
      Yfit[i, j] ~ dnorm(mu[i], tau2NSp)
    }
    mu[i] <- inprod(X[i, 1:p], beta[1:p]) + phi[i] * sigmaSp
  }
  phi[1:N] ~ dmnorm(mu_phi[1:N], Q_scaled[1:N, 1:N])
  for (i in 1:p) {
    beta[i] ~ dflat()
  }
  tau2NSp <- 1 / ((1 - rho) * sigma^2)
  tau2Sp <- 1 / (rho * sigma^2)
  sigmaSp <- sqrt(rho) * sigma
  sigma2Sp <- rho * sigma^2
  sigma2NSp <- (1 - rho) * sigma^2
  delta2 <- rho / (1 - rho)
  alpha <- delta2 / (1 + delta2)
  sigma <- 1 / sqrt(tau2)
  tau2 ~ dgamma(a0_tau2, b0_tau2)
  rho ~ dbeta(0.5, 0.5) # following Riebler (2016)
})
