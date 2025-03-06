names(phi_d) <- NULL
p <- length(beta)
X <- cbind(1, matrix(rnorm(N * (p - 1), 0, 1), ncol = p - 1))
mu <- X %*% beta
y <- rnorm(N, mean = as.vector(mu) + sqrt(sigma2Sp) * phi_d,
           sd = sqrt(sigma2NSp))
county_sf$y <- y
county_sf$phi_d <- phi_d
county_sf$phi <- phi
plot(county_sf)

