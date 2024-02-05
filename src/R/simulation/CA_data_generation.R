library(data.table)
library(sf)
library(spdep)
library(maps)
library(maptools)
library(magrittr)

# Import CA counties
county_poly <- maps::map("county","california", fill=TRUE, plot=FALSE)
county_state <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[1]]))
county_names <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[2]]))
sf_use_s2(FALSE)
county_sp <- maptools::map2SpatialPolygons(county_poly, IDs = county_names)
county_nbs <- poly2nb(county_sp)
W <- nb2mat(county_nbs, style="B")
D <- diag(rowSums(W))
alpha <- 0.99
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
scaling_factor <- exp(mean(log(diag(Sigma))))
Q_scaled <- Q * scaling_factor
Q_scaled_cholR <- Q_cholR * sqrt(scaling_factor)
Sigma_scaled <- Sigma / scaling_factor
N <- nrow(W)
# rho = proportion of overall variance attributed to spatial variance
# sigma2 = overall variance
sigma2 <- 4
rho <- 0.93
sigma2Sp <- rho * sigma2
sigma2NSp <- (1 - rho) * sigma2
# Generate random effects
z <- rnorm(N, 0, 1)
phi <- solve(Q_scaled_cholR, z)
beta <- c(2, 5)
p <- length(beta)
X <- cbind(1, matrix(rnorm(N * (p - 1), 0, 1), ncol = p - 1))
mu <- X %*% beta
# K replicates in each region
K <- 1
Y <- replicate(K, rnorm(N, mean = mu + sqrt(sigma2Sp) * phi, sd = sqrt(sigma2NSp)))
ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL
