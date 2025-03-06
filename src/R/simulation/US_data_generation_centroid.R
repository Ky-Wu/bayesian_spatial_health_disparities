# This file simulates data over Californian counties under a linear mixed
# effects model according to the BYM2 parameterization in Riebler (2016)

library(data.table)
library(sf)
library(spdep)
library(maps)
library(maptools)
library(magrittr)
library(stringr)
set.seed(122)
# Import US counties
county_poly <- maps::map("county", fill = TRUE, plot = FALSE)
county_state <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[1]]))
county_names <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[2]]))
sf_use_s2(FALSE)
county_sp <- maptools::map2SpatialPolygons(county_poly, IDs = county_poly$names)
county_nbs <- poly2nb(county_sp)
no_neighbors <- vapply(county_nbs, function(x) identical(x, 0L), logical(1))
# restrict to connected county map
county_sp <- county_sp[!no_neighbors,]
county_state <- county_state[!no_neighbors]
county_names <- county_names[!no_neighbors]
county_nbs <- poly2nb(county_sp)

county_cent <- st_centroid(st_as_sf(county_sp))
dist_matrix <- st_distance(county_cent)
Sigma <- Matern(dist_matrix)
Sigma_chol <- chol(Sigma)
Q <- chol2inv(Sigma_chol)
N <- nrow(Q)

sigma2Sp <- rho * sigma2
sigma2NSp <- (1 - rho) * sigma2
# Generate random effects
z <- rnorm(N, 0, 1)
phi <- Sigma_chol %*% z
beta <- c(2, 5)
p <- length(beta)
X <- cbind(1, matrix(rnorm(N * (p - 1), 0, 1), ncol = p - 1))
mu <- X %*% beta
# K replicates in each region
y <- rnorm(N, mean = mu + sqrt(sigma2Sp) * phi, sd = sqrt(sigma2NSp))

ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL

county_sf <- st_as_sf(county_sp)
rownames(county_sf) <- NULL
st_crs(county_sf) <- st_crs(st_as_sf(county_poly))
