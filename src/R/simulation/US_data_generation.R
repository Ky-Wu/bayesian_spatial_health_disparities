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
W <- nb2mat(county_nbs, style = "B")
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


# ARDP-DAGAR setup

## Adjacency matrix
county.nbs = poly2nb(county_sp)
n=length(county.nbs)
county.coords <- coordinates(county_sp)
Adj=sapply(county.nbs,function(x,n) {v=rep(0,n);v[x]=1;v},n)
#colnames(Adj)=county_id

num_edge <- sum(Adj)/2

## Reorder the map
ca.latrange=round(quantile(county.coords[,2],c(0.25,0.75)))
ca.albersproj=mapproject(county.coords[,1],county.coords[,2],projection = "albers",param=ca.latrange)

perm=order(ca.albersproj$x-ca.albersproj$y)
colnames(Adj)[perm]

Adj_new=Adj[perm,perm]

n=nrow(Adj_new)
ni=rowSums(Adj_new)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Adj_new[x,]==1))
#N(i): 2:n
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n
dni=sapply(dneighbors,length)
original_perm = 1:n
index2=c(1,which(dni==0)+1)

final_perm=c(original_perm[perm][index2],
             original_perm[perm][-index2])
final_perm[order(final_perm)]

Minc = Adj[final_perm,final_perm]
n=nrow(Minc)
ni=rowSums(Minc)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Minc[x,]==1))
#N(i): 2:n
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n

dni=sapply(dneighbors,length)
nmax=max(dni)
cni=cumsum(dni)
dneimat=sapply(dneighbors, function(nei,nmax,n) c(nei,rep(n+1,nmax+1-length(nei))),nmax,n)
udnei=unlist(dneighbors)

ni_wo = sapply(neighbors,length)
cni_wo = cumsum(ni_wo)
udnei_wo = unlist(neighbors)
cn = c(0, cni)
ns = dni

region = seq(1:n)
index = list()
for(i in 1:(n-2)){
  index[[i]] = region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 = unlist(index)
mns = max(dni) + 1
