library(maps)
library(spdep)
library(maptools)
library(mapproj)
library(data.table)
library(sf)
library(spdep)
library(magrittr)
library(stringr)
library(fields)
set.seed(122)

# Import US counties
ca.county = map("county","california", fill=TRUE, plot=FALSE)
county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords = coordinates(ca.poly)
n_county <- length(county.ID)

## Adjacency matrix
ca.neighbors = poly2nb(ca.poly)
n=length(ca.neighbors)

Adj=sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
#colnames(Adj)=county_id

num_edge <- sum(Adj)/2

## Reorder the map
ca.latrange=round(quantile(ca.coords[,2],c(0.25,0.75)))
ca.albersproj=mapproject(ca.coords[,1],ca.coords[,2],projection = "albers",param=ca.latrange)

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
original_perm = 1:58
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

# Construct CAR spatial kernel for analysis model
W <- Adj[final_perm,final_perm]
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

# data generation spatial variance: Matern covariance
county_sp <- ca.poly[final_perm,]
county_nbs <- poly2nb(county_sp)
county_cent <- st_centroid(st_as_sf(county_sp))
county_sf <- st_as_sf(county_sp)
sf_use_s2(TRUE)
rownames(county_sf) <- NULL
st_crs(county_sf) <- st_crs(st_as_sf(ca.county))
st_crs(county_cent) <- st_crs(st_as_sf(ca.county))
dist_matrix <- matrix(st_distance(county_cent), nrow = nrow(county_sf),
                      ncol = nrow(county_sf)) / 1000
#dist_matrix <- st_distance(county_cent[1,], county_cent[2,])
Sigma <- Matern(dist_matrix, range = 0.5 * 300, phi = 1, smoothness = 0.5, nu = 0.5)
Sigma_chol <- chol(Sigma)

# rho = proportion of overall variance attributed to spatial variance
# sigma2 = overall variance
beta <- c(2, 5)
sigma2 <- 5
rho <- 0.90
sigma2Sp <- rho * sigma2
sigma2NSp <- (1 - rho) * sigma2
# Generate random effects

ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL
# number of cuts on phi to create sharp boundaries
n_cut <- 5

z <- rnorm(N, 0, 1)
phi <- t(Sigma_chol) %*% z
#phi <- solve(Q_scaled_cholR, z)
#phi_d <- ceiling(phi)
phi_cut <- cut_number(phi, n_cut)
phi_means <- tapply(phi, phi_cut, mean)
phi_d <- phi_means[phi_cut]
true_diff <- vapply(seq_len(nrow(ij_list)), function(pair_indx) {
  i <- ij_list[pair_indx,]$i
  j <- ij_list[pair_indx,]$j
  abs(phi_d[i] - phi_d[j]) >= 0.001
}, logical(1))
ij_list$true_diff <- true_diff
print(paste0("Number of true difference boundaries: ", sum(true_diff)))
