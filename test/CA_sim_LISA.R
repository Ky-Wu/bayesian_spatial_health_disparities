library(data.table)
library(ggplot2)
library(ggpubr)
library(pROC)
library(rgeoda)
rm(list = ls())

outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim")

# Helper functions to compute posterior probabilities v_ij
source(file.path(getwd(), "src", "R", "simulation", "simulation_helper.R"))
source(file.path(getwd(), "src", "R", "kernel_matrix_construction.R"))

source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_setup.R"))

N <- nrow(W)
n_sim <- 100
current_i <- 1
w_rook <- rook_weights(county_sf)

all_ij_list <- data.frame()
for (sim_i in seq(current_i, n_sim)) {
  print(paste0("Running simulation ", sim_i, "/", n_sim, "..."))
  set.seed(sim_i)
  source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_data.R"))
  r <- resid(lm(y ~ X[,2]))
  lisa_rook <- local_moran(w_rook, as.data.frame(r), permutations = 50000,
                           significance_cutoff = 0.15)
  pvals <- lisa_pvalues(lisa_rook)
  pvals_adj <- p.adjust(pvals, method = "fdr")
  lisa_colors_rook <- lisa_colors(lisa_rook)
  lisa_labels_rook <- lisa_labels(lisa_rook)
  lisa_clusters_rook <- lisa_clusters(lisa_rook, cutoff = max(c(0, pvals[pvals_adj <= 0.15])))
  plot(st_geometry(county_sf),
       col=sapply(lisa_clusters_rook, function(x){return(lisa_colors_rook[[x+1]])}),
       border = "#333333", lwd=0.2)
  title(main = "LISA (Rook)")
  legend('bottomleft', legend = lisa_labels_rook, fill = lisa_colors_rook, border = "#eeeeee")
  region_labels <- vapply(lisa_clusters_rook,
                          function(x) lisa_labels_rook[x + 1],
                          character(1))
  ij_list_sim <- within(ij_list, {
    sim_i <- rep(sim_i, nrow(ij_list))
    pair_indx = seq_len(nrow(ij_list))
    i_label <- sapply(i, function(indx) region_labels[indx])
    j_label <- sapply(j, function(indx) region_labels[indx])
    LISA_diff_boundary <- (i_label %in% c("Low-High", "High-Low") |
                      j_label %in% c("Low-High", "High-Low")) & i_label != j_label
  })
  all_ij_list <- rbind(all_ij_list, ij_list_sim)
}

write.csv(all_ij_list, file = file.path(outputdir, "LISA_results.csv"))

### less smooth scenario: rho = 0.7 ###
rho <- 0.70
sigma2Sp <- rho * sigma2
sigma2NSp <- (1 - rho) * sigma2
all_ij_list2 <- data.frame()
outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim2")

for (sim_i in seq(current_i, n_sim)) {
  print(paste0("Running simulation ", sim_i, "/", n_sim, "..."))
  set.seed(sim_i)
  source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_data.R"))
  r <- resid(lm(y ~ X[,2]))
  lisa_rook <- local_moran(w_rook, as.data.frame(r), permutations = 50000,
                           significance_cutoff = 0.15)
  pvals <- lisa_pvalues(lisa_rook)
  pvals_adj <- p.adjust(pvals, method = "fdr")
  lisa_colors_rook <- lisa_colors(lisa_rook)
  lisa_labels_rook <- lisa_labels(lisa_rook)
  lisa_clusters_rook <- lisa_clusters(lisa_rook, cutoff = max(c(0, pvals[pvals_adj <= 0.15])))
  region_labels <- vapply(lisa_clusters_rook,
                          function(x) lisa_labels_rook[x + 1],
                          character(1))
  ij_list_sim <- within(ij_list, {
    sim_i <- rep(sim_i, nrow(ij_list))
    pair_indx = seq_len(nrow(ij_list))
    i_label <- sapply(i, function(indx) region_labels[indx])
    j_label <- sapply(j, function(indx) region_labels[indx])
    LISA_diff_boundary <- (i_label %in% c("Low-High", "High-Low") |
                             j_label %in% c("Low-High", "High-Low")) & i_label != j_label
  })
  all_ij_list2 <- rbind(all_ij_list2, ij_list_sim)
}

write.csv(all_ij_list2, file = file.path(outputdir, "LISA_results.csv"))
