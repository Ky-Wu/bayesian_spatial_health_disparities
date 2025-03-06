library(pROC)
library(data.table)
library(ggpubr)
rm(list = ls())

source(file.path(getwd(), "src", "R", "simulation", "CA_sharp_boundary_sim_setup.R"))
outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim")
all_vij_df <- fread(file.path(outputdir, "all_vij_df.csv"))
k <- nrow(all_vij_df[sim_i == 1])
T_edge <- seq_len(k)
sensSpec_df <- data.table()

computeSens <- function(vij, true_diff, target_T) {
  cutoff <- sort(vij, decreasing = TRUE)[target_T]
  pred_diff <- vij >= cutoff
  tp <- pred_diff & true_diff
  sum(tp) / sum(true_diff)
}
computeSpec <- function(vij, true_diff, target_T) {
  cutoff <- sort(vij, decreasing = TRUE)[target_T]
  pred_diff <- vij >= cutoff
  tn <- !pred_diff & !true_diff
  sum(tn) / sum(!true_diff)
}

for(target_t in T_edge) {
  new_row <- all_vij_df[, .(
    e_sens = computeSens(e_vij, true_diff, target_t),
    e_spec = computeSpec(e_vij, true_diff, target_t),
    ARDP_sens = computeSens(ARDP_vij, true_diff, target_t),
    ARDP_spec = computeSpec(ARDP_vij, true_diff, target_t)
  ), by = .(sim_i)][, .(
    n_edge = target_t,
    e_sens = mean(e_sens),
    e_spec = mean(e_spec),
    ARDP_sens = mean(ARDP_sens),
    ARDP_spec = mean(ARDP_spec)
  )]
  sensSpec_df <- rbind(sensSpec_df, new_row)
}

sens_df <- melt(sensSpec_df, id.vars = "n_edge",
                    measure.vars = c("e_sens", "ARDP_sens"),
                    variable.name = "Method", value.name = "sensitivity")
sens_df[, Method := ifelse(Method == "e_sens", "e-difference", "ARDP-DAGAR")]
spec_df <- melt(sensSpec_df, id.vars = "n_edge",
                measure.vars = c("e_spec", "ARDP_spec"),
                variable.name = "Method", value.name = "specificity")
spec_df[, Method := ifelse(Method == "e_spec", "e-difference", "ARDP-DAGAR")]
roc_df <- merge(sens_df, spec_df, by = c("n_edge", "Method"))

roc_plot <- ggplot(data = roc_df, aes(x = 1 - specificity, y = sensitivity,
                          group = Method, color = Method, linetype = Method)) +
  geom_line(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "grey") +
  theme_bw() +
  theme(legend.position = "bottom")

county_sf$`Spatial Effect` <- factor(round(phi_d, digits = 2), ordered = TRUE)
spatial_map <- ggplot(data = county_sf) +
  geom_sf(aes(fill = `Spatial Effect`)) +
  scale_fill_viridis_d() +
  coord_sf(crs = 5070) +
  theme_bw() +
  theme(legend.position = "bottom")

e_roc <- pROC::roc(sample(c(0, 1), k, replace = T), runif(k))
e_roc$sensitivities[1:k] <- rev(roc_df[Method == "e-difference",]$sensitivity)
e_roc$specificities[1:k] <- rev(roc_df[Method == "e-difference",]$specificity)
e_roc$auc <- auc(e_roc)
ARDP_roc <- pROC::roc(sample(c(0, 1), k, replace = T), runif(k))
ARDP_roc$sensitivities[1:k] <- rev(roc_df[Method == "ARDP-DAGAR",]$sensitivity)
ARDP_roc$specificities[1:k] <- rev(roc_df[Method == "ARDP-DAGAR",]$specificity)
ARDP_roc$auc <- auc(ARDP_roc)
roc.test(e_roc, ARDP_roc)

map_roc <- ggpubr::ggarrange(spatial_map, roc_plot, ncol = 2)
ggsave(file.path(outputdir, "map_roc.png"), map_roc,
       width = 10, height = 6, dpi = 500)
