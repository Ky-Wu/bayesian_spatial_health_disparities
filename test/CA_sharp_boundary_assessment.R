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

max_vijs <- all_vij_df[, .(ARDP_max = max(ARDP_vij), e_max = max(e_vij)),
                       keyby = .(sim_i)]
mean(max_vijs$ARDP_max <= 0.75)
mean(max_vijs$e_max <= 0.75)
FDR_estimate <- function(v, t, e = 0.5) {
  v_s <- v[v > t]
  sum(1 - v_s) / (length(v_s) + e)
}
computeNumPositives <- function(v, delta) {
  v_s <- sort(v, decreasing = FALSE)
  FDRs <- vapply(v_s, function(x) FDR_estimate(v, x, e = 0.01), numeric(1))
  tstar <- v_s[which.max(FDRs <= delta)]
  sum(v_s > tstar)
}
n_decisions <- all_vij_df[, .(ARDP_boundaries = computeNumPositives(ARDP_vij, delta = 0.15),
                              e_boundaries = computeNumPositives(e_vij, delta = 0.15)),
                          keyby = .(sim_i)]

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
  theme(legend.position = "bottom") +
  labs(x = "False Positive Rate", y = "True Positive Rate")

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
#roc.test(e_roc, ARDP_roc)

diff_prob <- melt(all_vij_df, measure.vars = c("e_vij", "ARDP_vij"))
diff_prob[, variable := ifelse(variable == "e_vij", "epsilon-difference", "ARDP-DAGAR")]
diff_prob_graph <- ggplot(data = diff_prob) +
  geom_violin(aes(x = variable, y = value), fill = "grey") +
  geom_jitter(aes(x = variable, y = value, color = true_diff),
              alpha = 0.15, width = 0.2, height = 0) +
  scale_color_discrete(name = "Boundary", labels = c("Not True Difference", "True Difference")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Difference Probability", x = "Difference Probability") +
  theme_bw() +
  theme(legend.position = "bottom")

decision_df <- melt(n_decisions, measure.vars = c("ARDP_boundaries", "e_boundaries"))
decision_graph <- ggplot(data = decision_df) +
  geom_histogram(aes(x =value, group = variable, fill = variable),
                 position = "identity", alpha = 0.7) +
  theme_bw() +
  labs(x = "Number of Reported Disparities", y = "Simulated Datasets") +
  scale_fill_discrete(name = "Difference Probability",
                      labels = c("ARDP-DAGAR", "epsilon-difference")) +
  theme(legend.position = "bottom")
decision_graph



### less smooth scenario: rho = 0.6 ###


outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim2")
all_vij_df <- fread(file.path(outputdir, "all_vij_df.csv"))
k <- nrow(all_vij_df[sim_i == 1])
T_edge <- seq_len(k)
sensSpec_df <- data.table()

max_vijs <- all_vij_df[, .(ARDP_max = max(ARDP_vij), e_max = max(e_vij)),
                       keyby = .(sim_i)]
mean(max_vijs$ARDP_max <= 0.75)
mean(max_vijs$e_max <= 0.75)
n_decisions <- all_vij_df[, .(ARDP_boundaries = computeNumPositives(ARDP_vij, delta = 0.15),
                              e_boundaries = computeNumPositives(e_vij, delta = 0.15)),
                          keyby = .(sim_i)]

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

roc_plot2 <- ggplot(data = roc_df, aes(x = 1 - specificity, y = sensitivity,
                                      group = Method, color = Method, linetype = Method)) +
  geom_line(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "grey") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "False Positive Rate", y = "True Positive Rate")

e_roc <- pROC::roc(sample(c(0, 1), k, replace = T), runif(k))
e_roc$sensitivities[1:k] <- rev(roc_df[Method == "e-difference",]$sensitivity)
e_roc$specificities[1:k] <- rev(roc_df[Method == "e-difference",]$specificity)
e_roc$auc <- auc(e_roc)
ARDP_roc <- pROC::roc(sample(c(0, 1), k, replace = T), runif(k))
ARDP_roc$sensitivities[1:k] <- rev(roc_df[Method == "ARDP-DAGAR",]$sensitivity)
ARDP_roc$specificities[1:k] <- rev(roc_df[Method == "ARDP-DAGAR",]$specificity)
ARDP_roc$auc <- auc(ARDP_roc)

diff_prob <- melt(all_vij_df, measure.vars = c("e_vij", "ARDP_vij"))
diff_prob[, variable := ifelse(variable == "e_vij", "epsilon-difference", "ARDP-DAGAR")]
diff_prob_graph2 <- ggplot(data = diff_prob) +
  geom_violin(aes(x = variable, y = value), fill = "grey") +
  geom_jitter(aes(x = variable, y = value, color = true_diff),
              alpha = 0.15, width = 0.2, height = 0) +
  scale_color_discrete(name = "Boundary", labels = c("Not True Difference", "True Difference")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Difference Probability", x = "Difference Probability") +
  theme_bw() +
  theme(legend.position = "bottom")

decision_df <- melt(n_decisions, measure.vars = c("ARDP_boundaries", "e_boundaries"))
decision_graph2 <- ggplot(data = decision_df) +
  geom_histogram(aes(x =value, group = variable, fill = variable),
                 position = "identity", alpha = 0.7) +
  theme_bw() +
  labs(x = "Number of Reported Disparities", y = "Simulated Datasets") +
  scale_fill_discrete(name = "Difference Probability",
                      labels = c("ARDP-DAGAR", "epsilon-difference")) +
  theme(legend.position = "bottom")


KLD <- function(rho, Lambda) {
  d <- -N / 2
  d <- d + sum(rho / Lambda + (1 - rho)) / 2
  d <- d - sum(log(rho / Lambda + (1 - rho))) / 2
}
rhoPCPrior <- function(rho, Lambda, l) {
  d <- sqrt(2 * KLD(rho, Lambda))
  l * exp(-l * d)
}
f <- function(xseq) vapply(xseq, function(x) rhoPCPrior(x, Lambda, chosen_lambda), numeric(1))
norm <- integrate(f, 0, 1)
Q_eigen <- eigen(Q_scaled)
Lambda <- Q_eigen$values
xseq <- seq(0, 1, length.out = 5000)
chosen_lambda <- 0.2
dseq <- vapply(xseq, function(x) rhoPCPrior(x, Lambda, chosen_lambda), numeric(1))

pc_prior <- ggplot() +
  geom_line(aes(x = xseq, y = dseq / norm$value), color = "dodgerblue") +
  coord_cartesian(ylim = c(0, 2.2)) +
  theme_bw() +
  labs(x = "rho", y = "Density")
setup_graph <- ggarrange(spatial_map, pc_prior, ncol = 2)

performance <- ggarrange(roc_plot, diff_prob_graph, decision_graph, ncol = 3)
performance2 <- ggarrange(roc_plot2, diff_prob_graph2, decision_graph2, ncol = 3)
roc_graphs <- ggarrange(roc_plot, roc_plot2, ncol = 2)
map_roc <- ggpubr::ggarrange(spatial_map, diff_prob_graph, roc_plot, ncol = 3)
map_roc

outputdir <- file.path(getwd(), "output", "CA_sharp_boundary_sim")
ggsave(file.path(outputdir, "setup_graph.png"), setup_graph,
       width = 8, height = 6, dpi = 500)
ggsave(file.path(file.path(getwd(), "output", "CA_sharp_boundary_sim"),
                 "performance.png"), performance,
       width = 12, height = 6, dpi = 500)
ggsave(file.path(file.path(getwd(), "output", "CA_sharp_boundary_sim2"),
                 "performance2.png"), performance2,
       width = 12, height = 6, dpi = 500)
p1 <- ggarrange(roc_plot, roc_plot2, ncol = 1, common.legend = TRUE)
p2 <- ggarrange(diff_prob_graph, diff_prob_graph2, ncol = 1, common.legend = TRUE)
p3 <- ggarrange(decision_graph, decision_graph2, ncol = 1, common.legend = TRUE)
all_performance <- ggarrange(p1, p2, p3, ncol = 3)
ggsave(file.path(file.path(getwd(), "output", "CA_sharp_boundary_sim"),
                 "all_performance.png"), all_performance,
       width = 12, height = 8, dpi = 500)
