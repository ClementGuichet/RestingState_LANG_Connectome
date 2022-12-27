##########################################################################################
# Script for Age-related analyses and visualization for the topologico-functional profile

# Written by CG
# 13-12-2022
##########################################################################################
library(dplyr)
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(rgl)
library(fitdistrplus)
library(sigclust2)


rm(list = ls())

source("_05_Clustering.R")
source("_radarplotting_function.R")

################################################################################
# Post-clustering analysis -----------------
################################################################################
data_post_clustering %>% 
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full")


data_post_clustering <- data_post_clustering %>%
  mutate(Age_group = ifelse(cluster == 2, "Young", ifelse(cluster == 3, "Old", 0)))

data_post_clustering %>%
  group_by(cluster) %>%
  count(Gender) %>% 
  mutate(n = prop.table(n))


a <- gghistogram(data_post_clustering %>% subset(Age_group == "Young"),
  x = "Age", y = "..density..", bins = 20,
  fill = "purple", add_density = TRUE
) + theme_pubclean()
b <- gghistogram(data_post_clustering %>% subset(Age_group == "Old"),
  x = "Age", y = "..density..", bins = 20,
  fill = "purple", add_density = TRUE
) + theme_pubclean()


Rmisc::multiplot(a, b)

# Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
data_cluster_selection <- function(cluster1, cluster2) {
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(Age_group == cluster1 | Age_group == cluster2)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID)
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  tmp_cluster_final <<- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, Age_group),
                             by = "Subj_ID"
  )
}

data_cluster_selection("Young", "Old")


################################################################################
# Box plot ---------------------------------------------------------------------
################################################################################

data_box <- data_post_clustering %>%
  pivot_longer(
    cols = !c("Subj_ID", "Gender", "Age", "cluster", "Age_group"),
    names_to = "Metrics",
    values_to = "Metric_value"
  )


data_box$Metrics <- factor(data_box$Metrics, levels = c(
  "Connector", "Satellite", "Provincial", "Peripheral",
  "Local_Bridge", "Global_Bridge", "Super_Bridge", "Not_a_Bridge", "Balance_eff"
))

# T-test after clustering is not very appropriate given that data is not drawn from a random
# distribution anymore, clustering has maximized the intra-cluster variability
# Furthermore, clustering already indicates that cluster profiles differ from one another
#
# data_box_sig <- data_box %>%
#   group_by(Metrics) %>%
#   rstatix::t_test(Metric_value ~ Median_Age, ref.group = "56") %>%
#   adjust_pvalue(method = "fdr") %>%
#   add_significance("p.adj") %>%
#   add_xy_position(x = "Metrics")

ggplot(data_box, aes(x = Metrics, y = Metric_value)) +
  geom_boxplot(aes(fill = Age_group)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_pubr()
# +
#   stat_pvalue_manual(data_box_sig,
#     label = "p.adj.signif",
#     tip.length = 0.01,
#     hide.ns = TRUE,
#     bracket.nudge.y = -5
#   )

data_box_eff_size <- data_box %>%
  group_by(Metrics) %>%
  rstatix::cohens_d(Metric_value ~ Age_group, comparison = list(c("Young", "Old")), paired = FALSE, hedges.correction = TRUE) %>%
  mutate(effsize = effsize * (-1))
# filter(magnitude != "negligible")

ggdotchart(
  data_box_eff_size,
  x = "Metrics", y = "effsize",
  ylab = "Cohen's d (hedge corrected)",
  palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "descending",
  rotate = TRUE, legend = "none"
)

################################################################################
# Topologico-functional profile across clusters  -------------------------------
################################################################################

# Log ratio of the percentage of each metric of a cluster to the geometric mean of all individuals
# equivalent to CLR-transform, preserves unit-sum constraint and removes value-range restriction

geometric_all <- data_post_clustering %>%
  summarize_at(vars(Connector:Super_Bridge), funs(compositions::geometricmean(.)))

Radar_functional_role_geometric_age <- data_post_clustering %>%
  dplyr::select(-Age) %>%
  group_by(Age_group) %>%
  summarize_at(vars(Connector:Super_Bridge), funs(compositions::geometricmean(.))) %>%
  mutate(Connector = log(Connector / geometric_all$Connector)) %>%
  mutate(Provincial = log(Provincial / geometric_all$Provincial)) %>%
  mutate(Satellite = log(Satellite / geometric_all$Satellite)) %>%
  mutate(Peripheral = log(Peripheral / geometric_all$Peripheral)) %>%
  mutate(Global_Bridge = log(Global_Bridge / geometric_all$Global_Bridge)) %>%
  mutate(Local_Bridge = log(Local_Bridge / geometric_all$Local_Bridge)) %>%
  mutate(Super_Bridge = log(Super_Bridge / geometric_all$Super_Bridge)) %>%
  mutate(Not_a_Bridge = log(Not_a_Bridge / geometric_all$Not_a_Bridge)) %>% 
  remove_rownames() %>%
  column_to_rownames(var = "Age_group")

radarplotting_overlap(Radar_functional_role_geometric_age, 1, -1, 1, 1,
  alpha = 0.5, label_size = 1,
  title_fill = "Topologico-functional profile of each cluster",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "topright", title = "Median age of each cluster\n Based on Model-based clustering\n with principal balances",
  legend = rownames(Radar_functional_role_geometric_age), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)

Radar_functional_role_age <- data_post_clustering %>%
  filter(Age_group != "Middle") %>%
  dplyr::select(-Age) %>%
  group_by(Age_group) %>%
  summarize_at(vars(Connector:Super_Bridge), funs(mean(.))) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Age_group")

radarplotting_overlap(Radar_functional_role_age, 50, 0, 1, 1,
  alpha = 0.05, label_size = 1,
  title_fill = "Topologico-functional profile of each cluster",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "topright", title = "Median age of each cluster\n Based on ILR-transformed Wald Hierarchichal clustering",
  legend = rownames(Radar_functional_role_geometric_age), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)


################################################################################
# QUALITY CHECK -------------
################################################################################
# Distribution of hubs across RSNs for each cluster for the individual hubs ----

# First averaging per subject then per clusters because grand mean is not equal to mean of means
# with unequal sample size i.e., subjects have a different total number of hubs
Radar_hub_RSN <- tmp_cluster_final %>%
  group_by(Subj_ID, Age_group, `1st_network`) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  group_by(Age_group, `1st_network`) %>%
  summarize_at(vars(freq), mean) %>%
  spread(`1st_network`, freq) %>%
  remove_rownames() %>%
  column_to_rownames("Age_group") %>%
  mutate_at(vars(everything()), funs(. * 100))

radarplotting_overlap(Radar_hub_RSN, 30, 0, 1, 1,
  alpha = 0.3, label_size = 1,
  title_fill = "Distribution of hubs regions across RSNs",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "topright", title = "Median age of each Age_group\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
  legend = rownames(Radar_hub_RSN), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)

# Distribution of hubs across communities for each Age_group for the individual hubs ----

Radar_hub_community <- tmp_cluster_final %>%
  group_by(Subj_ID, Age_group, Consensus_vector_0.15) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  group_by(Age_group, Consensus_vector_0.15) %>%
  summarize_at(vars(freq), mean) %>%
  spread(Consensus_vector_0.15, freq) %>%
  remove_rownames() %>%
  column_to_rownames("Age_group") %>%
  mutate_at(vars(everything()), funs(. * 100))

radarplotting_overlap(Radar_hub_community, 50, 0, 1, 1,
  alpha = 0.7, label_size = 1,
  title_fill = "Distribution of hubs regions across communities",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "topright", title = "Median age of each Age_group\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
  legend = rownames(Radar_hub_community), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)



