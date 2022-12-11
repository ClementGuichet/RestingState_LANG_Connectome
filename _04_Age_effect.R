##########################################################################################
# Script for Age-related analyses and visualization for the hubness profile

# Written by CG
# 26-11-2022
##########################################################################################
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(rgl)
library(FactoMineR)
library(factoextra)
library(FactoInvestigate)


rm(list = ls())
options(max.print = 99999)

source("_03_Hub_classification.R")
source("_radarplotting_function.R")

# What are the graph-based biomarkers of health aging ?
# Does the hubness profile, that is the proportion of functional role, changes throughout life ?
#
################################################################################
# Correlation between degree centrality and Age --------------------------------
# Inferential statistics
library(fitdistrplus)

data_stat_age <- data_functional_role %>%
  group_by(Subj_ID, CAB_NP_assign, Region, Age) %>%
  summarize_at(vars(degree), mean) %>%
  filter(Age != "NaN") %>%
  mutate(DC = degree / max(.$degree))


gghistogram(data_stat_age,
  x = "Age", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
gghistogram(data_stat_age,
  x = "DC", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()

fitdistrplus::descdist(data_stat_age$DC)

correlation_DC_age <- data_stat_age %>%
  group_by(CAB_NP_assign, Region) %>%
  group_split() %>%
  map_dfr(. %>%
    mutate(Estimate = cor.test(.$DC, .$Age, method = "kendall")$estimate) %>%
    mutate(p_value = cor.test(.$DC, .$Age, method = "kendall")$p.value)) %>%
  group_by(CAB_NP_assign, Region) %>%
  summarise_at(vars(Estimate, p_value), mean)

Radar_RSN_degree_age <- correlation_DC_age %>%
  group_by(CAB_NP_assign) %>%
  summarize_at(vars(Estimate), mean) %>%
  spread(CAB_NP_assign, Estimate)

radarplotting_unrounded_grad(Radar_RSN_degree_age, 0.2, -0.2, 1, 1,
  palette = "red", label_size = 1,
  alpha = 0.4, title_fill = "Mean correlation between degree centrality and Age for each RSN"
)

ggdotchart(
  correlation_DC_age %>% subset(p_value <= 0.05),
  x = "Region", y = "Estimate",
  ylab = "Mean correlation",
  palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "descending",
  rotate = TRUE, legend = "none", title = "Significant correlations between degree centrality and Age"
)


################################################################################
################################################################################
# ~~~~~~~~~~~ Hub Detection Procedure ~~~~~~~~~~~
################################################################################
# Method 1 ~ Detect top % regions for each metric ------------------------------

# Hubness profile computed with different hubs for 3 different age groups ------
top <- 131 * 0.2

Top_metric_Age_group_0 <- data_functional_role %>%
  filter(Age != "NaN") %>%
  mutate(Age_group = ifelse(Age < 39, "Young",
    ifelse(Age >= 39 & Age < 59, "Middle",
      ifelse(Age >= 59, "Old", "NaN")
    )
  )) %>%
  group_by(Age_group, Region) %>%
  summarize_at(vars(zK, Within_module_z_cons, zPC_cons, zBT, zFlow), mean) %>%
  # mutate(across(degree:PC, ~ rank(-.x), .names = "{.col}_rank")) %>%
  pivot_longer(
    cols = !c("Age_group", "Region"),
    names_to = "Metric_name",
    values_to = "Metric_value"
  ) %>%
  group_by(Age_group, Metric_name, .add = TRUE) %>%
  group_split() %>%
  map_dfr(. %>% slice_max(Metric_value, n = top) %>% distinct(Region, .keep_all = T) %>%
    mutate(rank = rep(seq(1:length(Region)))))

# Attach the corresponding RSN
get_RSN_label <- data_functional_role %>%
  filter(Region %in% Top_metric_Age_group_0$Region) %>%
  dplyr::select(CAB_NP_assign, Region) %>%
  distinct()

Top_metric_Age_group_1 <- merge(Top_metric_Age_group_0, get_RSN_label, by = "Region")

# Hubness profile
young_hub <- data_hub_selection_young %>%
  filter(Age < 39) %>%
  group_by(Subj_ID, Age, Hub_consensus) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Young") %>%
  plyr::rename(c("Hub_consensus" = "Functional_role"))

young_bridge <- data_hub_selection_young %>%
  filter(Age < 39) %>%
  group_by(Subj_ID, Age, Bridgeness) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Young") %>%
  plyr::rename(c("Bridgeness" = "Functional_role"))

middle_hub <- data_hub_selection_middle %>%
  filter(Age >= 39 & Age < 59) %>%
  group_by(Subj_ID, Age, Hub_consensus) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Middle") %>%
  plyr::rename(c("Hub_consensus" = "Functional_role"))

middle_bridge <- data_hub_selection_middle %>%
  filter(Age >= 39 & Age < 59) %>%
  group_by(Subj_ID, Age, Bridgeness) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Middle") %>%
  plyr::rename(c("Bridgeness" = "Functional_role"))

old_hub <- data_hub_selection_old %>%
  subset(Age >= 59) %>%
  group_by(Subj_ID, Age, Hub_consensus) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Old") %>%
  plyr::rename(c("Hub_consensus" = "Functional_role"))

old_bridge <- data_hub_selection_old %>%
  subset(Age >= 59) %>%
  group_by(Subj_ID, Age, Bridgeness) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Old") %>%
  plyr::rename(c("Bridgeness" = "Functional_role"))

data_hubness_profile_Age_group <- rbind(
  young_hub, young_bridge,
  middle_hub, middle_bridge,
  old_hub, old_bridge
)

data_hubness_profile_Age_group %>%
  subset(Functional_role != "None") %>%
  ggplot(aes(Age, freq, color = Functional_role)) +
  geom_point(size = 2, alpha = 0.2) +
  geom_jitter(height = 0.05, alpha = 0.2) +
  geom_smooth() +
  ggpubr::theme_pubr() +
  ggtitle("Evolution of functional roles across adult lifespan")

Radar_hubness_profile_Age_group <- data_hubness_profile_Age_group %>%
  subset(Functional_role != "None") %>%
  spread(Functional_role, freq) %>%
  dplyr::select(-c(Subj_ID, Age))

Radar_hubness_profile_Age_group$Age_group <- factor(Radar_hubness_profile_Age_group$Age_group, levels = c("Young", "Middle", "Old"))

Radar_hubness_profile_Age_group <- Radar_hubness_profile_Age_group %>%
  arrange(Age_group) %>%
  group_by(Age_group) %>%
  mutate_all(., ~ replace(., is.na(.), 0)) %>%
  summarize_at(vars(Connector:Super_Bridge), mean) %>%
  remove_rownames() %>%
  column_to_rownames("Age_group") %>%
  mutate_at(vars(everything()), funs(. * 100))

radarplotting_overlap(Radar_hubness_profile_Age_group, 50, 0, 1, 1,
  alpha = 0.2, label_size = 1,
  title_fill = "Distribution of functional roles for each Age group\n (after hub detection procedure)",
  palette = RColorBrewer::brewer.pal(3, "Dark2")
)

legend(
  x = "bottomleft", legend = rownames(Radar_hubness_profile_Age_group), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(3, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)

# Proportion of hub regions for each age group across RSN
data_hub_selection_young <- data_functional_role %>%
  filter(Region %in% (Top_metric_Age_group_0 %>% subset(Age_group == "Young"))$Region)
data_hub_selection_middle <- data_functional_role %>%
  filter(Region %in% (Top_metric_Age_group_0 %>% subset(Age_group == "Middle"))$Region)
data_hub_selection_old <- data_functional_role %>%
  filter(Region %in% (Top_metric_Age_group_0 %>% subset(Age_group == "Old"))$Region)


Radar_RSN_1 <- data_hub_selection_young %>%
  filter(Age < 39) %>%
  group_by(CAB_NP_assign) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Young")

Radar_RSN_2 <- data_hub_selection_middle %>%
  filter(Age >= 39 & Age < 59) %>%
  group_by(CAB_NP_assign) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Middle")

Radar_RSN_3 <- data_hub_selection_old %>%
  subset(Age >= 59) %>%
  group_by(CAB_NP_assign) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>%
  mutate(Age_group = "Old")

Radar_RSN_age <- rbind(Radar_RSN_1, Radar_RSN_2, Radar_RSN_3) %>%
  spread(CAB_NP_assign, freq)

Radar_RSN_age$Age_group <- factor(Radar_RSN_age$Age_group, levels = c("Young", "Middle", "Old"))

Radar_RSN_age <- Radar_RSN_age %>%
  arrange(Age_group) %>%
  remove_rownames() %>%
  column_to_rownames("Age_group") %>%
  mutate_at(vars(everything()), funs(. * 100))

radarplotting_overlap(Radar_RSN_age, 30, 0, 1, 1,
  alpha = 0.2, label_size = 1,
  title_fill = "Distribution of Hub regions for each Age group\n (after hub detection procedure)",
  palette = RColorBrewer::brewer.pal(3, "Dark2")
)

legend(
  x = "bottomleft", legend = rownames(Radar_RSN_age), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(3, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)


################################################################################
# Hubness profile with hub detection at the individual level -------------------

top <- 131 * 0.2

# LOG: Testing if clustering is stable when selecting another top% of regions
# Initial clustering at 20% yields 2 stable clusters and one who can be subdivided into 3 smaller components (young/old/middle with reduced var)
# Same results with 15% and 25%, ~ 70-80 regions for representativity
# When selecting top 1% of each GT metrics --> 6 regions in total per subject, --> maximizes effects for connector and provincial hubs between young and old

Top_metric_Age_ind <- data_functional_role %>%
  group_by(Subj_ID, Region) %>%
  summarize_at(vars(zK, Within_module_z_cons, zPC_cons, zBT, zFlow), mean) %>%
  # mutate(across(degree:PC, ~ rank(-.x), .names = "{.col}_rank")) %>%
  pivot_longer(
    cols = !c("Subj_ID", "Region"),
    names_to = "Metric_name",
    values_to = "Metric_value"
  ) %>%
  group_by(Subj_ID, Metric_name, .add = TRUE) %>%
  group_split() %>%
  map_dfr(. %>% slice_max(Metric_value, n = top) %>%
    mutate(rank = rep(seq(1:length(Region))))) %>%
  group_by(Subj_ID, .add = TRUE) %>%
  group_split()

Hub_selection <- list()
FR_list <- list()
for (i in 1:length(Top_metric_Age_ind)) {
  Hub_df <- rbindlist(Top_metric_Age_ind[i]) %>% distinct(Region, .keep_all = TRUE)
  # Here I subset the rows specific to each subject and their Hub regions
  tmp <- data_functional_role %>%
    filter(Region %in% Hub_df$Region) %>%
    filter(Subj_ID == i) %>%
    dplyr::select(Subj_ID, Region, Hub_consensus, Bridgeness)
  Hub_selection[[i]] <- tmp

  # Here I compute the proportion of functional roles regarding centrality and information flow
  # with the hubs specific to each subject
  FR_ind_hub <- tmp %>%
    group_by(Hub_consensus) %>%
    summarize(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Hub_consensus, freq)
  FR_ind_bridge <- tmp %>%
    group_by(Bridgeness) %>%
    summarize(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Bridgeness, freq)

  FR_ind <- cbind(FR_ind_hub, FR_ind_bridge)
  FR_list[[i]] <- FR_ind
}

# Adding a ratio score denoting propensity for segregation over integration
data_cluster_efficiency <- data_functional_role %>%
  dplyr::select(Subj_ID, Eglob, Eloc) %>%
  group_by(Subj_ID) %>%
  summarize_at(vars(Eglob, Eloc), mean) %>%
  mutate(Balance_eff = (Eloc - Eglob) / (Eloc + Eglob)) %>%
  dplyr::select(-c(Subj_ID, Eglob, Eloc)) %>%
  mutate_at(vars(Balance_eff), funs(. * 100))

# Putting everything together
data_hubness_profile_Age_ind <- cbind(
  rbindlist(FR_list, fill = TRUE) %>% dplyr::select(-None) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
  data_functional_role %>% group_by(Subj_ID) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID),
  Balance_eff = data_cluster_efficiency$Balance_eff
)

data_hubness_profile_Age_ind %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_point(size = 2, alpha = 0.2) +
  geom_jitter(height = 0.05, alpha = 0.2) +
  geom_smooth() +
  ggpubr::theme_pubr() +
  ggtitle("Evolution of functional roles across adult lifespan")

###############################################################################
# Getting ready for hierarchical clustering analysis
data_cluster <- data_hubness_profile_Age_ind %>%
  filter(Age != "NaN") %>%
  dplyr::select(-c(Subj_ID, Age))

# Clustering -----------------
performance::check_clusterstructure(data_cluster)
corrplot::corrplot(cor(data_cluster))

# Define linkage methods
# link <- c("average", "single", "complete", "ward")
# names(link) <- c("average", "single", "complete", "ward")
# # Function to compute the agglomerative coefficient (i.e., the strength of the clusters)
# ac <- function(x) {
#   cluster::agnes(data_cluster, method = x)$ac
# }
# sapply(link, ac)

# Defining optimal number of cluster with gap statistic
gap_stat <- cluster::clusGap(data_cluster,
  FUN = hcut,
  K.max = 5,
  B = 500,
  verbose = T
)

factoextra::fviz_gap_stat(gap_stat)

plot_cluster <- cluster::agnes(data_cluster, method = "ward", metric = "euclidian")
fviz_dend(plot_cluster,
  cex = 0.8,
  k = 4,
  palette = "jco",
  rect = T,
  color_labels_by_k = TRUE,
  main = "Ward Hierarchical clustering of Hub type proportion per age on regions defined as hubs previously\n (Connector, Provincial, Peripheral, Satellite)"
)


# Evaluate significance of hclust using 1000 Monte Carlo simulated null Gaussian
set.seed(2022)
shc_result <- sigclust2::shc(as.matrix(data_cluster),
  metric = "euclidian",
  linkage = "ward.D2",
  n_sim = 1000,
  alpha = 0.05
)

plot(shc_result, hang = .1)

# Compute distance matrix
d <- dist(data_cluster, method = "euclidean")
final_clust <- hclust(d, method = "ward.D2")
groups <- cutree(final_clust, k = 4)

data_post_clustering <- cbind(data_hubness_profile_Age_ind %>% filter(Age != "NaN"), cluster = groups)

################################################################################
# Post-clustering analysis -----------------
################################################################################
data_post_clustering %>%
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full")

# Young and old are more similar and stable than middle age

gghistogram(data_post_clustering %>% subset(cluster == "1"),
  x = "Age", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
gghistogram(data_post_clustering %>% subset(cluster == "2"),
  x = "Age", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
gghistogram(data_post_clustering %>% subset(cluster == "3"),
  x = "Age", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
gghistogram(data_post_clustering %>% subset(cluster == "4"),
  x = "Age", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()

ComplexHeatmap::Heatmap(scale(data_cluster), split = data_post_clustering$cluster)

data_post_clustering <- data_post_clustering %>%
  mutate(cluster = ifelse(cluster == "1", "23",
    ifelse(cluster == "2", "25",
      ifelse(cluster == "3", "35",
        ifelse(cluster == "4", "68", 0)
      )
    )
  ))

################################################################################
# PCA --------------------------------------------------------------------------
pc <- PCA(data_post_clustering %>% dplyr::select(-c(Subj_ID, Age, cluster)), axes = c(1, 2))

data_post_clustering$cluster <- factor(data_post_clustering$cluster)
fviz_pca_ind(pc,
  geom.ind = "point", pointshape = 21,
  axes = c(1, 2),
  pointsize = 2,
  fill.ind = data_post_clustering$cluster,
  col.ind = "black",
  palette = "jco",
  addEllipses = TRUE,
  label = "var",
  col.var = "black",
  repel = TRUE,
  legend.title = "Median age"
) +
  ggtitle("2D PCA-plot of functional roles") +
  theme(plot.title = element_text(hjust = 0.5))

################################################################################
# Hubness profile across clusters  ---------------------------------------------
################################################################################

Radar_functional_role_median_age <- data_post_clustering %>%
  dplyr::select(-Age) %>%
  group_by(cluster) %>%
  summarize_at(vars(Connector:Super_Bridge), mean) %>%
  remove_rownames() %>%
  column_to_rownames(var = "cluster")

radarplotting_overlap(Radar_functional_role_median_age, 60, 0, 1, 1,
  alpha = 0.05, label_size = 1,
  title_fill = "Adult Lifespan Hubness profile",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "bottomleft", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
  legend = rownames(Radar_functional_role_median_age), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)
################################################################################
# Box plot ---------------------------------------------------------------------
################################################################################

data_post_clustering %>%
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full") %>%
  arrange(cluster)


data_box <- data_post_clustering %>%
  pivot_longer(
    cols = !c("Subj_ID", "cluster", "Age"),
    names_to = "Metrics",
    values_to = "Metric_value"
  ) %>%
  plyr::rename(c("cluster" = "Median_Age"))


data_box$Metrics <- factor(data_box$Metrics, levels = c(
  "Connector", "Satellite", "Provincial", "Peripheral",
  "Local_Bridge", "Global_Bridge", "Super_Bridge", "Balance_eff"
))

data_box_sig <- data_box %>%
  group_by(Metrics) %>%
  rstatix::t_test(Metric_value ~ Median_Age, ref.group = "68") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Metrics")

ggplot(data_box, aes(x = Metrics, y = Metric_value)) +
  geom_boxplot(aes(fill = Median_Age)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_pubr() +
  stat_pvalue_manual(data_box_sig,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE,
    bracket.nudge.y = -5
  )

# Variabilité inter_sujet semble s'exprimer davantage tôt dans la vie

# Difference in the hubness profile --------------------------------------------
Radar_young_old <- data_post_clustering %>%
  dplyr::select(-Age) %>%
  group_by(cluster) %>%
  summarize_at(vars(Connector:Super_Bridge), mean) %>%
  filter(cluster == "23" | cluster == "68") %>%
  remove_rownames() %>%
  column_to_rownames(var = "cluster")

radarplotting_overlap(Radar_young_old, 60, 0, 1, 1,
  alpha = 0.05, label_size = 1,
  title_fill = "Adult Lifespan Hubness profile",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "bottomleft", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
  legend = rownames(Radar_young_old), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)

# Difference in the proportion of each functional role within each RSN for the two selected clusters
delta_hubness_profile <- function(cluster1, cluster2, max, min, max2, min2, alpha) {
  # Pick clusters
  a <- cluster1
  b <- cluster2

  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  # 1a
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(cluster == a | cluster == b)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID) %>%
    merge(., tmp_cluster_0 %>% dplyr::select(Subj_ID, cluster),
      by = "Subj_ID"
    ) %>%
    dplyr::select(Subj_ID, cluster, Region, `1st_network`)

  # 2a
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  # This df has the right number of subjects and regions, need to add the RSN networks
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )

  get_RSN_label <- data_functional_role %>%
    filter(Region %in% data_hub_selection_cluster$Region) %>%
    dplyr::select(`1st_network`, Region) %>%
    distinct()

  tmp_cluster_2 <- merge(data_hub_selection_cluster, get_RSN_label, by = "Region")

  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  tmp_cluster_final <- merge(tmp_cluster_2, tmp_cluster_0 %>%
    dplyr::select(Subj_ID, cluster),
  by = "Subj_ID"
  ) %>%
    mutate(Median_age = cluster)


  delta_proportion_a <- tmp_cluster_final %>%
    group_by(`1st_network`, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, cluster) %>%
    summarize_at(vars(Connector:Satellite), sum) %>%
    pivot_longer(cols = !c("1st_network", "cluster"), names_to = "Hub_consensus", values_to = "freq") %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(`1st_network`, Hub_consensus) %>%
    mutate(delta_freq = freq / lag(freq)) %>%
    dplyr::select(-cluster) %>%
    na.omit()

  delta_proportion_b <- tmp_cluster_final %>%
    group_by(`1st_network`, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-c(None, n)) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, cluster) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge), sum) %>%
    pivot_longer(cols = !c("1st_network", "cluster"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(`1st_network`, Bridgeness) %>%
    mutate(delta_freq = freq / lag(freq)) %>%
    dplyr::select(-cluster) %>%
    na.omit()

  Radar_functional_role_RSN_delta <-
    delta_proportion_a %>%
    dplyr::select(`1st_network`, Hub_consensus, delta_freq) %>%
    spread(`1st_network`, delta_freq) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    subset(Hub_consensus != "None") %>%
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus")


  radarplotting_overlap(Radar_functional_role_RSN_delta, max, min, 1, 1,
    alpha = alpha, label_size = 1,
    title_fill = paste("Delta of the hubness profile between cluster", cluster1, "and", cluster2, "\n Each value represent the ratio in proportion of functional role within each RSN. A positive ratio favors", cluster2),
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )

  Radar_functional_role_RSN_delta <-
    delta_proportion_b %>%
    dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
    spread(`1st_network`, delta_freq) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    subset(Bridgeness != "None") %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")

  radarplotting_overlap(Radar_functional_role_RSN_delta, max2, min2, 1, 1,
    alpha = alpha, label_size = 1,
    title_fill = paste("Delta of the hubness profile between cluster", cluster1, "and", cluster2, "\n Each value represent the ratio in proportion of functional role within each RSN. A positive ratio favors", cluster2),
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}

delta_hubness_profile("23", "68", 9, -3, 4, -2, 0.1)
