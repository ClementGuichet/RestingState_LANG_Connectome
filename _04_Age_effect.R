##########################################################################################
# Script for Age-related analyses

# Written by CG
# 26-11-2022
##########################################################################################
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(fitdistrplus)
library(sigclust2)
library(ggpubr)
library(cluster)
library(FactoMineR)
library(factoextra)
library(FactoInvestigate)

rm(list = ls())

source("_03_Hub_classification.R")
source("_radarplotting_function.R")

# What are the graph-based biomarkers of health aging ?
# Does the topological_functional profile, that is the proportion of the hubs' functional role, change throughout life ?

################################################################################
# Correlation between degree centrality and Age --------------------------------
# Inferential statistics

data_stat_age <- data_functional_role %>%
  group_by(Subj_ID, CAB_NP_assign, Region, Age, Consensus_vector_0.15, LANG_Net_assign) %>%
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
  group_by(CAB_NP_assign, Region, Consensus_vector_0.15, LANG_Net_assign) %>%
  group_split() %>%
  map_dfr(. %>%
    mutate(Estimate = cor.test(.$DC, .$Age, method = "kendall")$estimate) %>%
    mutate(p_value = cor.test(.$DC, .$Age, method = "kendall")$p.value)) %>%
  group_by(CAB_NP_assign, Region, Consensus_vector_0.15, LANG_Net_assign) %>%
  summarise_at(vars(Estimate, p_value), mean)

ggdotchart(
  correlation_DC_age %>% subset(p_value <= 0.05),
  x = "Region", y = "Estimate",
  ylab = "Mean correlation",
  palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "descending",
  rotate = TRUE, legend = "none", title = "Significant correlations between degree centrality and Age"
)

ComplexHeatmap::Heatmap(as.matrix(correlation_DC_age %>% subset(p_value <= 0.05) %>% arrange(desc(CAB_NP_assign)) %>% ungroup() %>% dplyr::select(Region, Estimate) %>% remove_rownames() %>% column_to_rownames("Region")), cluster_rows = FALSE, column_names_rot = TRUE, name = "Heatmap of significant correlations\n between degree centrality and Age")

################################################################################
# ~~~~~~~~~~~ Hub Detection Procedure ~~~~~~~~~~~
################################################################################
# Method ~ Detect top % regions for each metric ------------------------------

# Topologico-functional profile with hub detection at the individual level -------------------

top <- 131 * 0.2

# LOG: Testing if clustering is stable when selecting another top% of regions
# Initial clustering at 20% yields 3 stable clusters
#
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
    dplyr::select(Subj_ID, Region, `1st_network`, Consensus_vector_0.15, Hub_consensus, Bridgeness, zFlow)
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


# What are the most common hubs across subjects?
most_common_hubs <- rbindlist(Hub_selection) %>%
  count(Region, `1st_network`) %>%
  mutate(n = n / 72) %>%
  arrange(desc(n)) %>%
  filter(n > 0.8) %>%
  mutate_at(vars(n), funs(. * 100))

###############################################################################
# Getting ready for hierarchical clustering analysis

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
  rbindlist(FR_list, fill = TRUE) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
  data_functional_role %>% group_by(Subj_ID, Gender) %>% summarise_at(vars(Age), mean) %>% arrange(Subj_ID),
  Balance_eff = data_cluster_efficiency$Balance_eff
)

data_hubness_profile_Age_ind %>%
  pivot_longer(
    cols = !c("Subj_ID", "Age", "Gender"),
    names_to = "Functional_role",
    values_to = "Score"
  ) %>%
  ggplot(aes(Age, Score, color = Functional_role)) +
  geom_point(size = 2, alpha = 0.2) +
  geom_jitter(height = 0.05, alpha = 0.2) +
  geom_smooth() +
  ggpubr::theme_pubr() +
  ggtitle("Evolution of functional roles across adult lifespan")

data_pre_clustering <- data_hubness_profile_Age_ind %>%
  dplyr::select(-c(Subj_ID, Age, Gender, None, Not_a_Bridge))


################################################################################
# PCA --------------------------------------------------------------------------
# Keeping 5 components which explain 90% of variance
pc <- PCA(data_pre_clustering, ncp = 5, scale.unit = T, axes = c(1, 2))
pc$eig
dimdesc(pc, axes = 1:5)

data_cluster <- pc$ind$coord %>% as.data.frame()
# Factoshiny::Factoshiny(pc)

# Clustering -----------------
performance::check_clusterstructure(data_cluster)
corrplot::corrplot(cor(data_cluster))
ComplexHeatmap::Heatmap(scale(data_pre_clustering))


# Define linkage methods
link <- c("average", "single", "complete", "ward")
names(link) <- c("average", "single", "complete", "ward")
# Function to compute the agglomerative coefficient (i.e., the strength of the clusters)
ac <- function(x) {
  cluster::agnes(data_cluster, method = x)$ac
}
sapply(link, ac)

# Defining optimal number of cluster with gap statistic
gap_stat <- cluster::clusGap(data_cluster,
  FUN = hcut,
  K.max = 5,
  B = 1000,
  verbose = T
)

factoextra::fviz_gap_stat(gap_stat)

plot_cluster <- cluster::agnes(data_cluster, method = "ward", metric = "euclidian")
fviz_dend(plot_cluster,
  cex = 0.8,
  k = 3,
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
groups <- cutree(final_clust, k = 3)


data_post_clustering <- cbind(data_hubness_profile_Age_ind, cluster = groups) %>%
  dplyr::select(-c(None, Not_a_Bridge))
