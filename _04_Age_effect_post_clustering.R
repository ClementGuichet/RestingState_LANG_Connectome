##########################################################################################
# Script for Age-related analyses and visualization for the topologico-functional profile

# Written by CG
# 13-12-2022
##########################################################################################
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(rgl)
library(FactoMineR)
library(factoextra)
library(FactoInvestigate)
library(fitdistrplus)
library(sigclust2)


rm(list = ls())

source("_04_Age_effect.R")
source("_radarplotting_function.R")

################################################################################
# Post-clustering analysis -----------------
################################################################################
data_post_clustering %>%
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full")

data_post_clustering <- data_post_clustering %>%
  mutate(cluster = ifelse(cluster == "1", "25",
                          ifelse(cluster == "2", "50",
                                 ifelse(cluster == "3", "25_bis", 0)
                          )
  ))

data_post_clustering %>%
  group_by(cluster) %>% count(Gender) %>% mutate(n = prop.table(n))

# # Groups:   cluster [3]
# cluster Gender     n
# <int> <chr>  <dbl>
# 1 25      F      0.583 
# 2 25      M      0.417 
# 3 25_bis  F      0.571 
# 4 25_bis  M      0.4   
# 5 25_bis  NaN    0.0286
# 6 50      F      0.36  
# 7 50      M      0.6   
# 8 50      NaN    0.04 


a <- gghistogram(data_post_clustering %>% subset(cluster == "25"),
            x = "Age", y = "..density..",
            fill = "purple", add_density = TRUE
) + theme_pubclean()
b <- gghistogram(data_post_clustering %>% subset(cluster == "25_bis"),
            x = "Age", y = "..density..",
            fill = "purple", add_density = TRUE
) + theme_pubclean()
c <- gghistogram(data_post_clustering %>% subset(cluster == "50"),
            x = "Age", y = "..density..",
            fill = "purple", add_density = TRUE
) + theme_pubclean()

Rmisc::multiplot(a, b, c)



################################################################################
# PCA --------------------------------------------------------------------------
pc <- PCA(data_post_clustering %>% dplyr::select(-c(Subj_ID, Age, Gender, cluster)), axes = c(1, 2))
pc$var
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
# Box plot ---------------------------------------------------------------------
################################################################################

data_post_clustering %>%
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full") %>%
  arrange(cluster)


data_box <- data_post_clustering %>%
  pivot_longer(
    cols = !c("Subj_ID", "cluster", "Gender","Age"),
    names_to = "Metrics",
    values_to = "Metric_value"
  ) %>%
  plyr::rename(c("cluster" = "Median_Age"))


data_box$Metrics <- factor(data_box$Metrics, levels = c(
  "Connector", "Satellite", "Provincial", "Peripheral",
  "Local_Bridge", "Global_Bridge", "Super_Bridge", "Balance_eff"))

data_box_sig <- data_box %>%
  group_by(Metrics) %>%
  rstatix::t_test(Metric_value ~ Median_Age) %>%
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

data_box_eff_size <- data_box %>%
  group_by(Metrics) %>%
  rstatix::cohens_d(Metric_value ~ Median_Age, comparison = list(c("25", "50")), paired = FALSE, hedges.correction = TRUE) %>%
  mutate(effsize = effsize * (-1)) %>% filter(magnitude != "negligible")

ggdotchart(
  data_box_eff_size, x = "Metrics", y = "effsize",
  ylab = "Cohen's d (hedge corrected)",
  palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "descending",
  rotate = TRUE, legend = "none"
)

cor <- cor.test(data_post_clustering$Global_Bridge, data_post_clustering$Super_Bridge)
# Creating the plot
plot(data_post_clustering$Super_Bridge, data_post_clustering$Global_Bridge, pch = 19, col = "lightblue")
# Regression line
abline(lm(data_post_clustering$Global_Bridge ~ data_post_clustering$Super_Bridge), col = "red", lwd = 3)
# Pearson correlation
text(paste0("Correlation between Super & Global Bridge: ", round(cor$estimate, 2), "****"), x = 50, y = 30)

################################################################################
# Topologico-functional profile across clusters  -------------------------------
################################################################################

Radar_functional_role_median_age <- data_post_clustering %>%
  dplyr::select(-Age) %>%
  group_by(cluster) %>%
  summarize_at(vars(Connector:Super_Bridge), mean) %>%
  remove_rownames() %>%
  column_to_rownames(var = "cluster")

radarplotting_overlap(Radar_functional_role_median_age, 60, 0, 1, 1,
                      alpha = 0.05, label_size = 1,
                      title_fill = "Distribution of functional roles",
                      palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "topright", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
  legend = rownames(Radar_functional_role_median_age), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)

# Distribution of hubs across RSNs for each cluster for the individual hubs ----

# Get the associated Resting-state networks
tmp_cluster_1 <- data_functional_role %>% filter(Age != "NaN")

# Hub region specific to each subject yielded by hub detection procedure
data_hub_selection_per_subject <- rbindlist(Hub_selection)
# Removing rows from subjects with Age = NaN
data_hub_selection_cluster <- filter(
  data_hub_selection_per_subject,
  Subj_ID %in% tmp_cluster_1$Subj_ID
)

# Final dataframe with the subjects, their cluster assignment, their hub regions and the RSNs
tmp_cluster_final <- merge(data_hub_selection_cluster, data_post_clustering %>%
                             dplyr::select(Subj_ID, cluster),
                           by = "Subj_ID"
)

Radar_hub_RSN <- tmp_cluster_final %>% 
  group_by(cluster, `1st_network`) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::select(-n) %>% 
  spread(`1st_network`, freq) %>% 
  remove_rownames() %>% column_to_rownames("cluster") %>% 
  mutate_at(vars(everything()), funs(. * 100))

radarplotting_overlap(Radar_hub_RSN, 30, 0, 1, 1,
                      alpha = 0.05, label_size = 1,
                      title_fill = "Distribution of hubs regions across RSNs",
                      palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "topright", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
  legend = rownames(Radar_hub_RSN), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)


# What are the hubs that are most consistently assigned the same functional role across subjects? ----

# Radar plot of top % hub regions consistently labeled as a specific hub type across subjects
# This means that the funcitonal labeling is the most consistent for that hub across subjects

plot_functional_role_Region <- data_hub_selection_per_subject %>%
  group_by(Region, Hub_consensus) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Region, desc(freq))

# Select only the top n regions to be displayed
top <- 5

Radar_functional_role_Region <- plot_functional_role_Region %>%
  dplyr::select(Hub_consensus, freq) %>%
  ungroup() %>%
  group_by(Hub_consensus, .add = TRUE) %>%
  group_split() %>%
  map_dfr(. %>% slice_max(freq, n = top)) %>%
  spread(Region, freq) %>%
  subset(Hub_consensus != "None") %>%
  remove_rownames() %>%
  column_to_rownames(var = "Hub_consensus") %>%
  mutate_at(vars(everything()), funs(. * 100))

radarplotting(Radar_functional_role_Region, 100, 20, 2, 2,
              alpha = 0.3, label_size = 1,
              palette = RColorBrewer::brewer.pal(8, "Dark2")
)

plot_functional_role_Region <- data_hub_selection_per_subject %>%
  group_by(Region, Bridgeness) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Region, desc(freq))

# Select only the top n regions to be displayed
top <- 5

Radar_functional_role_Region <- plot_functional_role_Region %>%
  dplyr::select(Bridgeness, freq) %>%
  ungroup() %>%
  group_by(Bridgeness, .add = TRUE) %>%
  group_split() %>%
  map_dfr(. %>% slice_max(freq, n = top)) %>%
  spread(Region, freq) %>%
  subset(Bridgeness != "None") %>%
  remove_rownames() %>%
  column_to_rownames(var = "Bridgeness") %>%
  mutate_at(vars(everything()), funs(. * 100))

radarplotting(Radar_functional_role_Region, 100, 20, 2, 2,
              alpha = 0.3, label_size = 1,
              palette = RColorBrewer::brewer.pal(8, "Dark2")
)



##########################################################################
# Difference in the Topologico-functional profile ------------------------------
################################################################################

# Difference in the proportion of each functional role within each RSN for the two selected clusters

interaction_Age_FuncRole_RSN <- function(cluster1, cluster2, max, min, max2, min2, alpha) {
  # Pick clusters
  a <- cluster1
  b <- cluster2
  
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  # 1a
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(cluster == a | cluster == b)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID)
  # 2a
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  tmp_cluster_final <- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, cluster),
                             by = "Subj_ID"
  ) 
  
  # Distribution of functional roles for each cluster for the individual hubs
  Radar_functional_role <- data_post_clustering %>%
    dplyr::select(-Age) %>%
    group_by(cluster) %>%
    summarize_at(vars(Connector:Super_Bridge), mean) %>%
    filter(cluster == a | cluster == b) %>%
    remove_rownames() %>%
    column_to_rownames(var = "cluster")
  
  radarplotting_overlap(Radar_functional_role, 60, 0, 1, 1,
                        alpha = 0.05, label_size = 1,
                        title_fill = "Distribution of functional roles",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "topright", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
    legend = rownames(Radar_functional_role), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  
  # Distribution of hubs across RSNs for each cluster for the individual hubs
  Radar_hub_RSN <- tmp_cluster_final %>%
    group_by(cluster, `1st_network`) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(`1st_network`, freq) %>% 
    remove_rownames() %>% column_to_rownames("cluster") %>% 
    mutate_at(vars(everything()), funs(. * 100))
  
  radarplotting_overlap(Radar_hub_RSN, 25, 0, 1, 1,
                        alpha = 0.05, label_size = 1,
                        title_fill = "Distribution of hubs regions across RSNs",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "topright", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
    legend = rownames(Radar_hub_RSN), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
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
    dplyr::select(-n) %>%
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
    subset(Hub_consensus != "None") %>%
    mutate(`1st_network` = ifelse(delta_freq == "Inf"|delta_freq == 0, paste0(`1st_network`, "*"), `1st_network`)) %>% 
    mutate(delta_freq = ifelse(delta_freq == "Inf", 1, delta_freq)) %>%
    spread(`1st_network`, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus")
  
  
  radarplotting_overlap(Radar_functional_role_RSN_delta, max, min, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Ratio of the proportion of functional roles between cluster", cluster1, "and", cluster2, "\n A positive ratio means", cluster2, "has x times more {Connector} hubs within the RSN"),
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
    subset(Bridgeness != "None") %>%
    mutate(`1st_network` = ifelse(delta_freq == "Inf"|delta_freq == 0, paste0(`1st_network`, "*"), `1st_network`)) %>% 
    mutate(delta_freq = ifelse(delta_freq == "Inf", 1, delta_freq)) %>%
    spread(`1st_network`, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")
  
  radarplotting_overlap(Radar_functional_role_RSN_delta, max2, min2, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Ratio of the proportion of functional roles between cluster", cluster1, "and", cluster2, "\n A positive ratio means", cluster2, "has x times more {Global Bridge} hubs within the RSN"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", title = "* indicates there were no Global Bridge hubs in VMM for young subjects",
    legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}
interaction_Age_FuncRole_RSN("25", "50", 4, 0, 6, -2, 0.1)

# For each RSN, is there a difference in proportion of each functional role between clusters?
boxplot_sig_interaction_RSN <- function(cluster1, cluster2){
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
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  tmp_cluster_final <- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, cluster),
                             by = "Subj_ID"
  ) 
  
  Age_RSN_prop_1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Subj_ID, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    # Make sure comparisons with missing functional roles can be achieved
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% 
    group_by(`1st_network`, Subj_ID, cluster) %>%
    summarize_at(vars(Connector:Satellite), sum) %>% 
    pivot_longer(cols = !c("1st_network", "Subj_ID","cluster"), names_to = "Hub_consensus", values_to = "freq") 
  
  Age_RSN_prop_2 <- tmp_cluster_final %>%
    group_by(`1st_network`, Subj_ID, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    # Make sure comparisons with missing functional roles can be achieved
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, cluster) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge), sum) %>%
    pivot_longer(cols = !c("1st_network", "Subj_ID","cluster"), names_to = "Bridgeness", values_to = "freq")
  
  Age_RSN_prop_final <- bind_rows(Age_RSN_prop_1, Age_RSN_prop_2) %>% 
    mutate(Functional_role = ifelse(is.na(Hub_consensus) == TRUE, Bridgeness, Hub_consensus)) %>% subset(Functional_role != "None") 
  # Based on the ratio observed with the radar plot
  # filter(grepl("Auditory|Language|FPN|DMN|PMM", `1st_network`))
  
  data_box_sig_RSN <- Age_RSN_prop_final %>%
    group_by(`1st_network`, Functional_role) %>%
    rstatix::t_test(freq ~ cluster) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "1st_network")
  
  ggplot(Age_RSN_prop_final, aes(x = `1st_network`, y = freq)) +
    geom_boxplot(aes(fill = cluster)) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~Functional_role) +
    theme_pubr() +
    stat_pvalue_manual(data_box_sig_RSN,
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       hide.ns = TRUE
    ) 
}
boxplot_sig_interaction_RSN("25", "50")

# Difference in the proportion of each functional role within each community for the two selected clusters
interaction_Age_FuncRole_community <- function(cluster1, cluster2, max, min, max2, min2, alpha) {
  # Pick clusters
  a <- cluster1
  b <- cluster2
  
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  # 1a
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(cluster == a | cluster == b)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID) 
  # 2a
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  tmp_cluster_final <- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, cluster),
                             by = "Subj_ID"
  ) 
  
  # Distribution of hubs across RSNs for each cluster for the individual hubs
  Radar_hub_RSN <- tmp_cluster_final %>%
    group_by(cluster, Consensus_vector_0.15) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Consensus_vector_0.15, freq) %>% 
    remove_rownames() %>% column_to_rownames("cluster") %>% 
    mutate_at(vars(everything()), funs(. * 100))
  
  radarplotting_overlap(Radar_hub_RSN, 50, 0, 1, 1,
                        alpha = 0.05, label_size = 1,
                        title_fill = "Distribution of hubs regions across communities",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "topright", title = "Median age of each cluster\n Based on Ward Hierarchichal clustering (FWER corrected)\n (Kimes et al., 2017)",
    legend = rownames(Radar_hub_RSN), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  delta_proportion_a <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(Consensus_vector_0.15, cluster) %>%
    summarize_at(vars(Connector:Satellite), sum) %>%
    pivot_longer(cols = !c("Consensus_vector_0.15", "cluster"), names_to = "Hub_consensus", values_to = "freq") %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(Consensus_vector_0.15, Hub_consensus) %>%
    mutate(delta_freq = freq / lag(freq)) %>%
    dplyr::select(-cluster) %>%
    na.omit()
  
  delta_proportion_b <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(Consensus_vector_0.15, cluster) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge), sum) %>%
    pivot_longer(cols = !c("Consensus_vector_0.15", "cluster"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(Consensus_vector_0.15, Bridgeness) %>%
    mutate(delta_freq = freq / lag(freq)) %>% 
    dplyr::select(-cluster) %>%
    na.omit()
  
  Radar_functional_role_RSN_delta <-
    delta_proportion_a %>%
    dplyr::select(Consensus_vector_0.15, Hub_consensus, delta_freq) %>%
    subset(Hub_consensus != "None") %>%
    mutate(Consensus_vector_0.15 = ifelse(delta_freq == "Inf"|delta_freq == 0, paste0(Consensus_vector_0.15, "*"), Consensus_vector_0.15)) %>% 
    mutate(delta_freq = ifelse(delta_freq == "Inf", 1, delta_freq)) %>%
    spread(Consensus_vector_0.15, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus")
  
  
  radarplotting_overlap(Radar_functional_role_RSN_delta, max, min, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Ratio of the proportion of functional roles between cluster", cluster1, "and", cluster2, "\n A positive ratio means", cluster2, "has x times more {Connector} hubs within the community"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  Radar_functional_role_RSN_delta <-
    delta_proportion_b %>%
    dplyr::select(Consensus_vector_0.15, Bridgeness, delta_freq) %>%
    subset(Bridgeness != "None") %>%
    mutate(Consensus_vector_0.15 = ifelse(delta_freq == "Inf"|delta_freq == 0, paste0(Consensus_vector_0.15, "*"), Consensus_vector_0.15)) %>% 
    mutate(delta_freq = ifelse(delta_freq == "Inf", 1, delta_freq)) %>%
    spread(Consensus_vector_0.15, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")
  
  radarplotting_overlap(Radar_functional_role_RSN_delta, max2, min2, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Ratio of the proportion of functional roles between cluster", cluster1, "and", cluster2, "\n A positive ratio means", cluster2, "has x times more {Global Bridge} hubs within the community"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", title = "* indicates there were no Global Bridge hubs in VMM for young subjects",
    legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}
interaction_Age_FuncRole_community("25", "50", 4, 0, 4, 0, 0.1)

# For each community, is there a difference in proportion of each functional role between clusters?
boxplot_sig_interaction_community <- function(cluster1, cluster2){
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
    dplyr::select(Subj_ID, cluster, Region, Consensus_vector_0.15)
  
  # 2a
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  tmp_cluster_final <- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, cluster),
                             by = "Subj_ID"
  ) 
  
  Age_RSN_prop_1 <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, Subj_ID, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    # Make sure comparisons with missing functional roles can be achieved
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% 
    group_by(Consensus_vector_0.15, Subj_ID, cluster) %>%
    summarize_at(vars(Connector:Satellite), sum) %>% 
    pivot_longer(cols = !c("Consensus_vector_0.15", "Subj_ID","cluster"), names_to = "Hub_consensus", values_to = "freq") 
  
  Age_RSN_prop_2 <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, Subj_ID, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    # Make sure comparisons with missing functional roles can be achieved
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(Consensus_vector_0.15, Subj_ID, cluster) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge), sum) %>%
    pivot_longer(cols = !c("Consensus_vector_0.15", "Subj_ID","cluster"), names_to = "Bridgeness", values_to = "freq")
  
  Age_RSN_prop_final <- bind_rows(Age_RSN_prop_1, Age_RSN_prop_2) %>% 
    mutate(Functional_role = ifelse(is.na(Hub_consensus) == TRUE, Bridgeness, Hub_consensus)) %>% subset(Functional_role != "None") 
  # Based on the ratio observed with the radar plot
  # filter(grepl("Auditory|Language|FPN|DMN|PMM", Consensus_vector_0.15))
  
  data_box_sig_RSN <- Age_RSN_prop_final %>%
    group_by(Consensus_vector_0.15, Functional_role) %>%
    rstatix::t_test(freq ~ cluster) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Consensus_vector_0.15")
  
  ggplot(Age_RSN_prop_final, aes(x = Consensus_vector_0.15, y = freq)) +
    geom_boxplot(aes(fill = cluster)) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~Functional_role) +
    theme_pubr() +
    stat_pvalue_manual(data_box_sig_RSN,
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       hide.ns = TRUE
    ) 
}
boxplot_sig_interaction_community("25", "50")
