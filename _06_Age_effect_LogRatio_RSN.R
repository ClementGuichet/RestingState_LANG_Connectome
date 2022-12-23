
################################################################################
# LOG-RATIO of the Topologico-functional profile ------------------------------
################################################################################
source("_05_Age_effect_Kullback_Leibler.R")
# Difference in the proportion of each functional role within each RSN for the two selected clusters
cluster1 <- "23"
cluster2 <- "56"

interaction_Age_FuncRole_RSN <- function(cluster1, cluster2, max, min, max2, min2, alpha, divergent_RSN_modular, divergent_RSN_interareal) {
  # Pick clusters
  a <- cluster1
  b <- cluster2
  
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
  
  delta_proportion_a1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, cluster) %>%
    # Computing the PMFs for each region for each cluster for each functional role
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup()
  
    # Replace essential zeros at the regional and cluster level
  replacement_a1 <- multRepl(delta_proportion_a1 %>%
                               dplyr::select(Connector:Satellite), dl = rep(1, 5), label = 0, frac = 1e-5)
  
  delta_proportion_a <<- cbind(delta_proportion_a1 %>%
                                 dplyr::select(`1st_network`:cluster), replacement_a1) %>%
    dplyr::select(-Isolate) %>% 
    pivot_longer(cols = !c("1st_network", "Region", "cluster"), names_to = "Hub_consensus", values_to = "freq") %>%
    # Compute the geometric mean for each RSN
    group_by(`1st_network`, cluster, Hub_consensus) %>%
    summarise_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(`1st_network`, Hub_consensus) %>%
    # Take the log ratio of geometric mean between clusters
    # 
    # This is the same as taking the mean log ratio of all regions:
    # 
    # group_by(`1st_network`, Region, Hub_consensus) %>%
    #   mutate(delta_freq = log(freq / lag(freq))) %>% 
    #   dplyr::select(-cluster) %>%
    #   na.omit() %>% 
    #   group_by(`1st_network`, Hub_consensus) %>%
    #   summarize_at(vars(delta_freq), mean)
    mutate(delta_freq = log(freq / lag(freq))) %>% 
    dplyr::select(-c(freq, cluster)) %>%
    na.omit()
  
  
  delta_proportion_b1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, cluster) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup()
  
  replacement_b1 <- multRepl(delta_proportion_b1 %>%
                               dplyr::select(Global_Bridge:Super_Bridge), dl = rep(1, 4), label = 0, frac = 1e-5)
  
  delta_proportion_b <<- cbind(delta_proportion_b1 %>%
                                 dplyr::select(`1st_network`:cluster), replacement_b1) %>%
    pivot_longer(cols = !c("1st_network", "Region", "cluster"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the geometric mean for each RSN
    group_by(`1st_network`, cluster, Bridgeness) %>%
    summarise_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(`1st_network`, Bridgeness) %>%
    # Take the log ratio of geometric mean between clusters
    # 
    # This is the same as taking the mean log ratio of all regions:
    # 
    # group_by(`1st_network`, Region, Bridgeness) %>%
    #   mutate(delta_freq = log(freq / lag(freq))) %>% 
    #   dplyr::select(-cluster) %>%
    #   na.omit() %>% 
    #   group_by(`1st_network`, Bridgeness) %>%
    #   summarize_at(vars(delta_freq), mean)
    mutate(delta_freq = log(freq / lag(freq))) %>% 
    dplyr::select(-c(freq, cluster)) %>%
    na.omit()
  
  Radar_functional_role_RSN_delta_a <-
    delta_proportion_a %>%
    dplyr::select(`1st_network`, Hub_consensus, delta_freq) %>%
    filter(grepl(divergent_RSN_modular, `1st_network`)) %>%
    spread(`1st_network`, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus")
  
  
  radarplotting_overlap(Radar_functional_role_RSN_delta_a, max, min, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Log ratio of the clusters' geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older cluster"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_a), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  Radar_functional_role_RSN_delta_b <-
    delta_proportion_b %>%
    dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
    filter(grepl(divergent_RSN_interareal, `1st_network`)) %>%
    spread(`1st_network`, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")
  
  radarplotting_overlap(Radar_functional_role_RSN_delta_b, max2, min2, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Log ratios of geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older cluster"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}
interaction_Age_FuncRole_RSN(cluster1, cluster2, 6, -6, 6, -6, 0.1,
                             divergent_RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN",
                             divergent_RSN_interareal = "Auditory|CON|DMN|FPN|Language|PMM|SMN|Visual_2"
)

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
  
  delta_proportion_a1 <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, Region, Subj_ID, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% 
    group_by(Consensus_vector_0.15, Region, cluster) %>%
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup()
  
  replacement_a1 <- multRepl(delta_proportion_a1 %>%
                               dplyr::select(Connector:Satellite), dl = rep(1, 5), label = 0, frac = 1e-5)
  
  delta_proportion_a <<- cbind(delta_proportion_a1 %>%
                                 dplyr::select(Consensus_vector_0.15:cluster), replacement_a1) %>%
    dplyr::select(-Isolate) %>% 
    pivot_longer(cols = !c("Consensus_vector_0.15", "Region", "cluster"), names_to = "Hub_consensus", values_to = "freq") %>%
    # Compute the geometric mean for each community
    group_by(Consensus_vector_0.15, cluster, Hub_consensus) %>%
    summarise_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(Consensus_vector_0.15, Hub_consensus) %>%
    # Take the log ratio of geometric mean between clusters
    mutate(delta_freq = log(freq / lag(freq))) %>% 
    dplyr::select(-c(freq, cluster)) %>%
    na.omit()
  
  
  delta_proportion_b1 <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, Region, Subj_ID, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(Consensus_vector_0.15, Region, cluster) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup()
  
  replacement_b1 <- multRepl(delta_proportion_b1 %>%
                               dplyr::select(Global_Bridge:Super_Bridge), dl = rep(1, 4), label = 0, frac = 1e-5)
  
  delta_proportion_b <<- cbind(delta_proportion_b1 %>%
                                 dplyr::select(Consensus_vector_0.15:cluster), replacement_b1) %>%
    pivot_longer(cols = !c("Consensus_vector_0.15", "Region", "cluster"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the geometric mean for each community
    group_by(Consensus_vector_0.15, cluster, Bridgeness) %>%
    summarise_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(Consensus_vector_0.15, Bridgeness) %>%
    # Take the log ratio of geometric mean between clusters
    mutate(delta_freq = log(freq / lag(freq))) %>% 
    dplyr::select(-c(freq, cluster)) %>%
    na.omit()
  
  
  Radar_functional_role_RSN_delta_a <-
    delta_proportion_a %>%
    spread(Consensus_vector_0.15, delta_freq) %>% 
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus")
  
  radarplotting_overlap(Radar_functional_role_RSN_delta_a, max, min, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Log ratio of the clusters' geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older cluster"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_a), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  Radar_functional_role_RSN_delta_b <-
    delta_proportion_b %>%
    spread(Consensus_vector_0.15, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")
  
  radarplotting_overlap(Radar_functional_role_RSN_delta_b, max2, min2, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Log ratio of the clusters' geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older cluster"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}
interaction_Age_FuncRole_community(cluster1, cluster2, 4, -4, 6, -6, 0.1
)
# More global interfaces in RS-NET 2 & 4 but more Local & mixed interfaces in RS-NET 3
# Proportion in RS-NET 1 do not seem to change
# Topological reconfiguration is not equal across resting-state networks and the modules they form


################################################################################
# BOOTSTRAPPING
# ##############################################################################

# Bootstrapping procedure to compute confidence intervals of log-ratios
#
# The geometric mean of each behavior (%) in both clusters are calculated.
# The log-ratio of geometric means of both clusters is computed.
# First, 1000 virtual data sets are drawn with replacement from the source population and of the same size.
# That, for each combination of cluster and RSN, we take 1000 resample with replacement
# For each resample, the log-ratio of the geometric mean is calculated.
# The resulting distribution of 1000 log-ratios are averaged to calculate bootstrapped mean,
# and the 2.5th and 97.5th percentiles are selected as upper and lower limits of 95% confidence intervals of the bootstrapped mean.

# If the confidence interval contains the value ‘0’, no difference between the two groups for this particular cluster/RSN combination is concluded.
# Only combinations for which the intervals are outside 0 are considered responsible for the cluster differences.

bootstrap_ci <- function(n_boot, divergent_RSN_modular, divergent_RSN_interareal) {
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(cluster == "23" | cluster == "56")
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID)
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
  
  data_boot <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, cluster, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, cluster, Subj_ID) %>%
    summarize_at(vars(Connector:Satellite), sum) %>%
    ungroup() %>%
    filter(grepl(divergent_RSN_modular, `1st_network`)) %>% 
    # Draw samples with replacement for every region, RSN, and cluster
    group_by(cluster, `1st_network`, Region, .add = TRUE) %>%
    group_split()
  
  list_boot <- list()
  for (i in 1:length(data_boot)) {
    tmp <- rbindlist(data_boot[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    list_boot[[i]] <- rbindlist(list_boot_bis)
  }
  
  # 1000 samples for each cluster*1st_network*Region combination
  # So that for each region for each RSN, and within each cluster, 1000 resamples of the same size are taken
  resamples <- rbindlist(list_boot)
  
  delta_proportion_boot_a <- resamples %>%
    group_by(n_data_boot, n_data_split, `1st_network`, Region, cluster) %>%
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup()
  
  # Multiplicative replacement at the regional and cluster level
  replacement_a <- delta_proportion_boot_a %>%
    group_by(n_data_boot, .add = TRUE) %>%
    group_split()
  
  list_replacement <- list()
  for (i in 1:length(replacement_a)) {
    list_replacement[[i]] <- multRepl(rbindlist(replacement_a[i]) %>%
                                        dplyr::select(Connector:Satellite), dl = rep(1, 5), label = 0, frac = 1e-5)
  }
  
  delta_a_detail <<- cbind(
    delta_proportion_boot_a %>%
      dplyr::select(n_data_boot, n_data_split, `1st_network`:cluster),
    rbindlist(list_replacement)
  ) %>%
    dplyr::select(-Isolate) %>%
    pivot_longer(
      cols = !c("n_data_boot", "n_data_split", "1st_network", "Region", "cluster"),
      names_to = "Hub_consensus", values_to = "freq"
    ) %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(n_data_boot, `1st_network`, cluster, Hub_consensus) %>%
    summarize_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(n_data_boot, `1st_network`, Hub_consensus) %>%
    # Take the log ratios of geometric means
    mutate(delta_freq = log(freq / lag(freq))) %>%
    dplyr::select(-c(freq, cluster)) %>%
    na.omit()
  
  summary_bootstrap_a <<- delta_a_detail %>%
    group_by(`1st_network`, Hub_consensus) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Hub_consensus, mean, mean_025, mean_0975)
  
  data_boot <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, cluster, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    ungroup() %>%
    filter(grepl(divergent_RSN_interareal, `1st_network`)) %>% 
    group_by(cluster, `1st_network`, Region, .add = TRUE) %>%
    group_split()
  
  list_boot <- list()
  for (i in 1:length(data_boot)) {
    tmp <- rbindlist(data_boot[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    list_boot[[i]] <- rbindlist(list_boot_bis)
  }
  
  resamples <- rbindlist(list_boot)
  
  delta_proportion_boot_b <- resamples %>%
    group_by(n_data_boot, n_data_split, `1st_network`, Region, cluster) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup()
  
  # Multiplicative replacement
  replacement_b <- delta_proportion_boot_b %>%
    group_by(n_data_boot, .add = TRUE) %>%
    group_split()
  
  list_replacement <- list()
  for (i in 1:length(replacement_b)) {
    list_replacement[[i]] <- multRepl(rbindlist(replacement_b[i]) %>%
                                        dplyr::select(Global_Bridge:Super_Bridge), dl = rep(1, 4), label = 0, frac = 1e-5)
  }
  
  delta_b_detail <<- cbind(
    delta_proportion_boot_b %>%
      dplyr::select(n_data_boot, n_data_split, `1st_network`:cluster),
    rbindlist(list_replacement)
  ) %>%
    pivot_longer(
      cols = !c("n_data_boot", "n_data_split", "1st_network", "Region", "cluster"),
      names_to = "Bridgeness", values_to = "freq"
    ) %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(n_data_boot, `1st_network`, cluster, Bridgeness) %>%
    summarize_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(n_data_boot, `1st_network`, Bridgeness) %>%
    # Take the log ratios of geometric means
    mutate(delta_freq = log(freq / lag(freq))) %>%
    dplyr::select(-c(freq, cluster)) %>%
    na.omit()
  
  summary_bootstrap_b <<- delta_b_detail %>%
    group_by(`1st_network`, Bridgeness) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Bridgeness, mean, mean_025, mean_0975)
}

bootstrap_ci(n_boot = 1000, 
             divergent_RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN",
             divergent_RSN_interareal = "Auditory|CON|DMN|FPN|Language|PMM|SMN|Visual_2")

write.csv(summary_bootstrap_a, "summary_bootstrap_a.csv")
write.csv(summary_bootstrap_b, "summary_bootstrap_b.csv")

data_a <- read.csv("summary_bootstrap_a.csv")

library(ggpubr)

ggplot(data_a, aes(x = Hub_consensus, y = mean)) +
  geom_point(size = 2) + # Plot the individual points
  geom_errorbar(aes(ymin = mean_025, ymax = mean_0975), width = 0.2) + # Plot the confidence interval
  geom_line(aes(y = mean), color = "red", size = 1) + # Plot the mean
  labs(x = "Functional roles", y = "Log-ratio of the geometric means between clusters") + # Label the axes
  facet_wrap(~X1st_network, ncol = 3) +
  geom_hline(yintercept = 0, color = "red") +
  theme_pubclean()


data_b <- read.csv("summary_bootstrap_b.csv")

ggplot(data_b, aes(x = Bridgeness, y = mean)) +
  geom_point() + # Plot the individual points
  geom_errorbar(aes(ymin = mean_025, ymax = mean_0975), width = 0.2) + # Plot the confidence interval
  geom_line(aes(y = mean), color = "red", size = 1) + # Plot the mean
  labs(x = "Functional roles", y = "Log-ratio of the geometric means between clusters") + # Label the axes
  facet_wrap(~X1st_network, ncol = 3) +
  geom_hline(yintercept = 0, color = "red") +
  theme_pubclean()
