
################################################################################
# LOG-RATIO of the Topologico-functional profile ------------------------------
################################################################################
source("_07_Entropy.R")

# Difference in the proportion of each functional role within each RSN for the two selected clusters

interaction_Age_FuncRole_RSN <- function(cluster1, cluster2, max, min, max2, min2, alpha, divergent_RSN_modular, divergent_RSN_interareal) {
  
  data_cluster_selection(cluster1, cluster2)
  
  delta_proportion_a1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, Age_group) %>%
    # Computing the PMFs for each region for each Age_group for each functional role
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup()
  
    # Replace essential zeros at the regional and Age_group level
  replacement_a1 <- multRepl(delta_proportion_a1 %>%
                               dplyr::select(Connector:Satellite), dl = rep(1, 4), label = 0, frac = 1e-5)
  
  delta_proportion_a <<- cbind(delta_proportion_a1 %>%
                                 dplyr::select(`1st_network`:Age_group), replacement_a1) %>%
    pivot_longer(cols = !c("1st_network", "Region", "Age_group"), names_to = "Hub_consensus", values_to = "freq") %>%
    # Compute the geometric mean for each RSN
    group_by(`1st_network`, Age_group, Hub_consensus) %>%
    summarise_at(vars(freq), funs(geometricmean(.)))  %>% 
    group_by(`1st_network`, Hub_consensus) %>%
    # Take the log ratio of geometric mean between clusters
    # 
    # This is the same as taking the mean log ratio of all regions:
    # 
    # group_by(`1st_network`, Region, Hub_consensus) %>%
    #   mutate(delta_freq = log(freq / lag(freq))) %>% 
    #   dplyr::select(-Age_group) %>%
    #   na.omit() %>% 
    #   group_by(`1st_network`, Hub_consensus) %>%
    #   summarize_at(vars(delta_freq), mean)
    mutate(delta_freq = log(freq / lead(freq))) %>% 
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  
  delta_proportion_b1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, Age_group) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup()
  
  replacement_b1 <- multRepl(delta_proportion_b1 %>%
                               dplyr::select(Global_Bridge:Super_Bridge), dl = rep(1, 4), label = 0, frac = 1e-5)
  
  delta_proportion_b <<- cbind(delta_proportion_b1 %>%
                                 dplyr::select(`1st_network`:Age_group), replacement_b1) %>%
    pivot_longer(cols = !c("1st_network", "Region", "Age_group"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the geometric mean for each RSN
    group_by(`1st_network`, Age_group, Bridgeness) %>%
    summarise_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(`1st_network`, Bridgeness) %>%
    # Take the log ratio of geometric mean between clusters
    # 
    # This is the same as taking the mean log ratio of all regions:
    # 
    # group_by(`1st_network`, Region, Bridgeness) %>%
    #   mutate(delta_freq = log(freq / lag(freq))) %>% 
    #   dplyr::select(-Age_group) %>%
    #   na.omit() %>% 
    #   group_by(`1st_network`, Bridgeness) %>%
    #   summarize_at(vars(delta_freq), mean)
    mutate(delta_freq = log(freq / lead(freq))) %>% 
    dplyr::select(-c(freq, Age_group)) %>%
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
                        title_fill = paste("Log ratio of the clusters' geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older Age_group"),
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
                        title_fill = paste("Log ratios of geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older Age_group"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}
interaction_Age_FuncRole_RSN("Young", "Old", 2, -2, 2, -2, 0.1,
                             divergent_RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN",
                             divergent_RSN_interareal = "Auditory|CON|DMN|FPN|Language|SMN"
)

# Difference in the proportion of each functional role within each community for the two selected clusters
interaction_Age_FuncRole_community <- function(cluster1, cluster2, max, min, max2, min2, alpha, divergent_RSN_modular, divergent_RSN_interareal) {
  
  data_cluster_selection(cluster1, cluster2)
  
  delta_proportion_a1 <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(Consensus_vector_0.15, Region, Age_group) %>%
    # Computing the PMFs for each region for each Age_group for each functional role
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup()
  
  # Replace essential zeros at the regional and Age_group level
  replacement_a1 <- multRepl(delta_proportion_a1 %>%
                               dplyr::select(Connector:Satellite), dl = rep(1, 4), label = 0, frac = 1e-5)
  
  delta_proportion_a <<- cbind(delta_proportion_a1 %>%
                                 dplyr::select(Consensus_vector_0.15:Age_group), replacement_a1) %>%
    pivot_longer(cols = !c("Consensus_vector_0.15", "Region", "Age_group"), names_to = "Hub_consensus", values_to = "freq") %>%
    # Compute the geometric mean for each RSN
    group_by(Consensus_vector_0.15, Age_group, Hub_consensus) %>%
    summarise_at(vars(freq), funs(geometricmean(.)))  %>% 
    group_by(Consensus_vector_0.15, Hub_consensus) %>%
    # Take the log ratio of geometric mean between clusters
    # 
    # This is the same as taking the mean log ratio of all regions:
    # 
    # group_by(Consensus_vector_0.15, Region, Hub_consensus) %>%
    #   mutate(delta_freq = log(freq / lag(freq))) %>% 
    #   dplyr::select(-Age_group) %>%
    #   na.omit() %>% 
    #   group_by(Consensus_vector_0.15, Hub_consensus) %>%
    #   summarize_at(vars(delta_freq), mean)
    mutate(delta_freq = log(freq / lead(freq))) %>% 
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  
  delta_proportion_b1 <- tmp_cluster_final %>%
    group_by(Consensus_vector_0.15, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(Consensus_vector_0.15, Region, Age_group) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup()
  
  replacement_b1 <- multRepl(delta_proportion_b1 %>%
                               dplyr::select(Global_Bridge:Super_Bridge), dl = rep(1, 4), label = 0, frac = 1e-5)
  
  delta_proportion_b <<- cbind(delta_proportion_b1 %>%
                                 dplyr::select(Consensus_vector_0.15:Age_group), replacement_b1) %>%
    pivot_longer(cols = !c("Consensus_vector_0.15", "Region", "Age_group"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the geometric mean for each RSN
    group_by(Consensus_vector_0.15, Age_group, Bridgeness) %>%
    summarise_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(Consensus_vector_0.15, Bridgeness) %>%
    # Take the log ratio of geometric mean between clusters
    # 
    # This is the same as taking the mean log ratio of all regions:
    # 
    # group_by(Consensus_vector_0.15, Region, Bridgeness) %>%
    #   mutate(delta_freq = log(freq / lag(freq))) %>% 
    #   dplyr::select(-Age_group) %>%
    #   na.omit() %>% 
    #   group_by(Consensus_vector_0.15, Bridgeness) %>%
    #   summarize_at(vars(delta_freq), mean)
    mutate(delta_freq = log(freq / lead(freq))) %>% 
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  Radar_functional_role_RSN_delta_a <-
    delta_proportion_a %>%
    dplyr::select(Consensus_vector_0.15, Hub_consensus, delta_freq) %>%
    spread(Consensus_vector_0.15, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus")
  
  
  radarplotting_overlap(Radar_functional_role_RSN_delta_a, max, min, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Log ratio of the clusters' geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older Age_group"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_a), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  Radar_functional_role_RSN_delta_b <-
    delta_proportion_b %>%
    dplyr::select(Consensus_vector_0.15, Bridgeness, delta_freq) %>%
    spread(Consensus_vector_0.15, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")
  
  radarplotting_overlap(Radar_functional_role_RSN_delta_b, max2, min2, 1, 1,
                        alpha = alpha, label_size = 1,
                        title_fill = paste("Log ratios of geometric mean proportion of modular roles\n A positive ratio reflects a higher proportion for the older Age_group"),
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}
interaction_Age_FuncRole_community("Young", "Old", 1, -1, 1, -1, 0.1)


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

bootstrap_ci <- function(n_boot, cluster1, cluster2, divergent_RSN_modular, divergent_RSN_interareal) {
  
  data_cluster_selection(cluster1, cluster2)
  
  data_boot <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, Age_group, Subj_ID) %>%
    summarize_at(vars(Connector:Satellite), sum) %>%
    ungroup() %>%
    filter(grepl(divergent_RSN_modular, `1st_network`)) %>% 
    # Draw samples with replacement for every region, RSN, and Age_group
    group_by(Age_group, `1st_network`, Region, .add = TRUE) %>%
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
  
  # 1000 samples for each Age_group*1st_network*Region combination
  # So that for each region for each RSN, and within each Age_group, 1000 resamples of the same size are taken
  resamples <- rbindlist(list_boot)
  
  delta_proportion_boot_a <- resamples %>%
    group_by(n_data_boot, n_data_split, `1st_network`, Region, Age_group) %>%
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup()
  
  # Multiplicative replacement at the regional and Age_group level
  replacement_a <- delta_proportion_boot_a %>%
    group_by(n_data_boot, .add = TRUE) %>%
    group_split()
  
  list_replacement <- list()
  for (i in 1:length(replacement_a)) {
    list_replacement[[i]] <- multRepl(rbindlist(replacement_a[i]) %>%
                                        dplyr::select(Connector:Satellite), dl = rep(1, 4), label = 0, frac = 1e-5)
  }
  
  delta_a_detail <<- cbind(
    delta_proportion_boot_a %>%
      dplyr::select(n_data_boot, n_data_split, `1st_network`:Age_group),
    rbindlist(list_replacement)
  ) %>%
    pivot_longer(
      cols = !c("n_data_boot", "n_data_split", "1st_network", "Region", "Age_group"),
      names_to = "Hub_consensus", values_to = "freq"
    ) %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(n_data_boot, `1st_network`, Age_group, Hub_consensus) %>%
    summarize_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(n_data_boot, `1st_network`, Hub_consensus) %>%
    # Take the log ratios of geometric means
    mutate(delta_freq = log(freq / lead(freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  summary_bootstrap_a <<- delta_a_detail %>%
    group_by(`1st_network`, Hub_consensus) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Hub_consensus, mean, mean_025, mean_0975)
  
  data_boot <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    # Make sure comparisons with missing functional roles can be achieved
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    ungroup() %>%
    filter(grepl(divergent_RSN_interareal, `1st_network`)) %>% 
    group_by(Age_group, `1st_network`, Region, .add = TRUE) %>%
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
    group_by(n_data_boot, n_data_split, `1st_network`, Region, Age_group) %>%
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
      dplyr::select(n_data_boot, n_data_split, `1st_network`:Age_group),
    rbindlist(list_replacement)
  ) %>%
    pivot_longer(
      cols = !c("n_data_boot", "n_data_split", "1st_network", "Region", "Age_group"),
      names_to = "Bridgeness", values_to = "freq"
    ) %>%
    # Compute the difference in proportion of a given functional role within each RSN
    group_by(n_data_boot, `1st_network`, Age_group, Bridgeness) %>%
    summarize_at(vars(freq), funs(geometricmean(.))) %>% 
    group_by(n_data_boot, `1st_network`, Bridgeness) %>%
    # Take the log ratios of geometric means
    mutate(delta_freq = log(freq / lead(freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()
  
  summary_bootstrap_b <<- delta_b_detail %>%
    group_by(`1st_network`, Bridgeness) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Bridgeness, mean, mean_025, mean_0975)
}

bootstrap_ci(n_boot = 100, "Young", "Old",
             divergent_RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN",
             divergent_RSN_interareal = "Auditory|CON|DMN|FPN|Language|SMN")

write.csv(summary_bootstrap_a, "summary_bootstrap_modular.csv")
write.csv(summary_bootstrap_b, "summary_bootstrap_interareal.csv")


data_a <- read.csv("summary_bootstrap_modular.csv")

library(ggpubr)

ggplot(data_a, aes(x = Hub_consensus, y = mean)) +
  geom_point(size = 2) + # Plot the individual points
  geom_errorbar(aes(ymin = mean_025, ymax = mean_0975), width = 0.2) + # Plot the confidence interval
  geom_line(aes(y = mean), color = "red", size = 1) + # Plot the mean
  labs(x = "Modular functional roles", y = "Log-ratio of the geometric means between clusters") + # Label the axes
  facet_wrap(~X1st_network, ncol = 3) +
  geom_hline(yintercept = 0, color = "red") +
  theme_pubclean()


data_b <- read.csv("summary_bootstrap_interareal.csv")

ggplot(data_b, aes(x = Bridgeness, y = mean)) +
  geom_point() + # Plot the individual points
  geom_errorbar(aes(ymin = mean_025, ymax = mean_0975), width = 0.2) + # Plot the confidence interval
  geom_line(aes(y = mean), color = "red", size = 1) + # Plot the mean
  labs(x = "Interareal functional roles", y = "Log-ratio of the geometric means between clusters") + # Label the axes
  facet_wrap(~X1st_network, ncol = 3) +
  geom_hline(yintercept = 0, color = "red") +
  theme_pubclean()
