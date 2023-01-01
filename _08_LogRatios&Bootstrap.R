################################################################################
# LOG-RATIO of the Topologico-functional profile at the RSN level --------------
################################################################################
source("_07_Entropy.R")
source("_geometricmeanCruz.R")

# Ratio of the proportion of each functional role within each RSN between the two selected clusters
interaction_Age_FuncRole_RSN <- function(cluster1, cluster2, max, min, max2, min2, alpha, RSN_modular, RSN_interareal) {
  data_cluster_selection(cluster1, cluster2)

  ##############################################################################
  # For modular roles
  ##############################################################################

  # Proportion of each functional role within each RSN for each subject
  delta_proportion_a <<- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup() %>%
    pivot_longer(
      cols = !c("1st_network", "Subj_ID", "Age_group"),
      names_to = "Hub_consensus", values_to = "freq"
    ) %>%
    group_by(`1st_network`, Age_group, Hub_consensus) %>%
    summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = 1e-3))) %>%
    group_by(`1st_network`, Hub_consensus) %>%
    # Compute the log ratios
    mutate(delta_freq = log(freq / lead(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()

  if (RSN_modular == "All") {
    Radar_functional_role_RSN_delta_a <<-
      delta_proportion_a %>%
      dplyr::select(`1st_network`, Hub_consensus, delta_freq) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Hub_consensus")
  } else {
    Radar_functional_role_RSN_delta_a <<-
      delta_proportion_a %>%
      dplyr::select(`1st_network`, Hub_consensus, delta_freq) %>%
      filter(grepl(RSN_modular, `1st_network`)) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Hub_consensus")
  }

  radarplotting_overlap(Radar_functional_role_RSN_delta_a, max, min, 1, 1,
    alpha = alpha, label_size = 1,
    title_fill = paste("Log ratio of geometric mean proportions of modular roles"),
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_a), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )



  ##############################################################################
  # For interareal roles
  ##############################################################################

  delta_proportion_b <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup() %>%
    pivot_longer(cols = !c("1st_network", "Subj_ID", "Age_group"), names_to = "Bridgeness", values_to = "freq") %>%
    group_by(`1st_network`, Age_group, Bridgeness) %>%
    summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = 1e-3))) %>%
    group_by(`1st_network`, Bridgeness) %>%
    # Compute the log ratios
    mutate(delta_freq = log(freq / lead(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()

  if (RSN_interareal == "All") {
    Radar_functional_role_RSN_delta_b <<-
      delta_proportion_b %>%
      dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Bridgeness")
  } else {
    Radar_functional_role_RSN_delta_b <<-
      delta_proportion_b %>%
      dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
      filter(grepl(RSN_interareal, `1st_network`)) %>%
      spread(`1st_network`, delta_freq) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Bridgeness")
  }


  radarplotting_overlap(Radar_functional_role_RSN_delta_b, max2, min2, 1, 1,
    alpha = alpha, label_size = 1,
    title_fill = paste("Log ratios of geometric mean proportions of interareal roles"),
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", legend = rownames(Radar_functional_role_RSN_delta_b), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}


# Choose the clusters
# Max and min value on radar plot grid
# mfrow settings
# alpha level
# RSN to be displayed -- either "All" or specify the RSN with | for separation

interaction_Age_FuncRole_RSN("Young", "Old", 2, -2, 2, -2, 0.1,
  RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN",
  RSN_interareal = "Auditory|CON|DMN|FPN|Language|SMN"
)




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
set.seed(123)
bootstrap_ci <- function(n_boot, cluster1, cluster2, RSN_modular, RSN_interareal) {
  data_cluster_selection(cluster1, cluster2)

  data_boot <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Hub_consensus, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Connector:Satellite), mean) %>%
    ungroup() %>%
    filter(grepl(RSN_modular, `1st_network`)) %>%
    # Draw samples with replacement for every RSN, and Age_group
    group_by(Age_group, `1st_network`, .add = TRUE) %>%
    group_split()

  list_boot <- list()
  for (i in 1:length(data_boot)) {
    tmp <- rbindlist(data_boot[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    delta_boot_a <- rbindlist(list_boot_bis) %>%
      as.data.frame() %>%
      pivot_longer(cols = !c("n_data_boot", "n_data_split", "1st_network", "Subj_ID", "Age_group"), names_to = "Hub_consensus", values_to = "freq") %>%
      group_by(n_data_boot, n_data_split, `1st_network`, Age_group, Hub_consensus) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = 1e-3)))

    list_boot[[i]] <- delta_boot_a
  }

  # 1000 samples for each Age_group*1st_network*Region combination
  # So that for each <- <- <- <-  region for each RSN, and within each Age_group, 1000 resamples of the same size are taken

  resamples <- rbindlist(list_boot) %>% as.data.frame()

  log_ratio_boot_a <- resamples %>%
    group_by(n_data_boot, `1st_network`, Hub_consensus) %>%
    # Compute the log ratios
    mutate(delta_freq = log(freq / lead(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()

  summary_bootstrap_a <<- log_ratio_boot_a %>%
    group_by(`1st_network`, Hub_consensus) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Hub_consensus, mean, mean_025, mean_0975)

  write.csv(summary_bootstrap_a, "summary_bootstrap_modular.csv")

  ##############################################################################
  ##############################################################################

  data_boot_bis <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, Subj_ID, Age_group, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Subj_ID, Age_group) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>%
    ungroup() %>%
    filter(grepl(RSN_interareal, `1st_network`)) %>%
    # Draw samples with replacement for every RSN, and Age_group
    group_by(Age_group, `1st_network`, .add = TRUE) %>%
    group_split()

  list_boot <- list()
  for (i in 1:length(data_boot_bis)) {
    tmp <- rbindlist(data_boot_bis[i]) %>% mutate(n_data_split = rep(i, times = length(nrow(.))))
    tmp_bis <- lapply(1:n_boot, function(i) tmp[sample(nrow(tmp), replace = TRUE), ])
    list_boot_bis <- list()
    for (j in 1:length(tmp_bis)) {
      list_boot_bis[[j]] <- rbindlist(tmp_bis[j]) %>% mutate(n_data_boot = rep(j, times = length(nrow(.))))
    }
    delta_boot_b <<- rbindlist(list_boot_bis) %>%
      as.data.frame() %>%
      pivot_longer(
        cols = !c("n_data_boot", "n_data_split", "1st_network", "Subj_ID", "Age_group"),
        names_to = "Bridgeness",
        values_to = "freq"
      ) %>%
      group_by(n_data_boot, n_data_split, `1st_network`, Age_group, Bridgeness) %>%
      summarise_at(vars(freq), funs(geomMeanExtension(., epsilon = 1e-3)))

    list_boot[[i]] <- delta_boot_b
  }

  # 1000 samples for each Age_group*1st_network*Region combination
  # So that for each <- <- <- <-  region for each RSN, and within each Age_group, 1000 resamples of the same size are taken

  resamples <- rbindlist(list_boot) %>% as.data.frame()

  log_ratio_boot_b <- resamples %>%
    group_by(n_data_boot, `1st_network`, Bridgeness) %>%
    mutate(delta_freq = log(freq / lead(freq))) %>%
    filter(!(grepl("NaN|-Inf", delta_freq))) %>%
    dplyr::select(-c(freq, Age_group)) %>%
    na.omit()

  summary_bootstrap_b <<- log_ratio_boot_b %>%
    group_by(`1st_network`, Bridgeness) %>%
    get_summary_stats(delta_freq, type = "full") %>%
    mutate(mean_025 = mean - ci) %>%
    mutate(mean_0975 = mean + ci) %>%
    dplyr::select(`1st_network`, Bridgeness, mean, mean_025, mean_0975)

  write.csv(summary_bootstrap_b, "summary_bootstrap_interareal.csv")
  ##############################################################################
  ##############################################################################
}

bootstrap_ci(
  n_boot = 1000, "Young", "Old",
  RSN_modular = "Auditory|CON|DMN|FPN|Language|SMN",
  RSN_interareal = "Auditory|CON|DMN|FPN|Language|SMN"
)


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
