rm(list=ls())

source("_05_Age_effect_post_clustering.R")

library(philentropy)

################################################################################
# For resting state networks
################################################################################

# Evolution of within-RSN variability between clusters
gJSD_heuristic <- function(cluster1, cluster2) {
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(cluster == cluster1 | cluster == cluster2)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID)
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  tmp_cluster_final <- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, cluster),
                             by = "Subj_ID"
  )
  
  modular_gJSD_2 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, cluster, Subj_ID, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Hub_consensus, freq) %>%
    mutate_all(., ~replace(., is.na(.), 0)) %>% 
    group_by(`1st_network`, Region, cluster) %>%
    summarise_at(vars(Connector:Satellite), mean) %>% 
    ungroup() %>%
    arrange(cluster, `1st_network`) %>% 
  
  # modular_gJSD_2 <- cbind(
  #   modular_gJSD_1 %>% dplyr::select(`1st_network`:cluster),
  #   multRepl(modular_gJSD_1 %>% dplyr::select(Connector:Satellite),
  #            label = 0, dl = rep(1, 5), frac = 1e-5
  #   )
  # ) %>%
  #   ungroup() %>%
    pivot_longer(
      cols = !c("1st_network", "Region", "cluster"),
      names_to = "Hub_consensus", values_to = "freq"
    )
  
  modular_gJSD_split <- modular_gJSD_2 %>%
    ungroup() %>%
    group_by(cluster, `1st_network`, .add = TRUE) %>%
    group_split()
  
  gJSD_list <- list()
  for (i in 1:length(modular_gJSD_split)) {
    tmp <- rbindlist(modular_gJSD_split[i]) %>%
      group_by(cluster,`1st_network`, Region, Hub_consensus) %>%
      summarize_at(vars(freq), mean) %>%
      spread(Hub_consensus, freq) %>%
      remove_rownames() %>%
      column_to_rownames("Region") %>% 
      dplyr::select(-c(`1st_network`, cluster))
    gJSD <- philentropy::gJSD(tmp %>% as.matrix())
    gJSD_list[[i]] <- gJSD
  }
  
  modular_gJSD_3 <<- t(rbindlist(list(gJSD_list))) %>% 
    as.data.frame() %>% 
    plyr::rename(c("V1" = "Generalized_Jensen_Shannon_divergence")) %>% 
    cbind(., modular_gJSD_2 %>% group_by(cluster, `1st_network`) %>% 
            summarize_at(vars(freq), mean))
  
  interareal_gJSD_1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, cluster, Subj_ID, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Bridgeness, freq) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, cluster) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>% 
    ungroup() %>%
    arrange(cluster, `1st_network`)
  
  interareal_gJSD_2 <- cbind(
    interareal_gJSD_1 %>% dplyr::select(`1st_network`:cluster),
    multRepl(interareal_gJSD_1 %>% dplyr::select(Global_Bridge:Super_Bridge),
             label = 0, dl = rep(1, 4), frac = 1e-5
    )
  ) %>%
    ungroup() %>%
    pivot_longer(
      cols = !c("1st_network", "Region", "cluster"),
      names_to = "Bridgeness", values_to = "freq"
    ) 
  
  interareal_gJSD_split <- interareal_gJSD_2 %>%
    ungroup() %>%
    group_by(cluster, `1st_network`, .add = TRUE) %>%
    group_split()
  
  gJSD_list_bis <- list()
  for (i in 1:length(interareal_gJSD_split)) {
    tmp <- rbindlist(interareal_gJSD_split[i]) %>%
      group_by(cluster, `1st_network`, Region, Bridgeness) %>%
      summarize_at(vars(freq), mean) %>%
      spread(Bridgeness, freq) %>%
      remove_rownames() %>%
      column_to_rownames("Region") %>% 
      dplyr::select(-c(`1st_network`, cluster))
    gJSD <- philentropy::gJSD(tmp %>% as.matrix())
    gJSD_list_bis[[i]] <- gJSD
  }
  
  interareal_gJSD_3 <<- t(rbindlist(list(gJSD_list_bis))) %>%
    as.data.frame() %>%
    plyr::rename(c("V1" = "Generalized_Jensen_Shannon_divergence")) %>% 
    cbind(., interareal_gJSD_2 %>% group_by(cluster, `1st_network`) %>% 
            summarize_at(vars(freq), mean))
  
  Radar_RSN_modular_gJSD <- modular_gJSD_3 %>% dplyr::select(-freq) %>% 
    spread(`1st_network`, Generalized_Jensen_Shannon_divergence) %>% 
    remove_rownames() %>% column_to_rownames("cluster")
  
  radarplotting_overlap(Radar_RSN_modular_gJSD, 1, 0, 1, 1,
                        alpha = 0.3, label_size = 1, title_fill = "Relative entropy of the modular role probability mass functions (PMF) between regions within each RSN\n (Generalized Jensen-Shannon divergence)",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "topright", title = "Type of functional roles",
    legend = rownames(Radar_RSN_modular_gJSD), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  Radar_RSN_interareal_gJSD <- interareal_gJSD_3 %>% dplyr::select(-freq) %>%
    spread(`1st_network`, Generalized_Jensen_Shannon_divergence) %>% 
    remove_rownames() %>% column_to_rownames("cluster")
  
  radarplotting_overlap(Radar_RSN_interareal_gJSD, 1, 0, 1, 1,
                        alpha = 0.3, label_size = 1, title_fill = "Relative entropy of the interareal role probability mass functions (PMF) between regions within each RSN\n (Generalized Jensen-Shannon divergence)",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "topright", title = "Type of functional roles",
    legend = rownames(Radar_RSN_interareal_gJSD), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  final_heuristic <<- modular_gJSD_3 %>% group_by(`1st_network`) %>% 
    mutate(between_cluster_var = Generalized_Jensen_Shannon_divergence / lag(Generalized_Jensen_Shannon_divergence)) %>% 
    na.omit()
  
  final_heuristic_bis <<- interareal_gJSD_3 %>% group_by(`1st_network`) %>% 
    mutate(between_cluster_var = Generalized_Jensen_Shannon_divergence / lag(Generalized_Jensen_Shannon_divergence)) %>% 
    na.omit()
}
gJSD_heuristic("23.5", "56")

# Strength of between-cluster reconfiguration
KL_heuristic <- function(cluster1, cluster2, norm = NULL) {
  
  # Final dataframe with only the subjects of chosen clusters, their hub regions and the RSNs -----
  
  # Retain only the rows specific of the two clusters
  tmp_cluster_0 <- data_post_clustering %>% subset(cluster == cluster1 | cluster == cluster2)
  # Get the associated Resting-state networks
  tmp_cluster_1 <- filter(data_functional_role, Subj_ID %in% tmp_cluster_0$Subj_ID)
  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject <- rbindlist(Hub_selection)
  # Select the subjects from the clusters
  data_hub_selection_cluster <- filter(
    data_hub_selection_per_subject,
    Subj_ID %in% tmp_cluster_1$Subj_ID
  )
  tmp_cluster_final <- merge(data_hub_selection_cluster, tmp_cluster_0 %>%
                               dplyr::select(Subj_ID, cluster),
                             by = "Subj_ID"
  )
  # For modular-level functional roles (i.e., Connector, Satellite ...)
  modular_KL_1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, cluster, Subj_ID, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>% 
    spread(Hub_consensus, freq) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% 
    group_by(`1st_network`, Region, cluster) %>%
    summarize_at(vars(Connector:Satellite), mean) %>% 
    ungroup() %>% 
    arrange(Region)
  
  # modular_KL_2 <<- cbind(
  #   modular_KL_1 %>% dplyr::select(`1st_network`:cluster),
  #   # Replace essential zeros with 1e-5
  #   multRepl(modular_KL_1 %>% dplyr::select(Connector:Satellite),
  #            label = 0, dl = rep(1, 5), frac = 1e-5
  #   )
  # )  %>% 
  modular_KL_2 <<- modular_KL_1 %>% 
    pivot_longer(
      cols = !c("1st_network", "Region", "cluster"),
      names_to = "Hub_consensus", values_to = "freq"
    )
  
  # Computing the Kullback-Leibler divergence for each region between clusters
  modular_KL_split <- modular_KL_2 %>%
    ungroup() %>%
    group_by(Region, .add = TRUE) %>%
    group_split()
  
  KL_list <- list()
  for (i in 1:length(modular_KL_split)) {
    tmp <- rbindlist(modular_KL_split[i]) %>%
      dplyr::select(cluster, Hub_consensus, freq) %>%
      spread(Hub_consensus, freq) %>%
      remove_rownames() %>%
      column_to_rownames("cluster")
    KL_list[[i]] <- philentropy::KL(tmp %>% as.matrix())
  }
  
  modular_KL_3 <- t(rbindlist(list(KL_list))) %>%
    as.data.frame() %>%
    # First make sure Regions are in the same order as the KL measures
    cbind(., Region = unique(modular_KL_1$Region)) %>% 
    # Add the KL measure per region back to the original dataframe
    merge(., modular_KL_1 %>% dplyr::select(`1st_network`, Region), by = "Region") %>% distinct(Region, .keep_all = TRUE) %>% 
    plyr::rename(c("V1" = "Kullback_Leibler_divergence"))
  
  
  # For interareal-level functional roles (i.e., Global Bridge, Local Bridge ...)
  interareal_KL_1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, cluster, Subj_ID, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Bridgeness, freq) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Region, cluster) %>%
    summarize_at(vars(Global_Bridge:Super_Bridge), mean) %>% 
    ungroup() %>%
    arrange(Region)
  
  interareal_KL_2 <<- cbind(
    interareal_KL_1 %>% dplyr::select(`1st_network`:cluster),
    multRepl(interareal_KL_1 %>% dplyr::select(Global_Bridge:Super_Bridge),
             label = 0, dl = rep(1, 4), frac = 1e-5
    )
  ) %>%
    pivot_longer(
      cols = !c("1st_network", "Region", "cluster"),
      names_to = "Bridgeness", values_to = "freq"
    )
  
  interareal_KL_split <- interareal_KL_2 %>%
    ungroup() %>%
    group_by(Region, .add = TRUE) %>%
    group_split()
  
  KL_list_bis <- list()
  for (i in 1:length(interareal_KL_split)) {
    tmp <- rbindlist(interareal_KL_split[i]) %>%
      dplyr::select(cluster, Bridgeness, freq) %>%
      spread(Bridgeness, freq) %>%
      remove_rownames() %>%
      column_to_rownames("cluster")
    KL_list_bis[[i]] <- philentropy::KL(tmp %>% as.matrix())
  }
  
  interareal_KL_3 <- t(rbindlist(list(KL_list_bis))) %>%
    as.data.frame() %>%
    cbind(., Region = unique(interareal_KL_1$Region)) %>%
    merge(., interareal_KL_1 %>% dplyr::select(`1st_network`, Region), by = "Region") %>% distinct(Region, .keep_all = TRUE) %>% 
    plyr::rename(c("V1" = "Kullback_Leibler_divergence"))
  
  # Two type of compensatory mechanisms
  # Reconfiguration of proportion across modular functional roles
  # Reconfiguration of proportion within each connectomic functional role, meaning that overall proportion
  # within the functional role is preserved, only the regions changed assignment
  
  # To quantify how topological reconfiguration is taking place within each RSN between clusters
  if (norm == TRUE) {
    Radar_RSN_modular <<- modular_KL_3 %>%
    merge(., final_heuristic %>% dplyr::select(`1st_network`, between_cluster_var), by = "1st_network") %>% 
    mutate(norm_KL = Kullback_Leibler_divergence/between_cluster_var) %>% 
    group_by(`1st_network`) %>%
    summarize_at(vars(norm_KL), mean) %>%
    spread(`1st_network`, norm_KL)
  
  Radar_RSN_interareal <<- interareal_KL_3 %>%
    merge(., final_heuristic_bis %>% dplyr::select(`1st_network`, between_cluster_var), by = "1st_network") %>% 
    mutate(norm_KL = Kullback_Leibler_divergence/between_cluster_var) %>% 
    group_by(`1st_network`) %>%
    summarize_at(vars(norm_KL), mean) %>%
    spread(`1st_network`, norm_KL)
  
  Radar_KL <- rbind(Radar_RSN_modular, Radar_RSN_interareal)
  rownames(Radar_KL) <- c("Modular-level\n PMF", 
                          "Interareal-level\n PMF")
  
  radarplotting_overlap(Radar_KL, 3, 0, 1, 1,
                        alpha = 0.3, label_size = 1, title_fill = "Normalized relative entropy of the functional role probability mass functions (PMF) between the two clusters\n Computed for the PMF of each region between clusters and averaged at the RSN level",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "topright", title = "Type of functional roles",
    legend = rownames(Radar_KL), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
  
  } else {
    Radar_RSN_modular <<- modular_KL_3 %>%
      group_by(`1st_network`) %>%
      summarize_at(vars(Kullback_Leibler_divergence), mean) %>%
      spread(`1st_network`, Kullback_Leibler_divergence)
    
    Radar_RSN_interareal <<- interareal_KL_3 %>%
      group_by(`1st_network`) %>%
      summarize_at(vars(Kullback_Leibler_divergence), mean) %>%
      spread(`1st_network`, Kullback_Leibler_divergence)
    
    Radar_KL <- rbind(Radar_RSN_modular, Radar_RSN_interareal)
    rownames(Radar_KL) <- c("Modular-level\n PMF", 
                            "Interareal-level\n PMF")
    
    radarplotting_overlap(Radar_KL, 2, 0, 3, 1,
                          alpha = 0.3, label_size = 1, title_fill = "Unnormalized relative entropy of the functional role probability mass functions (PMF) between the two clusters\n Computed for the PMF of each region between clusters and averaged at the RSN level",
                          palette = RColorBrewer::brewer.pal(8, "Dark2")
    )
    
    legend(
      x = "topright", title = "Type of functional roles",
      legend = rownames(Radar_KL), horiz = TRUE,
      bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
      text.col = "black", cex = 1, pt.cex = 2
    )
  }
}

# If normalized, then the radar plot displays the RSN with the most amount of topological reconfiguration
# combined with high specificity: all regions within the RSN reconfigure in unison towards a specific functional role
KL_heuristic("23.5", "56", norm = TRUE)

# Arbitrary threshold to pick the RSN to investigate the difference in proportions
median(as.numeric(Radar_RSN_modular))
median(as.numeric(Radar_RSN_interareal))


# Consistency computation ----
# A RSN is considered to have a consistent reconfiguration if regions with the same functional role have the same most probable trajectory

# Specificity necessarily entails consistency but the reciprocal is not true:
# If all regions reconfigure towards a same functional role as a whole, this is necessarily true at for every functional role
# Conversely, regions with a same functional role can be very consistent in their reconfiguration but not very specific, that is in unison, with the other consistent trajectories


consistency <- function(type_func_df, cluster1, cluster2) {
  adjacency_to_2col <- function(dataframe) {
    crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
    crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
    crossdatamat <- t(crossdatatmp)
    colnames(crossdatamat) <- c(cluster1, cluster2, "Value")
    crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
    crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
    return(crossdatadf)
  }
  
  # For each region, compute the probability of every possible trajectory
  # 4*4 functional roles = 16 possible trajectories for each region
  if (colnames(type_func_df[4]) == "Hub_consensus") {
    Cond_PMF <- type_func_df %>%
      spread(cluster, freq) %>%
      arrange(Hub_consensus) %>% 
      group_by(Region, .add = TRUE) %>%
      group_split()
  } else {
    Cond_PMF <- type_func_df %>%
      spread(cluster, freq) %>%
      arrange(Bridgeness) %>% 
      group_by(Region, .add = TRUE) %>%
      group_split()
  }
  
  pairwise_prob_list <- list()
  for (i in 1:length(Cond_PMF)) {
    tmp <- rbindlist(Cond_PMF[i])
    data <- outer(tmp[,4] %>% as.matrix(), tmp[,5] %>% as.matrix()) %>% as.data.frame()
    if (colnames(type_func_df[4]) == "Hub_consensus") {
      colnames(data) <- c("Connector", "Peripheral", "Provincial", "Satellite")
      rownames(data) <- c("Connector", "Peripheral", "Provincial", "Satellite")
    } else {
      colnames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
      rownames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
    }
    crossdatadf <- adjacency_to_2col(data) 
    pairwise_prob_list[[i]] <- cbind(tmp %>% slice(1:nrow(crossdatadf)) %>% dplyr::select(`1st_network`, Region), crossdatadf)
  }
  
  
  
  # Remove non-existent trajectories
  pairwise_prob <- rbindlist(pairwise_prob_list) %>% 
    filter(Value != 0)
  
  consistency_1 <- pairwise_prob %>% unite(., trajectory, c(colnames(.[,3]), colnames(.[,4])), remove = FALSE) %>% 
    arrange(`1st_network`, Region, trajectory)
  
# Other approach with gJSD

# Compute each possible trajectories and their probability
# Normalize the probabilities for each RSN for each starting functional role
# Compare region, functional role avec gJSD
# 
  consistency_gJSD_split <- consistency_1 %>%
    ungroup() %>%
    arrange(`1st_network`, Region, trajectory) %>% 
    group_by(`1st_network`, .[,4], .add = TRUE) %>%
    group_split()
  
  gJSD_list <- list()
  for (i in 1:length(consistency_gJSD_split)) {
    tmp <- rbindlist(consistency_gJSD_split[i]) %>%
      group_by(Region) %>% 
      mutate(norm_prob = Value/sum(Value)) %>% 
      dplyr::select(`1st_network`, Region, trajectory, colnames(.[,4]), norm_prob)
    
    tmp_bis <- tmp %>% 
      spread(trajectory, norm_prob) %>% ungroup()
    
    tmp_ter <- tmp_bis %>% 
      dplyr::select(-c(`1st_network`, colnames(.[,3]))) %>%
      remove_rownames() %>%
      column_to_rownames("Region") %>% 
      mutate_all(., ~replace(., is.na(.), 0))
    
    gJSD <- philentropy::gJSD(tmp_ter %>% as.matrix())
    gJSD_list[[i]] <- cbind(tmp_bis %>% slice(1) %>% dplyr::select(`1st_network`, colnames(.[,3])), gJSD)
  }
  
  consistency_gJSD <<- rbindlist(gJSD_list) %>% as.data.frame() %>% 
    plyr::rename(c("gJSD" = "Generalized_Jensen_Shannon_divergence"))
  
  Radar_RSN_consistency_gJSD <- consistency_gJSD %>%  
    group_by(.[1]) %>% summarize_at(vars(everything()), mean) %>% 
    dplyr::select(-colnames(.[,2])) %>% 
    spread(`1st_network`, Generalized_Jensen_Shannon_divergence) 
  
  radarplotting_overlap(Radar_RSN_consistency_gJSD, 1, 0, 1, 1, inverse_grad = TRUE,
                        alpha = 0.5, label_size = 1, 
                        title_fill = "gJSD Consistency",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
}

consistency(modular_KL_2, "23.5", "56")

  
  # # Count the number of similar trajectories
  # consistency_2 <- consistency_1 %>% 
  #   group_by(`1st_network`, trajectory, .[,4]) %>% 
  #   summarize(n = n()) %>% 
  #   mutate(freq = n/sum(n))
  # 
  # consistency_2_bis <- consistency_1 %>% distinct(Region, .keep_all = TRUE) %>%  
  #   group_by(`1st_network`) %>% 
  #   summarize(n_region = n())
  # 
  # consistency_2_merge <- merge(consistency_2, consistency_2_bis, by = "1st_network") %>% 
  #   mutate(freq = n/n_region)
  # 
  # # Normalize this count by the mean probability that this trajectory occurs
  # consistency_3 <- consistency_1 %>%  
  #   group_by(`1st_network`, trajectory, .[,4]) %>% 
  #   summarize_at(vars(Value), mean) %>% 
  #   # normalize the probability of trajectory occurence by each functional role
  #   group_by(`1st_network`, .[,3]) %>% 
  #   mutate(norm_prob = Value/sum(Value)) %>% 
  #   cbind(., freq = consistency_2_merge$freq) %>% 
  #   # Frequency of the trajectory within the RSN multiplied by the probability
  #   mutate(consistency_score = freq*norm_prob) %>% 
  #   # Keep only the highest consistency score for each functionak role
  #   group_by(`1st_network`, .[,3]) %>% slice_max(consistency_score) %>% 
  #   distinct(consistency_score, .keep_all = TRUE)
  #         
  # consistency_final <- consistency_3 %>% 
  #   # Normalize for each starting functional role
  #   group_by(`1st_network`) %>% 
  #   summarise_at(vars(consistency_score), mean)
  # 
  # Radar_RSN_consistency <- consistency_final %>%
  #   spread(`1st_network`, consistency_score) %>% 
  #   mutate_at(vars(everything()), funs(. * 10))
  # 
  # radarplotting_overlap(Radar_RSN_consistency, 5, 1, 1, 1,
  #                       alpha = 0.3, label_size = 1, 
  #                       title_fill = "Consistency score",
  #                       palette = RColorBrewer::brewer.pal(8, "Dark2")
  # )

  



  









