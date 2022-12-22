rm(list=ls())

source("_05_Age_effect_post_clustering.R")
library(philentropy)

################################################################################
# For resting state networks
################################################################################



################################################################################
# This looks at the variability between regions within each RSN
# For each cluster and each RSN, the generalized Jensen-Shannon divergence is computed which reflects
# the amount of variability between regions penalized by within-cluster variability
# High gJSD indicates high between region variability and low within-cluster variability

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
  
  modular_gJSD_1 <- tmp_cluster_final %>%
    group_by(`1st_network`, Region, cluster, Subj_ID, Hub_consensus) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Hub_consensus, freq) %>%
    mutate_all(., ~replace(., is.na(.), 0)) %>% 
    group_by(`1st_network`, Region, cluster) %>%
    summarise_at(vars(Connector:Satellite), mean) %>% 
    ungroup() %>%
    arrange(cluster, `1st_network`)
  
  modular_gJSD_2 <- cbind(
    modular_gJSD_1 %>% dplyr::select(`1st_network`:cluster),
    multRepl(modular_gJSD_1 %>% dplyr::select(Connector:Satellite),
             label = 0, dl = rep(1, 5), frac = 1e-5
    )
  ) %>%
    ungroup() %>%
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
}

gJSD_heuristic("23", "56")


final_heuristic <<- modular_gJSD_3 %>% group_by(`1st_network`) %>% 
  mutate(between_cluster_var = Generalized_Jensen_Shannon_divergence / lag(Generalized_Jensen_Shannon_divergence)) %>% 
  na.omit()

final_heuristic_bis <<- interareal_gJSD_3 %>% group_by(`1st_network`) %>% 
  mutate(between_cluster_var = Generalized_Jensen_Shannon_divergence / lag(Generalized_Jensen_Shannon_divergence)) %>% 
  na.omit()

# This function returns the median relative relative entropy of the PMF between clusters
# High KLD indicates high between cluster variability and low within-cluster variability

# The function also returns two dataframes, 
# one for each PMF that is used later when looking at the topological trajectory

# The function also returns a boolean vector for each PMF indicating the RSN to be kept when computing the log-ratios
# That is, the RSN with the most amount of reconfiguration

KL_heuristic <- function(norm = NULL) {
  
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
  
  modular_KL_2 <<- cbind(
    modular_KL_1 %>% dplyr::select(`1st_network`:cluster),
    # Replace essential zeros with 1e-5
    multRepl(modular_KL_1 %>% dplyr::select(Connector:Satellite),
             label = 0, dl = rep(1, 5), frac = 1e-5
    )
  )  %>% 
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
    Radar_RSN_modular <- modular_KL_3 %>%
    merge(., final_heuristic %>% dplyr::select(`1st_network`, between_cluster_var), by = "1st_network") %>% 
    mutate(norm_KL = Kullback_Leibler_divergence/between_cluster_var) %>% 
    group_by(`1st_network`) %>%
    summarize_at(vars(norm_KL), mean) %>%
    spread(`1st_network`, norm_KL)
  
  Radar_RSN_interareal <- interareal_KL_3 %>%
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
  
  # Auditory   CON  DAN   DMN   FPN Language   PMM  SMN Visual_1 Visual_2   VMM
  choice_modular <<- scale(as.numeric(Radar_RSN_modular)) >= 0
  choice_interareal <<- scale(as.numeric(Radar_RSN_interareal)) >= 0
  } else {
    Radar_RSN_modular <- modular_KL_3 %>%
      group_by(`1st_network`) %>%
      summarize_at(vars(Kullback_Leibler_divergence), mean) %>%
      spread(`1st_network`, Kullback_Leibler_divergence)
    
    Radar_RSN_interareal <- interareal_KL_3 %>%
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
    
    # Auditory   CON  DAN   DMN   FPN Language   PMM  SMN Visual_1 Visual_2   VMM
    choice_modular <<- scale(as.numeric(Radar_RSN_modular)) >= 0
    choice_interareal <<- scale(as.numeric(Radar_RSN_interareal)) >= 0
  }
}

# Should the KLD for each region be normalized by the difference between clusters in the between-region variability for each RSN?
KL_heuristic(norm = TRUE)

choice_modular
choice_interareal

# PDF Visualisation, it's like the log-ratio radar plot but with added information of whether individuals do agree

# agreement_pdf <- agreement_KL %>%
#   arrange(Region, cluster, Hub_consensus) %>%
#   filter(Hub_consensus != "Isolate") %>%
#   filter(grepl("DMN|FPN|Language", `1st_network`))
# filter(grepl("Auditory|SMN", `1st_network`))
# 
# # ridge plot
# ggplot(agreement_pdf, aes(x = freq, y = forcats::fct_rev(`1st_network`), fill = cluster, height = ..density..)) +
#   ggridges::geom_density_ridges(scale = 1, show.legend = TRUE, alpha = 0.7, stat = "density") +
#   scale_x_continuous(name = "Relative proportion within the topologico-functional profile", limits = c(0, 1), labels = scales::percent) +
#   # control space at top and bottom of plot
#   scale_y_discrete(name = "") +
#   scale_fill_viridis_d(option = "D") + # colourblind-safe colours
#   facet_wrap(~Hub_consensus) +
#   theme_pubclean() +
#   ggtitle("Within and between cluster divergence as a function of the relative proportion of each functional role")

# The variance of the density functions represents the degree of within-cluster divergence
# the assigned probability is the proportion of functional role within the composition
# A PDF with high kurtosis that individuals within the cluster agree a lot of the proportion of the functional role
# For example, individuals within cluster 56 tend to agree a lot that Provincial hubs represent 0 to 10% of their TFP
# whether there's is much higher disagreement within cluster 23






