rm(list = ls())
source("_07_Entropy.R")

# ******************************************************************************
# Most likely trajector with Conditional PMFs
# ******************************************************************************

trajectory <- function(type_func_df, cluster1, cluster2, list_RSN) {
  # Convert an adjacency dataframe to a 2-column dataframe
  adjacency_to_2col <- function(dataframe) {
    crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
    crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
    crossdatamat <- t(crossdatatmp)
    colnames(crossdatamat) <- c(cluster1, cluster2, "Value")
    crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
    crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
    return(crossdatadf)
  }

  # For each region, find the most probable trajectory
  if (colnames(type_func_df[4]) == "Hub_consensus") {
    Cond_PMF <- type_func_df %>%
      spread(Age_group, freq) %>%
      arrange(Hub_consensus) %>%
      group_by(Region, .add = TRUE) %>%
      group_split()
  } else {
    Cond_PMF <- type_func_df %>%
      spread(Age_group, freq) %>%
      arrange(Bridgeness) %>%
      group_by(Region, .add = TRUE) %>%
      group_split()
  }


  Cond_PMF_list <- list()
  for (i in 1:length(Cond_PMF)) {
    tmp <- rbindlist(Cond_PMF[i])
    # Beware of the order here, first young then old
    data <- outer(tmp[, 5] %>% as.matrix(), tmp[, 4] %>% as.matrix()) %>% as.data.frame()
    if (colnames(type_func_df[4]) == "Hub_consensus") {
      colnames(data) <- c("Connector", "Peripheral", "Provincial", "Satellite")
      rownames(data) <- c("Connector", "Peripheral", "Provincial", "Satellite")
    } else {
      colnames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
      rownames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
    }
    crossdatadf <- adjacency_to_2col(data)
    data_2 <- crossdatadf %>% dplyr::slice_max(Value)
    Cond_PMF_list[[i]] <- cbind(tmp %>% slice(1:nrow(data_2)) %>% dplyr::select(`1st_network`, Region), data_2)
  }

  Cond_PMF_final <<- rbindlist(Cond_PMF_list) %>%
    # Identify each region with a unique label
    mutate(helper_vector = rep(seq(nrow(.)))) %>%
    unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
    # Transform  back to an alluvial format
    pivot_longer(
      cols = c(cluster1, cluster2),
      names_to = "Age_group",
      values_to = colnames(type_func_df[4])
    )

  library(ggalluvial)
  if (list_RSN != "All") {
    if (colnames(type_func_df[4]) == "Hub_consensus") {
      display_cluster <- Cond_PMF_final %>%
        filter(grepl(list_RSN, `1st_network`)) %>%
        group_by(Age_group, Hub_consensus) %>%
        summarize(s = n()) %>%
        arrange(desc(Age_group), desc(Hub_consensus)) %>%
        .$Hub_consensus
      
      alluvial_cluster <- ggplot(
        Cond_PMF_final %>%
          filter(grepl(list_RSN, `1st_network`)),
        aes(x = forcats::fct_rev(Age_group), stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
      ) +
        geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
        geom_stratum(alpha = .8) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
        geom_text(stat = "stratum", label = display_cluster) +
        scale_fill_brewer(palette = "Oranges", direction = 1) +
        labs(title = "Most probable topological reconfiguration trajectory for each region between age clusters") +
        labs(x = "Age group") +
        theme_pubclean()
      
      alluvial_cluster
    } else {
      display_cluster <- Cond_PMF_final %>%
        filter(grepl(list_RSN, `1st_network`)) %>%
        group_by(Age_group, Bridgeness) %>%
        summarize(s = n()) %>%
        arrange(desc(Age_group), desc(Bridgeness)) %>%
        .$Bridgeness
      
      alluvial_cluster <- ggplot(
        Cond_PMF_final %>%
          filter(grepl(list_RSN, `1st_network`)),
        aes(x = forcats::fct_rev(Age_group), stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
      ) +
        geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
        geom_stratum(alpha = .8) +
        scale_x_discrete(expand = c(.1, .1)) +
        # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
        geom_text(stat = "stratum", label = display_cluster) +
        scale_fill_brewer(palette = "Oranges", direction = 1) +
        labs(title = "Most probable topological reconfiguration trajectory for each region between age clusters") +
        labs(x = "Age group") +
        theme_pubclean()
      
      alluvial_cluster
    }
  } else {
      if (colnames(type_func_df[4]) == "Hub_consensus") {
        display_cluster <- Cond_PMF_final %>%
          group_by(Age_group, Hub_consensus) %>%
          summarize(s = n()) %>%
          arrange(desc(Age_group), desc(Hub_consensus)) %>%
          .$Hub_consensus
        
        alluvial_cluster <- ggplot(
          Cond_PMF_final,
          aes(x = forcats::fct_rev(Age_group), stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
        ) +
          geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
          geom_stratum(alpha = .8) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
          geom_text(stat = "stratum", label = display_cluster) +
          scale_fill_brewer(palette = "Oranges", direction = 1) +
          labs(title = "Most probable topological reconfiguration trajectory for each region between age clusters") +
          labs(x = "Age group") +
          theme_pubclean()
        
        alluvial_cluster
      } else {
        display_cluster <- Cond_PMF_final %>%
          group_by(Age_group, Bridgeness) %>%
          summarize(s = n()) %>%
          arrange(desc(Age_group), desc(Bridgeness)) %>%
          .$Bridgeness
        
        alluvial_cluster <- ggplot(
          Cond_PMF_final,
          aes(x = forcats::fct_rev(Age_group), stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
        ) +
          geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
          geom_stratum(alpha = .8) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
          geom_text(stat = "stratum", label = display_cluster) +
          scale_fill_brewer(palette = "Oranges", direction = 1) +
          labs(title = "Most probable topological reconfiguration trajectory for each region between age clusters") +
          labs(x = "Age group") +
          theme_pubclean()
        
        alluvial_cluster
      }
    }
}

# Pick the appropriate dataframe - modular_KL_2 or interareal_KL_2
# Pick the desired RSN - to be specified in a grepl format e.g., "DMN|FPN|Language"

trajectory(interareal_KL_2, "Young", "Old", "All")
