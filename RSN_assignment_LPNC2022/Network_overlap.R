# Script for formatting output from network overlap calculator
# CG - November 2022

################################################################################
################################################################################

library(tidyverse)
library(jsonlite)
library(plyr)
library(data.table)
library(janitor)

rm(list = ls())

# Read file from JSON
overlap <- as.data.frame(fromJSON("Overlap.json")) # JSON for the liberal CAB-NP
overlap_tight <- as.data.frame(fromJSON("Overlap_tight.json")) # JSON for the tight CAB-NP


# Import Power labels and InLang data
Label <- read.delim("path/to/Power_264.txt", header = F) %>%
  plyr::rename(c("V1" = "Index.Power264")) %>%
  plyr::rename(c("V2" = "Region"))

# DO NOT MODIFY - ORIGINAL DATA PROVIDED ONLY (Roger et al., 2022)
InLang <- readxl::read_xlsx("path/to/LANG_131_FC.xlsx", col_names = T) %>%
  janitor::row_to_names(1) %>%
  select(Index.Power264, CAB_NP_assign, LANG_Net_assign)


# Data manipulation------------------------------------------------------------
# Change accordingly
version <- overlap


# Rbind every two columns------------------------------------------------------
lst <- split.default(version, cumsum(rep(c(TRUE, FALSE), ncol(overlap) / 2)))
data <- rbindlist(setNames(lst, seq_along(lst)), idcol = "Index.Power264", use.names = FALSE) %>%
  mutate_at("Index.Power264", funs(as.numeric(as.character(.)))) %>%
  plyr::rename(c("V1" = "Network")) %>%
  plyr::rename(c("V2" = "voxel_proportion"))

data_wider <- data %>%
  # Create a unique identifier row to avoid identification problem w/ pivot_wider
  group_by(Index.Power264, Network) %>%
  dplyr::mutate(row = row_number()) %>%
  pivot_wider(
    names_from = "Network",
    values_from = "voxel_proportion"
  ) %>%
  subset(row == "1") %>%
  dplyr::select(-c(row, `0`)) %>%
  replace(is.na(.), 0) %>%
  plyr::rename(c("1" = "Visual_1")) %>%
  plyr::rename(c("2" = "Visual_2")) %>%
  plyr::rename(c("3" = "SMN")) %>%
  plyr::rename(c("4" = "CON")) %>%
  plyr::rename(c("5" = "DAN")) %>%
  plyr::rename(c("6" = "Language")) %>%
  plyr::rename(c("7" = "FPN")) %>%
  plyr::rename(c("8" = "Auditory")) %>%
  plyr::rename(c("9" = "DMN")) %>%
  plyr::rename(c("10" = "PMM")) %>%
  plyr::rename(c("11" = "VMM")) %>%
  plyr::rename(c("12" = "Orbito-affective"))

# Rearrange data in terms of Primary to last network by region by voxel_proportion
list_dfs <- list()

for (i in 1:nrow(data_wider)) {
  df <- data.frame(data_wider[i, ])
  # Get the network label with the highest value for each region (rowise)
  df$`1st_network` <- colnames(df[, 2:12])[max.col(df[, 2:12])]
  df$`1st_value` <- apply(df[, 2:12], 1, max)
  df <- df %>%
    select(-(unique(df$`1st_network`)))

  df$`2nd_network` <- colnames(df[, 2:11])[max.col(df[, 2:11])]
  df$`2nd_value` <- apply(df[, 2:11], 1, max)
  df <- df %>%
    select(-(unique(df$`2nd_network`)))

  df$`3rd_network` <- colnames(df[, 2:10])[max.col(df[, 2:10])]
  df$`3rd_value` <- apply(df[, 2:10], 1, max)
  df <- df %>%
    select(-(unique(df$`3rd_network`)))

  df$`4th_network` <- colnames(df[, 2:9])[max.col(df[, 2:9])]
  df$`4th_value` <- apply(df[, 2:9], 1, max)
  df <- df %>%
    select(-(unique(df$`4th_network`)))

  df$`5th_network` <- colnames(df[, 2:8])[max.col(df[, 2:8])]
  df$`5th_value` <- apply(df[, 2:8], 1, max)
  df <- df %>%
    select(-(unique(df$`5th_network`)))

  df$`6th_network` <- colnames(df[, 2:7])[max.col(df[, 2:7])]
  df$`6th_value` <- apply(df[, 2:7], 1, max)
  df <- df %>%
    select(-(unique(df$`6th_network`)))

  df$`7th_network` <- colnames(df[, 2:7])[max.col(df[, 2:7])]
  df$`7th_value` <- apply(df[, 2:7], 1, max)
  df <- df %>%
    select(-(unique(df$`7th_network`)))

  df$`8th_network` <- colnames(df[, 2:6])[max.col(df[, 2:6])]
  df$`8th_value` <- apply(df[, 2:6], 1, max)
  df <- df %>%
    select(-(unique(df$`8th_network`)))

  df$`9th_network` <- colnames(df[, 2:5])[max.col(df[, 2:5])]
  df$`9th_value` <- apply(df[, 2:5], 1, max)
  df <- df %>%
    select(-(unique(df$`9th_network`)))

  df$`10th_network` <- colnames(df[, 2:4])[max.col(df[, 2:4])]
  df$`10th_value` <- apply(df[, 2:4], 1, max)
  df <- df %>%
    select(-(unique(df$`10th_network`)))

  df$`11th_network` <- colnames(df[, 2:3])[max.col(df[, 2:3])]
  df$`11th_value` <- apply(df[, 2:3], 1, max)
  df <- df %>%
    select(-(unique(df$`11th_network`)))

  df$`12th_network` <- colnames(df)[2]
  colnames(df)[2] <- "12th_value"
  df <- df %>%
    relocate(`12th_value`, .after = "12th_network")


  list_dfs[[i]] <- df
}

Network_overlap <- rbindlist(list_dfs, use.names = FALSE)
Network_overlap <- merge(Label, Network_overlap, by.y = "Index.Power264")

Network_overlap_InLang <- merge(InLang, Network_overlap, by.y = "Index.Power264") %>%
  arrange(as.numeric(Index.Power264)) %>%
  relocate(Region, .after = Index.Power264)

# library(writexl)
# write_xlsx(Network_overlap_InLang, "Network_overlap_InLang.xlsx")
