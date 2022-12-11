library(tidyverse)
library(janitor)
library(plyr)
library(readxl)
library(car)
library(Rmisc)
library(tidylog)
library(plotly)
library(jsonlite)


rm(list = ls())
options(max.print = 99999)

# Import processed data---------------------------------------------------------
source("_01_DataManip.R")

################################################################################
# Get mean r per threshold per subject with absolute threshold
# Mean_r <- as.data.frame((fromJSON("Mean_r.json"))) %>%
#   plyr::rename(c("V1" = "0.05")) %>%
#   plyr::rename(c("V2" = "0.07")) %>%
#   plyr::rename(c("V3" = "0.1")) %>%
#   plyr::rename(c("V4" = "0.15")) %>%
#   plyr::rename(c("V5" = "0.2")) %>%
#   plyr::rename(c("V6" = "0.25")) %>%
#   pivot_longer(
#     cols = c("0.05", "0.07", "0.1", "0.15", "0.2", "0.25"),
#     names_to = "threshold",
#     values_to = "Mean_r"
#   ) %>%
#   mutate(Subj_ID = rep(seq_len(72), each = 6))
#
# data_per_subject_combined_bu <- merge(data_per_subject_combined_bu, Mean_r, by = c("Subj_ID", "threshold")) %>%
#   arrange(Subj_ID, threshold)
#
# test_mean_r <- data_per_subject_combined_bu %>% subset(threshold == "0.2")
# cor.test(test_mean_r$Mean_r, test_mean_r$Age)
# # Get proportion of LLC per threshold
# # Connected components generate with consensus based approach, niter = 1000, tau = .5
# LLC <- as.data.frame((fromJSON("consensus_LLC.json"))) %>%
#   pivot_longer(
#     cols = c("V1", "V2", "V3", "V4", "V5", "V6"),
#     names_to = "threshold"
#   ) %>%
#   count(threshold, value) %>%
#   group_by(threshold) %>%
#   mutate(prop = prop.table(n))
#
# threshold <- c(.05, .07, .1, .15, .2, .25)
# LLC_coeff <- c(0.442, 0.557, 1, 1, 1, 1) # Input the largest connected component
# LLC_plot <- cbind(threshold, LLC_coeff)
# data_per_subject_combined_bu <- merge(data_per_subject_combined_bu, LLC_plot, by = "threshold") %>% arrange(Subj_ID, threshold)

# Evolution of Global metrics - what is the optimal threshold?
evo <- data_full_per_subject %>%
  group_by(threshold) %>%
  summarise_at(vars(Eglob, Clustering_coeff_glob, Eloc, Q_consensus), funs(mean)) %>%
  pivot_longer(
    cols = !c("threshold"),
    names_to = "Metrics"
  )


# Choosing threshold = 20% for subsequent analyses
# Number of communities is 3 throughout the thresholds

plotly::ggplotly(
  evo %>%
    ggplot(aes(threshold, value, color = Metrics)) +
    geom_line() +
    geom_point(size = 2) +
    geom_jitter(height = 0.05, alpha = 0.2) +
    ggtitle("Evolution of global metrics as a function of the threshold") +
    xlab("Threshold") +
    ylab("") +
    scale_x_continuous(breaks = c(0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    coord_cartesian(ylim = c(0.2, 0.8)) +
    ggpubr::theme_pubclean()
)


# AMI for all threshold consensus clustering
# source("_NMI&AMI_functions.R")
# # Make sure this is indexed on 131 observations only
# data_AMI <- data_per_region_combined_bu %>% subset(threshold == "0.2")
#
# a <- AMI_func(factor(data_AMI$Community_vector_0.05), factor(data_AMI$Community_vector_0.2))
# b <- AMI_func(factor(data_AMI$Community_vector_0.07), factor(data_AMI$Community_vector_0.2))
# c <- AMI_func(factor(data_AMI$Community_vector_0.1), factor(data_AMI$Community_vector_0.2))
# d <- AMI_func(factor(data_AMI$Community_vector_0.15), factor(data_AMI$Community_vector_0.2))
# e <- AMI_func(factor(data_AMI$Community_vector_0.25), factor(data_AMI$Community_vector_0.2))
# # Mean AMI with partitions at other thresholds
# (a + b + c + d + e)/5
