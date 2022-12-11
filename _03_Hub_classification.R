##########################################################################################
# Script for statistical analyses with the data processed from _DataManip_01_bin.R

# Original hubness score
# Written by CG
# 26-11-2022
##########################################################################################
library(tidyverse)
library(janitor)
library(plyr)
library(readxl)
library(lmerTest)
library(car)
library(Rmisc)
library(tidylog)
library(lme4)
library(ggpubr)
library(broom)
library(purrr)
library(effectsize)
library(jsonlite) # for working with json files
library(emmeans) # for post-hoc tests
library(factoextra)
library(data.table) # for working with lists
library(FactoMineR) # for PCA

rm(list = ls())
options(max.print = 99999)

################################################################################
# Import processed data---------------------------------------------------------
source("_01_DataManip.R")
# Import local functions
source("_radarplotting_function.R")
source("_NMI&AMI_functions.R")

################################################################################
# ~~~~~~~~~~~ Hub classification ~~~~~~~~~~~ -----------------------------------

# High Betweenness centrality = global bridge
# High Flow centrality = local bridge
# High Participation coefficient (based on consensus group-level modular decomposition after 1000 iterations)
# High Within-z (based on consensus group-level modular decomposition after 1000 iterations)

# High_zPC/High z = connector
# High_zPC/low_z = satellite
# Low_zPC/High_z = provincial
# Low_zPC/Low_z = peripheral

# 72 Subjects, 131 Regions, one threshold = 0.15

# !!!!!!!!!! PC & Wz have been generated using GraphVar consensus affiliation vector !!!!!!!!!!!!!!!

PC_consensus <- as.data.frame(fromJSON("Participation_coefficient_consensus.json")) %>%
  mutate(Subj_ID = rep(seq_len(72))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "PC_cons"
  )

Within_module_z_consensus <- as.data.frame(fromJSON("Within_module_z_score_consensus.json")) %>%
  mutate(Subj_ID = rep(seq_len(72))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "Within_module_z_cons"
  )

nodal_metrics_cons <- cbind(PC_consensus, Within_module_z_cons = Within_module_z_consensus$Within_module_z_cons) %>%
  dplyr::select(-Region)

# Make sure dataframe is ordered identically to nodal_metrics
data_full_thresholded <- data_full %>%
  subset(threshold == "0.15") %>%
  arrange(Subj_ID, Region)

data_bind_PC_Wz <- cbind(data_full_thresholded,
  PC_cons = nodal_metrics_cons$PC_cons, Within_module_z_cons = nodal_metrics_cons$Within_module_z_cons
)



data_functional_role <- data_bind_PC_Wz %>%
  # Normalizing at the connectomic level
  mutate(zK = as.numeric(scale(degree))) %>%
  mutate(zBT = as.numeric(scale(Betweenness))) %>%
  mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
  mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
    ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
      ifelse(zBT > 0 & zFlow > 0, "Super_Bridge", "None")
    )
  )) %>%
  # Normalizing at the community level with the affiliation vector from consensus clustering
  group_by(Consensus_vector_0.15) %>%
  mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
  mutate(Hub_consensus = ifelse(zPC_cons > 0 & Within_module_z_cons > 0, "Connector",
    ifelse(zPC_cons > 0 & Within_module_z_cons < 0, "Satellite",
      ifelse(zPC_cons < 0 & Within_module_z_cons > 0, "Provincial",
        ifelse(zPC_cons < 0 & Within_module_z_cons < 0, "Peripheral",
          "None"
        )
      )
    )
  )) %>%
  arrange(Subj_ID, Region) %>%
  ungroup()


################################################################################
# # Method 2 - For PCA-weighted mean ---------------------------
# cov <- data_hub_1_thresholded %>%
#   group_by(Region) %>%
#   summarize_at(vars(zK, Within_module_z, zPC, zBT, zFlow), mean) %>%
#   dplyr::select(-Region)
#
# # with the square cosine from the PCA for the PCA-weighted mean
# contrib <- (PCA(cov)$var$contrib) / 100
# dim1 <- (PCA(cov)$eig / 100)[1, 3]
# dim2 <- (PCA(cov)$eig / 100)[2, 3]
# dim3 <- (PCA(cov)$eig / 100)[3, 3]
#
# # Hubness score for each subject for each region
# optim_weight_ind <- data_hub_1_thresholded %>%
#   dplyr::select(Subj_ID, Region, zK, Within_module_z, zPC, zBT, zFlow) %>%
#   group_by(Subj_ID, Region) %>%
#   mutate(Hubness_score = (zK * (dim1 * contrib[1, 1] + dim2 * contrib[1, 2] + dim3 * contrib[1, 3]) +
#                             Within_module_z * (dim1 * contrib[2, 1] + dim2 * contrib[2, 2] + dim3 * contrib[2, 3]) +
#                             zPC * (dim1 * contrib[3, 1] + dim2 * contrib[3, 2] + dim3 * contrib[3, 3]) +
#                             zBT * (dim1 * contrib[4, 1] + dim2 * contrib[4, 2] + dim3 * contrib[4, 3]) +
#                             zFlow * (dim1 * contrib[5, 1] + dim2 * contrib[5, 2] + dim3 * contrib[5, 3])) / 5) %>%
#   ungroup() %>%
#   dplyr::select(-Community_vector_0.2) %>%
#   group_by(Subj_ID, .add = TRUE) %>%
#   group_split() %>%
#   map_dfr(. %>%
#             slice_max(Hubness_score, n = 131) %>% arrange(desc(Hubness_score)) %>%
#             mutate(rank = rep(seq(1:length(Region)))))
#
#
# data_hub_selection_131 <- merge(data_hub_1_thresholded, optim_weight_ind %>%
#                                   dplyr::select(Subj_ID, Region, Hubness_score, rank), by = c("Subj_ID", "Region"))
#
#
# # Selection top 20% regions as hubs
# top <- 131 * 0.3
#
# optim_weight_group <- data_hub_1_thresholded %>%
#   group_by(Region) %>%
#   summarize_at(vars(zK, Within_module_z, zPC, zBT, zFlow), mean) %>%
#   mutate(Hubness_score = (zK * (dim1 * contrib[1, 1] + dim2 * contrib[1, 2] + dim3 * contrib[1, 3]) +
#                             Within_module_z * (dim1 * contrib[2, 1] + dim2 * contrib[2, 2] + dim3 * contrib[2, 3]) +
#                             zPC * (dim1 * contrib[3, 1] + dim2 * contrib[3, 2] + dim3 * contrib[3, 3]) +
#                             zBT * (dim1 * contrib[4, 1] + dim2 * contrib[4, 2] + dim3 * contrib[4, 3]) +
#                             zFlow * (dim1 * contrib[5, 1] + dim2 * contrib[5, 2] + dim3 * contrib[5, 3])) / 5) %>%
#   slice_max(Hubness_score, n = top) %>%
#   arrange(desc(Hubness_score)) %>%
#   mutate(rank = rep(seq(1:length(Region))))
#
# data_hub_selection <- data_hub_selection_131 %>% filter(Region %in% optim_weight_group$Region)
#
# test <- data_hub_selection %>%
#   group_by(Subj_ID, Region) %>%
#   summarize_at(vars(zK, Within_module_z, zPC, zBT, zFlow, Hubness_score), mean)
#
# cor.test(test$Hubness_score, test$zK)$estimate
# cor.test(test$Hubness_score, test$Within_module_z)$estimate
# cor.test(test$Hubness_score, test$zPC)$estimate
# cor.test(test$Hubness_score, test$zBT)$estimate
# cor.test(test$Hubness_score, test$zFlow)$estimate
#
# data_mod <- data_hub_selection %>%
#   group_by(Subj_ID, Age) %>%
#   summarize_at(vars(zK, Within_module_z, zPC, zBT, zFlow, Hubness_score), mean)
#
# ggplot(data_mod, aes(Age, zK)) +
#   geom_smooth() +
#   theme_pubclean()
#
#
# mod <- lm(zK ~ Age, data_mod)
# summary(mod)
# # performance::check_model(mod)
# performance::check_normality(mod)
# performance::check_outliers(mod, threshold = list("zscore" = 3), method = "zscore")
#
# library(flexplot)
# visualize(mod)
#
# library(lmerTest)
# mod2 <- lmer(Hubness_score~Gender + (1|Subj_ID) + (1|Region), data_hub_selection_131)
# summary(mod2)
################################################################################
