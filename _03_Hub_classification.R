##########################################################################################
# Script for functional hub classification with the data processed

# Original definition of Bridgeness for topological-functional profiles

# Written by CG
# 26-11-2022
##########################################################################################
library(tidyverse)
library(janitor)
library(plyr)
library(car)
library(Rmisc)
library(tidylog)
library(jsonlite) # for working with json files
library(emmeans) # for post-hoc tests
library(data.table) # for working with lists


rm(list = ls())
################################################################################
# Import processed data---------------------------------------------------------
source("_01_DataManip.R")

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


# Wz is mathematicaaly classified as zero when belonging to its own module

data_functional_role <- data_bind_PC_Wz %>%
  # Normalizing at the connectomic level
  mutate(zK = as.numeric(scale(degree))) %>%
  mutate(zBT = as.numeric(scale(Betweenness))) %>%
  mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
  mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
    ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
      ifelse(zBT > 0 & zFlow > 0, "Super_Bridge", 
             ifelse(zBT < 0 & zFlow < 0, "Not_a_Bridge", 0))
    )
  )) %>%
  # Normalizing at the community level with the affiliation vector from consensus clustering
  group_by(Consensus_vector_0.15) %>%
  mutate(zPC_cons = as.numeric(scale(PC_cons))) %>%
  mutate(Hub_consensus = ifelse(zPC_cons > 0 & Within_module_z_cons > 0, "Connector",
    ifelse(zPC_cons > 0 & Within_module_z_cons < 0, "Satellite",
      ifelse(zPC_cons < 0 & Within_module_z_cons > 0, "Provincial",
        ifelse(zPC_cons < 0 & Within_module_z_cons < 0, "Peripheral", "Isolate")
      )
    )
  )) %>%
  arrange(Subj_ID, Region) %>%
  ungroup()

