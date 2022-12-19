##########################################################################################
# Script for Age-related analyses

# Written by CG
# 26-11-2022
##########################################################################################
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(fitdistrplus)
library(ggpubr)

rm(list = ls())

source("_03_Hub_detection.R")

# What are the graph-based biomarkers of health aging ?

################################################################################
# Correlation between degree centrality and Age --------------------------------
# Inferential statistics

data_stat_age <- data_functional_role %>%
  group_by(Subj_ID, CAB_NP_assign, Region, Age, Consensus_vector_0.15, LANG_Net_assign) %>%
  summarize_at(vars(degree), mean) %>%
  filter(Age != "NaN") %>%
  mutate(DC = degree / max(.$degree))


gghistogram(data_stat_age,
  x = "Age", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
gghistogram(data_stat_age,
  x = "DC", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()

fitdistrplus::descdist(data_stat_age$DC)

correlation_DC_age <- data_stat_age %>%
  group_by(CAB_NP_assign, Region, Consensus_vector_0.15, LANG_Net_assign) %>%
  group_split() %>%
  map_dfr(. %>%
    mutate(Estimate = cor.test(.$DC, .$Age, method = "kendall")$estimate) %>%
    mutate(p_value = cor.test(.$DC, .$Age, method = "kendall")$p.value)) %>%
  group_by(CAB_NP_assign, Region, Consensus_vector_0.15, LANG_Net_assign) %>%
  summarise_at(vars(Estimate, p_value), mean)

ggdotchart(
  correlation_DC_age %>% subset(p_value <= 0.05),
  x = "Region", y = "Estimate",
  ylab = "Mean correlation",
  palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "descending",
  rotate = TRUE, legend = "none", title = "Significant correlations between degree centrality and Age"
)

ComplexHeatmap::Heatmap(as.matrix(correlation_DC_age %>% subset(p_value <= 0.05) %>% arrange(desc(CAB_NP_assign)) %>% ungroup() %>% dplyr::select(Region, Estimate) %>% remove_rownames() %>% column_to_rownames("Region")), cluster_rows = FALSE, column_names_rot = TRUE, name = "Heatmap of significant correlations\n between degree centrality and Age")
