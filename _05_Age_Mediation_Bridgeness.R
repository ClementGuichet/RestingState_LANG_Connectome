##########################################################################################
# Script for the cross-sectional mediation analysis with SEM

# Written by CG
# 14-12-2022
##########################################################################################
library(tidyverse)
library(lavaan)
library(semPlot)


source("_04_Age_effect_post_clustering.R")

# Reconfiguration of global into super bridges with age ---------------
# Do the hub regions previously considered global bridge become super bridge later in life?

cor <- cor.test(data_post_clustering$Global_Bridge, data_post_clustering$Super_Bridge)
# Creating the plot
plot(data_post_clustering$Super_Bridge, data_post_clustering$Global_Bridge, pch = 19, col = "lightblue")
# Regression line
abline(lm(data_post_clustering$Global_Bridge ~ data_post_clustering$Super_Bridge), col = "red", lwd = 3)
# Pearson correlation
text(paste0("Correlation between Super & Global Bridge: ", round(cor$estimate, 2), "****"), x = 50, y = 30)


# Testing informatin flow reconfig
data_reconfig <- data_functional_role %>%
  group_by(Subj_ID, CAB_NP_assign, `1st_network`, Region, Age, Consensus_vector_0.15, LANG_Net_assign) %>%
  summarize_at(vars(zFlow, zBT), mean)

gghistogram(data_reconfig,
  x = "zFlow", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()

cor.test(data_reconfig$zFlow, data_reconfig$Age, method = "kendall")


gghistogram(data_reconfig,
  x = "zBT", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()

cor.test(data_reconfig$zBT, data_reconfig$Age, method = "kendall")

# Mediation

data_reconfig_mediation_zFlow <- tmp_cluster_final %>%
  filter(cluster == "25" | cluster == "50") %>%
  group_by(Region, cluster) %>%
  summarise_at(vars(zFlow), mean) %>%
  spread(cluster, zFlow) %>%
  mutate(zFlow_diff = `50` - `25`)

data_reconfig_mediation <- tmp_cluster_final %>%
  filter(cluster == "25" | cluster == "50") %>%
  group_by(Region, cluster, Bridgeness) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  mutate(reconfig = ifelse(cluster == "25" & Bridgeness == "Global_Bridge", "GB_25",
    ifelse(cluster == "25" & Bridgeness == "Super_Bridge", "SP_25",
      ifelse(cluster == "50" & Bridgeness == "Super_Bridge", "SP_50",
        ifelse(cluster == "25" & Bridgeness == "Not_a_Bridge", "NB_25",
          ifelse(cluster == "50" & Bridgeness == "Local_Bridge", "LB_50", 0)
        )
      )
    )
  )) %>%
  filter(reconfig != 0) %>%
  ungroup() %>%
  dplyr::select(Region, freq, reconfig) %>%
  spread(reconfig, freq) %>%
  mutate_all(., ~ replace(., is.na(.), 0))

mediation <- cbind(data_reconfig_mediation, zFlow = data_reconfig_mediation_zFlow$zFlow_diff) %>%
  dplyr::select(-Region) %>%
  scale(.)


mod <- "
  zFlow ~ ind1*GB_25
  # + ind2*NB_25
  SP_50 ~ direct1*GB_25 + ind3*zFlow
  # LB_50 ~ direct2*NB_25 + ind4*zFlow
  dir_up := direct1
  # dir_down := direct2
  mediation_up := ind1*ind3
  # mediation_down := ind2*ind4
  tot_up := (ind1*ind3) + direct1
  # tot_down := (ind2 + ind4) + direct2
  # total := (ind1*ind3) + direct1 + (ind2 + ind4) + direct2
"

fit <- sem(mod, data = mediation)
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE, ci = TRUE)
parameterEstimates(fit, level = 0.95, boot.ci.type = "bca.simple", standardized = TRUE)
semPaths(fit, "est", layout = "tree2", edge.label.cex = 1.25, fade = FALSE)
