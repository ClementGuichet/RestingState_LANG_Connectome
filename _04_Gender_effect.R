################################################################################
# Script for Gender-related analyses

# Written by CG
# 26-11-2022
################################################################################
library(ComplexHeatmap)
library(ggalluvial)
library(data.table)
library(RColorBrewer)
library(rstatix)
library(FactoMineR)
library(factoextra)
library(brunnermunzel)
library(jsonlite)

rm(list = ls())

source("_03_Hub_classification.R")
source("_radarplotting_function.R")
source("_NMI&AMI_functions.R")

# ~~~~~~~~~~~ Gender effect ~~~~~~~~~~~
by(data_full_per_subject %>% subset(threshold == "0.15"), factor((data_full_per_subject %>% subset(threshold == "0.15"))$Gender), summary)

data_stat_gender <- data_functional_role %>% 
  group_by(Subj_ID, CAB_NP_assign, Region, Gender) %>%
  summarize_at(vars(degree), mean) %>%
  filter(Gender != "NaN") %>%
  mutate(DC = degree / 59) %>%
  mutate_at(vars(Gender), funs(factor(.)))

data_stat_gender %>%
  group_by(Region, Gender) %>%
  get_summary_stats(DC, type = "median_iqr")

M <- gghistogram(data_stat_gender %>% subset(Gender == "M"),
  x = "DC", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
F <- gghistogram(data_stat_gender %>% subset(Gender == "F"),
  x = "DC", y = "..density..",
  fill = "purple", add_density = TRUE
) + theme_pubclean()
Rmisc::multiplot(M, F)
# Is there any significant difference between men and women degree centrality in each region?
# Mann-Withney U for right-skewed data
# Brunner-Munzel test for heteroskedastic-robustness

t_test_gender <- data_stat_gender %>%
  group_by(CAB_NP_assign, Region) %>%
  group_split() %>% # split per region and perform a t_test
  map_dfr(. %>%
    mutate(statistic = wilcox.test(DC ~ Gender, alternative = "two.sided", paired = FALSE, exact = FALSE, digits.rank = 7)$statistic) %>%
    mutate(p_value = wilcox.test(DC ~ Gender, alternative = "two.sided", paired = FALSE, exact = FALSE, digits.rank = 7)$p.value) %>%
    mutate(Mean_male = t.test(DC ~ Gender, alternative = "two.sided", paired = FALSE, var.equal = FALSE)$estimate[[2]]) %>%
    mutate(Mean_female = t.test(DC ~ Gender, alternative = "two.sided", paired = FALSE, var.equal = FALSE)$estimate[[1]]) %>%
    mutate(statistic_BM = brunnermunzel.test(DC ~ Gender, alternative = "two.sided")$statistic) %>%
    mutate(p_value_BM = brunnermunzel.test(DC ~ Gender, alternative = "two.sided")$p.value) %>%
    mutate(Estimate_BM = brunnermunzel.test(DC ~ Gender, alternative = "two.sided")$estimate)) %>%
  group_by(CAB_NP_assign, Region) %>%
  summarize_at(vars(statistic, p_value, Mean_male, Mean_female, statistic_BM, p_value_BM, Estimate_BM), mean)
# mutate(p_adjusted = p.adjust(p_value, method = "fdr"))

ggdotchart(
  data_gender_selection <- t_test_gender %>% subset(p_value <= 0.05),
  x = "Region", y = "statistic_BM",
  ylab = "Mann-Withney U (Brunner-Munzel Test if unequal variances)\n Positive statistic suggests higher degree centralities for men",
  palette = "jco",
  add = "segment", position = position_dodge(0.3),
  sorting = "descending",
  rotate = TRUE, legend = "none"
)

################################################################################
# Modular composition between Men & Women --------------------------------------
################################################################################

data_RS <- data_per_region_combined_bu %>% subset(threshold == "0.15")

# NMI & AMI
NMI_func(factor(data_RS$Consensus_vector_male), factor(data_RS$Consensus_vector_female))
AMI_func(factor(data_RS$Consensus_vector_male), factor(data_RS$Consensus_vector_female))

# Contingency table
addmargins(table(data_RS$Consensus_vector_male, data_RS$Consensus_vector_female))

# Create alluvial diagram between the two community structures

data_alluvial_community <- data_RS %>%
  dplyr::select(Consensus_vector_male, Consensus_vector_female, Region) %>%
  plyr::rename(c("Consensus_vector_male" = "Community_structure_men")) %>%
  plyr::rename(c("Consensus_vector_female" = "Community_structure_women")) %>%
  pivot_longer(
    cols = c("Community_structure_men", "Community_structure_women"),
    names_to = "Community_structure",
    values_to = "Communities"
  )

display_percentage <- data_alluvial_community %>%
  group_by(Community_structure, Communities) %>%
  summarize(s = n()) %>%
  group_by(Community_structure) %>%
  mutate(s = scales::percent(s / sum(s), accuracy = 0.1)) %>%
  arrange(Community_structure, desc(Communities)) %>%
  .$s

display_communities <- data_alluvial_community %>%
  group_by(Community_structure, Communities) %>%
  summarize(s = n()) %>%
  arrange(Community_structure, desc(Communities)) %>%
  .$Communities


alluvial_community <- ggplot(
  data_alluvial_community,
  aes(x = Community_structure, stratum = Communities, alluvium = Region, fill = Communities)
) +
  geom_flow(alpha = .7, curve_type = "sigmoid", width = .2, na.rm = TRUE) +
  geom_stratum(alpha = .8) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
  geom_text(stat = "stratum", label = display_communities) +
  scale_fill_brewer(palette = "Oranges", direction = 1) +
  labs(title = "Community structure differences between Men and Women") +
  theme_pubclean()

alluvial_community

################################################################################
# RSN Composition
################################################################################

# For men ----------------------------------------------------------------------
# With InLang method for network assignment
Net_proportion_male <- data_RS %>%
  # filter(!grepl("NaN", CAB_NP_assign)) %>%
  group_by(Consensus_vector_male, CAB_NP_assign) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Consensus_vector_male, desc(freq))

# Now using our method for network assignement
Net_proportion_male_our_method <- data_RS %>%
  mutate_at(vars(ends_with("value")), funs(as.numeric(as.character(.)))) %>%
  group_by(Consensus_vector_male, `1st_network`) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Consensus_vector_male, desc(freq))

Weighted_RSN_vector_male_our_method <- data_RS %>%
  mutate_at(vars(ends_with("value")), funs(as.numeric(as.character(.)))) %>%
  # 1 = High certainty, 0 = Low certainty
  mutate(certainty_factor = `1st_value` / 100) %>%
  group_by(Consensus_vector_male, `1st_network`) %>%
  # Mean proportion by RSNs when a ROI is assigned with that RSN
  summarise_at(vars(certainty_factor), mean)

Net_proportion_male_weighted <- merge(Net_proportion_male_our_method,
  Weighted_RSN_vector_male_our_method,
  by = c("Consensus_vector_male", "1st_network")
) %>%
  group_by(Consensus_vector_male) %>%
  mutate(n_adjusted = n * certainty_factor) %>%
  mutate(adjusted_freq = n_adjusted / sum(n_adjusted)) %>%
  arrange(Consensus_vector_male, desc(adjusted_freq))



p3 <- ggplot(Net_proportion_male, aes(
  # the group argument allows to stack according to the increasing values instead of the labels
  x = forcats::fct_rev(Consensus_vector_male), y = freq, group = factor(freq), fill = CAB_NP_assign
)) +
  geom_bar(stat = "identity", position = position_stack(), width = 0.5) +
  geom_text(
    data = Net_proportion_male %>% filter(freq > 0.12),
    aes(
      x = Consensus_vector_male, y = freq,
      label = scales::percent(freq, accuracy = .1)
    ),
    position = position_stack(vjust = .25)
  ) +
  geom_text(
    data = Net_proportion_male %>% filter(freq > 0.12),
    aes(
      x = Consensus_vector_male, y = freq,
      label = CAB_NP_assign
    ),
    position = position_stack(vjust = .7)
  ) +
  guides(fill = guide_legend(title = "CAB-NP Networks")) +
  # y implements how certain, on average, LANG nets correspond to each network
  labs(x = "Community structure for men", y = "Proportion using the InLang method for network assignement") +
  coord_flip() +
  scale_fill_manual(values = custom_palette) +
  theme_pubr(legend = "none")

p4 <- ggplot(Net_proportion_male_weighted, aes(
  # the group argument allows to stack according to the increasing values instead of the labels
  x = forcats::fct_rev(Consensus_vector_male), y = adjusted_freq, group = factor(adjusted_freq), fill = `1st_network`
)) +
  geom_bar(stat = "identity", position = position_stack(), width = 0.5) +
  geom_text(
    data = Net_proportion_male_weighted %>% filter(adjusted_freq > 0.12),
    aes(
      x = Consensus_vector_male, y = adjusted_freq,
      label = scales::percent(adjusted_freq, accuracy = .1)
    ),
    position = position_stack(vjust = .25)
  ) +
  geom_text(
    data = Net_proportion_male_weighted %>% filter(adjusted_freq > 0.12),
    aes(
      x = Consensus_vector_male, y = adjusted_freq,
      label = `1st_network`
    ),
    position = position_stack(vjust = .7)
  ) +
  guides(fill = guide_legend(title = "CAB-NP Networks")) +
  # y implements how certain, on average, LANG nets correspond to each network
  labs(x = "Community structure for men", y = "Proportion using our method for network assignment\n Weighted by the certainty of the overlap") +
  coord_flip() +
  scale_fill_manual(values = custom_palette) +
  theme_pubr(legend = "none")

# Compare outputs
gridExtra::grid.arrange(p3, p4)

# Radar plot of RSN per Community
Radar_RSN_community_male <- Net_proportion_male_weighted %>%
  dplyr::select(`1st_network`, adjusted_freq) %>%
  spread(`1st_network`, adjusted_freq) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Consensus_vector_male") %>%
  mutate_at(vars(everything()), funs(. * 100))


radarplotting_overlap(Radar_RSN_community_male, 100, 0, 1, 1,
  alpha = 0.2, label_size = 1,
  title = "Composition of each community for men\n weighted by the certainty factor of overlap",
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#FF99FF")
)
legend(
  x = "bottomleft", legend = rownames(Radar_RSN_community_male), horiz = TRUE,
  bty = "n", pch = 20, col = c("#00AFBB", "#E7B800", "#FC4E07", "#FF99FF"),
  text.col = "black", cex = 1, pt.cex = 2
)

# For women --------------------------------------------------------------------
# With InLang method for network assignment
Net_proportion_female <- data_RS %>%
  # filter(!grepl("NaN", CAB_NP_assign)) %>%
  group_by(Consensus_vector_female, CAB_NP_assign) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Consensus_vector_female, desc(freq))

# Now using our method for network assignement
Net_proportion_female_our_method <- data_RS %>%
  mutate_at(vars(ends_with("value")), funs(as.numeric(as.character(.)))) %>%
  group_by(Consensus_vector_female, `1st_network`) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(Consensus_vector_female, desc(freq))

Weighted_RSN_vector_female_our_method <- data_RS %>%
  mutate_at(vars(ends_with("value")), funs(as.numeric(as.character(.)))) %>%
  # 1 = High certainty, 0 = Low certainty
  mutate(certainty_factor = `1st_value` / 100) %>%
  group_by(Consensus_vector_female, `1st_network`) %>%
  # Mean proportion by RSNs when a ROI is assigned with that RSN
  summarise_at(vars(certainty_factor), mean)

Net_proportion_female_weighted <- merge(Net_proportion_female_our_method,
  Weighted_RSN_vector_female_our_method,
  by = c("Consensus_vector_female", "1st_network")
) %>%
  group_by(Consensus_vector_female) %>%
  mutate(n_adjusted = n * certainty_factor) %>%
  mutate(adjusted_freq = n_adjusted / sum(n_adjusted)) %>%
  arrange(Consensus_vector_female, desc(adjusted_freq))


p3 <- ggplot(Net_proportion_female, aes(
  # the group argument allows to stack according to the increasing values instead of the labels
  x = forcats::fct_rev(Consensus_vector_female), y = freq, group = factor(freq), fill = CAB_NP_assign
)) +
  geom_bar(stat = "identity", position = position_stack(), width = 0.5) +
  geom_text(
    data = Net_proportion_female %>% filter(freq > 0.12),
    aes(
      x = Consensus_vector_female, y = freq,
      label = scales::percent(freq, accuracy = .1)
    ),
    position = position_stack(vjust = .25)
  ) +
  geom_text(
    data = Net_proportion_female %>% filter(freq > 0.12),
    aes(
      x = Consensus_vector_female, y = freq,
      label = CAB_NP_assign
    ),
    position = position_stack(vjust = .7)
  ) +
  guides(fill = guide_legend(title = "CAB-NP Networks")) +
  # y implements how certain, on average, LANG nets correspond to each network
  labs(x = "Community structure for women", y = "Proportion using the InLang method for network assignement") +
  coord_flip() +
  scale_fill_manual(values = custom_palette) +
  theme_pubr(legend = "none")

p4 <- ggplot(Net_proportion_female_weighted, aes(
  # the group argument allows to stack according to the increasing values instead of the labels
  x = forcats::fct_rev(Consensus_vector_female), y = adjusted_freq, group = factor(adjusted_freq), fill = `1st_network`
)) +
  geom_bar(stat = "identity", position = position_stack(), width = 0.5) +
  geom_text(
    data = Net_proportion_female_weighted %>% filter(adjusted_freq > 0.12),
    aes(
      x = Consensus_vector_female, y = adjusted_freq,
      label = scales::percent(adjusted_freq, accuracy = .1)
    ),
    position = position_stack(vjust = .25)
  ) +
  geom_text(
    data = Net_proportion_female_weighted %>% filter(adjusted_freq > 0.12),
    aes(
      x = Consensus_vector_female, y = adjusted_freq,
      label = `1st_network`
    ),
    position = position_stack(vjust = .7)
  ) +
  guides(fill = guide_legend(title = "CAB-NP Networks")) +
  # y implements how certain, on average, LANG nets correspond to each network
  labs(x = "Community structure for women", y = "Proportion using our method for network assignment\n Weighted by the certainty of the overlap") +
  coord_flip() +
  scale_fill_manual(values = custom_palette) +
  theme_pubr(legend = "none")

# Compare outputs
gridExtra::grid.arrange(p3, p4)

# Radar plot of RSN per Community
Radar_RSN_community_female <- Net_proportion_female_weighted %>%
  dplyr::select(`1st_network`, adjusted_freq) %>%
  spread(`1st_network`, adjusted_freq) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Consensus_vector_female") %>%
  mutate_at(vars(everything()), funs(. * 100))


radarplotting_overlap(Radar_RSN_community_female, 100, 0, 1, 1,
  alpha = 0.4, label_size = 1,
  title = "Composition of each community for women\n weighted by the certainty factor of overlap",
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "#FF99FF")
)
legend(
  x = "bottomleft", legend = rownames(Radar_RSN_community_female), horiz = TRUE,
  bty = "n", pch = 20, col = c("#00AFBB", "#E7B800", "#FC4E07", "#FF99FF"),
  text.col = "black", cex = 1, pt.cex = 2
)



################################################################################
# Hubness profile
################################################################################

# Normalize PC & Wz for each gender and merge it back into one dataframe -------

# Import data ------------------------------------------------------------------
PC_consensus_M <- as.data.frame(fromJSON("Participation_coefficient_consensus_male.json")) %>%
  mutate(Subj_ID = rep(seq_len(34))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "PC_cons_M"
  )

Within_module_z_consensus_M <- as.data.frame(fromJSON("Within_module_z_score_consensus_male.json")) %>%
  mutate(Subj_ID = rep(seq_len(34))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "Within_module_z_cons_M"
  )

nodal_metrics_cons_M <- cbind(PC_consensus_M, Within_module_z_cons_M = Within_module_z_consensus_M$Within_module_z_cons_M)

PC_consensus_F <- as.data.frame(fromJSON("Participation_coefficient_consensus_female.json")) %>%
  mutate(Subj_ID = rep(seq_len(36))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "PC_cons_F"
  )

Within_module_z_consensus_F <- as.data.frame(fromJSON("Within_module_z_score_consensus_female.json")) %>%
  mutate(Subj_ID = rep(seq_len(36))) %>%
  pivot_longer(
    cols = !c("Subj_ID"),
    names_to = "Region",
    values_to = "Within_module_z_cons_F"
  )

nodal_metrics_cons_F <- cbind(PC_consensus_F, Within_module_z_cons_F = Within_module_z_consensus_F$Within_module_z_cons_F)

nodal_metrics_gender <- bind_rows(nodal_metrics_cons_M, nodal_metrics_cons_F)


tmp_gender <- data_full %>%
  # Retaining the 70 self-gendered subjects
  filter(Gender != "NaN") %>%
  subset(threshold == "0.15") %>%
  arrange(desc(Gender), Subj_ID, Region)


data_gender_0 <- cbind(tmp_gender,
  PC_cons_M = nodal_metrics_gender$PC_cons_M, Within_module_z_cons_M = nodal_metrics_gender$Within_module_z_cons_M,
  PC_cons_F = nodal_metrics_gender$PC_cons_F, Within_module_z_cons_F = nodal_metrics_gender$Within_module_z_cons_F
)

# Normalizing each gender independently with their own affiliation vector
data_gender_M <- data_gender_0 %>%
  filter(is.na(PC_cons_F)) %>%
  # Normalizing at the community level with the affiliation vector from male consensus clustering
  group_by(Consensus_vector_male) %>%
  mutate(zPC_cons_gender = as.numeric(scale(PC_cons_M))) %>%
  plyr::rename(c("Within_module_z_cons_M" = "Within_module_z_cons_gender")) %>%
  mutate(Hub_consensus_gender = ifelse(zPC_cons_gender > 0 & Within_module_z_cons_gender > 0, "Connector",
    ifelse(zPC_cons_gender > 0 & Within_module_z_cons_gender < 0, "Satellite",
      ifelse(zPC_cons_gender < 0 & Within_module_z_cons_gender > 0, "Provincial",
        ifelse(zPC_cons_gender < 0 & Within_module_z_cons_gender < 0, "Peripheral",
          "None"
        )
      )
    )
  )) %>%
  # Normalizing at the connectomic level
  mutate(zK = as.numeric(scale(degree))) %>%
  mutate(zBT = as.numeric(scale(Betweenness))) %>%
  mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
  mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
    ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
      ifelse(zBT > 0 & zFlow > 0, "Super_Bridge", "None")
    )
  )) %>%
  arrange(Subj_ID, Region)

data_gender_F <- data_gender_0 %>%
  filter(is.na(PC_cons_M)) %>%
  # Normalizing at the community level with the affiliation vector from male consensus clustering
  group_by(Consensus_vector_female) %>%
  mutate(zPC_cons_gender = as.numeric(scale(PC_cons_F))) %>%
  plyr::rename(c("Within_module_z_cons_F" = "Within_module_z_cons_gender")) %>%
  mutate(Hub_consensus_gender = ifelse(zPC_cons_gender > 0 & Within_module_z_cons_gender > 0, "Connector",
    ifelse(zPC_cons_gender > 0 & Within_module_z_cons_gender < 0, "Satellite",
      ifelse(zPC_cons_gender < 0 & Within_module_z_cons_gender > 0, "Provincial",
        ifelse(zPC_cons_gender < 0 & Within_module_z_cons_gender < 0, "Peripheral",
          "None"
        )
      )
    )
  )) %>%
  # Normalizing at the connectomic level
  mutate(zK = as.numeric(scale(degree))) %>%
  mutate(zBT = as.numeric(scale(Betweenness))) %>%
  mutate(zFlow = as.numeric(scale(Flow_coeff))) %>%
  mutate(Bridgeness = ifelse(zBT > 0 & zFlow < 0, "Global_Bridge",
    ifelse(zFlow > 0 & zBT < 0, "Local_Bridge",
      ifelse(zBT > 0 & zFlow > 0, "Super_Bridge", "None")
    )
  )) %>%
  arrange(Subj_ID, Region)

# Merge everything back into one df
data_gender_1 <- bind_rows(data_gender_M, data_gender_F) %>% ungroup()


# Apply hub detection at 20% and output the hubness profile --------------------

top <- 131 * 0.2

# Recreating a sequential vector to iterate through the list later
# When filtering the dataframe for each subject, we must subset the rows corresponding to the current subject which may differ from
# the ith position in the list if subjects labels' are not sequential e.g., 1, 3, 4
# This is because I have removed the two NaN gender subjects before normalizing PC & Wz
data_gender_1_helper_vector <- data_gender_1 %>%
  ungroup() %>%
  mutate(helper_vector = rep(seq(70), each = 131))

Top_metric_Gender_ind <- data_gender_1_helper_vector %>%
  group_by(helper_vector, Region) %>%
  summarize_at(vars(zK, Within_module_z_cons_gender, zPC_cons_gender, zBT, zFlow), mean) %>%
  # mutate(across(degree:PC, ~ rank(-.x), .names = "{.col}_rank")) %>%
  pivot_longer(
    cols = !c("helper_vector", "Region"),
    names_to = "Metric_name",
    values_to = "Metric_value"
  ) %>%
  group_by(helper_vector, Metric_name, .add = TRUE) %>%
  group_split() %>%
  map_dfr(. %>% slice_max(Metric_value, n = top) %>%
    mutate(rank = rep(seq(1:length(Region))))) %>%
  group_by(helper_vector, .add = TRUE) %>%
  group_split()

Hub_selection_gender <- list()
FR_list_gender <- list()
for (i in 1:length(Top_metric_Gender_ind)) {
  Hub_df <- rbindlist(Top_metric_Gender_ind[i]) %>% distinct(Region, .keep_all = TRUE)
  # Here I subset the rows specific to each subject and their Hub regions
  tmp <- data_gender_1_helper_vector %>%
    filter(Region %in% Hub_df$Region) %>%
    filter(helper_vector == i) %>%
    dplyr::select(helper_vector, Region, Hub_consensus_gender, Bridgeness)
  Hub_selection_gender[[i]] <- tmp

  FR_ind_hub <- tmp %>%
    group_by(Hub_consensus_gender) %>%
    summarize(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Hub_consensus_gender, freq)
  FR_ind_bridge <- tmp %>%
    group_by(Bridgeness) %>%
    summarize(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>%
    spread(Bridgeness, freq)

  FR_ind <- cbind(FR_ind_hub, FR_ind_bridge)
  FR_list_gender[[i]] <- FR_ind
}

data_cluster_efficiency_gender <- data_gender_1 %>%
  dplyr::select(Subj_ID, Eglob, Eloc) %>%
  group_by(Subj_ID) %>%
  summarize_at(vars(Eglob, Eloc), mean) %>%
  mutate(Balance_eff = (Eloc - Eglob) / (Eloc + Eglob)) %>%
  dplyr::select(-c(Subj_ID, Eglob, Eloc)) %>%
  mutate_at(vars(Balance_eff), funs(. * 100))

data_FR_Gender_ind <- cbind(
  rbindlist(FR_list_gender, fill = TRUE) %>% dplyr::select(-None) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% mutate_at(vars(everything()), funs(. * 100)),
  Balance_eff = data_cluster_efficiency_gender$Balance_eff,
  data_full_per_subject %>% filter(Gender != "NaN") %>%
    subset(threshold == "0.15") %>%
    dplyr::select(Subj_ID, Gender)
)


Radar_functional_role_gender <- data_FR_Gender_ind %>%
  dplyr::select(-c(Subj_ID)) %>%
  group_by(Gender) %>%
  summarize_at(vars(Connector:Super_Bridge), mean) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Gender")

radarplotting_overlap(Radar_functional_role_gender, 60, 0, 1, 1,
  alpha = 0.05, label_size = 1,
  title_fill = "Gender Hubness profile",
  palette = RColorBrewer::brewer.pal(8, "Dark2")
)

legend(
  x = "bottomleft", title = "Gender",
  legend = rownames(Radar_functional_role_gender), horiz = TRUE,
  bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
  text.col = "black", cex = 1, pt.cex = 2
)

# Hubness profile within each RSN
delta_hubness_profile <- function(cluster1, cluster2, alpha) {
  # Retain only the regions yielded by hub detection procedure across the two genders

  # Hub region specific to each subject yielded by hub detection procedure
  data_hub_selection_per_subject_gender <- rbindlist(Hub_selection_gender)

  # Get the associated Resting-state networks
  get_RSN_label <- data_functional_role %>%
    filter(Region %in% data_hub_selection_per_subject_gender$Region) %>%
    dplyr::select(`1st_network`, Region) %>%
    distinct()

  data_gender_final <- merge(data_hub_selection_per_subject_gender, get_RSN_label, by = "Region") %>%
    merge(., data_gender_1_helper_vector %>% dplyr::select(helper_vector, Subj_ID, Gender, Region),
      by = c("helper_vector", "Region")
    )
  
  Radar_hub_RSN <- data_gender_final %>% 
    group_by(Gender, `1st_network`) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    dplyr::select(-n) %>% 
    spread(`1st_network`, freq) %>% 
    remove_rownames() %>% column_to_rownames("Gender") %>% 
    mutate_at(vars(everything()), funs(. * 100))
  
  radarplotting_overlap(Radar_hub_RSN, 25, 0, 1, 1,
                        alpha = 0.05, label_size = 1,
                        title_fill = "Distribution of hubs regions across RSNs",
                        palette = RColorBrewer::brewer.pal(8, "Dark2")
  )
  
  legend(
    x = "bottomleft", title = "Gender",
    legend = rownames(Radar_hub_RSN), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )

  delta_proportion_a <- data_gender_final %>%
    group_by(`1st_network`, Gender, Hub_consensus_gender) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    # Make sure comparisons with missing functional roles can be achieved
    spread(Hub_consensus_gender, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>% 
    group_by(`1st_network`, Gender) %>%
    summarize_at(vars(Connector:Satellite), sum) %>% 
    pivot_longer(cols = !c("1st_network", "Gender"), names_to = "Hub_consensus_gender", values_to = "freq") %>%
    # Compute the difference in proportion of a given functional role within each RSN
<<<<<<< HEAD
    arrange(Gender) %>% 
=======
    arrange(Gender) %>%
>>>>>>> 81460f0e160007f36784da9dbdba95e2ccbea91a
    group_by(`1st_network`, Hub_consensus_gender) %>%
    mutate(delta_freq = freq / lag(freq)) 
    dplyr::select(-Gender) %>%
    na.omit()

  delta_proportion_b <- data_gender_final %>%
    group_by(`1st_network`, Gender, Bridgeness) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    # Make sure comparisons with missing functional roles can be achieved
    spread(Bridgeness, freq) %>%
    dplyr::select(-n) %>%
    mutate_all(., ~ replace(., is.na(.), 0)) %>%
    group_by(`1st_network`, Gender) %>%
    summarize_at(vars(Global_Bridge, Local_Bridge, Super_Bridge), sum) %>%
    pivot_longer(cols = !c("1st_network", "Gender"), names_to = "Bridgeness", values_to = "freq") %>%
    # Compute the difference in proportion of a given functional role within each RSN
    arrange(Gender) %>%
    group_by(`1st_network`, Bridgeness) %>%
<<<<<<< HEAD
    mutate(delta_freq = freq / lag(freq)) %>% 
=======
    mutate(delta_freq = freq / lag(freq))
    filter(delta_freq != "Inf") %>% 
>>>>>>> 81460f0e160007f36784da9dbdba95e2ccbea91a
    dplyr::select(-Gender) %>%
    na.omit()

  Radar_functional_role_RSN_delta <-
    delta_proportion_a %>%
<<<<<<< HEAD
    dplyr::select(`1st_network`, Hub_consensus_gender, delta_freq) %>% 
=======
    dplyr::select(`1st_network`, Hub_consensus_gender, delta_freq) %>%
    spread(`1st_network`, delta_freq) %>%
    # mutate_all(., ~ replace(., is.na(.), 0)) %>%
>>>>>>> 81460f0e160007f36784da9dbdba95e2ccbea91a
    subset(Hub_consensus_gender != "None") %>%
    mutate(`1st_network` = ifelse(delta_freq == "Inf"|delta_freq == 0, paste0(`1st_network`, "*"), `1st_network`)) %>% 
    mutate(delta_freq = ifelse(delta_freq == "Inf", 1, delta_freq)) %>%
    spread(`1st_network`, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Hub_consensus_gender")


<<<<<<< HEAD
  radarplotting_overlap(Radar_functional_role_RSN_delta, 3, -1, 1, 1,
=======
  radarplotting_overlap(Radar_functional_role_RSN_delta, 4, 0, 1, 1,
>>>>>>> 81460f0e160007f36784da9dbdba95e2ccbea91a
    alpha = alpha, label_size = 1,
    title_fill = "Ratio of the proportion of functional roles between genders. A positive ratio favors men",
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", title = "\n * indicates there were no Provincial hubs in Visual 2 for men",
    legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )

  Radar_functional_role_RSN_delta <-
    delta_proportion_b %>%
    dplyr::select(`1st_network`, Bridgeness, delta_freq) %>%
<<<<<<< HEAD
=======
    spread(`1st_network`, delta_freq) %>%
    # mutate_all(., ~ replace(., is.na(.), 0)) %>%
>>>>>>> 81460f0e160007f36784da9dbdba95e2ccbea91a
    subset(Bridgeness != "None") %>%
    mutate(`1st_network` = ifelse(delta_freq == "Inf"|delta_freq == 0, paste0(`1st_network`, "*"), `1st_network`)) %>% 
    mutate(delta_freq = ifelse(delta_freq == "Inf", 1, delta_freq)) %>%
    spread(`1st_network`, delta_freq) %>%
    remove_rownames() %>%
    column_to_rownames(var = "Bridgeness")

  radarplotting_overlap(Radar_functional_role_RSN_delta, 2, 0, 1, 1,
    alpha = alpha, label_size = 1,
    title_fill = "Ratio of the proportion of functional roles between genders. A positive ratio favors men",
    palette = RColorBrewer::brewer.pal(8, "Dark2")
  )

  legend(
    x = "bottomleft", title = "* indicates there were no Global Bridge hubs in VMM for women",
    legend = rownames(Radar_functional_role_RSN_delta), horiz = TRUE,
    bty = "n", pch = 20, col = RColorBrewer::brewer.pal(8, "Dark2"),
    text.col = "black", cex = 1, pt.cex = 2
  )
}

delta_hubness_profile("M", "F", 0.2)
