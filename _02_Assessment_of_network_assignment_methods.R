################################################################################
# Written by CG
# 23-11-2022

# Assessment of network assignment method between InLang and Network overlap calc
################################################################################
library(tidyverse)
library(ggpubr)
library(ggalluvial)

rm(list = ls())

##########################################################################################
# Import processed data------------------------------------------------------------------
source("_01_DataManip.R")
source("_NMI&AMI_functions.R")
source("_radarplotting_function.R")

# Helpers
# Define normalizing function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Define palette for visualization
custom_palette <- c(
  "Auditory" = "#FF99FF",
  "Language" = "#FF6600",
  "CON" = "#9900CC",
  "DMN" = "#FF0033",
  "FPN" = "#FFCC33",
  "SMN" = "#0099CC",
  "DAN" = "#33FF66",
  "Visual_1" = "#CCCCCC",
  "Visual_2" = "#CCCCCC",
  "PMM" = "#CC0033",
  "VMM" = "#CC0033",
  "NaN" = "white",
  "Multi" = "white",
  "Multi/SM" = "white"
)


################################################################################
# Look at the contingency table between the two network assignment methods for the 131 LANG ROIs
data_contingency <- data_full_per_region %>%
  subset(threshold == "0.15") %>%
  dplyr::select(Region, CAB_NP_assign, `1st_network`) %>%
  plyr::rename(c("CAB_NP_assign" = "InLang method")) %>%
  plyr::rename(c("1st_network" = "Our method"))

addmargins(table(data_contingency$`Our method`, data_contingency$`InLang method`, useNA = "always"))
AMI_func(factor(data_contingency$`InLang method`), factor(data_contingency$`Our method`))

# Create alluvial diagram between the two network assignment methods
data_alluvial <- data_contingency %>%
  pivot_longer(
    cols = c("InLang method", "Our method"),
    names_to = "Labeling",
    values_to = "Networks"
  )

display_percentage <- data_alluvial %>%
  group_by(Labeling, Networks) %>%
  summarize(s = n()) %>%
  group_by(Labeling) %>%
  mutate(s = scales::percent(s / sum(s), accuracy = 0.1)) %>%
  arrange(Labeling, desc(Networks)) %>%
  .$s

display_Networks <- data_alluvial %>%
  group_by(Labeling, Networks) %>%
  summarize(s = n()) %>%
  arrange(Labeling, desc(Networks)) %>%
  .$Networks

display_Regions <- data_alluvial %>%
  subset(Labeling == "InLang method") %>%
  arrange(Networks) %>%
  .$Region

alluvial <- ggplot(
  data_alluvial,
  aes(x = Labeling, stratum = Networks, alluvium = Region, fill = Networks)
) +
  geom_flow(alpha = .3, curve_type = "arctangent", width = .2, na.rm = TRUE) +
  geom_stratum(alpha = .8) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
  geom_text(stat = "stratum", label = display_Networks, nudge_x = 0.05) +
  scale_fill_manual(values = custom_palette) +
  labs(title = "Difference between network assignment methods across the 131 LANG ROIs") +
  theme(
    axis.text.y = element_blank(),
    panel.background = element_blank()
  ) +
  # See contingency table
  annotate("text",
    fontface = "bold",
    x = 1.5, y = 128, label = "30%",
    color = "black", size = 4
  ) +
  annotate("text",
    fontface = "bold",
    x = 1.5, y = 120, label = "36.8%",
    color = "black", size = 4
  ) +
  annotate("text",
    fontface = "bold",
    x = 1.5, y = 91, label = "38.7%",
    color = "black", size = 4
  ) +
  annotate("text",
    fontface = "bold",
    x = 1.5, y = 62, label = "30.7%",
    color = "black", size = 4
  ) +
  annotate("text",
    fontface = "bold",
    x = 1.5, y = 42, label = "28.6%",
    color = "black", size = 4
  ) +
  annotate("text",
    fontface = "bold",
    x = 1.5, y = 10.75, label = "64.7%",
    color = "black", size = 4
  ) +
  annotate("text",
    fontface = "bold",
    x = 1.5, y = -5, label = "% of regions preserving their network assigment (forward)\n NMI = 0.27",
    color = "black", size = 4
  )


alluvial
