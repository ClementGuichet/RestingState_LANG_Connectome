##########################################################################################
# Script for compositional data analysis

# Written by CG
# 26-11-2022
##########################################################################################
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(sigclust2)
library(ggpubr)
library(cluster)
library(tidyverse)
library(compositions)
library(zCompositions)
library(robCompositions)
library(mvoutlier)
library(factoextra)
library(mclust)


source("_04_Hub_detection.R")

###############################################################################
# Getting ready for hierarchical clustering analysis
###############################################################################

data_coda_modular <- TFP_General %>%
  filter(Gender != "NaN" & Age != "NaN") %>%
  dplyr::select(Connector, Satellite, Provincial, Peripheral) %>%
  acomp(.)

data_coda_interareal <- TFP_General %>%
  filter(Gender != "NaN" & Age != "NaN") %>%
  dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>%
  acomp(.) %>%
  # Bayesian non-parametric multiplicative replacement that preserves the ratios between non-zero components
  cmultRepl(., output = "prop")



################################################################################
# The compositional data set is expressed in isometric logratio coordinates.
# Then, robust principal component analysis is performed with a robust covariance MCD-estimator.

# ILR - moving the D-part compositional data from the simplex to a D-1 dimensional real space
# Orthonormal basis recommended by Egozcue et al. (2003)
ilrV <- function(x) {
  x.ilr <- matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
  for (i in 1:ncol(x.ilr)) {
    x.ilr[, i] <- sqrt((i) / (i + 1)) * log(((apply(as.matrix(x[, 1:i]), 1, prod))^(1 / i)) / (x[, i + 1]))
  }
  return(x.ilr)
}
data_ilr <- cbind(ilrV(data_coda_modular), ilrV(data_coda_interareal))

library(vegan)

vegan::protest(TFP_General %>% filter(Gender != "NaN" & Age != "NaN") %>% 
                 dplyr::select(Connector, Satellite, Provincial, Peripheral,
                               Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge) %>% 
                 as.matrix(),
               data_ilr)
# Correlation in a symmetric Procrustes rotation: 0.8901 

# PCA of ILR-transformed data because a non-singular covariance matrix is needed/ robust covariance estimation need a full-rank matrix
# CLR removes the value-range restriction but not the unit-sum constraint which makes PCA sensitive to outliers ----
set.seed(4)
cv <- robustbase::covMcd(data_ilr, nsamp = "deterministic")
pcaIlr <- princomp(data_ilr, covmat = cv)
pcaIlr$scores
biplot(pcaIlr, scale = 0)


# CLR-backtransform for interpretability of compositional variability and compositional biplot
data_CODA <- cbind(data_coda_modular, data_coda_interareal) %>% as.data.frame()

# Back-transform to clr fo rinterpretability of compositional variability
robCODA <- robCompositions::pcaCoDa(data_CODA,
  method = "robust", solve = "eigen",
  mult_comp = list(
    c(1, 2, 3, 4), c(5, 6, 7, 8)
  )
)

summary(robCODA)
biplot(robCODA, scale = 0)
plot(robCODA, type = "l")
robCODA$scores
robCODA$loadings

# Outlier detection
par(mfrow = c(2, 1))
res <- mvoutlier.CoDa(data_CODA)
plot(res, which = "parallel", onlyout = TRUE)
plot(outCoDa(data_CODA))

# Keeping log contrasts which explain 95% of variance
data_cluster <- pcaIlr$scores %>%
  as.data.frame() %>%
  dplyr::select(Comp.1:Comp.5)
# Adding Age and Gender features
# cbind(., Age = (TFP_General %>% filter(Gender != "NaN" & Age != "NaN"))$Age)

##########################################################################
# Finite Gaussian mixture model clustering ----
################################################################################

par(mfrow = c(1, 1))
mc <- Mclust(data_cluster, G = 3:9, verbose = TRUE, scale = FALSE)
summary(mc)
# Mclust EEI (diagonal, equal volume and shape) model with 3 components: 
#   
#   log-likelihood  n df       BIC       ICL
# -149.3084 69 22 -391.7671 -409.1081
# 
# Clustering table:
#   1  2  3 
# 35 18 16
mc$parameters

library(fpc)
fpc::calinhara(data_cluster, mc$classification)
fpc::cluster.stats(dist(data_cluster), mc$classification)

######
# BIC values used for choosing the number of clusters
fviz_mclust(mc, "BIC", palette = "jco") + theme_pubclean()
# Classification: plot showing the clustering
fviz_mclust(mc, "classification",
  geom = "point",
  pointsize = 1.5, palette = "jco"
) + theme_pubclean()
# Classification uncertainty
fviz_mclust(mc, "uncertainty", palette = "jco") + theme_pubclean()


data_post_clustering <- cbind(TFP_General %>% filter(Gender != "NaN" & Age != "NaN"),
  cluster = mc$classification, mc$z
) %>%
  mutate(filter = ifelse(V1 > .8, 1,
    ifelse(V2 > .8, 1,
      ifelse(V3 > .8, 1, 0)
    )
  )) %>%
  filter(filter == 1) %>%
  .[, 1:13]

data_post_clustering %>%
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full")

################################################################################
# Finite Gaussian and Binary mixture modelling ----
################################################################################
# library(flexmix)
#
# set.seed(4)
# flex_mod <- flexmix(data_cluster %>% as.matrix() ~ 1, k = 3, model = FLXMCmvnorm(diagonal = TRUE))
# parameters(flex_mod)
# summary(flex_mod)
#
# plot(flex_mod)
#
# data_flex <- flex_mod@posterior %>% as.data.frame() %>%
#   .[,1:3]
#
# data_post_clustering <- cbind(TFP_General %>%  filter(Gender != "NaN" & Age != "NaN"),
#                               cluster = flex_mod@cluster, data_flex) %>%
#   mutate(filter = ifelse(scaled.1 > .8, 1,
#                        ifelse(scaled.2 > .8, 1,
#                               ifelse(scaled.3 > .8, 1, 0)))) %>%
#   filter(filter == 1) %>%
#   .[,1:13]
#
# data_post_clustering %>%
#   group_by(cluster) %>%
#   get_summary_stats(Age, type = "full")

################################################################################
# Hierarchical clustering ----
################################################################################
#
# # Clustering tendency
# performance::check_clusterstructure(data_cluster)
# corrplot::corrplot(cor(data_cluster))
# ComplexHeatmap::Heatmap(scale(data_cluster))
# #
# # # Define linkage methods
# # link <- c("average", "single", "complete", "ward")
# # names(link) <- c("average", "single", "complete", "ward")
# # # Function to compute the agglomerative coefficient (i.e., the strength of the clusters)
# # ac <- function(x) {
# #   cluster::agnes(data_cluster, method = x)$ac
# # }
# # sapply(link, ac)
#
# # Defining optimal number of cluster with gap statistic
# gap_stat <- cluster::clusGap(data_cluster,
#   FUN = hcut,
#   K.max = 5,
#   B = 1000,
#   verbose = T
# )
#
# factoextra::fviz_gap_stat(gap_stat)
#
# plot_cluster <- cluster::agnes(data_cluster, method = "ward", metric = "gower")
# fviz_dend(plot_cluster,
#   cex = 0.8,
#   k = 4,
#   palette = "jco",
#   rect = T,
#   color_labels_by_k = TRUE,
#   main = "Ward Hierarchical clustering of compositional topologico-functional profiles"
# )
#
#
# # Evaluate significance of hclust using 1000 Monte Carlo simulated null Gaussian
# shc_result <- sigclust2::shc(as.matrix(data_cluster),
#   metric = "euclidian",
#   linkage = "ward.D2",
#   n_sim = 1000,
#   alpha = 0.05
# )
#
# base::plot(shc_result, hang = .1)
#
# library(NbClust)
# output <- NbClust::NbClust(data_cluster, distance = "euclidean", method = "ward.D2")
# output$Best.partition
# # # Compute distance matrix
# d <- dist(data_cluster, method = "euclidean")
# final_clust <- hclust(d, method = "ward.D2")
# groups <- cutree(final_clust, k = 2)
#
#
# data_post_clustering <- cbind(TFP_General %>%  filter(Gender != "NaN" & Age != "NaN"),
#                               cluster = groups)
#
# data_post_clustering %>%
#   group_by(cluster) %>%
#   get_summary_stats(Age, type = "full")
#
# data_post_clustering %>%
#   group_by(cluster) %>%
#   count(Gender)

################################################################################
# Multiple Factor Analysis ----
################################################################################

# data_mfa <- data_ilr %>%
#   as.data.frame() %>%
#   # Adding Age and Gender features
#   cbind(., Age = (TFP_General %>% filter(Gender != "NaN" & Age != "NaN"))$Age,
#         Gender = (TFP_General %>% filter(Gender != "NaN" & Age != "NaN"))$Gender)
#
#
# mfa <- FactoMineR::MFA(data_mfa,
#                 group = c(3, 3, 1, 1),
#                 type = c("c", "c", "s", "n"),
#                 name.group = c("modular", "interareal", "age", "gender")
#                 )
# summary(mfa)
#
# get_eigenvalue(mfa)
# fviz_screeplot(mfa)
# fviz_mfa_var(mfa, "group")
#
# fviz_mfa_var(mfa, "quanti.var", palette = "jco",
#              col.var = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE)
#
# fviz_contrib(mfa, choice = "quanti.var", axes = 1, top = 20,
#              palette = "jco")
# fviz_contrib(mfa, choice = "quanti.var", axes = 2, top = 20,
#              palette = "jco")
#
# fviz_mfa_ind(mfa, col.ind = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE)
#
# data_cluster <- mfa$ind$coord %>% as.data.frame()
