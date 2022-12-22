##########################################################################################
# Script for compositional data analysis

# Written by CG
# 26-11-2022
##########################################################################################
library(data.table)
library(RColorBrewer)
library(rstatix)
library(Rmisc)
library(fitdistrplus)
library(sigclust2)
library(ggpubr)
library(cluster)
library(tidyverse)
library(compositions)
library(zCompositions)
library(robCompositions)


rm(list = ls())

set.seed(2022)
source("_03_Hub_detection.R")

# What are the graph-based biomarkers of health aging ?
# Does the topological_functional profile, that is the proportion of the hubs' functional role, change throughout life ?

###############################################################################
# Getting ready for hierarchical clustering analysis
###############################################################################

data_compositional <- data_hubness_profile_Age_ind %>%
  dplyr::select(Connector:Satellite)

data_compositional_bis <- data_hubness_profile_Age_ind %>%
  dplyr::select(Global_Bridge, Local_Bridge, Super_Bridge, Not_a_Bridge)

# Bayesian non-parametric multiplicative replacement that preserves the ratios between non-zero components ----
imputated_data <- cmultRepl(data_compositional, output = "prop")
imputated_data_bis <- cmultRepl(data_compositional_bis, output = "prop")

data_ilr_bis <- acomp(imputated_data) # Declare the dataset to be compositional
# and use relative geometry
plot(data_ilr_bis) # plot.acomp : ternary diagram
ellipses(mean(data_ilr_bis), var(data_ilr_bis), r = 2, col = "red") # Simplex 2sigma predictive region
pr_bis <- princomp(data_ilr_bis)
straight(mean(data_ilr_bis), pr_bis$Loadings)


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

one <- ilrV(imputated_data)
two <- ilrV(imputated_data_bis)
data_ilr <- cbind(one, two)

# PCA of ILR-transformed data because a non-singular covariance matrix is needed/ robust covariance estimation need a full-rank matrix
# CLR removes the value-range restriction but not the unit-sum constraint which makes PCA sensitive to outliers ----

cv <- robustbase::covMcd(data_ilr, cor = FALSE)
pcaIlr <- princomp(data_ilr, covmat = cv, cor = FALSE)
pcaIlr$loadings
# Ilr-robust biplot
biplot(pcaIlr, scale = 0)

# Keeping 5 log contrasts which explain 95% of variance
data_cluster <- pcaIlr$scores %>%
  as.data.frame() %>%
  dplyr::select(Comp.1:Comp.4)

# CLR-backtransform for interpretability of compositional variability and compositional biplot
data_CODA <- cbind(imputated_data, imputated_data_bis) %>% as.data.frame()


# Back-transform to clr fo rinterpretability of compositional variability
robCODA <- robCompositions::pcaCoDa(data_CODA,
  method = "robust", solve = "eigen",
  mult_comp = list(
    c(1, 2, 3, 4, 5), c(6, 7, 8, 9)
  )
)
# Clr-robust biplot
biplot(robCODA, scale = 0)
robCODA$loadings
summary(robCODA)
# out <- outCoDa(data_cluster, coda = FALSE)
# out$limit
# out$outlierIndex
################################################################################
# Clustering on Ilr-transformed princiapl balances ----
################################################################################

performance::check_clusterstructure(data_cluster)
corrplot::corrplot(cor(data_cluster))
ComplexHeatmap::Heatmap(scale(data_cluster))

# Define linkage methods
link <- c("average", "single", "complete", "ward")
names(link) <- c("average", "single", "complete", "ward")
# Function to compute the agglomerative coefficient (i.e., the strength of the clusters)
ac <- function(x) {
  cluster::agnes(data_cluster, method = x)$ac
}
sapply(link, ac)

# Defining optimal number of cluster with gap statistic
gap_stat <- cluster::clusGap(data_cluster,
  FUN = hcut,
  K.max = 5,
  B = 1000,
  verbose = T
)

factoextra::fviz_gap_stat(gap_stat)

plot_cluster <- cluster::agnes(data_cluster, method = "ward", metric = "euclidian")
fviz_dend(plot_cluster,
  cex = 0.8,
  k = 4,
  palette = "jco",
  rect = T,
  color_labels_by_k = TRUE,
  main = "Ward Hierarchical clustering of compositional topologico-functional profiles"
)


# Evaluate significance of hclust using 1000 Monte Carlo simulated null Gaussian
shc_result <- sigclust2::shc(as.matrix(data_cluster),
  metric = "euclidian",
  linkage = "ward.D2",
  n_sim = 1000,
  alpha = 0.05
)

base::plot(shc_result, hang = .1)

# Compute distance matrix
d <- dist(data_cluster, method = "euclidian")
final_clust <- hclust(d, method = "ward.D2")
groups <- cutree(final_clust, k = 4)


data_post_clustering <- cbind(data_hubness_profile_Age_ind, cluster = groups)


data_post_clustering %>%
  group_by(cluster) %>%
  get_summary_stats(Age, type = "full")
