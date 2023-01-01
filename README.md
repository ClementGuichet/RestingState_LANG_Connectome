<h1 align="center">Analysis of the topological changes in the LANG Connectome across the adult lifespan</h1>

Method Overview

- Preprocessing of resting-state data using SPM12
- Relying on 6mm-264-regions from Power atlas (Power et al., 2011)
- Connectivity matrices obtained via CONN Toolbox processing
- Graph-theory measures with GraphVar and BCT Toolbox on the previously defined 131 ROIs of the LANG Connectome (Roger et al., 2022)
- Hub classification, hub detection, and computation of the individual topologico-functional profile (TFP)
- Robust PCA on ILR-transformed data
- Clustering analysis with principal log-contrasts

- Investigation of the three-way interaction between AGE, TFP, & RSNs using log-ratios of the (extended) geometric mean
- Multilevel bootstrapping to compute the confidence intervals of log-ratios

- Entropy-based measures to explore the modes of topological reconfiguration across RSNs
- Investigation of the most probable topological reconfiguration trajectory

Features

- More accurate method for RSN assignment (with voxel-level precision) to ROIs 
- Original topological and functional hub classification according to modular and interareal nodal centrality metrics
(i.e., degree, Participation coefficient, Within-module_z-score, Betweenness centrality, Flow coefficient)
- Compositional data analysis
- Shiny app to explore the most probable topological reconfiguration trajectory across AGE: https://clementguichet.shinyapps.io/shinyapp/
 
## Author
**Clément Guichet**
- Clement.Guichet@univ-grenoble-alpes.fr
- [ResearchGate](https://www.researchgate.net/profile/Clement-Guichet)

## Lab
- CNRS UMR 5105 - Laboratoire Psychologie et NeuroCognition (LPNC)
- https://lpnc.univ-grenoble-alpes.fr/laboratoire-0
