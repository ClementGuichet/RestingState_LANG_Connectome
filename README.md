<h1 align="center">RestingState_LANG_Connectome</h1>

Method Overview

- Preprocessing of resting-state data using SPM12
- Relying on 6mm-264-regions from Power atlas (Power et al., 2011)
- Connectivity matrices obtained via CONN Toolbox processing
- Graph-theory measures with GraphVar and BCT Toolbox on the previously defined 131 ROIs of the LANG Connectome (Roger et al., 2022)
- Hub classification, hub detection, and computation of the individual topologico-functional profile
- Robust PCA on ILR-transformed compositional data
- Wald hierarchical clustering analysis with log-contrasts
- Kullback-Leibler & Generalized Jensen-Shannon divergence to investigate the three-way interaction between AGE, Topologico-functional profile, Resting-state-networks or community structure
- Multilevel bootstrapping to compute the confidence intervals of log-ratios between the young and old cluster
- Investigation of the most probable topological reconfiguration trajectory

Features

- Original methodology for accurate RSN assignment (with voxel-level precision) to ROIs 
- New topological and functional hub classification according to modular and interareal nodal centrality measures
(i.e., degree, Participation coefficient, Within-module_z-score, Betweenness centrality, Flow coefficient)
 
## Author
**Clément Guichet**
- Clement.Guichet@univ-grenoble-alpes.fr
- [ResearchGate](https://www.researchgate.net/profile/Clement-Guichet)

## Lab
- CNRS UMR 5105 - Laboratoire Psychologie et NeuroCognition (LPNC)
- https://lpnc.univ-grenoble-alpes.fr/laboratoire-0
