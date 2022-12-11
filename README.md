"# RestingState_LANG_Connectome" 

Method Overview
1/ Preprocessing of resting-state data using SPM12
2/ Relying on 6mm-264-regions from Power atlas (Power et al., 2011)
3/ Connectivity matrices obtained via CONN Toolbox processing
4/ Graph-theory measures with GraphVar and BCT Toolbox on the previously defined 131 ROIs of the LANG Connectome (Roger et al., 2022)

Features:
- Original nodal functional role classification according to nodal centrality and information flow
(i.e., degree, Participation coefficient, Within-module_z-score, Betweenness centrality, Flow coefficient)
- Individual hub detection procedure
- Data-driven approach for Age and Gender-related effects
- Original hubness profile visualization
 
