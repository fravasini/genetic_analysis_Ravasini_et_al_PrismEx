# genetic_analysis_Ravasini_et_al_PrismEx

Code to replicate IDW interpolation maps of MDS values and phylogenetic tree starting from EIGENSTRAT file as in Ravasini et al.


`Compute_f2_and_make_simple_tree.R` --> To compute f2-statistics from EIGENSTRAT values with [admixtools](https://uqrmaie1.github.io/admixtools/index.html).

`IDW_Interpolation_map.R` --> To produce the interpolation maps based on the first dimension of the MDS, as in figure 4 of the paper.

`metameta_for_interpolation.txt` --> Example of the metadata file required for the interplation maps.

`plot_f2_NJ_tree.R` --> To produce the tree with the collapsed branches if all the samples downstream of the node belong to the same culture (ATU 2.1 in the paper), as in figure 5 of the paper. Please note that the tree in figure 5 was additionally manually curated for a better representation.

`metadata_for_tree_figure.txt` --> Example of the metadata file required for the tree figure.
