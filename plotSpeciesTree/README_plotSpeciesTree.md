### Prune and plot a species tree

__WRT species tree inference:__  
I selected 'core genes' in Fol4287 (defined as genes located on core chromosomes, where core is defined wrt to _F. verticillioides_: hence excluding genes located on chr 1b, 2b, 3, 6, 14 and 15). I did a megaBLAST against genomic sequences (many genoems in our dataset are not annotated), and used quite arbitrary length and similarity cutoffs to include or exclude genes. All genes that had a single hit in all genomes in the dataset according to these cutoffs were included. For these genes and their BLAST hits, I inferred multi sequence alignments (MSAs) using Clustal Omega, concatenated the MSAs using a Python script and inferred a tree from the concatenated alignment. See van Dam et al 2016 (doi: 10.1111/1462-2920.13445).  

I rooted the tree using _F. verticillioides_ as an outgroup.

I also excluded Focub_B2 from our analyses.

*species_tree_outgroup_removed_root.newick*  
As it's name suggests, this file contains the tree described above, with the outgroup pruned away. It is in newick format. 

  
*data4speciestreeplot.txt*
a tab separated file where it is described per species how it needs to plotted (color, shape).  

*exclude.txt*
In Figure 1, we don't show all genomes included in our analyses because of size and readability constraints. 
  
Run   
`python PATH_TO_CLONED_REPO/species_tree/prune_and_plot_species_tree_bootstrap.py -tree species_tree_outgroup_removed_root.newick -datafile data4speciestreeplot.txt -exclude exclude.txt -figfilename species_tree__pruned.svg -treefilename species_tree__pruned.nw`  
to obtain tree that is depicted in Figure 1 (tree is polished in Adobe Illustrator (spheres are removed, lines thickened etc.))    