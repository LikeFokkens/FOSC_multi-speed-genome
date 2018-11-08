# FOSC_multi-speed-genome

## Figure 1  
####1A   
see `LikeFokkens/whole_genome_alignments` for pipelines to align genomes using nucmer, create coordinate files and sort these for plotting.  

See `LikeFokkens/genome-wide_plots` for scripts used for plotting presence-absence plots and `LikeFokkens/species_tree` for scripts used to plot the species tree.

A species tree was inferred as described in van Dam et al. (doi: 10.1111/1462-2920.13445) and saved in newick format. 
  
Plotting the species tree is described in `plotSpeciesTree`, plotting presence-absence, color-coded according to sequence similarity or synteny (length of the alignment) is described in `presence-absencePlots`. Figure 1 was made by cropping the png generated as described in `presence-absencePlots` and positioning this adjacent to the tree. The order of rows in a presence-absence plot can be adjusted so that it matches the species tree.  
  
####1B  
see `map_reads_to_reference.py` for mapping reads to the recipient genome, extracting unmapped reads and mapping those to the genome sequence of the donor strain. Bamfiles with extracted unmapped reads mapped tot he donor genome, with putative PCR duplicates removed can be downloaded from Zenodo (10.5281/zenodo.1479943), with duplications removed. Read densities were plotted on the genome using `plot_read_density_withGnuplot.py` (see `--help` for more detail on how to use this).


## Figure 2
Similar to Figure 1A.

## Figure 3
Up- and downregulated genes, see `notebooks/DEGS.ipnb`
imports `pandas`, `scipy`, `matplotlib`.  
  








