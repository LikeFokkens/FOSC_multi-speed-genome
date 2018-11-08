To make these plots you need  `genome-wide_plots/plot_presence_absence_wrt_referenceGenome_inTree.py`.   
  
For example:
`plot_presence_absence_wrt_referenceGenome_inTree.py -nameReferenceGenome Fol4287broad -queryGenomes Fol4287broad Fol007illumina Fol026 Fol018 Fol029 Fol074 Fol014 Fom001 FOSC-3a Forc016 Fol069 Fo47 Fol051 FolMN25 Fol077 Forl_CL57 Fov_NRRL25433 Fon020 Foq015 Foq037 Foq001 Focub_N2 Fo5176 Focon_PHW808 Fop_HDV247 Fon037 Fon015blob2c Fom009pilon Fom010pilon Foq021 For_PHW815 Fom005pilon Fom006pilon Foq013 Foq011 Focub_II5  -reference_fasta Fol4287broad.fasta -inDir presence-absencePlots/coords/ -outDir presence-absencePlots/ -minSim 90 -minLength 1000 -colorBy length  -palette 7,5,15 -includePreset Fol4287_4speed -name Fol4287_4speed_byLength`
