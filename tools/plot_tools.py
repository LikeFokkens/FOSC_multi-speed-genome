import sys, os, math, argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import fasta_tools


preset_keywords = ['Fol4287_3speed', 'Fol4287_core', 'Fol4287_fast-core', 'Fol4287_accessory', 'Fol4287_LS', 'Fol4287_pathogenicity', 'Fol4287_4speed', 'Fol4287_perChr', 'Fol4287_perChr_noGAP']
chrs = ['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
for chr in chrs:
	preset_keywords.append('Fol4287_'+chr)


def argument_parser_for_genome_wide_plots():
	parser = argparse.ArgumentParser()

	parser.add_argument('-gnupath', dest='gnupath', default = '', help="Path to gnuplot in case it is not invokable by just typing 'gnuplot'")

	group_density = parser.add_argument_group('Arguments for density plots (sliding windows)')
	group_density.add_argument('-windowsize', dest='windowSize', type = int, default = 10000, help = 'Plot average read density per window: provide the window size here.')
	group_density.add_argument('-stepsize', dest='stepSize', type = int, default = 0, help = 'If you want a sliding window, you can provide the slide size here. If not provided, windows will be non-overlapping.')
	
	group_input = parser.add_argument_group('Generic input genome-wide plots')

	group_input.add_argument("-referenceFasta", dest='Rfasta', type = str, help='Name of the reference genome fastafile')
	group_input.add_argument("-includePreset", dest='presetKeyword', choices=preset_keywords, default = None, help='Presets for plotting the Fol4287 genome.\
		To plot the Fol4287 genome with gaps between the core, fast-core and accessory genome, type "Fol4287_3speed".\n\
		To plot the Fol4287 genome with gaps between the core, fast-core, lineage-specific and pathogenicity genome, type "Fol4287_4speed".\n\
		For Plot with just the core, fast-core, accessory (lineage-specific + pathogenicity), lineage-specific, or pathogenicity chromosomes, \
		type Fol4287_core, Fol4287_fast-core, Fol4287_accessory, Fol4287_LS or Fol4287_pathogenicity respectively. Fol4287_chr01 will plot only chr01, etc.')

	group_input.add_argument("-include", dest='included', nargs='*', default = [], help='List of contigs in the reference genome, in the order in which they should appear on the x-axis.\
		You can add gaps by adding a fake "gap" contig: __GAP#<gapsize>.\n')
	group_input.add_argument('-exclude', dest='excluded', nargs='+', default = [], help="In case you want to exclude part of the genome you can provide a list of supercontigs here")
	group_input.add_argument("-inDir", dest='inDir', default = None, help='Name of the directory with the inputfiles.')
	group_input.add_argument("-inFiles", dest='inFiles', nargs='*', default = None, help='Names of the inputfiles, full path or when -inDir is specified, path from inDir')
	group_input.add_argument("-postfix", dest='postfix', default = '.bed', help='Postfix of input files, e.g. .bam, .bed, .nucmer_maxmatch.coords')
	group_input.add_argument('-start', dest='start', type = int, default = 0, help="In case you only want to plot part of the genome you can define a region by providing a supercontig with -include and a the start of the region here and the end of the region with -end. If only an end is provided, start is set at 0")
	group_input.add_argument('-end', dest='end', default = 'end', help="In case you only want to plot part of the genome you can define a region by providing a supercontig with -include and a the start of the region with -start and the end of the region here. If only a start is provided, this is set at the length of the contig in -include")
	
	group_output = parser.add_argument_group('Generic output genome-wide plots')
	group_output.add_argument("-outDir", dest='outDir', default = os.path.abspath(os.path.dirname(sys.argv[0])), type = str, \
		help='Name of the directory where the output will be saved.	Default is directory of this code:'+os.path.abspath(os.path.dirname(sys.argv[0])))
	group_output.add_argument("-name", dest='name', type=str, default = '', help='Name for the plot.')
	group_output.add_argument('-outfmt', dest='outfmt', choices = ['eps', 'svg','png'], default = 'png', help = 'Output format of the density plot. Png is usually recommended.')
	group_output.add_argument('-overwrite', dest='overwrite', default = False, action = 'store_true', help='Force overwriting existing files')
	group_output.add_argument('-sizeSorted', dest='sizeSorted', default = False, action = 'store_true', help='Sort contigs according to size (longest first).')
	group_output.add_argument('-minContigSize', dest='minContigSize', default = 0, help='Minimum size of contigs to be included in the plot.')
	#Figure settings genome-wide plots
	group_figfmt = parser.add_argument_group('Generic figure formatting genome-wide plots')
	group_figfmt.add_argument("-palette", dest='palette', type = str, default = '33,13,10', help='Specification of the colormap that is used: see Gnuplot documentation for more detail,\
		e.g. http://gnuplot.sourceforge.net/demo/pm3dcolors.html')
	group_figfmt.add_argument('-colors', dest='colors', nargs='+', default = ['#00008B', '#008B8B', '#B22222', '#006400', '#9400D3', '#000000'], help='(List of) Hex colorcode(s), defaults are darkblue, darkcyan, firebrick, darkgreen, darkviolet and black')
	group_figfmt.add_argument("-figwidth", dest='figwidth', type = int, default = None, help='Preferred width of the resulting png. If you set this to very large, please also adjust the font (-font). If not set, the width is adjusted to the length of the sequence that is plotted.')
	group_figfmt.add_argument("-figheight", dest='figheight', type = int, default = None, help='Preferred width of the resulting png. If you set this to very large, please also adjust the font (-font).  If not set, the height is adjusted to the number of query genomes that are plotted.')
	group_figfmt.add_argument("-font", dest='font', type = str, default = 'Arial,12', help='Preferred font: consists of FontName,size. Check Gnuplot documentation for available fonts.')
	group_figfmt.add_argument("-bptics", dest='tics', type = str, default = None, help="Set tics at regular intervals, put interval size as number of bp. Tics will be marked at top of the plot, \
		unless only a single contig is plotted, then they will be at the bottom.")

	group = parser.add_argument_group('Generic options genome-wide plots')
	group.add_argument("-includeUnpScaff", dest='inclUnPos', default = False, help="Whether to include unpositioned scaffolds in the plot. They will be size-sorted.", action="store_true")

	# args = parser.parse_args()
	# ## run checks on input:
	# if args.Rfasta == None:
	# 	print 'You must provide a fastafile with the genome sequence with -reference_fasta'
	# 	sys.exit()

	# if (args.start != 0 or args.end != 'end') and args.included == []:
	# 	print "You must provide the name of the supercontig that you want to plot. If you want to plot everything in the reference_fasta, don't give a value for -start and -end."
	# 	sys.exit()

	return parser

#names		: list of names that have to put on the axis (e.g. scaffold ids)
#name2axis_pos	: dictionary of where each name should go on the axis
#axisrange	    : length of the axis
#plotwidth: width of the plot in pixels
#x_or_y		    : 'x' or 'y' (or 'z' even) what is the names of the axis
#readable	    : if True, it will skip some names if they are to close by, to avoid them printed on top of eachother making the axis unreadable
def get_gnuplot_tic_lines(names, name2axis_pos, x_or_y, axis_range, rotate = False, trim_name = True, name2strand = None, readable = False, plotwidth = None, fontsize = 10, middle = False, scalar = 1000):
	
	#what should be put on the x-axis:
	ticlines = 'set '+x_or_y+'tics '
	if rotate: ticlines += 'rotate '
	ticlines += '( \\\n'
	for i, name in enumerate(names):
			
		pos_on_axis = name2axis_pos[name]
		
		if trim_name:
			name = name.split('contig')[-1]
			name = name.split('#')[0]
			name = name.replace('scf718000000', '').replace('Supercontig_', '').replace('|quiver', '').replace('_', '').replace('_', '')
		
		if name2strand != None and name in name2strand.keys():
			name += name2strand[name]
	
		if not readable or i == len(names)-1:
			ticlines += ' "'+name+'" '+str(pos_on_axis)+', \\\n'

		else:	#do not put the name if the distance to the next name is too small, it will not be visible anyway...

			dist_unit        = axis_range/float(plotwidth)              # how many bp are 'covered' in one pixel (approximately)
			dist2next        = name2axis_pos[names[i+1]] - pos_on_axis  # how many bp to next tic
			dist2next_pixels = dist2next / dist_unit                    # how many pixels to next tic
			#print dist2next_pixels * scalar, fontsize
			if fontsize > scalar * dist2next_pixels:                    # is this enough space given the fontsize?
				ticlines += ' "" '+str(pos_on_axis)+', \\\n'
			else:	
				ticlines += ' "'+name+'" '+str(pos_on_axis)+', \\\n'
			
	if axis_range != None:		
		ticlines += ' "" '+str(axis_range)+')\n\n\n'
		ticlines += 'set '+x_or_y+'range[0:'+str(axis_range)+']\n'

	else:	ticlines = ticlines[:-4]+')\n\n\n'
	
	ticlines += '\n'


	return ticlines


def get_default_gnuplot_xtics_lines(idlist, id2xstart, maxX, plotwidth, font, ticdist = None, Ntics = 10):

	# set readable axis labels:
	# X-axis
	gnu = ''
	if len(idlist) > 1:
		fontsize = int(font.split(',')[-1])
		gnu += get_gnuplot_tic_lines(idlist, id2xstart, 'x', maxX, rotate = True, trim_name = True, name2strand = None, readable = True, plotwidth = int(plotwidth), fontsize = fontsize)	
		
		# put secondary tics on top:
		if ticdist != None:
			gnu += 'set link x2\n'
			x2tic   = ticdist
			gnu += 'set x2tics scale 0.5 ("" '+str(x2tic)
			x2tic += ticdist
			while x2tic < maxX:
				gnu += ', "" '+str(x2tic)
				x2tic += ticdist
			gnu += ')\n'
	else:
		# if no distance between xtics is given, put default from gnuplot
		if ticdist != None: gnu += "set xtics "+str(ticdist)+"\n"
	
	gnu += '\n'
	gnu += 'set xrange[0:'+str(maxX)+']\n'
	
	return gnu + '\n'



def get_gnuline_relativesize_png(max_length, length, align = 'left', axis = 'x'):
	'''
	These lines can used to generate a multiplot, where individual subplots are scaled to 'length' 
	(e.g. the length of individual chromosomes). To achieve this, we need to adjust the right margin
	(to align subplots at the left) or the left margin (to align subplots at the right)
	'''

	whitespace = (max_length - length)/float(max_length)
	if axis == 'x':
		if align == 'left':
			return 'set rmargin at screen '+str(1-whitespace)+'\n'

		return 'set lmargin at screen '+str(whitespace)+'\n'
	else:
		return 'set tmargin at screen '+str(1-whitespace)+'\n'




def keyword_to_contigslist(keyword, prefix = 'Supercontig_2.', includeUnposScaffolds = False, gapsize = 500000):

	
	# if keyword does not exist, assume it is not a keyword but just the name of a single contig
	if keyword not in preset_keywords:
		return [keyword]
	else:
		
		chromosomes = ['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', \
		'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
		# get list of scaffols per category
		chr_categories = ['core', 'fast-core', 'accessory', 'lineage-specific', 'pathogenicity']
		cat2chrs = {}

		cat2chrs['core']             = ['chr01_C','chr02_C','chr04','chr05','chr07','chr08','chr09','chr10',]
		cat2chrs['fast-core']        = ['chr11','chr12','chr13']
		cat2chrs['accessory']        = ['chr01_LS', 'chr02_LS', 'chr03P', 'chr03A', 'chr06A','chr06P','chr14','chr15']
		cat2chrs['lineage-specific'] = ['chr01_LS', 'chr02_LS', 'chr03A', 'chr06A','chr15']
		cat2chrs['pathogenicity']    = ['chr14', 'chr03P', 'chr06P']
  
		chr2scaffolds = {}
		chr2scaffolds['chr01']    = [14, 1, 27]
		chr2scaffolds['chr01_C']  = [14, 1]
		chr2scaffolds['chr01_LS'] = [27]
		chr2scaffolds['chr02']    = [6,10, 31]
		chr2scaffolds['chr02_C']  = [6, 10]
		chr2scaffolds['chr02_LS'] = [31]
		chr2scaffolds['chr03']    = [47,18,32,7,25]
		chr2scaffolds['chr03A']   = [32, 7, 25]
		chr2scaffolds['chr03P']   = [47, 18]
		chr2scaffolds['chr04']    = [8, 4]
		chr2scaffolds['chr05']    = [26, 2]
		chr2scaffolds['chr06']    = [9, 33, 41, 21, 53, 42]
		chr2scaffolds['chr06A']   = [9, 33]
		chr2scaffolds['chr06P']   = [41, 21, 53, 42]
		chr2scaffolds['chr07']    = [5, 13]
		chr2scaffolds['chr08']    = [3, 29]
		chr2scaffolds['chr09']    = [11, 17]
		chr2scaffolds['chr10']    = [20, 15, 45]
		chr2scaffolds['chr11']    = [35, 12]
		chr2scaffolds['chr12']    = [19, 23]
		chr2scaffolds['chr13']    = [16, 39]
		chr2scaffolds['chr14']    = [22, 43, 51, 36]
		chr2scaffolds['chr15']    = [37, 38, 24, 28]

		positioned_scaffolds_list = []
		for chr in chromosomes:
			positioned_scaffolds_list += chr2scaffolds[chr]
		positioned_scaffolds = set(positioned_scaffolds_list)

		#size-sorted unpositionaed scaffolds
		unpositioned_scaffolds = list(set(range(1, 115)).difference(positioned_scaffolds))
		unpositioned_scaffolds.sort()

		# get idlists for plotting
		idlist = []
		if keyword == 'Fol4287_core':
			for chr in cat2chrs['core']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			

		elif keyword == 'Fol4287_fast-core':
			for chr in cat2chrs['fast-core']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			

		elif keyword == 'Fol4287_accessory':
			for chr in cat2chrs['accessory']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			

		elif keyword == 'Fol4287_LS':
			for chr in cat2chrs['lineage-specific']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			

		elif keyword == 'Fol4287_pathogenicity':
			for chr in cat2chrs['pathogenicity']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			

		elif keyword == 'Fol4287_3speed':
			for chr in cat2chrs['core']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP1#'+str(gapsize))

			for chr in cat2chrs['fast-core']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP2#'+str(gapsize))

			for chr in cat2chrs['accessory']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP3#'+str(gapsize))
			

		elif keyword == 'Fol4287_4speed':
			for chr in cat2chrs['core']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP1#'+str(gapsize))

			for chr in cat2chrs['fast-core']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP2#'+str(gapsize))

			for chr in cat2chrs['lineage-specific']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP3#'+str(gapsize))

			for chr in cat2chrs['pathogenicity']:
				for scaffold in chr2scaffolds[chr]:
					idlist.append(prefix+str(scaffold))
			idlist.append('__GAP4#'+str(gapsize))

			

		elif keyword == 'Fol4287_perChr':
			# use chromosomes as they are, without e.g. separating chr01b and chr02b.
			
			for chr in chromosomes:
			    for s in chr2scaffolds[chr]:
			        idlist.append(prefix+str(s))
			    idlist.append('__GAP'+chr+'#'+str(gapsize)) 
			
		elif keyword == 'Fol4287_perChr_noGAP':
			# use chromosomes as they are, without separating chr01b and chr02b.
			for chr in chromosomes:
			    for s in chr2scaffolds[chr]:
			        idlist.append(prefix+str(s)) 
			return idlist

		elif keyword[:11] == 'Fol4287_chr':
			for s in chr2scaffolds[keyword.split('_')[-1]]:
			    idlist.append(prefix+str(s))
		
		if includeUnposScaffolds:
			for sc in unpositioned_scaffolds:
				idlist.append('Supercontig_2.'+str(sc))

		return idlist



def get_idlist_id2length_and_id2xstart(fasta_fname, idlist = [], size_sorted = True, complete_idlist = True, min_contig_size = 0, splitheader = True, exclude = set([])):

	'''
	Based on the fastafile (to get the length of sequences), an idlist with the preferred sequence of scaffold (optional), a set of contig ids that needs to be excluded (optional),
	whether only contigs of a certain minimum size should be included and size sorted, this method returns:

	The ordered list of contig IDs as they will appear on the axis (left-to-right) in case of X, bottom-to-top in case of Y)
	A dictionary contig_id -> length of the contig
	A dictionary contig_id -> axis position, the position where the contig will start on the X or Y-axis.

	If there is an idlist provided with a the preferred sequence of contig_ids as they should appear ont e axis, and complete_idlist = True, the rest of the contigs in the fasta will 
	be plotted as welll, after the contig_ids in idlist. If size_sorted is True, these will be sorted accoridng to size, largest first.
	Gaps can be added by ading '__GAPn__#gapsize' to the idlist, where n is a number to distinguish one gap from the other

	'''
	
	id2seq, idlist_fastaorder = fasta_tools.fasta2dict_and_genelist(open(fasta_fname))

	# get dictionary contig_id -> length
	# and a dictionary length -> contig_ids (a list, because we can have contigs with the same length)
	# the latter will be used to size-sort contigs if necessary
	id2length  = {}
	length2ids = {}
	for cid in id2seq.keys(): 	
		seqlength      = len(id2seq[cid])
		id2length[cid] = seqlength
		
		if seqlength in length2ids: 
			length2ids[seqlength].append(cid)
		else:	length2ids[seqlength]=[cid]



	# gaps can be added by ading '__GAPn__#gapsize' to the idlist
	# we need to add these to the id2length, and at the same time we can check whether all ids that are provided in idlist are actually valid (i.e. either a gap or a sequence that exists in the fastaheader)
	for cid in idlist:

		if cid in exclude: 		
			print('\n\n*** WARNING\n'+cid+' also occurs in the set of contigs to exclude. I will include anyway.')
		
		if cid not in id2seq:
			if cid[:5] == '__GAP':
				id2length[cid] = int(cid.split('#')[-1])
			else:
				print('\n\n*** WARNING\n'+cid+' does not occur in the fasta, please check whether it occurs '\
					+'in the fasta '+fasta_fname+' or not and whether you split the header (at first whitespace, as BLAST does as well) or not')


	# Complete provided idlist with the rest of the contigs in the fastafile, or obtain a list of ids
	if idlist == [] or complete_idlist:
		idset = set(idlist) # these ids are in the list that is provided and don't need to be added to the idlist
							# if idlist is empyt this will jsut be an empty set and have no effect

		if size_sorted:
			lengths = list(length2ids.keys())
			lengths.sort()
			lengths = lengths[::-1]

			for l in lengths:
				if l > min_contig_size:
					for cid in length2ids[l]:
						if cid not in idset and cid not in exclude:
							idlist.append(cid)
		else:
			for cid in idlist_fastaorder: #if not sorted accoriding to size, we will just plot the contigs in the same order as tehy occur in the fasta file
				if cid not in idset and cid not in exclude:
					if cid in id2length.keys() and id2length[cid] > min_contig_size:
						idlist.append(cid)		

	
	id2xstart = {}
	id2xstart[idlist[0]] = 1
	xstart = 0
	for index, id in enumerate(idlist[1:]):
		xstart += id2length[idlist[index]] # add length from previous supercontig
		id2xstart[id] = xstart
		
	return idlist, id2length, id2xstart



if __name__ == "__main__":

	print('Test plot for Fol4287 chromosome sizes')

	parser = argument_parser_for_genome_wide_plots()
	args   = parser.parse_args()


	id2seq, idlist = fasta_tools.fasta2dicts(open(args.Rfasta))
	
	chr_lengths = []
	for chr in chrs:
		length =  0
		for cid in chr2scaffolds[chr]:
			length += len(id2seq['Supercontig_2.'+str(cid)])

		chr_lengths.append(length)

	gnu_fname = args.outDir+'/test_rightAlign.gnu'
	gnufile   = open(gnu_fname, 'w')

	gnu = 'set terminal pngcairo size '+str(args.figwidth)+', '+str(args.figheight)+' font "Arial, 8"\n'
	gnu += 'set output "'+args.outDir+'/test_rightAlign.png\n\n'

	gnu += 'set multiplot layout 16, 1\n'
	gnu += 'f(x) = 1\nset yrange[0:2]\n'

	max_length = chr_lengths[0]
	for length in chr_lengths:
		gnu += get_gnuline_relativesize_png(max_length, length)
		gnu+= 'plot f(x) with lines\n\n'

	gnufile.write(gnu)
	gnufile.close()

	os.system('gnuplot '+gnu_fname)



	
