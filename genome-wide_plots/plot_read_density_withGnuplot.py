import argparse
import os
import sys
import pysam
import glob
import numpy

#sys.path.append(os.getcwd())
#sys.path.append(os.path.dirname(os.getcwd())+'/tools/')

parent_dir  = os.path.dirname(os.path.abspath(os.path.dirname(sys.argv[0])))
PATH_PYTHON_TOOLS = parent_dir+'/tools/'

sys.path.append(PATH_PYTHON_TOOLS)

import fasta_tools, plot_tools



def parse_cmndline_arguments():
	general_parser = plot_tools.argument_parser_for_genome_wide_plots()
	parser = argparse.ArgumentParser(description="Plot read density from one or multiple bamfiles on (part of) a genome, using gnuplot.\n\
		If you want semi-transparent pngcairo output check whether your version of gnuplot supports this.\n\
		You may need to install a newer version of gnuplot and/or it's depencies like libcairo2-dev or libpango1.0-dev \
		(e.g. brew install gnuplot -with-cairo)", parents = [general_parser], conflict_handler='resolve')
	
	parser.add_argument('-bamInSinglePlot', dest='bamInSinglePlot', default = False, action = 'store_true', help='Flag whether, in case of multiple bamfiles, these should be plotted in a single or in separate plots.')
	parser.add_argument('-normalize', dest='normalize', default = False, action = 'store_true', help='if set normalize read density wrt RPM (number of reads per million bases), to compare experiments with different sequencing depths')
	#planned:
	#parser.add_argument('-scaled', dest='scaled', default = False, action = 'store_true', help='If separate plots are generated for a set of bamfiles, scale the plots such that the y-axis range is the same for all plots and read depth can be directly compared.')
	parser.add_argument('-maxDepth', dest='maxDepth', type = int, default = None, help='Set maximum depth: maximum value for y-axis. Comes in handy if you want to directly compare read depth of different plots.')
	parser.add_argument('-sortBam', dest='sortBam', default = False, action = 'store_true', help='Set this flag is bamfiles are not sorted')
	parser.add_argument('-v', '--verbose', dest='verbose', default = False, action = 'store_true', help='Print additional messages')
	args = parser.parse_args()


	
	if args.inFiles == None and args.inDir == None:
		print('You must provide a set of bamfiles or a directory with bamfiles with either -inFiles or -inDir')
		sys.exit()
	
	if args.bamInSinglePlot and len(args.name) == 0:
		print("You must provide the name of the plot if you plot multiple bamfiles in single plot.")
		sys.exit()

	report = '\n\n##################################\n#\n#   SETTINGS\n#'
	argsdict = vars(args)
	for var in argsdict.keys():
		report += '#\t'+var+'\t'+str(argsdict[var])+'\n'
	report += '#\n##################################\n\n'

	os.system('mkdir -p '+args.outDir)

	print(report)

	

	return args, report




def init_genome(args):
	
	id2seq, idlist_fastaorder = fasta_tools.fasta2dict_and_genelist(open(args.Rfasta))

	idlist = args.included
	if args.presetKeyword != None:
		idlist = plot_tools.keyword_to_contigslist(args.presetKeyword)
		
	idlist, id2length, id2xstart = plot_tools.get_idlist_id2length_and_id2xstart(args.Rfasta, idlist = idlist, complete_idlist = args.inclUnPos, min_contig_size = args.minContigSize, exclude = set(args.excluded))

	return idlist, id2length, id2xstart, id2seq


#write read density to a datfile
#idlist     : list of scaffolds in the order you wish them to be not necessarily the order in the fasta file used for mapping)
#id2xstart  : dictionary: scaffoldid --> position x-axis: x-start of scaffoldX is the sum of the lengths of all the scaffolds to the left of scaffoldX
#id2length  : dictionary scaffoldid --> length
#bamfilename: the mapping file in bam format
#datfilename: the name of the output file
def write_datfile(idlist, id2xstart, id2length, bamfilename, datfilename, normalize = False, overwrite = False, verbose = False):
	max_depth = 0

	if not os.path.exists(datfilename) or overwrite:
		datfile = open(datfilename, 'w')
		bamfile = pysam.AlignmentFile(bamfilename, 'rb')
		norm = 1
		if normalize != 0:
			norm = 1000000.0/bamfile.mapped

		for id in idlist:
			if id in bamfile.references:  #when using an idlist that include gaps (e.g. to separate individual chromosomes in a plot), these are of course not in the bamfile
				if verbose:
					print(id)
				prev_pos = 0
				xstart = id2xstart[id]
				iter = bamfile.pileup(id, 0, id2length[id])
				for column in iter:
					for i in range(prev_pos, column.pos):
						datfile.write(str(xstart+i)+'\t0\n')

					datfile.write(str(xstart+column.pos)+'\t'+str(column.n)+'\n')
					if column.n > max_depth: max_depth = column.n
					
					prev_pos = column.pos + 1

				for i in range(prev_pos+1, id2length[id]):
					datfile.write(str(xstart+i)+'\t0\n')
			else:
				datfile.write('\n')

		datfile.close()
	else:
		print(datfilename+" exists. Skipping...")
		return datfilename, None

	return max_depth




def write_datfile_region(cid, start, end, bamfilename, datfilename, normalize = False, overwrite = False, verbose = False):
	max_depth = 0

	if not os.path.exists(datfilename) or overwrite:
		datfile = open(datfilename, 'w')
		bamfile = pysam.AlignmentFile(bamfilename, 'rb')
		norm = 1
		if normalize != 0:
			norm = 1000000.0/bamfile.mapped

		prev_pos = 0
		iter = bamfile.pileup(cid, start, end)
		for column in iter:

			#include positions with zero read density:
			for i in range(prev_pos, column.pos):
				datfile.write(str(i)+'\t0\n')

			datfile.write(str(column.pos)+'\t'+str(column.n*norm)+'\n')
			if column.n*norm > max_depth: max_depth = column.n*norm
			
			prev_pos = column.pos + 1

		for i in range(prev_pos+1, end):
			datfile.write(str(i)+'\t0\n')

		datfile.close()


	else:
		print(datfilename+" exists. Skipping...")
		return datfilename, None

	return max_depth


#Writes average read density per sliding window. Writes it as a bedfile, only the last two columns will be used to plot the data
# output fileformat: 
# chromosome (or contig) \t windowstart \t window end \t windowend mapped to x-axis \t mean read density

#idlist     : list of scaffolds in the order you wish them to be (not necessarily the order in the fasta file used for mapping)
#id2xstart  : dictionary: scaffoldid --> position x-axis: x-start of scaffoldX is the sum of the lengths of all the scaffolds to the left of scaffoldX
#id2length  : dictionary scaffoldid --> length
#bamfilename: the mapping file in bam format
#bedfilename: the name of the output file
def write_bedfile_averagePerSlidingWindow(idlist, id2xstart, id2length, bamfilename, bedfilename, normalize = False, windowsize = 10000, slidesize = 1000, overwrite = False, verbose = False):
	max_depth = 0
	if not os.path.exists(bedfilename) or overwrite:
		bedfile = open(bedfilename, 'w')
		bamfile = pysam.AlignmentFile(bamfilename, 'rb')
		norm = 1.0
		if normalize != 0:
			norm = 1000000.0/bamfile.mapped

		for cid in idlist:
			if cid in bamfile.references:  #when using an idlist that include gaps (e.g. to separate individual chromosomes in a plot), these are of cours enot int he bamfile
				if verbose: print(cid)
				prev_pos = 0
				xstart   = id2xstart[cid] # start position on x-axis
				cidstart = 0 			  # start position in chromosome/contig

				iter = bamfile.pileup(cid, 0, id2length[cid])

				wsize       = min([windowsize,id2length[cid]])  # added on 2015_10_18 to avoid short contigs to have no density at all
				window      = []
				for column in iter:
					for i in range(prev_pos, column.pos):
						if len(window) == wsize:
							meandepth = numpy.mean(window)*norm
							bedfile.write(cid+'\t'+str(cidstart)+'\t'+str(i)+'\t'+str(xstart+i)+'\t'+str(meandepth)+'\n')
							cidstart = column.pos

							if meandepth > max_depth: max_depth = meandepth

							if slidesize == 0:
								window = []
							else: window = window[slidesize:]

						window.append(0)
						
					if len(window) == wsize:
						meandepth = numpy.mean(window)*norm
						bedfile.write(cid+'\t'+str(cidstart)+'\t'+str(column.pos)+'\t'+str(xstart+column.pos)+'\t'+str(meandepth)+'\n')
						cidstart = column.pos

						if meandepth > max_depth: max_depth = meandepth

						if slidesize == 0:
							window = []
						else: window = window[slidesize:]

					window.append(column.n)

					prev_pos = column.pos + 1

				for i in range(prev_pos+1, id2length[cid]):
					if len(window) == wsize:
						meandepth = numpy.mean(window)*norm
						bedfile.write(cid+'\t'+str(cidstart)+'\t'+str(i)+'\t'+str(xstart+i)+'\t'+str(meandepth)+'\n')
						cidstart = column.pos

						if meandepth > max_depth: max_depth = meandepth

						if slidesize == 0:
							window = []
						else: window = window[slidesize:]
					window.append(0)

				# last one, may be a smaller window...
				meandepth = numpy.mean(window)*norm
				bedfile.write(cid+'\t'+str(cidstart)+'\t'+str(id2length[cid])+'\t'+str(xstart+id2length[cid]-1)+'\t'+str(meandepth)+'\n')

				if meandepth > max_depth: max_depth = meandepth


		bedfile.close()
	else:
		print(bedfilename+' exists, will not overwrite. Skipping....')

	return max_depth


#Writes average read density per sliding window of a region. Writes it as a bedfile, only the last two columns will be used to plot the data
# output fileformat: 
# chromosome (or contig) \t windowstart \t window end \t windowend mapped to x-axis \t mean read density

#idlist     : list of scaffolds in the order you wish them to be (not necessarily the order in the fasta file used for mapping)
#id2xstart  : dictionary: scaffoldid --> position x-axis: x-start of scaffoldX is the sum of the lengths of all the scaffolds to the left of scaffoldX
#id2length  : dictionary scaffoldid --> length
#bamfilename: the mapping file in bam format
#bedfilename: the name of the output file
def write_bedfile_region_averagePerSlidingWindow(cid, start, end, bamfilename, bedfilename, normalize = False, windowsize = 10000, slidesize = 1000, overwrite = False, verbose = False):
	max_depth = 0

	if not os.path.exists(bedfilename) or overwrite:
		bedfile = open(bedfilename, 'w')
		bamfile = pysam.AlignmentFile(bamfilename, 'rb')
		norm = 1.0
		if normalize != 0:
			norm = 1000000.0/bamfile.mapped

		
		prev_pos = 1
		iter = bamfile.pileup(cid, start, end)

		wsize       = min([windowsize,id2length[cid]])  # added on 2015_10_18 to avoid short contigs to have no density at all
		window      = []
		for column in iter:
			for i in range(prev_pos, column.pos):
				if len(window) == wsize:
					meandepth = numpy.mean(window)*norm
					bedfile.write(str(i)+'\t'+str(meandepth)+'\n')
					if meandepth > max_depth: max_depth = meandepth
					window = window[slidesize:]
				window.append(0)
				
			if len(window) == wsize:
				meandepth = numpy.mean(window)*norm
				bedfile.write(str(column.pos)+'\t'+str(meandepth)+'\n')
				if meandepth > max_depth: max_depth = meandepth

				window = window[slidesize:]

			window.append(column.n)

			prev_pos = column.pos + 1

		for i in range(prev_pos+1, end):
			if len(window) == wsize:
				meandepth = numpy.mean(window)*norm
				bedfile.write(str(i)+'\t'+str(meandepth)+'\n')
				if meandepth > max_depth: max_depth = meandepth
				window = window[slidesize:]
			window.append(0)

		meandepth = numpy.mean(window)*norm
		bedfile.write(str(end)+'\t'+str(meandepth)+'\n')
		if meandepth > max_depth: max_depth = meandepth


		bedfile.close()
	else:
		print(bedfilename+' exists, will not overwrite. Skipping....')

	return max_depth




def plot_read_density(datfilenames, colors, plotname, idlist, id2xstart, id2length, plotdir, gnudir, out_fmt,  plotwidth=None, plotheight=None, font = 'Arial,12', xtics = 1000000, maxDepth = None):
	
	outfile_base = os.path.basename(datfilenames[0]).split('.bed')[0]
	print(idlist)
	last = idlist[-1]
	maxX = id2xstart[last]+id2length[last]

	gnu = ''
	if out_fmt == 'svg':
		gnu += "set terminal svg font \""+font+"\"\n"
		gnu += 'set size 4,1\n'
	elif out_fmt == 'eps':
		gnu += "set terminal postscript eps enhanced color font \""+font+"\"\n"
		gnu += 'set size 4,1\n'
	elif out_fmt == 'png':
		# gnu = 'set terminal pngcairo size '
		if plotwidth == None: 
			plotwidth  = (0.0005 * maxX)
		if plotheight == None:
			plotheight = 1000

		gnu = 'set terminal pngcairo size '+str(plotwidth)+','+str(plotheight)+' font "'+font+'"\n'
	

	fig_fname = plotdir+'/'+outfile_base+"."+out_fmt
	gnu += "set output '"+fig_fname+"'\n\n"
	gnu += 'set grid\n'
	gnu += 'unset key\n'
	
	
	gnu += plot_tools.get_default_gnuplot_xtics_lines(idlist, id2xstart, maxX, plotwidth, font, ticdist = xtics)

	gnu += 'set xlabel "Position"\n'
	gnu += 'set ylabel "Read depth"\n'
	if maxDepth != None: gnu += 'set yrange[0:'+str(maxDepth)+']\n'
	for line_index, color in enumerate(colors):
		gnu += 'set style line '+str(line_index+1)+' lt 1 lw 2 pt 7 ps 0.2 lc rgb "'+color+'"\n'

	gnu += "plot "
	for line_index, datfilename in enumerate(datfilenames):
		if datfilename.split('.')[-1]=='dat':
			gnu += "'"+datfilename+"' using 1:2 w l ls "+str(line_index+1)+", " 
		else:
			#bedfile with windows
			gnu += "'"+datfilename+"' using 4:5 w l ls "+str(line_index+1)+", " 
	gnu = gnu[:-2] #trim of last ', '

	gnufilename = gnudir+'/'+outfile_base+'.gnu'
	gnufile = open(gnufilename, 'w')
	gnufile.write(gnu)
	gnufile.close()
	
	return gnufilename, fig_fname
	

if __name__ == "__main__":


	args, settings_report  = parse_cmndline_arguments()
	idlist, id2length, id2xstart, id2seq = init_genome(args)
	print(idlist)
	
	bamfiles = args.inFiles
	print(bamfiles)
	if bamfiles == None:
		bamfiles = glob.glob(args.inDir+'/*.bam')

	if len(args.colors) < len(bamfiles):
		print('*** WARNING *** Number of bamfiles exceeds number of colors: some will be plotted in the same color!')

	####################
	#
	#  get datfiles
	#
	####################
	datdir = None
	if args.windowSize > 1:
		datdir  = args.outDir+'/BEDFILES'
	else:
		datdir  = args.outDir+'/DATFILES'
	os.system('mkdir -p '+datdir)

	datfilenames = []
	max_depths   = [] #keep maximum depths even if it is not used yet (normally I normalize to compare or manually set the max depth)

	for bfname in bamfiles:
		if not os.path.exists(bfname+'.bai'):
			if args.sortBam:
				sorted_dir    = os.path.dirname(bfname)+'/sorted/'
				os.system('mkdir -p '+sorted_dir)
				sorted_bfname = sorted_dir+os.path.basename(bfname).replace('.bam', '.sorted')
			
				cmnd = 'samtools sort '+bfname+' '+sorted_bfname
				if args.verbose:
					print(cmnd)
				if os.system(cmnd) == 0:
					bfname = sorted_bfname+'.bam'
			
			cmnd = 'samtools index '+bfname
			if args.verbose: print(cmnd)
			os.system(cmnd)
					
					
		max_depth_bf = 0

		datfilename = datdir + '/'
		if len(args.name) > 0:
			datfilename += args.name+'__'
		datfilename += os.path.basename(bfname).replace('.bam', '')

		if args.normalize != 0:
			datfilename += '.normalizedRPM'

		if args.start != 0 or args.end != 'end':
			cid = args.included[0]
			datfilename += '.'+cid+'_'+str(args.start)+'-'+str(args.end)
			if args.end == 'end': args.end = len(id2seq[cid])

			if args.windowSize > 1:
				datfilename += '.w'+str(args.windowSize)
				if args.stepSize > 0:
					datfilename += '-'+str(args.stepSize)
				datfilename += '.bed'
				max_depth_bf = write_bedfile_region_averagePerSlidingWindow(cid, start, end, bfname, datfilename, normalize=args.normalize, windowsize = args.windowSize, slidesize = args.stepSize, overwrite = args.overwrite, verbose = args.verbose)

			else:
				datfilename += '.dat'
				max_depth_bf = write_datfile_region(cid, start, end, bfname, datfilename, normalize=args.normalize, overwrite = args.overwrite, verbose = args.verbose)


		elif args.windowSize > 1:
			datfilename += '.w'+str(args.windowSize)
			if args.stepSize > 0:
				datfilename += '-'+str(args.stepSize)
			datfilename += '.bed'
			
			max_depth_bf = write_bedfile_averagePerSlidingWindow(idlist, id2xstart, id2length, bfname, datfilename, normalize=args.normalize, windowsize = args.windowSize, slidesize = args.stepSize, overwrite = args.overwrite, verbose = args.verbose)
			
		else:
			datfilename += '.dat'
			max_depth_bf = write_datfile(idlist, id2xstart, id2length, bfname, datfilename, normalize=args.normalize, overwrite = args.overwrite, verbose = args.verbose)

		if args.verbose: print(datfilename+' '+str(max_depth_bf))

		datfilenames.append(datfilename)
		max_depths.append(max_depth_bf)



	####################
	#
	#  make plots
	#
	####################
	gnudir  = args.outDir+'/GNUFILES'
	os.system('mkdir -p '+gnudir)
	plotdir  = args.outDir+'/PLOTS'
	os.system('mkdir -p '+plotdir)

	#if not os.path.exists(args.outDir+'/README'):
	readme = open(args.outDir+'/README', 'a')

	if args.bamInSinglePlot:
		gnufilename, fig_fname = plot_read_density(datfilenames, args.colors, args.name, idlist, id2xstart, id2length, plotdir, gnudir, args.outfmt, args.figwidth, args.figheight, maxDepth = args.maxDepth)
		os.system('gnuplot '+gnufilename)
		log = fig_fname+' generated by '+gnufilename+' based on:\n'
		for dat_fname in datfilenames:
			log += dat_fname+'\n'
		log += '\ngenerated by '+__file__+'\n'
		log += 'SETTINGS:\n'+settings_report+'\n'
		readme.write(log)
		
	else:

		for index, datfname in enumerate(datfilenames):	
			gnufilename, fig_fname = plot_read_density([datfname], [args.colors[min([len(args.colors)-1, index])]], args.name, idlist, id2xstart, id2length, plotdir, gnudir, args.outfmt, args.figwidth, args.figheight, maxDepth = args.maxDepth)
			os.system('gnuplot '+gnufilename)

			log = fig_fname+' generated by '+gnufilename+' based on:\n'
			for dat_fname in datfilenames:
				log += dat_fname+'\n'
			
			log += '\ngenerated by '+__file__+'\n'
			log += 'SETTINGS:\n'+settings_report+'\n'
			readme.write(log)

	readme.close()




	
	
	
