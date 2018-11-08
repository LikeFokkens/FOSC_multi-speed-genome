import sys, os, argparse, numpy
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import fasta_tools, plot_tools


###########################
#
#  parse options
#
###########################
def parse_args():

	parser = argparse.ArgumentParser(description='Calculate average densities, or gene-based values (e.g. expression levels or dN/dS) per sliding window on a genome. Implements bedtools and/or samtools.')
	parser.add_argument('-windowsize', dest='windowSize', type = int, default = 10000, help = 'Plot average read density per window: provide the window size here.')
	parser.add_argument('-stepsize', dest='stepSize', type = int, default = 0, help = 'If you want a sliding window, you can provide the slide size here. If not provided, windows will be non-overlapping.')

	input_options1 = parser.add_argument_group('Input options for the reference genome. You can provide a fasta- or genome file for the reference genome, then a bedfile with windows will be generated.\n\
You can also provide a bedfile with windows directly.')
	input_options1.add_argument("-fasta", dest='fasta_fname', type = str, default = None, help='Name of the reference genome fastafile.')
	input_options1.add_argument("-genome", dest='genome_fname', type = str, default = None, help='Name of the reference genome file (= tabulated file: Supercontig{tab}length, used in bedtools often).')
	input_options1.add_argument("-windows_bed", dest='windows_fname', type = str, default = None, help='Name of the bedfile with sliding windows on the reference genome.')


	input_options2 = parser.add_argument_group('Input options for what you need to calculate averages per window for.\n\
Possible input: the bam- or bedfile to determine read density per window, or amount of overlap per window, or provide a tabulated file with values per gene.')
	input_options2.add_argument("-bam", dest='bam_fname', type = str, default = None, help='Name of the bamfile containing reads that are mapped to the genome: will return a file with average read density per window.')
	input_options2.add_argument("-bed", dest='bed_fname', type = str, default = None, help='Name of the bedfile containing regions: will return a file with average read density per window.')	
	input_options2.add_argument("-gene2value", dest='gene2value_fname', type = str, default = None, help='Name of the tabulated file containing values per gene: will return a file with mean, min, max, .25 and .75 of gene values per window. \
		A bedfile with the positions of each gene on the genome is required.')
	input_options2.add_argument("-genebed", dest='gene2pos_fname', type = str, default = None, help='Name of the bedfile with the positions of each gene on the genome (required when calculating average values per gene.')

	output_options = parser.add_argument_group('Output options')
	#output_options.add_argument("-outfilename", dest='out_fname', type = str, default = None, help='Name of the file with the output.')
	output_options.add_argument("-name", dest='name', type = str, default = None, help='Name of the analysis, will be used in output files.')

	args = parser.parse_args()

	#check arguments
	if args.fasta_fname == None and args.genome_fname == None and args.windows_fname == None: 
		print('Please provide a reference genome, with either -fasta, -genome or -windows_bed. See help.')
		sys.exit()

	if args.bed_fname == None and args.bam_fname == None and args.gene2value_fname == None: 
		print('Please provide a file with something to calculate the average per window of, with either -bam, -bed or -gene2value combined with -genebed. See help.')
		sys.exit()

	if args.gene2value_fname != None and args.gen2pos_fname == None:
		print('Please provide bedfile with positions of the genes on the genome. See help.')
		sys.exit()

	return args


	

#####################
#
#  run
#
#####################
def genome2windows(genome_fname, windowsize, stepsize):
	# make windowfile
	wname = windowsize
	if stepsize != "0": wname += '-'+stepsize

	windows_fname = genome_fname.replace("contig2size.tab", wname+'_windows.bed')
	result = 0
	if not os.path.exists(windows_fname):
		cmnd = cmnd = 'bedtools makewindows -g '+genome_fname+' -w '+windowsize+' > '+windows_fname
		if stepsize != '0':
			cmnd = 'bedtools makewindows -g '+genome_fname+' -w '+windowsize+' -s '+stepsize+' > '+windows_fname
		print('Make sliding windows on genome:'+windows_fname)
		result = os.system(cmnd)
		print(cmnd, result)


#cmnd = "sed -i 's/Fol4287broad_contig_/Supercontig_/g' "+windows_fname
#print cmnd, os.system(cmnd)

# calculate density wrt bedfile (can also be a gtf-file)
def density_bedfile(windows_fname, wname, bed_fname):

	density_fname = bed_fname.replace(bed_fname.split('.')[-1], 'density_'+wname+'_windows.bedtools.coverage.out')
	
	cmnd = ''
	if not os.path.exists(density_fname):
		cmnd = 'bedtools coverage -a '+windows_fname+' -b '+bed_fname+" | awk -v OFS='\\t' '{print $1,$2,$3,$7}' > "+density_fname
		print('Coverage within sliding windows', density_fname)
		print(cmnd, os.system(cmnd))


	else:
		print('WARNING!', density_fname, 'exists, will not overwrite')
		
	return density_fname
		


#calculate density wrt bamfile
def density_bamfile(windows_fname, wname, bam_fname):
	if not os.path.exists(bam_fname+'.bai'):
		os.system('samtools index '+bam_fname)
		
	density_fname = bam_fname.replace('.bam', '.density'+wname+'_windows.samtools_bedcov.avgPerWindow.out')
	if not os.path.exists(density_fname):
		''' samtools bedcov calculates the accumulative coverage: sum of all the per-base coverage.
		We use awk to divide this sum ($4) by the windowsize (($3-$2)) to obtain the average rather than the sum.
		'''
		cmnd = 'samtools bedcov '+windows_fname+' '+bam_fname+" | awk -v OFS='\\t' '{print $1,$2,$3,$4/($3-$2)}'  > "+density_fname
		print('bedcov', density_fname)
		result = os.system(cmnd)
		print(cmnd, result)

		if result != 0:
			print('Error in', cmnd,':', result)
			return None
	else:
			print('WARNING!', density_fname, 'exists, will not overwrite')

	return density_fname


# put gene locations in a dictionary per scaffold
def genebed2dict(gene2pos_fname):
	
	scaffold2start2gene = {}
	for line in open(gene2pos_fname).readlines():
		scaffold, start, end, gene = line.strip().split()

		if scaffold not in scaffold2start2gene:  scaffold2start2gene[scaffold] = {}

		start = int(start)
		end   = int(end)
		gene_start = start
		if start > end:
			gene_start = end

		if gene_start in scaffold2start2gene[scaffold]:
			scaffold2start2gene[scaffold][gene_start].add(gene.split('T')[0])  # just in case two genes happen to start at the same site,
																		  # but make it a set because we don't want genes with multiple transcripts to bias the average
		else:
			scaffold2start2gene[scaffold][gene_start] = set([gene.split('T')[0]])

	return scaffold2start2gene

# make a dictionary: per window, the genes that are located in that window (based on their start sites)
def get_window2genes(windows_fname, scaffold2start2gene):
	
	window2genes     = {}
	current_scaffold = None
	gene_starts      = []
	index            = 0
	processed_scaffolds = set([])
	for line in open(windows_fname).readlines():
		window = line.strip()
		window2genes[window] = []

		scaffold, start, end = window.split()
		wstart = int(start)
		wend   = int(end)

		# if we start at a new scaffold, get the positions of all the genes that lie on the scaffold, and order them for efficiency
		if scaffold != current_scaffold:
			if scaffold in processed_scaffolds:
				print("ERROR window2genes: Please provide a windows-file ordered per Supercontig! Use sort -k1,1 -k2,2n!")
				sys.exit()
			else:
				processed_scaffolds.add(scaffold)

			if scaffold in scaffold2start2gene:
				gene_starts = list(scaffold2start2gene[scaffold].keys())
				gene_starts.sort()
				index = 0
			else:
				gene_starts = []

			current_scaffold = scaffold

		if len(gene_starts) > 0:			
			for gs in gene_starts[index:]:  #[index:] : discard genes that lie in previous windows
				#print scaffold, wstart, wend, gs
				if gs >= wstart:
					if gs < wend:
						for gene_id in scaffold2start2gene[scaffold][gs]:
							window2genes[window].append(gene_id)
						index += 1 #move to next gene_start

					else: break #gene falls into next window 

				else:	
					index += 1
					print('Error with gene indices in window2genes!') # this should not occur, given keeping of the indices

	return window2genes




def average_gene_values_per_window(windows_fname, wname, gene2value_fname, window2genes):

	# put values per gene in a dictionary
	gene2value = {}
	for line in open(gene2value_fname).readlines():
		tabs = line.strip().split()
		gene2value[tabs[0]] = float(tabs[1].replace(',', '.'))

	# determine for each window, which genes lie in it (based on start sites) and what is the average (median...etc.) value for the genes that lie in the window
	out_fname = gene2value_fname.replace('.tab', '.'+wname+'_windows')
	if not os.path.exists(out_fname):
		outfile   = open(out_fname, 'w')
		
		for window in window2genes.keys():
			values = []

			for gene in window2genes[window]:
				if gene in gene2value: values.append(gene2value[gene])

			
			if len(values) > 0:
				outfile.write(window+'\t'+str(numpy.mean(values))+'\t'+str(min(values))+'\t'+str(max(values))+'\t'+str(numpy.percentile(values, 25))+'\t'+str(numpy.percentile(values, 75))+'\n')
			else:
				outfile.write(window+'\t\t\t\t\t\n')

		outfile.close()
	else:
			print('WARNING!', out_fname, 'exists, will not overwrite')
			
	return out_fname





if __name__ == "__main__":

	args = parse_args()

	# make 'genome-file' that contains the size for each scaffold:	
	if args.fasta_fname != None and args.genome_fname == None and args.windows_fname == None:
		genome_fname = fasta_fname.replace('.fasta', '.contig2size.tab')
		if not os.path.exists(genome_fname):
			genome_fname = fasta_tools.fasta2genomefile(fasta_fname, genome_fname = genome_fname, outdir = None)

	if args.windows_fname == None:
		args.windows_fname = genome2windows(args.genome_fname, str(args.windowSize), str(args.stepSize))

	wname = str(args.windowSize)
	if args.name != None:
		wname = args.name

	if args.bed_fname != None:
		density_bedfile(args.windows_fname, args.wname, args.bed_fname)

	if args.bam_fname != None:
		density_bamfile(args.windows_fname, args.wname, args.bam_fname)

	if args.gene2value_fname != None:

		scaffold2start2gene = genebed2dict(args.gene2pos_fname)
		window2genes 		= get_window2genes(windows_fname, scaffold2start2gene)
		outfname     		= average_gene_values_per_window(args.windows_fname, wname, args.gene2value_fname, window2genes)


