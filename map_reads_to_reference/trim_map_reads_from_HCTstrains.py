import os, sys, glob, argparse

def init():
	parser = argparse.ArgumentParser(description='Trim and map sequencing reads from putative Horizontal Chromsome Transfer strains to the acceptor and donor genomes to identify which genomic regions have been transferred')
	parser.add_argument('task', type=str, choices = ['all','trim', 'map2acceptor', 'mapall2donor', 'extractUnmappedReads','map2donor'], default = 'all', help='which step has to be performed')
	parser.add_argument('-v', dest='verbose', action = 'store_true', default = True, help='print messages')
	parser.add_argument('-readDir', dest='readDir', type=str, help='name of the directory with fastq or fastq.gz files')
	parser.add_argument('-outDir', dest='outDir', type=str, help='name of the output directory, output will be put in outDir/01.Trimming/, outDir/02.MappingToAcceptor/, outDir/03.UnmappedReads/, outDir/04.MappingToDonor/')
	parser.add_argument('-adapterFasta', dest='adapterFasta', type=str, help='fastafile with adapter sequences')
	parser.add_argument('-donorFasta', dest='donorFasta', type=str, help='fastafile with the genome of the donor strain')
	parser.add_argument('-acceptorFasta', dest='acceptorFasta', type=str, help='fastafile with the genome of the acceptor strain')
	parser.add_argument('-nCPU', dest='nCPU', type=str, default = '10', help='number of CPUs to use during mapping')
	parser.add_argument('-I', dest='minInsertSize', type=str, default = '200', help='lower bound of insert size of read library (-I in bowtie2)')
	parser.add_argument('-X', dest='maxInsertSize', type=str, default = '700', help='upper bound of insert size of read library (-X in bowtie2)')
	parser.add_argument('-mp', dest='isMatePair', action = 'store_true', help='if this flag is set, it is assumed that data is mate pair, hence bowtie2 will run with -rf (default it runs with -fr)')
	parser.add_argument('-path2fastq-mcf', dest='path2fastqmcf', type=str, default = '/Applications/ExpressionAnalysis-ea-utils-27a4809/clipper/', help='path to fastq-mcf')
	args = parser.parse_args()

	print '\n\n##################################'
	print '#'
	print '#   SETTINGS'
	print '#'
	if args.verbose:
		argsdict = vars(args)
		for var in argsdict.keys():
			print '# ', var, argsdict[var]
	print '#'
	print '##################################\n\n'
	return args


def trim(trimmer, adapterfasta, dirname_read_data, outdirname, verbose):

	#trimmer = '/Applications/ExpressionAnalysis-ea-utils-27a4809/clipper/fastq-mcf -q 20 '
	cmnd = ''
	failInTotal = 0
	cmnds       = []
	for fastq_R1 in glob.glob(dirname_read_data+'/*_R1_*.fastq') + glob.glob(dirname_read_data+'/*_R1_*.fq'):

		fastq_R2   = fastq_R1.replace('_R1_', '_R2_')
		trimmed_R1 = outdirname+fastq_R1.split('/')[-1].replace('.fastq', '.trimmed_fastq-mcf_q20.fastq').replace('.fq', '.trimmed_fastq-mcf_q20.fq')
		trimmed_R2 = outdirname+fastq_R2.split('/')[-1].replace('.fastq', '.trimmed_fastq-mcf_q20.fastq').replace('.fq', '.trimmed_fastq-mcf_q20.fq')

		if verbose:
			print fastq_R1
			print fastq_R2
			print 'paired data'
			print 'trimmed reads in', outdirname+':'
			print trimmed_R1
			print trimmed_R2

		if not os.path.exists(trimmed_R1) or not os.path.exists(trimmed_R2):
			cmnd = trimmer+'fastq-mcf -q 20 -o '+trimmed_R1+' -o'+trimmed_R2+' '+adapterfasta+' '+fastq_R1+' '+fastq_R2 +' >& '+trimmed_R1.replace('_R1_', '_')+'.TRIM.LOG'
			cmnds.append(cmnd)
			print cmnd
			fail = os.system(cmnd)
			print fail
			if fail != 0:
				failInTotal += 1
		else:
			print
			print trimmed_R1, 'and', trimmed_R2, 'exist, skipping...'

	print len(cmnds), 'commands, of which', failInTotal,'failed'
	return failInTotal


def mapReads(fasta, dirname_read_data, outdirname, I, X, matePair, nCPU, verbose):

	base_cmnd = 'bowtie2 -x '+fasta+' --end-to-end'
	if matePair:
		base_cmnd += ' --rf'
	else: base_cmnd += ' --fr'

	base_cmnd += ' -p '+nCPU
	base_cmnd += ' --reorder -I '+I+' -X '+X

	failInTotal = 0
	cmnds       = []
	for fastq_R1 in glob.glob(dirname_read_data+'/*_R1_*.fastq') + glob.glob(dirname_read_data+'/*_R1_*.fq'):
		fastq_R2 = fastq_R1.replace('_R1_', '_R2_')

		if verbose:
			print "Will map reads from", fastq_R1, 'and', fastq_R2, 'to', fasta

		outfile_base = outdirname+fastq_R1.split('/')[-1].replace('_R1_', '_')+'__mapped2__'+fasta.split('/')[-1].split('.fa')[0]+'.bowtie2_I'+I+'X'+X
		if matePair: outfile_base+='-rf'

		if not os.path.exists(outfile_base+'.bam'):
			cmnd = base_cmnd +' -1 '+fastq_R1 +' -2 '+fastq_R2+' -S '+outfile_base+'.sam >& '+outfile_base+'.MAPPING.LOG'
			cmnds.append(cmnd)
			print cmnd
			fail = os.system(cmnd)
			print fail
			if fail == 0:
				cmnd = 'samtools view -S -b '+outfile_base+'.sam > '+outfile_base+'.bam'
				if os.system(cmnd) == 0:
					os.system('rm '+outfile_base+'.sam')
				else:
					failInTotal += 1
			else:
				failInTotal += 1
			

		else:
			print
			print outfile_base+'.bam exists, skipping...'

	print len(cmnds), 'commands, of which', failInTotal,'failed'
	return failInTotal


def extractUnmappedReads(dirname_read_data, outdirname, verbose):
	failInTotal = 0
	for bam_fname in glob.glob(dirname_read_data+'/*.bam'):
		fail   = 0
		unmapped_bams = []
		unmapped_bam = outdirname+bam_fname.split('/')[-1].replace('.bam', '.unmapped_matemapped.bam')
		cmnd   = 'samtools view -bf 4 -F264 '+bam_fname+' > '+ unmapped_bam

		if verbose:
			print '\nExtract unmapped reads for '+bam_fname+':'
			print cmnd
		fail = os.system(cmnd)
		if fail != 0:
			print 'Command failed, skipping....'
		else:
			unmapped_bams.append(unmapped_bam)
			unmapped_bam = outdirname+bam_fname.split('/')[-1].replace('.bam', '.mapped_mateunmapped.bam') 
			cmnd = 'samtools view -bf 8 -F260 '+bam_fname+' > '+unmapped_bam
			if verbose:
				print cmnd
			fail = os.system(cmnd)
			if fail != 0:
				print 'Command failed, skipping....'
			else:
				unmapped_bams.append(unmapped_bam)
				unmapped_bam = outdirname+bam_fname.split('/')[-1].replace('.bam', '.both_mates_unmapped')
				cmnd = 'samtools view -bf 12 -F256 '+bam_fname+' > '+unmapped_bam
				if verbose: 
					print cmnd
				fail = os.system(cmnd)
				if fail != 0:
					print 'Command failed, skipping....'
				else:
					unmapped_bams.append(unmapped_bam)

		if fail != 0:
			failInTotal += 1
		else:
			print '\n\nSort 3 classes of unmapped reads and merge into a single file...'
			for unmapped_bam in unmapped_bams:
				cmnd = 'samtools sort -n '+unmapped_bam+' '+unmapped_bam.replace('.bam', '.sorted')
				fail += os.system(cmnd)
			if fail != 0:
				failInTotal += 1
				print cmnd, ' Failed'
			else:
				merged_unmapped_bam_fname = outdirname+bam_fname.split('/')[-1].replace('.bam', '.UNMAPPED.bam')
				cmnd = 'samtools merge -f -n '+ merged_unmapped_bam_fname
				for unmapped_bam in unmapped_bams:
					cmnd += ' '+unmapped_bam.replace('.bam', '.sorted.bam')
				cmnd += ' >& '+merged_unmapped_bam_fname.replace('.UNMAPPED.bam', '.MERGE.LOG')

				if verbose: 
					print cmnd
				fail = os.system(cmnd)
				if fail != 0:
					failInTotal += 1
					print cmnd, ' Failed'
				else:
					cmnd = "bedtools bamtofastq -i "+merged_unmapped_bam_fname
					cmnd +=' -fq '+merged_unmapped_bam_fname.replace('.bam', '_R1_.fastq')
					cmnd +=' -fq2 '+merged_unmapped_bam_fname.replace('.bam', '_R2_.fastq')
					cmnd +=' >& '+merged_unmapped_bam_fname.replace('.UNMAPPED.bam', '.BAM2FASTQ.LOG')
					if verbose: 
						print cmnd
					fail = os.system(cmnd)
					if fail != 0:
						failInTotal += 1
						print cmnd, ' Failed'

	print 'Processed ', len(glob.glob(dirname_read_data+'/*.bam')), 'files, of which', failInTotal, 'failed'
	return failInTotal





if __name__ == "__main__":

	args   = init()
	readDir = args.readDir 
	can_savely_continue = True
	trimOutDir = args.outDir+'/01.Trimming/'

	if args.task == 'trim' or args.task == 'all':
		if args.verbose:
			print 'Will trim reads in', readDir
		
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+trimOutDir)

		if trim(args.path2fastqmcf, args.adapterFasta, readDir, trimOutDir, verbose = args.verbose) == 0:
			readDir = trimOutDir
		else: can_savely_continue = False

	if can_savely_continue and (args.task == 'map2acceptor' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill map reads in', readDir, 'to', args.acceptorFasta

		mapOutDir = args.outDir+'/02.MappingToAcceptor/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+mapOutDir)
		if len(glob.glob(args.acceptorFasta+'*.bt2')) == 0:
			os.system('bowtie2-build '+args.acceptorFasta+' '+args.acceptorFasta)

		if mapReads(args.acceptorFasta, readDir, mapOutDir, args.minInsertSize, args.maxInsertSize, args.isMatePair, args.nCPU, args.verbose) == 0:
			readDir = mapOutDir
		else: can_savely_continue = False

	if can_savely_continue and (args.task == 'mapall2donor' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill map reads in', trimOutDir, 'to', args.donorFasta

		if len(glob.glob(args.donorFasta+'*.bt2')) == 0:
			cmnd = 'bowtie2-build '+args.donorFasta+' '+args.donorFasta
			print cmnd, os.system(cmnd)

		mapOutDir = args.outDir+'/02a.MappingAllReadsToDonor/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+mapOutDir)

		print mapReads(args.donorFasta, trimOutDir, mapOutDir, args.minInsertSize, args.maxInsertSize, args.isMatePair, args.nCPU, args.verbose)

	if can_savely_continue and (args.task == 'extractUnmappedReads' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill extract unmappedreads from bamfiles in ', readDir

		unmappedOutDir = args.outDir+'/03.UnmappedReads/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+unmappedOutDir)

		if extractUnmappedReads(readDir, unmappedOutDir, args.verbose) == 0:
			readDir = unmappedOutDir
		else: can_savely_continue = False

	if can_savely_continue and (args.task == 'map2donor' or args.task == 'all'):
		if args.verbose:
			print '\n\nWill map reads in', readDir, 'to', args.acceptorFasta
		mapOutDir = args.outDir+'/04.MappingToDonor/'
		os.system('mkdir -p '+args.outDir)
		os.system('mkdir -p '+mapOutDir)

		if len(glob.glob(args.donorFasta+'*.bt2')) == 0:
			os.system('bowtie2-build '+args.donorFasta+' '+args.donorFasta)

		print mapReads(args.donorFasta, readDir, mapOutDir, args.minInsertSize, args.maxInsertSize, args.isMatePair, args.nCPU, args.verbose)
			