import sys, os , glob, argparse

# tasks to add: filter, reorder, sizesort

def init():
	parser = argparse.ArgumentParser(description='Suite of tools to manipulate and extract data from fastafiles')
	parser.add_argument("-task", dest='task', choices = ['fa2phy', 'printseqlen', 'fasta2genomefile', 'csv2fasta', 'printseq', 'toComplement', 'concatenateMSA', 'split_fasta', 'bed2fasta'],\
	 help='task: convert fasta 2 sequential phylip [fa2phy], print length of a sequence in the fasta [printseqlen <contig(s)> <start> <end>],\
	  convert a csv file (seqid, sequence) to a fasta (removing non-DNA sequences, removing comments (split sequence at first whitespace)) [csv2fasta], \
	  print (sub)sequence in the fasta [printseq <contig(s)> <start> <end>], create a file with the length of each sequence in the fasta (often required when using samtools) [fasta2genomefile]\
	  create a new file *.RC.fasta containing the reverse complement of each sequence in the fastafiles provided [toComplement], concatenate a multiple sequence alignment [concatenateMSA], \
	  split a large fasta with a lot of sequences in several smaller ones [splitfasta ,maxNseqs> or <Nfiles>, useful when running signalP for example.' )
	parser.add_argument("-fasta", dest='fasta', help='name of the fastafile on which to perform the task')
	parser.add_argument("-inDir", dest='inDir', help='dirname containing the fastafiles (end with .fasta) on which to perform the task')
	parser.add_argument("-outDir", dest='outDir', help='outdirname if output should not be saved in same dir as fasta or inDir')
	parser.add_argument("-contig", nargs='*', dest='contigs', help='(list of) contigname(s) for which you want to retrieve the sequence (length) or that specifies how to reorder the fastafile')
	parser.add_argument("-start", dest='start', type=int, default = None, help='start of the region for which you want to retrieve the sequence (length), default is 0, starting at the start of the contig. WARNING: indexing in sequence is 0-based!')
	parser.add_argument("-end", dest='end', type=int, default = None, help='end of the region for which you want to retrieve the sequence (length), default is -1, will retunr sequence until the end of the contig')
	parser.add_argument("-splitHeaderAt", dest='split_char', default = '\s', help='Symbol where to split the fasta header, default is \s (any whitespace).')
	parser.add_argument("-csvfile", dest='csv', help='name of the csv file that needs to be converted to a fasta')
	parser.add_argument("-bed", dest='bed', help='name of the bed file that needs to be converted to a fasta')
	parser.add_argument("-sep", dest='sep', default = ';', choices = ['tab', ';'], help='symbol that seprates columns in the csv')
	parser.add_argument("-maxNseqs", dest='maxNseqs',type=int, help='if task is splitfasta: maximum number of sequences per fastafile')
	parser.add_argument("-Nfiles", dest='Nfiles', type = int, help='if task is splitfasta: number of files to split this fasta into')
	
	parser.add_argument("-v", "--verbose", default = False, help="print progress messages", action="store_true")
	

	########################
	#
	#  Catch possible errors here
	#
	########################
	args = parser.parse_args()
	if not args.task == 'csv2fasta' and args.inDir == None and args.fasta == None:
		print('Please provide the name of the fastafile (-fasta) or of the directory that contains the fastafile (-inDir) which I need to work on.')
		sys.exit()

	if args.task == 'printseq' and args.contigs == None:
		print("Please provide name of the contig for which I need to print the sequence")
		sys.exit()

	if args.task == 'concatenateMSA' and args.inDir == None:
		print('Please provide an input directory with fastafiles to concatenate')
		sys.exit() 

	if args.task == 'csv2fasta' and args.csv == None:
		print('Please provide a csv file that should be converted to a fasta with -csvfile')
		sys.exit()

	if args.task == 'bed2fasta' and args.bed == None:
		print('Please provide a bed file that should be converted to a fasta with -bed')
		sys.exit()

	if args.task =='split_fasta' and args.maxNseqs == None and args.Nfiles == None:
		print('You need to specify either the maximum number of sequences per file (with -maxNseqs) \
			or the number of files (-Nfiles) you want to split the fasta(s) into')

	return args






cs = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G','a': 't', 'g': 'c', 't': 'a', 'c': 'g'}
def to_complement(seq, verbose = False):
	cseq = ''
	for p in seq:
		try: cseq+=cs[p]
		except: 
			cseq += p
			if verbose: print('Encountered non-DNA character: '+p)

	return cseq	



	
def is_DNA(seq):
	if len(set(seq).difference(set(['A','T','G','C','a','t','g','c', 'Y', 'y', 'X', 'x', 'N', 'n']))) == 0: return True
	else: 
		return False




	

	
def fasta2dicts(fastafile, splitheader = True):
	
	id2seq    = {}
	id2header = {}
	line      = fastafile.readline()
	if len(line) < 2:
		print('*** WARNING! this fastafile is empty!')
		return {},{}

	currentid = line.strip()[1:]
	if splitheader:
		currentid = line.strip().split()[0][1:]
	id2header[currentid] = line.strip()
	
	seq=''
	line = fastafile.readline()
	while len(line)>0:
		if line[0]=='>': #new header
			if currentid in id2seq and id2seq[currentid] != seq: 
				print('WARNING! key', currentid, "already exists, with a different sequence! I will overwrite it.")
			
			if seq[-1] == '*': seq = seq[:-1]
			id2seq[currentid]=seq.upper()
		
			#intialize new
			if splitheader: currentid = line.strip().split()[0][1:]
			else:           currentid = line.strip()[1:]
			id2header[currentid] = line.strip()
			
			seq=''
		else:
			#add this line to sequence
			seq += line.strip()
			
		line = fastafile.readline()
	
	#the last one
	if seq[-1] == '*': seq = seq[:-1]
	id2seq[currentid]=seq.upper()
	
	return id2seq, id2header




def fasta2dict_and_genelist(geneseqfile, keep_whole_header = False):
	
	line        = geneseqfile.readline()
	currentgene = ''
	if keep_whole_header: currentgene = line[1:]
	else:                 currentgene = line.split()[0][1:] #locusname
	seq=''
	gene_seq = {}
	line     = geneseqfile.readline()
	genelist = []
	while len(line)>0:
		if line[0]=='>':
			if currentgene in gene_seq and gene_seq[currentgene] != seq: 
				print('key '+currentgene+" already exists, with a different sequence! I will overwrite it. \
					I suggest you will set keep_whole_header to 'True'")
			
			gene_seq[currentgene]=seq
			genelist.append(currentgene)
			#intialize new
			if keep_whole_header: currentgene = line[1:]
			else:                 currentgene = line.split()[0][1:] #locusname
			seq=''
		else:
			#add this line to sequence
			seq += line.strip()
		line = geneseqfile.readline()
	
	#don't forget the last one :-)
	gene_seq[currentgene]=seq
	genelist.append(currentgene)
	
	return gene_seq, genelist


	
def fasta2genomefile(fasta_fname, genome_fname = None, outdir = None):

	if genome_fname == None:
		genome_fname = fasta_fname.replace('.fasta', '.contig2size.tab')
		if genome_fname == fasta_fname: 
			genome_fname = fasta_fname.replace('.fa', '.contig2size.tab')
		if genome_fname == fasta_fname: 
			genome_fname = fasta_fname.replace('.fna', '.contig2size.tab')
		if genome_fname == fasta_fname: 
			genome_fname = fasta_fname+'.contig2size.tab'
		if outdir != None:
			genome_fname = outdir + genome_fname.split('/')[-1]


	id2seq, id2header = fasta2dicts(open(fasta_fname))
	total = 0
	outfile = open(genome_fname, 'w')
	for id in id2seq.keys():
		seqlen = len(id2seq[id])
		total += seqlen
		outfile.write(id+'\t'+str(seqlen)+'\n')

	outfile.close()
	return genome_fname, total

	

def get_sequence_from_file(fastafile, header, start=None, end=None):
	infile = open(fastafile)
	line   = infile.readline()
	seq = ''
	while len(line)>0:
		if line[0] == '>':
			if line[1:].strip().split()[0] == header:
				line = infile.readline()
				while line[0] != '>':
					seq  += line.strip().upper()
					line = infile.readline()
				
		line   = infile.readline()
	infile.close()
	
	
	if start == None:
		return '>'+header, seq
	elif end ==None:
		newheader = '>'+header+':'+str(start)+'-'+str(len(seq))
		return newheader, seq[start:]
	else:
		if start < end:
			newheader = '>'+header+':'+str(max([0, start]))+'-'+str(min([end+1, len(seq)]))
			return newheader, seq[max([0, start]):min([end+1, len(seq)])]
		else:
			newheader = '>'+header+':'+str(min([start+1, len(seq)]))+'-'+str(max([0, end]))
			return newheader, to_complement(seq[max([0, end]):min([start+1, len(seq)])])[::-1]
	



def cut_region_from_sequence(refstart, refend, seqstart, queryseq, targetseq):
	
	i     = seqstart
	iseq  = 0
	t     = 0
	#'walk' to start of region shared by all HSPs which is where we want to start 'recording' the sequence of the subject:
	while i < max([refstart, seqstart]):
		if queryseq[iseq] != '-' : i+=1
		iseq+=1
	
	outseq = ''
	t = len(targetseq[:iseq].replace('-', ''))
	
	while i < refend:
		outseq+=targetseq[iseq]
		if queryseq[iseq]!='-': i+=1
		
		iseq+=1
		
	outseq = outseq.replace('-', '')
	
	
	return outseq, t
	

	

def concatenateAll(list_of_fastafnames, outfilename):

	empty_fastas = []
	species2concatenatedseq = {}
	species_set = set([])
	for fasta in list_of_fastafnames:
		print(fasta)
		id2seq, id2header = fasta2dicts(open(fasta))
		if id2seq == {}: empty_fastas.append(fasta)
		lengths = set([])

		species_set_f = set(id2seq.keys())
		if len(species_set) > 0: 
			if species_set_f != species_set:
				print('Irregular set of species!')
				extra = species_set_f.difference(species_set)
				if len(extra) > 0: print(extra+' in this set but not other fasta')
				missing = species_set.difference(species_set_f)
				if len(missing) > 0: print(missing+' missing in this set')

		else: species_set = species_set_f
		for species in id2seq.keys():
			seq = id2seq[species]
			if species2concatenatedseq.has_key(species): 
				species2concatenatedseq[species] += seq.strip()
			else: species2concatenatedseq[species]  = seq.strip()
			lengths.add(len(seq.strip()))

		if len(lengths) > 1: 
			print('ERROR! not all species have sequences of the same length '+fasta)
		
	outfile = open(outfilename, 'w')
	lengths = set([])
	for species in species2concatenatedseq.keys():
		seq = species2concatenatedseq[species]
		outfile.write('>'+species+'\n'+seq+'\n')
		lengths.add(len(seq))
		print(species+': '+str(len(seq)))
	
	if len(lengths) > 1: 
		print('ERROR! not all species have sequences of the same length '+lengths)
	
	return empty_fastas




def get_base_id(long_id, split_char):
	return long_id.split(split_char)[0]




def overlap2fasta(file1, file2, outfile1, outfile2, split_char):

	id2seq1, id2header1 = fasta2dicts(open(file1))
	id2seq2, id2header2 = fasta2dicts(open(file1))
	
	ids1 = id2seq1.keys()
	ids2 = id2seq2.keys()
	base_ids1 = map(get_base_id, ids1)
	base_ids2 = map(get_base_id, ids2)

	shared = set(base_ids1).intersection(set(base_ids2))

	for index, id1 in enumerate(base_ids1):
		if id1 in shared:
			outfile.write('>'+id2header1[ids1[index]]+'\n'+id2seq1[ids1[index]]+'\n')

	for index, id2 in enumerate(base_ids2):
		if id2 in shared:
			outfile.write('>'+id2header2[ids2[index]]+'\n'+id2seq2[ids2[index]]+'\n')	
	


	
def fasta2phylip(fasta_fname, phylip_fname):
	
	id2seq, id2header = fasta2dicts(open(fasta_fname))
		
	# first line contains number and length of the sequences
	out=str(len(id2seq.keys()))+' '+str(len(id2seq[id2seq.keys()[0]]))+'\n'
	for sid in id2seq.keys():
		out+=sid+' '+id2seq[sid]+'\n'
	
	phylipfile = open(phylip_fname, 'w')
	phylipfile.write(out)

	phylipfile.close()
		



def csv2fasta(csv_filename, out_fname, separator = '\t'):
	csv_file = open(csv_filename)
	lines = csv_file.read().split('\r')
	
	print(lines)
	outfile = open(out_fname, 'w')
	Nseqs = 0
	for line in lines:
		
		data =line.strip().split(separator)

		if len(data) == 2:
			id  = data [0]
			seq = data [1]
			seq = seq.replace('\s', '')  # remove any whitespaces
			seq = seq.replace("5'-", '') # remove that comment
			seq = seq.replace("3'-", '') # remove that comment
			seq = seq.replace("5'", '')  # remove that comment
			seq = seq.replace("3'", '')  # remove that comment

			if is_DNA(seq) and len(seq) > 0:
				outfile.write('>'+id+'\n'+seq+'\n')
				Nseqs += 1
	print('wrote '+str(Nseqs)+' sequences')
	outfile.close()



def split_fasta(fasta, max_n_seqs = None, Nfiles = None):

	id2seq, genelist = fasta2dict_and_genelist(open(fasta), keep_whole_header = True)

	if max_n_seqs == None:
		if Nfiles == None:
			print("I can't split a fasta if I don't know into how many pieces....")
			sys.exit()
			
		else:
			max_n_seqs = round(id2seq/float(Nfiles), 0)

	start = 0
	outfile_base = fasta.split('.fa')[0]
	outfile      = None
	for i, header in enumerate(genelist):
		if i%max_n_seqs == 0:
			if i > 0: outfile.close()
			start       = i
			outfilename = outfile_base+'_'+str(start)+'-'+str(min([len(id2seq.keys()), start+max_n_seqs]))+'.fasta'
			outfile     = open(outfilename, 'w')
			

		if i==0: print(outfilename+' '+header)
		outfile.write('>'+header.strip()+'\n'+id2seq[header]+'\n')




def bed2fasta(bedfname, fastafname):
	bedfile   = open(bedfname)
	fastafile = open(fastafname)

	contig2seq, contig2header = fasta2dicts(fastafile)
	
	outfile   = open(bedfname.replace('.bed', '.fasta'), 'w') 
	for line in bedfile.readlines():
		contig, start, end = line.strip().split()[:3]

		seq   = contig2seq[contig]
		start = int(start)
		end   = int(end)
		if start < end:
			newheader = '>'+contig+':'+str(max([0, start]))+'-'+str(min([end+1, len(seq)]))
			outfile.write(newheader+'\n'+seq[max([0, start]):min([end+1, len(seq)])]+'\n')
		else:
			newheader = '>'+contig+':'+str(min([start+1, len(seq)]))+'-'+str(max([0, end]))
			outfile.write(newheader+'\n'+to_complement(seq[max([0, end]):min([start+1, len(seq)])])[::-1]+'\n')

	outfile.close()



##############################
#
# TEST
#
#############################


	
def test_cut_region_from_sequence():
	queryseq  = 'HEE-BRAV-E---HOND'
	targetseq = 'HEELBR--RENGEHOND'

	refstart = 7
	refend = 13
	seqstart = 2
	seq, t = cut_region_from_sequence(refstart, refend, seqstart, queryseq, targetseq)
	print(seq+' '+str(t))
	
	refstart = 1
	refend = 13
	seqstart = 1
	seq, t = cut_region_from_sequence(refstart, refend, seqstart, queryseq, targetseq)
	print(seq+' '+str(t))

#############################












if __name__ == "__main__":

	args = init()

	if args.task == 'csv2fasta':
		out_fname = args.csv.replace('.csv', '.fasta')
		if out_fname == args.csv: out_fname = args.csv+'.fasta'
		if args.outDir != None:
			out_fname = outDir + out_fname.split('/')[-1]
		if args.sep == 'tab':
			csv2fasta(args.csv, out_fname)
		else:
			csv2fasta(args.csv, out_fname, separator = args.sep)
	
	else:
		fastafiles = []
		if args.fasta == None:
			fastafiles = glob.glob(args.inDir+'*.fasta')
		else:
			fastafiles = [args.fasta]

		if args.task == 'fa2phy':
			for fastafile in fastafiles:
				phylipfile = fastafile.replace('.fasta', '.phylips')
				if phylipfile == fastafile:
					phylipfile = fastafile.replace('.fa', '.phylips')
				if phylipfile == fastafile:
					phylipfile = fastafile+'.phylips'

			if args.outDir != None:
				pylipfile = outDir + phylipfile.split('/')[-1]

			if args.verbose:
				print('Converting '+fastafile+' to sequential phylip format: '+phylipfile)

			fasta2phylip(fastafile, phylipfile)


		elif args.task == 'printseqlen':
			for fastafile in fastafiles:
				print(fastafile)
				id2seq, id2header = fasta2dicts(open(fastafile))

				contigs = id2seq.keys()
				if args.contigs != None: contigs = args.contigs
				
				total = 0
				for cid in contigs:
					seqlen = len(id2seq[cid])
					print(cid+' '+str(seqlen)+' '+str(round(seqlen/1000000.0, 2))+' Mb')

					total += seqlen
				print('\nTOTAL\n'+str(total)+' '+str(round(total/1000000.0, 2))+' Mb')
				print('\n\n\n')
		

		elif args.task == 'printseq':
			for fastafile in fastafiles:
				header, seq = get_sequence_from_file(fastafile, args.contigs[0], start=args.start, end=args.end)
				if len(seq) == 0:
					print('This header is not found in the fasta, check again.')
					sys.exit()
				print(header)
				print(seq)

				
		elif args.task == 'fasta2genomefile':
			for fastafile in fastafiles:
				outfilename = None
				#if args.outDir !=None:
				#	outfilename = args.outDir + fastafile.split('/')[-1].split('.fa')[0]+'contig2size.tab'

				outfilename = fasta2genomefile(fastafile, genome_fname = outfilename, outdir = args.outDir)
				if args.verbose: print(fastafilename+' --> '+outfilename)


		elif args.task == 'toComplement':
			for fastafile in fastafiles:
				outfilename = fastafile.split('.fa')[0]+'.RC.fasta'
				if args.outDir != None:
					outfilename = outDir+fastafile.split('/')[-1].split('.fa')[0]+'.RC.fasta'

				outfile = open(outfilename)
				id2seq, id2header = fasta2dicts(open(fastafile))
				for cid in id2seq.keys():
					outfile.write('>'+cid+' Reverse Complement\n'+to_complement(id2seq[cid])[::-1]+'\n')
				outfile.close()


		elif args.task == 'concatenateMSA':
			outfilename = ''
			if args.outDir == None:
				outfilename = args.outDir + 'concatenated_aligment.fasta'
			else: outfilename = args.inDir + 'concatenated_aligment.fasta'

			concatenateAll(fastafiles, outfilename)


		elif args.task == 'split_fasta':
			for fastafile in fastafiles:
				split_fasta(fastafile, max_n_seqs = args.maxNseqs, Nfiles = args.Nfiles)

	
		elif args.task == 'bed2fasta':
			bed2fasta(args.bed, fastafiles[0])

	

		
		

