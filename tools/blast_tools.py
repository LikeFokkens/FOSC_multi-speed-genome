import cluster, glob, re

debug = False #True # 
debug2 = False
debug3 = False #True
#print 'importing scripts/tools/blast_tools.py'

# format blast output
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen


cs = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
def to_complement(seq):
	cseq = ''
	for p in seq:
		try: cseq+=cs[p]
		except: cseq += p
			
	return cseq
	
	
#returns the absolute overlap of two regions
#can only calculate overlap for +strand regions
def overlap(start1, end1, start2, end2):
	if end1<=start1 or end2<=start2: 
		print("incorrect format for regions, region of size 0 or negative size")
		return -1
		
	#complete overlap
	if start2<=start1 and end2>=end1: return end1-start1 #region1 falls completely in region2  
	if start1<=start2 and end1>=end2: return (end2-start2) #region2 falls completely in region1 
	
	#partial overlap:
	#region 2 starts somewhere in region 1
	if start2>=start1 and start2<=end1: return end1-start2
	#region 2 ends somewhere in region 1
	if end2>=start1 and end2<=end1: return end2-start1
	
	return -1
	
#returns the absolute overlap of two regions
# def overlap2((start1, end1), (start2, end2)):
	# if end1<=start1 and end2<=start2: 
		# end    = 0+start1
		# start1 = 0+end1
		# end1   = 0+end
		# 
		# end    = 0+start2
		# start2 = 0+end2
		# end2   = 0+end
		# 
	# elif end1<=start1 or end2<=start2: 
		# return -1
		# 
	# #complete overlap
	# if start2<=start1 and end2>=end1: return end1-start1 #region1 falls completely in region2  
	# if start1<=start2 and end1>=end2: return (end2-start2) #region2 falls completely in region1 
	# 
	# #partial overlap:
	# #region 2 starts somewhere in region 1
	# if start2>=start1 and start2<=end1: return end1-start2
	# #region 2 ends somewhere in region 1
	# if end2>=start1 and end2<=end1: return end2-start1
	# 
	# return -2
	
def test_overlap():
	errors=""
	if overlap((90,10), (10,90)) != -1: errors+="overlap((90,10), (10,90))\n"
	if overlap((10,90), (10,90)) != 80: errors+="overlap((10,90), (10,90))\n"
	if overlap((6,100),(200,300)) != -1: errors+="overlap((6,100),(200,300))\n"
	if overlap((6,100),(20,300)) != 80: errors+="overlap((6,100),(20,300))\n"
	if overlap((20,300),(6,100)) != 80: errors+="overlap((20,300),(6,100))\n"
	if overlap((0,300),(6,100)) != 94: errors+="overlap((0,300),(6,100))\n"
	if overlap((0,100),(100,200)) != 0: errors+="overlap((0,100),(100,200))\n"
	
	if len(errors)>0:
		print("ERRORS:")
		print(errors)
	else: print("fine")
	
	
# merge overlapping regions
def merge_overlapping_regions(list_of_start_end_tuples, min_overlap = 0, minusstrand = False):
	if len(list_of_start_end_tuples)==1: return list_of_start_end_tuples
	
	if minusstrand:
		tmp = [] + list_of_start_end_tuples
		list_of_start_end_tuples = []
		for (s,e) in tmp:
			list_of_start_end_tuples.append((e,s))
			
	region_overlap={}
	Noverlap=0
	for x in range(len(list_of_start_end_tuples)-1):
		if x not in region_overlap.keys(): region_overlap[x]=[]
		for y in range(x+1, len(list_of_start_end_tuples)):
			if overlap(list_of_start_end_tuples[x][0],list_of_start_end_tuples[x][1], list_of_start_end_tuples[y][0], list_of_start_end_tuples[y][1])>=min_overlap:
				Noverlap=1
				region_overlap[x].append(y)
				
				if y in region_overlap.keys():region_overlap[y].append(x)
				else:  region_overlap[y]=[x]
				
	if len(list_of_start_end_tuples)-1 not in region_overlap.keys(): region_overlap[len(list_of_start_end_tuples)-1]=[]  #bugfix 05-07-2013
	
	#print region_overlap
	if Noverlap==0: return list_of_start_end_tuples
	
	groups_of_overlapping_regions=cluster.single_linkage(region_overlap)
	
	regions=[]
	for g in groups_of_overlapping_regions:
		#get min of start sites
		g=list(g)
		s=list_of_start_end_tuples[g[0]][0]
		e=list_of_start_end_tuples[g[0]][1]
		for region_index in g[1:]:
			region=list_of_start_end_tuples[region_index]
			if region[0]<s: s=region[0]
			if region[1]>e: e=region[1]
		
		regions.append((s,e))
	
	if minusstrand:
		tmp = []+regions
		regions = []
		for (s,e) in tmp:
			regions.append((e,s))
			
	return regions
	

	
	
#returns the absolute distance or overlap of two regions
def distance(start1, end1, start2, end2):
	if start1 == end1 or start2==end2:
		print('WARNING! At least one of the regions has size 0!', (start1, end1), (start2, end2))
		
	#check whether both regions are in the same direction
	if (end1>start1 and end2<start2) or (end1<start1 and end2>start2): 
		if debug2: print('Can not merge regions that are not in the same direction!', (start1, end1), (start2, end2))
		return('F', abs(end1-start1), abs(end2-start2))
	
	if type(start1) == str or type(end1) == str or type(start2) == str or type(end2) == str:
		print('ERROR! I am comparing strings! Please convert to int or float. It will make so much more sense!')
		return('F')

	#if regions are in the opposite direction, change direction
	if end1<=start1 and end2<=start2: 
		end    = 0+start1
		start1 = 0+end1
		end1   = 0+end
		
		end    = 0+start2
		start2 = 0+end2
		end2   = 0+end
		
		
	#check for OVERLAP: return overlap as negative distance
	#complete overlap
	if start2<=start1 and end2>=end1: return -1*(end1-start1), abs(end1-start1), abs(end2-start2) #region1 falls completely in region2  
	if start1<=start2 and end1>=end2: return -1*(end2-start2), abs(end1-start1), abs(end2-start2) #region2 falls completely in region1 
	
	#partial overlap:
	#region 2 starts somewhere in region 1
	if start2>=start1 and start2<=end1: return -1*(end1-start2), abs(end1-start1), abs(end2-start2)
	#region 2 ends somewhere in region 1
	if end2>=start1 and end2<=end1: return -1*(end2-start1), abs(end1-start1), abs(end2-start2)
	
	#If not overlapping, check their DISTANCE
	if start1 >= end2: return start1-end2, abs(end1-start1), abs(end2-start2) #region1 lies downstream of region2
	elif start2 >= end1: return start2-end1, abs(end1-start1), abs(end2-start2) #region1 lies upstream of region2
	
	print('could not determine distance for', (start1, end1), (start2, end2))
	return('F', abs(end1-start1), abs(end2-start2))
	
	
	
# merge overlapping regions
# if max_distance < 0: minimum overlap, e.g. max_distance -100: regions should at least overlap for 100bp/aa.
def merge_nearby_regions(list_of_start_end_tuples, max_distance = 0, min_overlap = 0):
	if len(list_of_start_end_tuples)==1: return list_of_start_end_tuples
	
	region_overlap = {}
	Noverlap       = 0
	for x in range(len(list_of_start_end_tuples)-1):
		if x not in region_overlap.keys(): region_overlap[x]=[]
		for y in range(x+1, len(list_of_start_end_tuples)):
			dist, sizex, sizey = distance(list_of_start_end_tuples[x],list_of_start_end_tuples[y])
			if debug3: print(list_of_start_end_tuples[x], list_of_start_end_tuples[y], dist)
			if dist != 'F':
				if dist <= max_distance:
					if min_overlap <= 0 or (dist < 0 and abs(dist) > min_overlap*min([sizex, sizey])):
						Noverlap=1
						region_overlap[x].append(y)
						
						if y in region_overlap.keys():region_overlap[y].append(x)
						else:  region_overlap[y]=[x]
			
	if len(list_of_start_end_tuples)-1 not in region_overlap.keys(): 
		region_overlap[len(list_of_start_end_tuples)-1]=[]  #bugfix 05-07-2013
	
	#print region_overlap
	if Noverlap==0: return list_of_start_end_tuples
	
	groups_of_overlapping_regions=cluster.single_linkage(region_overlap)
	#print region_overlap
	regions=[]
	for g in groups_of_overlapping_regions:
		#get min of start sites
		g=list(g)
		s=list_of_start_end_tuples[g[0]][0]
		e=list_of_start_end_tuples[g[0]][1]
		for region_index in g[1:]:
			region=list_of_start_end_tuples[region_index]
			if region[0]<s: s=region[0]
			if region[1]>e: e=region[1]
		
		regions.append((s,e))
	
	
	return regions
	
	
#merge coupled regions
#the two lists are coupled, list_of_start_end_tuples1[x] is coupled to list_of_start_end_tuples2[x]
#e.g. query and subject alignment regions in BLAST output

#if you want to merge overlapping regions only (that overlap both in the query as well as in the subject":
#merge_nearby_coupled_regions(list_of_start_end_tuples1, list_of_start_end_tuples2, max_overlap1 = -1, max_distance1 = 0, max_overlap2 = -1, max_distance2 = 0)

#if you want to merge regions that do not overlap and are no more than 100 bp apart:
#merge_nearby_coupled_regions(list_of_start_end_tuples1, list_of_start_end_tuples2, max_overlap1 = 0, max_distance1 = 100, max_overlap2 = 0, max_distance2 = 100)

#if you want to merge regions that are about equal distance (the one not more than twice as much as the other), that is maximum 1kb,  apart:
#merge_nearby_coupled_regions(list_of_start_end_tuples1, list_of_start_end_tuples2, max_overlap1 = -1, max_distance1 = 1000, max_overlap2 = -1, max_distance2 = 1000, distance_difference = 0.5)

#if you want to merge regions that are about equal distance (not more than 200bp difference), that is maximum 1kb,  apart:
#merge_nearby_coupled_regions(list_of_start_end_tuples1, list_of_start_end_tuples2, max_overlap1 = -1, max_distance1 = 1000, max_overlap2 = -1, max_distance2 = 1000, distance_difference = 200)

def merge_nearby_coupled_regions(list_of_start_end_tuples1, list_of_start_end_tuples2, max_overlap1 = -1, max_distance1 = 0, max_overlap2 = -1, max_distance2 = 0, distance_difference = None):
	
	if len(list_of_start_end_tuples1)==1: return list_of_start_end_tuples1, list_of_start_end_tuples2, {}
	
	#make a dictionary regionIndex --> indices of regions that overlap
	region2linkedregions     = {}
	NregionsWithLinkedRegion = 0
	
	for x in range(len(list_of_start_end_tuples1)-1):
		if x not in region2linkedregions.keys(): region2linkedregions[x]=[]
		for y in range(x+1, len(list_of_start_end_tuples1)):
			
			dist1, sizex1, sizey1 = distance(list_of_start_end_tuples1[x],list_of_start_end_tuples1[y]) #e.g. distance in query
			dist2, sizex2, sizey2 = distance(list_of_start_end_tuples2[x],list_of_start_end_tuples2[y]) #e.g. distance in subject
			
			if debug:
				print(x, y)
				print(list_of_start_end_tuples1[x],list_of_start_end_tuples1[y], dist1)
				print(list_of_start_end_tuples2[x],list_of_start_end_tuples2[y], dist2)
				
			if dist1 != 'F' and dist2 != 'F':
				#check whether the distance between the regions is not to big:
				if dist1 <= max_distance1 and dist2 <= max_distance2:
					if debug:	print('both distances within maximum distance')
					if distance_difference != None:
						if distance_difference > 0 and distance_difference < 1: #relative maximum difference in distance 
							if dist1 > dist2 and dist2/float(dist1) > distance_difference:	
								region2linkedregions[x].append(y)
								
							elif dist2 > dist1 and dist1/float(dist2) > distance_difference:	
								region2linkedregions[x].append(y)
						elif abs(dist1 - dist2) <= distance_difference:	
							region2linkedregions[x].append(y)
							
					else:
						if debug:	print('checking maximum overlap')
						#first check whether they do not overlap to much
						#check if regions1 fullfill requirements to be linked
						if max_overlap1 < 0 or (max_overlap1 > 0 and dist1 < 0 and abs(dist1) <= max_overlap1) or dist1 > 0:
							#check if regions2 also fullfill requirements to be linked
							if max_overlap2 < 0 or (max_overlap2 > 0 and dist2 < 0 and abs(dist2) <= max_overlap2) or dist2 > 0:
								
								if debug: print('linking regions', list_of_start_end_tuples1[x],list_of_start_end_tuples1[y],\
								 'and', list_of_start_end_tuples2[x],list_of_start_end_tuples2[y])
								NregionsWithLinkedRegion=1
								region2linkedregions[x].append(y)
								
								
							
	if len(list_of_start_end_tuples1)-1 not in region2linkedregions.keys(): region2linkedregions[len(list_of_start_end_tuples1)-1]=[]  
	
	#if no regions can be linked, return original lists
	if NregionsWithLinkedRegion==0: return list_of_start_end_tuples1, list_of_start_end_tuples2, region2linkedregions
	
	
	#if there are regions that can be linked, cluster linked regions
	groups_of_linked_regions = cluster.single_linkage2(region2linkedregions.copy())
	
	regions1=[]
	regions2=[]
	for g in groups_of_linked_regions:
		
		#merge linked regions 
		#REGIONS1
		#check direction
		g=list(g)
		s=list_of_start_end_tuples1[g[0]][0]
		e=list_of_start_end_tuples1[g[0]][1]
		while s == e and len(g) > 0:
			g=g[1:]
			s=list_of_start_end_tuples1[g[0]][0]
			e=list_of_start_end_tuples1[g[0]][1]
			
		if s<e: #+strand
			#get min of start sites and max of end sites for set of linked regions
			for region_index in g[1:]:
				region=list_of_start_end_tuples1[region_index]
				if region[0]<s: s=region[0]
				if region[1]>e: e=region[1]
		elif s>e: #-strand
			#get min of end sites and max of start sites for set of linked regions
			for region_index in g[1:]:
				region=list_of_start_end_tuples1[region_index]
				if region[0]>s: s=region[0]
				if region[1]<e: e=region[1]
			
		
		regions1.append((s,e))
		
		#REGIONS2
		#check direction
		g=list(g)
		s=list_of_start_end_tuples2[g[0]][0]
		e=list_of_start_end_tuples2[g[0]][1]
		while s == e and len(g) > 0:
			g=g[1:]
			s=list_of_start_end_tuples1[g[0]][0]
			e=list_of_start_end_tuples1[g[0]][1]
			
		if s<e: #+strand
			#get min of start sites and max of end sites for set of linked regions
			for region_index in g[1:]:
				region=list_of_start_end_tuples2[region_index]
				if region[0]<s: s=region[0]
				if region[1]>e: e=region[1]
		elif s>e: #-strand
			#get min of end sites and max of start sites for set of linked regions
			for region_index in g[1:]:
				region=list_of_start_end_tuples2[region_index]
				if region[0]>s: s=region[0]
				if region[1]<e: e=region[1]
		
		regions2.append((s,e))
		
	return regions1, regions2, region2linkedregions
	
	
	
def test_merge_nearby_coupled_regions():
	
	#list_of_start_end_tuples1
	#list_of_start_end_tuples2
	#max_overlap1 = -1
	#max_distance1 = 0
	#max_overlap2 = -1
	#max_distance2 = 0
	#distance_difference = None
	
	regions1 = [(10,50),(0,40),(60,100),(100,200),(300,350),(150,350)]
	regions2 = [(60,10),(80,50), (10000, 10050), (12000, 120100), (6050, 6000), (6055, 6260)]
	print(merge_nearby_coupled_regions(regions1, regions2))
	print(merge_nearby_coupled_regions(regions1, regions2, max_distance1 = 100, max_distance2 = 2000))
	print(merge_nearby_coupled_regions(regions1, regions2, max_overlap1 = 10, max_distance1 = 100, max_overlap2 = 20, max_distance2 = 2000))
	
	regions1 = [(10,50),(0,40),(60,100),(100,200),(300,350),(150,350),(150,350)]
	regions2 = [(60,10),(80,50),(100, 150), (12000, 120100), (6050, 6000), (6260, 6055), (6360, 6155)]
	print(merge_nearby_coupled_regions(regions1, regions2))
	print(merge_nearby_coupled_regions(regions1, regions2, max_distance1 = 100, max_distance2 = 2000))
	print(merge_nearby_coupled_regions(regions1, regions2, max_overlap1 = 200, max_distance1 = 100, max_overlap2 = 200, max_distance2 = 2000))
	
	
	
	
	
def sort_regions(regions):
	start_end=dict(regions)
	
	starts=start_end.keys()
	starts.sort()
	
	sorted_regions=[]
	for s in starts:
		sorted_regions.append((s,start_end[s]))
	
	return sorted_regions
	


def get_total_length_regions(regions):
	length = 0
	for (s,e) in regions:
		if e>s: length += (e-s)
		else:	length += (s-e)
		
	return length


def get_total_length_and_longest_region(regions):

	total_length = 0
	max_length   = 0
	longest_region = None
	for (s,e) in regions:
		length = 0
		if e>s: length=(e-s)
		else:	length=(s-e)
		
		total_length += length
		if length > max_length:
			max_length = length
			longest_region = (s,e)

	return total_length, longest_region
	


def test_merge_overlapping_regions():
	regions=[(10,50),(0,40),(60,100),(100,200),(300,350),(150,350)]
	print(merge_overlapping_regions(regions))
	
	regions=[(435, 935), (1, 221)]
	print(merge_overlapping_regions(regions))
	
	regions=[(1, 1307), (1408, 2064), (1408, 2064)] 
	print(merge_overlapping_regions(regions))
	
	regions=[(917, 1921), (1, 626)] 
	print(merge_overlapping_regions(regions))
	
	regions=[(100,200), (150, 250), (917, 1921), (1, 626)] 
	print(merge_overlapping_regions(regions))
	
	regions=[(2100,2200), (2150, 2250), (917, 1921), (1, 626)] 
	print(merge_overlapping_regions(regions))
	



def merge_blast_output(blastfiles, cutoff=None, filename2species=None):
	
	query_target_hsps = {}
	for bf in blastfiles:
		infile = open(bf)
		line = infile.readline()
		while len(line)>0:
			tabs = line.split('\t')
			
			if cutoff==None or float(tabs[10]) < cutoff:
				q = tabs[0]
				t = tabs[1]
				if filename2species != None: t = filename2species[bf]+'_'+t
				if q in query_target_hsps.keys(): 	
					
					if t in query_target_hsps[q].keys(): query_target_hsps[q][t].append(tabs)
					else: query_target_hsps[q][t] = [tabs]
				else: 
					query_target_hsps[q] = {}
					query_target_hsps[q][t] = [tabs]
				
			line = infile.readline()
	return query_target_hsps
	
	
	
def merge_blast_output_in_new_file_unsorted(blastfiles, outfilename, cutoff=None):
	
	query_target_hsps = {}
	for bf in blastfiles:
		infile = open(bf)
		line = infile.readline()
		while len(line)>0:
			tabs = line.split('\t')
			
			if cutoff==None or float(tabs[10]) < cutoff:
				q = tabs[0]
				t = tabs[1]
				
				tabs.append(bf.split('/')[-1])
				if q in query_target_hsps.keys(): 	
					
					if t in query_target_hsps[q].keys(): query_target_hsps[q][t].append(tabs)
					else: query_target_hsps[q][t] = [tabs]
				else: 
					query_target_hsps[q] = {}
					query_target_hsps[q][t] = [tabs]
				
			line = infile.readline()
	
	outfile=open(outfilename, 'w')
	for q in query_target_hsps.keys():
		print(q)
		for t in query_target_hsps[q]:
			for hsp in query_target_hsps[q][t]:
				out=''
				for d in hsp:
					out+=d.strip()+'\t'
				outfile.write(out.strip()+'\n')
			


#parse blastoutput, put everything in a dictionary so all data is gathered per query and then per target
def blast2hash(blastfile, max_Evalue = 10, min_perc_identity = 0, minoverlap_query = 0, minoverlap_subject = 0, minlength_HSP = 0):
	
	# qseqid:0 sseqid:1 pident:2 length:3 qstart:4 qend:5 sstart:6 send:7 evalue:8 bitscore:9 qlen:10 slen:11 qseq:12 sseq:13
	query2subject2data = {}
	line = blastfile.readline()
	tabs = line.strip().split('\t')
	while len(tabs) > 1:
		#print len(tabs)
		Evalue    = float(tabs[8])
		percIdent = float(tabs[2])
		if Evalue <= max_Evalue and percIdent >= min_perc_identity:
			#print tabs
			qstart  = int(tabs[4])
			qend    = int(tabs[5])
			qlenHSP = qend - qstart
			
			sstart  = int(tabs[6])
			send    = int(tabs[7])
			slenHSP = abs(send - sstart)
			
			if qlenHSP >= minlength_HSP and slenHSP >= minlength_HSP and qlenHSP >= minoverlap_query*float(tabs[10]) and slenHSP >= minoverlap_subject*float(tabs[11]):
				q = tabs[0].strip()
				s = tabs[1].strip()
				
				if q not in query2subject2data.keys(): query2subject2data[q] = {}
				if s not in query2subject2data[q].keys(): query2subject2data[q][s] = []
				
				query2subject2data[q][s].append((qstart, qend, sstart, send, percIdent, Evalue, float(tabs[9]), int(tabs[10]), int(tabs[11])))
		line = blastfile.readline()
		tabs = line.strip().split('\t')
				
	return query2subject2data

		

#parse blastoutput, put everything in a dictionary so all data is gathered per query and then per target
def blast2besthit(blastfile, max_Evalue = 10, min_perc_identity = 0, minoverlap_query = 0, minoverlap_subject = 0, minlength_HSP = 0, species_naming_function = None):
	
	# qseqid:0 sseqid:1 pident:2 length:3 qstart:4 qend:5 sstart:6 send:7 evalue:8 bitscore:9 qlen:10 slen:11 qseq:12 sseq:13
	query2species2dataBH = {}
	line = blastfile.readline()
	tabs = line.strip().split('\t')
	while len(tabs) > 1:
		#print len(tabs)
		Evalue    = float(tabs[8])
		percIdent = float(tabs[2])
		if Evalue <= max_Evalue and percIdent >= min_perc_identity:
			#print tabs
			qstart  = int(tabs[4])
			qend    = int(tabs[5])
			qlenHSP = qend - qstart
			
			sstart  = int(tabs[6])
			send    = int(tabs[7])
			slenHSP = abs(qend - qstart)
			
			if qlenHSP >= minlength_HSP and slenHSP >= minlength_HSP and qlenHSP >= minoverlap_query*float(tabs[10]) and slenHSP >= minoverlap_subject*float(tabs[11]):
				q = tabs[0].strip()
				s = tabs[1].strip()
				sp = tabs[1].strip().split('_')[0]
				if species_naming_function != None:
					sp = species_naming_function(s)

				if q not in query2species2dataBH.keys():    query2species2dataBH[q] = {}
				if sp not in query2species2dataBH[q].keys(): 
					query2species2dataBH[q][sp] = []
					query2species2dataBH[q][sp].append((s, qstart, qend, sstart, send, percIdent, Evalue, float(tabs[9]), int(tabs[10]), int(tabs[11])))

		line = blastfile.readline()
		tabs = line.strip().split('\t')
				
	return query2species2dataBH


#parse blastoutput, put everything in a dictionary so all data is gathered per query and then per target
def blast2species2hits(blastfile, max_Evalue = 10, min_perc_identity = 0, minoverlap_query = 0, minoverlap_subject = 0, minlength_HSP = 0, species_naming_function = None):
	
	# qseqid:0 sseqid:1 pident:2 length:3 qstart:4 qend:5 sstart:6 send:7 evalue:8 bitscore:9 qlen:10 slen:11 qseq:12 sseq:13
	query2species2hits = {}
	line = blastfile.readline()
	#print line
	tabs = line.strip().split('\t')
	while len(tabs) > 1:
		#print tabs
		#print len(tabs)
		Evalue    = float(tabs[8])
		percIdent = float(tabs[2])
		if Evalue <= max_Evalue and percIdent >= min_perc_identity:
			#print tabs
			qstart  = int(tabs[4])
			qend    = int(tabs[5])
			qlenHSP = qend - qstart
			
			sstart  = int(tabs[6])
			send    = int(tabs[7])
			slenHSP = abs(qend - qstart)
			
			if qlenHSP >= minlength_HSP and slenHSP >= minlength_HSP and qlenHSP >= minoverlap_query*float(tabs[10]) and slenHSP >= minoverlap_subject*float(tabs[11]):
				q = tabs[0].strip()
				s = tabs[1].strip()
				sp = tabs[1].strip().split('_')[0]
				if species_naming_function != None:
					sp = species_naming_function(s)

				if q not in query2species2hits.keys():    query2species2hits[q] = {}
				if sp not in query2species2hits[q].keys(): 
					query2species2hits[q][sp] = []
				
				query2species2hits[q][sp].append((s, qstart, qend, sstart, send, percIdent, Evalue, float(tabs[9]), int(tabs[10]), int(tabs[11])))

		line = blastfile.readline()
		tabs = line.strip().split('\t')
				
	return query2species2hits


# # Inlcude introns into the alignment
# def get_subject2similarityScores_perQuery(blastfile, max_Evalue = 10, min_perc_identity = 0, minoverlap_query = 0, minoverlap_subject = 0, max_gap = 1000):
	
# 	# data = qstart[0], qend[1], sstart[2], send[3], percIdent[4], Evalue[5], qlen[6], slen[7]

# 	query2subject2data = blast2hash(blastfile, max_Evalue = max_Evalue, min_perc_identity = min_perc_identity, minoverlap_query = 0, minoverlap_subject = 0, minlength_HSP = 0)
# 	query2subject2regions_and_similarityScores = {}
# 	for query in query2subject2data.keys():
# 		for subject in query2subject2data[query].keys():
# 			query_regions   = []
# 			subject_regions = []
# 			minEvalue       = query2subject2data[query][subject][0][5]
# 			maxPercIdent    = query2subject2data[query][subject][0][4]
# 			for tabs in query2subject2data[query][subject]:
# 				query_regions.append((data[0], data[1]))
# 				subject_regions.append((data[2], data[3]))
# 				if tabs[5] < minEvalue:    minEvalue = tabs[5]
# 				if tabs[4] > maxPercIdent: maxPercIdent = tabs[4]

# 			nr_query_regions, nr_subject_regions, region2linkedregions = merge_nearby_coupled_regions(list_of_start_end_tuples1, list_of_start_end_tuples2, max_overlap1 = -1, max_distance1 = 100, max_overlap2 = -1, max_distance2 = 2000, distance_difference = None)
# 			length_hit_query   = get_total_length_regions(nr_query_regions)
# 			length_hit_subject, (sstart, send) = get_total_length_and_longest_region(nr_subject_regions)

# 			qlen = float(query2subject2data[query][subject][0][6])
# 			slen = float(query2subject2data[query][subject][0][7])

# 			if length_hit_query/qlen > minoverlap_query and length_hit_subject/slen) > minoverlap_subject:
# 				if not query2subject2similarityScores.has_key(query): query2subject2similarityScores[query] = {}

# 				query2subject2regions_and_similarityScores[query].has_key(subject): query2subject2regions_and_similarityScores[query][subject] = [nr_subject_regions, minEvalue, maxPercIdent, length_hit_query/qlen, length_hit_subject/slen]

# 	return query2subject2regions_and_similarityScores





if __name__ == "__main__":
	
	test_merge_nearby_coupled_regions()