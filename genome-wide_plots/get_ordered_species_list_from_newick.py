import sys

newick_fname = sys.argv[1]
nw = ''
lines = open(newick_fname).readlines()
for line in lines:
	nw += line.strip()

nw = nw.replace('(', ',').replace(')', ',').replace(':', ',').replace(';', ',')

leaflist = []
for e in nw.split(','):
	if len(e) > 0:
		try: 
			fl = float(e)
		except:
			leaflist.append(e)

print(leaflist)

leafstr = ''
for l in leaflist[::-1]:
	leafstr+=l+' '

print(leafstr[:-1])