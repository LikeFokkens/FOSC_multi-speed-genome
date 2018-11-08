import sys, argparse
from ete3 import PhyloTree, TextFace, NodeStyle, TreeStyle

def init():
	parser = argparse.ArgumentParser(description='Plot a species tree using ete3')
	parser.add_argument("-tree", dest='tree_fname', default = None, help='Name of a file containing a phylogenetic tree in newick format')
	parser.add_argument("-datafile", dest='data_fname', help='Name of a tabulated file containing information on how to plot individual nodes: nodename{tab}nodenameInFigure{tab}shape{tab}color')
	parser.add_argument("-include", dest='include_fname', default = None, help='Name of a file containing node ids to include in the tree (others will be removed))')
	parser.add_argument("-exclude", dest='exclude_fname', default = None, help='Name of a file containing node ids to exclude in the tree (others will be kept))')
	
	parser.add_argument("-figfilename", dest='fname_fig', default = None, help='Name of a file in which the figure will be saved')
	parser.add_argument("-treefilename", dest='fname_tree', default = None, help='Name of a file in which the pruned tree will be saved (as newick)')
	parser.add_argument("-minBootstrap", dest='minB',type=float, default = 0.0, help='Minimum bootstrap support for a node to be kept')
	
	parser.add_argument("-show_branch_support", dest='show_support', default = False, help="Show branch support in the tree", action="store_true")
	

	########################
	#
	#  Catch possible errors here
	#
	########################
	args = parser.parse_args()
	if args.data_fname == None:
		print('WARNING: no data file, will plot in default style')

	return args


def prune_tree(tree, minB, include, exclude):

	for node in tree.traverse():
		if not node.is_leaf():
			if node.support < minB: 
				node.delete(prevent_nondicotomic=True, preserve_branch_length=True)

		else:
			if len(exclude) > 0:
				if node.name in exclude:
					node.delete(prevent_nondicotomic=True, preserve_branch_length=True)

			elif len(include) > 0:
				if node.name not in include:
					node.delete(prevent_nondicotomic=True, preserve_branch_length=True)

	return tree


def layout(node):
	
	print(node.name, leafname2data)
	nst_internal = NodeStyle()   
	nst_internal['size'] = 1
	nst_internal['fgcolor'] = 'black'
	nst_internal["vt_line_color"] = "black"
	nst_internal["hz_line_color"] = "black"
	nst_internal["vt_line_width"] = 1
	nst_internal["hz_line_width"] = 1
	nst_internal["vt_line_type"] = 0 #0 solid, 1 dashed, 2 dotted
	nst_internal["hz_line_type"] = 0
	
	node.set_style(nst_internal) 

	if leafname2data != None and node.is_leaf():  
		if not node.name[0] == '_':
			print(node.name)
			nst = NodeStyle()

			node_data  = leafname2data[node.name]
			node_name  = node.name
			node_shape = None
			node_color = None
			if len(node_data) >= 3:
				node_name  = node_data[0]
				node_shape = node_data[1]
				node_color = '#'+node_data[2].lower()
			else:
				node_shape = node_data[0]
				node_color = '#'+node_data[1].lower()

			nst['size']    = 15
			nst['shape']   = node_shape
			nst['fgcolor'] = node_color
			
			node.set_style(nst)
			node.add_face(TextFace(node_name, ftype='Verdana', fsize=12), column = 0, position = "aligned")



def plot_tree(tree, leafname2data = None , show_support = False, figname = None, outfname_tree = None):

	
	def layout(node):
		
		nst_internal = NodeStyle()   
		nst_internal['size'] = 1
		nst_internal['fgcolor'] = 'black'
		nst_internal["vt_line_color"] = "black"
		nst_internal["hz_line_color"] = "black"
		nst_internal["vt_line_width"] = 1
		nst_internal["hz_line_width"] = 1
		nst_internal["vt_line_type"] = 0 #0 solid, 1 dashed, 2 dotted
		nst_internal["hz_line_type"] = 0
		
		node.set_style(nst_internal) 

		if leafname2data != None and node.is_leaf():  
			if not node.name[0] == '_':
				print(node.name)
				nst = NodeStyle()

				node_data  = leafname2data[node.name]
				node_name  = node.name
				node_shape = None
				node_color = None
				if len(node_data) >= 3:
					node_name  = node_data[0]
					node_shape = node_data[1]
					node_color = '#'+node_data[2].lower()
				else:
					node_shape = node_data[0]
					node_color = '#'+node_data[1].lower()

				nst['size']    = 15
				nst['shape']   = node_shape
				nst['fgcolor'] = node_color
				
				node.set_style(nst)
				node.add_face(TextFace(node_name, ftype='Verdana', fsize=12), column = 0, position = "aligned")

	ts = TreeStyle()
	ts.layout_fn = layout
	ts.show_leaf_name = False
	if show_support: 
		ts.show_branch_support = True

	ts.branch_vertical_margin  = 10

	if outfname_tree != None:
		tree.write(outfile = outfname_tree)

	if figname != None:
		tree.render(figname, layout=layout, w=None, h=None, tree_style=ts, units='px', dpi=300)
	else:
		tree.show(tree_style=ts)	



if __name__ == "__main__":

	args = init()

	lines  = open(args.tree_fname).readlines()
	newick = ''
	for line in lines: newick += line.strip()
	tree = PhyloTree(newick, sp_naming_function=None)

	#prune_tree
	include = set([])
	if args.include_fname != None:
		for line in open(args.include_fname).readlines():
			include.add(line.strip())

	exclude = set([])
	if args.exclude_fname != None:
		for line in open(args.exclude_fname).readlines():
			exclude.add(line.strip())
	
	pruned_tree = prune_tree(tree, args.minB, include, exclude)


	# plot tree
	leafname2data = None
	print(args.data_fname)
	if args.data_fname != None:
		leafname2data = {}
		lines = open(args.data_fname).readlines()
		for line in lines:
			tabs = line.strip().split('\t')
			leafname2data[tabs[0]] = tabs[1:]

	print(leafname2data)

	plot_tree(pruned_tree, leafname2data = leafname2data, show_support = args.show_support, figname = args.fname_fig, outfname_tree = args.fname_tree)




