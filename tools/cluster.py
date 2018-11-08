debug=False

#print ('Importing cluster.py'

def single_linkage(node_partners):

    clusters=[]
    
    nodes=list(node_partners.keys())
    while len(nodes)>0:
        if debug:
            print ('new cluster, seed node:', nodes[0])
            print ('# keys in hash:', len(nodes))
        cluster=update_cluster(node_partners, [nodes[0]], set([]))
        clusters.append(cluster)
        nodes=list(node_partners.keys())
        if debug:
            print ('finished cluster:', cluster)
            print ('# keys in hash:', len(nodes))
        
    return clusters


def update_cluster(node_partners, todo, cluster=set([])):
    if debug: print (todo, cluster)
    new_todo=set([])
    for t in todo:
        partners=node_partners[t]
        new_todo=new_todo.union(set(partners))
        cluster.add(t)
        del node_partners[t]

    new_todo=new_todo.difference(cluster)
    if len(new_todo)>0: return update_cluster(node_partners, new_todo, cluster)
    else: return cluster
    
    
# # for key -> values hashes in which not all values are a key themselves
# def single_linkage2(node_partners):

#     clusters=[]
    
#     nodes=node_partners.keys()
#     while len(nodes)>0:
#         if debug:
#             print ('new cluster, seed node:', nodes[0]
#             print ('# keys in hash:', len(nodes)
#         cluster=update_cluster2(node_partners, [nodes[0]], set([]))
#         clusters.append(cluster)
#         nodes=node_partners.keys()
#         if debug:
#             print ('finished cluster:', cluster
#             print ('# keys in hash:', len(nodes)
        
#     return clusters


# def update_cluster2(node_partners, todo, cluster=set([])):
#     if debug: print (todo, cluster 
#     new_todo=set([])
#     for t in todo:
#     	if node_partners.has_key(t):
#     		partners=node_partners[t]
#     		new_todo=new_todo.union(set(partners))
#         	del node_partners[t]
#         cluster.add(t)

#     new_todo=new_todo.difference(cluster)
#     if len(new_todo)>0: return update_cluster2(node_partners, new_todo, cluster)
#     else: return cluster    
    
    
def test():
    links=[(0,1), (1,2), (2,3), (3,4), (1,4), (2,4), (5,6), (6,8), (6,9), (8,9), (10,11), (10,12), (11,13), (5,7)]

    testhash={}
    for (x,y) in links:
        try:    testhash[x].add(y)
        except: testhash[x]=set([y])

        try:    testhash[y].add(x)
        except: testhash[y]=set([x])

    if single_linkage(testhash) != [set([0, 1, 2, 3, 4]), set([8, 9, 5, 6, 7]), set([10, 11, 12, 13])]:
    	print ("failed to cluster:", links) 
    	
    	
    links=[(0,1), (2, 1), (2,3), (3,4), (4, 1), (2,4), (6,5), (6,8), (6,9), (8,9), (10,11), (10,12), (11,13), (7,5)]

    # testhash={}
    # for (x,y) in links:
    #     try:    testhash[x].add(y)
    #     except: testhash[x]=set([y])

    # print ('hash', testhash
    
    # if single_linkage2(testhash) != [set([0, 1, 2, 3, 4]), set([8, 9, 5, 6, 7]), set([10, 11, 12, 13])]:
    # 	print ("failed to cluster:", links 
    	

if __name__ == "__main__":
    test()
