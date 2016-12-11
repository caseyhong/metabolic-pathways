from igraph import *

#USED AS A TESTER 
# x = Graph(directed=True)
# y = Graph(directed=True)
# x.add_vertices(['a','b','c'])
# y.add_vertices(['a','b','c'])
# x.add_edge('a','b',weight=5)
# x.add_edge('a','c',weight=10)
# y.add_edge('a','b',weight=2)

def difference(x,y): 
	diff_vertex_set = set()
	edges = {} #key is tuple and value is weight 

	for xEdge in x.es: 
		s = x.vs[xEdge.source]['name']
		t = x.vs[xEdge.target]['name']
		weight = x[s,t]
		edges[(s,t)] = weight 
		diff_vertex_set.update([s,t])#adds both to set 

	for yEdge in y.es: 
		s = y.vs[yEdge.source]['name']
		t = y.vs[yEdge.target]['name']
		weight = y[s,t]
		diff_vertex_set.update([s,t])#adds both to set 

		#performs difference 
		if (s,t) in edges.keys(): 
			edges[(s,t)] -= weight 
		else: 
			edges[(s,t)] = weight 

	diffG = Graph(directed=True) 
	diffG.add_vertices(list(diff_vertex_set))
	for e in edges: 
		diffG.add_edge(e[0],e[1],weight=edges[e])

	return diffG 

def createNetwork(metabolites,rxnMap,rxnExpVals): 
	network = Graph(directed=True)
	network.add_vertices(setofallmetabolites)

	for rxn_id in rxnMap.keys(): 
		sources = rxnMap[rxn_id][0]
		targets = rxnMap[rxn_id][1]

		# Need to add the reaction id prefix for each metabolite edge 
		edges = [(rxn_id+"_"+s,rxn_id+"_"+st) for s in sources for t in targets]

		for (source,target) in edges: 
			network.add_edge(source,target,weight=rxnExpVals[rxn_id]) 
	return network

if __name__ == '__main__':
	print 'diff_paths'
	#createNetwork(metabolites,rxnMap,rxnExpVals)


# ## DO NOT REMOVE 
# # >>> x = Graph(directed=True)
# # >>> y = Graph(directed=True)
# # >>> x.add_vertices(['a','b','c'])
# # Each graph represents the network for each patient for each tumor type 
# # N = number of vertices 
# # edges = pairs of integers or names of two endpoints  
# # directed = True 
# # graph_attrs
# # vertex_attrs
# # edge_attrs 
