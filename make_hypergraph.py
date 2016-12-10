from halp.directed_hypergraph import DirectedHypergraph
import gpr_mapping as GPR
from halp.utilities.directed_graph_transformations import to_networkx_digraph
import networkx as nx
import igraph as ig
import numpy as np

H = DirectedHypergraph()

met_map = GPR.get_metabolite_associations('./RECON1.json')
for reaction in met_map:
	source = met_map[reaction][0]
	target = met_map[reaction][1]
	H.add_hyperedge(set(source), set(target), weight=1)

g = to_networkx_digraph(H)
graph = ig.Graph(len(g), zip(*zip(*nx.to_edgelist(g))[:2]))
# print nx.to_edgelist(g)
# graph = ig.Graph.Adjacency((nx.to_numpy_matrix(g) > 0).tolist())
