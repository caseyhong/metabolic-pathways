# import diff_paths as dp
import gpr_mapping as GPR
# from igraph import *
import numpy as np

met_map = GPR.get_metabolite_associations('./RECON1.json')
src_set = set()
for m in met_map:
	sources = met_map[m][0]
	for i in sources:
		if i not in src_set:
			src_set.add(i)
print len(src_set)

## which function do we call to get the diff graph?
## diff_graph is an iGraph.Graph object
# diff_graph