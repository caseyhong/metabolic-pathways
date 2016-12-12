import gpr_mapping as GPR

def path_to_genes():
	path_to_genes = {}
	with open('c2.cp.kegg.v5.2.entrez.gmt.txt', 'r') as file:
		for line in file:
			items = line.split()
			pathway = items[0]
			genes = items[2:]
			path_to_genes[pathway] = genes
	return path_to_genes

def get_src_tar_genes(path_to_genes):
	candidates = []
	for path in path_to_genes:
		genes = path_to_genes[path]
		candidates.append((genes[0], genes[-1]))
	return candidates

def get_reactions_of_interest(candidates, gene_to_reaction):
	candidate_reactions = []
	for candidate in candidates:
		src = candidate[0]
		tar = candidate[1]
		if src in gene_to_reaction and tar in gene_to_reaction:
			s = gene_to_reaction[src]
			t = gene_to_reaction[tar]
			candidate_reactions.append((s,t))
	return candidate_reactions
	# src_rxn = []
	# tar_rxn = []
	# for s in src_genes:
	# 	if s in gene_to_reaction:
	# 		rxn = gene_to_reaction[s]
	# 		for i in rxn:
	# 			if i not in src_rxn:
	# 				src_rxn.append(i)
	# for t in tar_genes:
	# 	if t in gene_to_reaction:
	# 		rxn = gene_to_reaction[t]
	# 		for i in rxn:
	# 			if i not in tar_rxn:
	# 				tar_rxn.append(i)
	# return src_rxn, tar_rxn

def get_candidate_nodes(candidate_reactions, met_associations):
	# print "candidate reactions..."
	# print candidate_reactions
	# print "met_associations..."
	# print met_associations
	candidate_pairs = []
	for c in candidate_reactions:
		src = c[0] ## reactions associated with src candidate
		tar = c[1] ## reactions associated with tar candidate
		for s in src:
			met_src = met_associations[s][0]
		for t in tar:
			met_tar = met_associations[t][1]
		for m1 in met_src:
			for m2 in met_tar:
				candidate_pairs.append((m1, m2))
	# for s in src_rxn:
	# 	srcs = met_associations[s][0]
	# 	for src in srcs:
	# 		if src not in candidate_src:
	# 			candidate_src.append(src)
	# for t in tar_rxn:
	# 	tars = met_associations[t][0]
	# 	for tar in tars:
	# 		if tar not in candidate_tar:
	# 			candidate_tar.append(tar)
	return candidate_pairs


def parse(): 
	ptg = path_to_genes()
	candidates = get_src_tar_genes(ptg)
	candidate_reactions = get_reactions_of_interest(candidates, GPR.gpr_genes())
	candidate_pairs = get_candidate_nodes(candidate_reactions, GPR.get_metabolite_associations('./RECON1.json'))
	#print len(candidate_pairs)
	#print candidate_pairs
	src_tar = {}
	for pair in candidate_pairs:
		src = pair[0]
		tar = pair[1]
		if src not in src_tar:
			src_tar[src] = [tar]
		else:
			src_tar[src].append(tar)
	#c = 0
	return src_tar
	# for s in src_tar:
	# 	if len(src_tar[s]) >= 3:
	# 		#c += 1
	# 		print 'src: ' + s
	# 		print 'targets: ' + str(src_tar[s])
	#print c


if __name__ == '__main__':
	src_tar = parse()
	print src_tar














	