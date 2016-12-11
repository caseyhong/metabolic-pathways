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
	src = set()
	tar = set()
	for path in path_to_genes:
		genes = path_to_genes[path]
		src.add(genes[0])
		tar.add(genes[-1])
	return src,tar

def get_reactions_of_interest(src_genes, tar_genes, gene_to_reaction):
	src_rxn = []
	tar_rxn = []
	for s in src_genes:
		rxn = gene_to_reaction[s]
		for i in rxn:
			if i not in src_rxn:
				src_rxn.append(i)
	for t in tar_genes:
		rxn = gene_to_reaction[t]
		for i in rxn:
			if i not in tar_rxn:
				tar_rxn.append(i)
	return src_rxn, tar_rxn

def get_candidate_nodes(src_rxn, tar_rxn, met_associations):
	candidate_src = []
	candidate_tar = []
	for s in src_rxn:
		srcs = met_associations[s][0]
		for src in srcs:
			if src not in candidate_src:
				candidate_src.append(src)
	for t in tar_rxn:
		tars = met_associations[t][0]
		for tar in tars:
			if tar not in candidate_tar:
				candidate_tar.append(tar)
	return candidate_src, candidate_tar

if __name__ == '__main__':
	ptg = path_to_genes()
	src_genes = get_src_tar_genes(ptg)[0]
	tar_genes = get_src_tar_genes(ptg)[1]
	src_rxn, tar_rxn = get_reactions_of_interest(src_genes, tar_genes, GPR.gpr_genes())
	candidate_src, candidate_tar = get_candidate_nodes(src_rxn, tar_rxn, GPR.get_metabolite_associations('./RECON1.json'))
	print len(candidate_src)
	print len(candidate_tar)