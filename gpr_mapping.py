import json

# TODO get metabolites and stoich constants
def parse_recon(filename):
	relevant_genes = set()
	filtered_probes = set() 
	bools = ['or', 'and']
	with open(filename) as json_file:
		data = json.load(json_file)
		for rxn in data["reactions"]:
			if len(rxn["gene_reaction_rule"]) > 0:
				rule = rxn["gene_reaction_rule"]
				words = rule.split()
				genes = [w for w in words if w not in bools]
				for g in genes:
					if g not in relevant_genes:
						g = g.replace('(','')
						g = g.replace(')','')
						relevant_genes.add(g)
						filtered_g = g.split('_AT')[0]+"_at" #removes numbers, lowercases 
						filtered_probes.add(filtered_g) 

	return relevant_genes, filtered_probes

if __name__ == '__main__':
	genes_of_interest,filtered_probes = parse_recon('./RECON1.json')
	open('./recon1_genes.txt', 'w').close() #clears file 
	output = open('./recon1_genes.txt','r+')
	fg = open('./recon1_filteredGenes.txt','w')

	for gene in filtered_probes: 
		fg.write(gene+",")

	for gene in genes_of_interest:
		output.write(gene+",")

	print 'finished writing!'

