import json

# TODO get metabolites and stoich constants
def parse_recon(filename):
	relevant_genes = set()
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
						g=g.replace('(','')
						g=g.replace(')','')
						relevant_genes.add(g)
	return relevant_genes

if __name__ == '__main__':
	genes_of_interest = parse_recon('./RECON1.json')
	open('./recon1_genes.txt', 'w').close() #clears file 
	output = open('./recon1_genes.txt','r+')

	for gene in genes_of_interest:
		output.write(gene+",")
	print 'finished writing!'

