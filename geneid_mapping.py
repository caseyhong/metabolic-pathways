#Generates dictionary for entrez ids in Recon1 
def prepareEntrez(): 
	entrezDict = {}
	entrez = open('recon1_entrez_locus.txt','r+')
	for eid in entrez:
		eid = eid.replace('\n','')
		if eid not in entrezDict.keys(): 
			eid,locus = eid.split('.')
			biggID = eid+"_AT"+locus 
			entrezDict[biggID] = True 
	return entrezDict

#Returns the dictionary that maps ensembl keys to entrez keys found in Recon 1 
def ensemblToEntrez(reconEntrez): 

	#Generates the entrez w/o loci -> ensembl 
	entrezEnsemblDict = {}
	ee = open('recon1_entrez_ensembl.txt','r+')
	for line in ee: 
		line = line.replace('\n','')
		ensembl,entrez = line.split(',')
		entrezEnsemblDict[entrez] = ensembl

	#Generates the ensembl -> entrez w/ loci mapping 
	#Even without loci from Biomart, loci are added to ensembl mapping 
	#where numbers only mean different loci for same gene 
	ensemblEntrezDict = {}
	for entrezLociID in reconEntrez:
		entrezOnlyID = entrezLociID.split('_AT')[0]
		if entrezOnlyID in entrezEnsemblDict: 
			ensemblKey = entrezEnsemblDict[entrezOnlyID]
			if ensemblKey in ensemblEntrezDict: 
				if entrezLociID not in ensemblEntrezDict[ensemblKey]: 
					ensemblEntrezDict[ensemblKey].append(entrezLociID)
			else: 
				ensemblEntrezDict[ensemblKey] = [entrezLociID]

	return ensemblEntrezDict

# FOR TESTING PURPOSES ONLY
# To confirm recon1 json and recon1 entrez map
def mapToEntrez(entrezDict): 
	idMap = {} # key: BIGGID, value: entrezid 
	reconGenes = open('recon1_json_genes.txt','r+')
	reconGenes = reconGenes.readline().split(',')
	for rGene in reconGenes:
		rGene = rGene.upper()
		if rGene in entrezDict.keys() and rGene not in idMap: 
			idMap[rGene] = True 
	return idMap

if __name__ == '__main__':
	entrezDict = prepareEntrez()
	ensemblEntrezDict = ensemblToEntrez(entrezDict)
	print ensemblEntrezDict
