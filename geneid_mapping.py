def prepareEntrez(): 
	entrezDict = {}
	entrez = open('recon1_entrez.txt','r+')
	for eid in entrez:
		eid = eid.replace('\n','')
		if eid not in entrezDict.keys(): 
			entrezDict[eid] = True 
	return entrezDict

def mapToEntrez(entrezDict): 
	probeMap = {} # key: probeid, value: entrezid 
	reconGenes = open('recon1_json_genes.txt','r+')
	reconGenes = reconGenes.readline().split(',')
	for rGene in reconGenes:
		rGene = rGene.lower()
		mappableGeneID = rGene.split('_at')[0]
		if mappableGeneID in entrezDict.keys() and mappableGeneID not in probeMap: 
			probeMap[mappableGeneID+"_at"] = mappableGeneID

	return probeMap

if __name__ == '__main__':
	entrezDict = prepareEntrez()
	probe_entrez = mapToEntrez(entrezDict)
	print probe_entrez
	print len(probe_entrez)

#12/08/2016 133_plus2.csv only maps to 4 in recon as well... 
#12/08/2016 tried mapping 133a_2.tsv , only got 5-7 mappings from Recon1


