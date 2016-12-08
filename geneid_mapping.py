def prepareEnsembl(): 
	ensemblDict = {}
	mart = open('mart_export.txt','r+')#Gene ID,Affy HG U133-PLUS-2 probeset
	for line in mart: 
		geneid,affy = line.split(',')
		affy = affy.replace('\n','')
		ensemblDict[affy] = geneid
	print len(ensemblDict)
	return ensemblDict

def mapToEnsembl(ensembleDict): 
	finalDict = {}
	reconGenes = open('recon1_genes.txt','r+')
	reconGenes = reconGenes.readline().split(',')
	for rGene in reconGenes:
		rGene = rGene.lower()
		rGene = rGene.split('_at')[0]+'_at'
		if rGene in ensembleDict.keys(): 
			finalDict[ensembleDict[rGene]] = rGene
	return finalDict

if __name__ == '__main__':
	ensembleDict = prepareEnsembl() 
	#print 'ENSEMBL DICT: ',ensembleDict
	geneID_Affy = mapToEnsembl(ensembleDict)
	print geneID_Affy
	print 'done!'
