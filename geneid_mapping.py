def prepareEnsembl():
	entrezDict = {}
	affymetrix = open('133_plus2.csv','r+')#Human Genom U133 2.0 Array 
	for line in affymetrix: 
		entrezid,probeid= line.split('\t')
		probeid = probeid.replace('\n','')
		for pid in probeid.split(','): 
			if '_s_at' not in pid and '_x_at' not in pid: 
				if pid in entrezDict.keys(): 
					entrezDict[pid].append(entrezid)
				else: 
					entrezDict[pid] = [entrezid]
	return entrezDict

def mapToEnsembl(ensemblDict): 
	finalDict = {}
	reconGenes = open('recon1_genes.txt','r+')
	reconGenes = reconGenes.readline().split(',')
	for rGene in reconGenes:
		rGene = rGene.lower()
		rGene = rGene.split('_at')[0]+'_at'
		if rGene in ensemblDict.keys() and rGene not in finalDict.keys():
			finalDict[rGene] = ensemblDict[rGene]
	return finalDict

if __name__ == '__main__':
	entrezDict = prepareEnsembl() 
	geneID_Affy = mapToEnsembl(entrezDict)
	print geneID_Affy
	print 'done!'

#12/08/2016 133_plus2.csv only maps to 4 in recon as well... 


#12/08/2016 tried mapping, only got 5-7 mappings from Recon1
# def prepareEnsembl():
# 	entrezDict = {} 
# 	ensemblDict = {}
# 	reverseEnsemblDict = {}
# 	affymetrix = open('133a_2.tsv','r+')#Human Genom U133 2.0 Array 
# 	for line in affymetrix: 
# 		probeid, gene,entrezid,ensembl= line.split('\t')
# 		if '_s_at' not in probeid and '_x_at' not in probeid: 
# 			ensembl = ensembl.replace('\n','')	

# 			if '---' != entrezid: 
# 				for eid in entrezid.split(' /// '): 
# 					entrezDict[eid] = probeid 

# 			if '---' != ensembl: 
# 				for enid in ensembl.split(' /// '):
# 					if 'ENSG' in enid: 
# 						ensemblDict[enid] = probeid

# 						if probeid in reverseEnsemblDict.keys(): 
# 							reverseEnsemblDict[probeid].append(enid)
# 						else: 
# 							reverseEnsemblDict[probeid] = [enid]

# 	return ensemblDict,reverseEnsemblDict,entrezDict

