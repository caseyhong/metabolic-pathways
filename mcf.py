import json 
from gpr_mapping import * #aggregate, parse_gpr_mapping 
from geneid_mapping import * #

def createGeneIdMapping(): 
	entrezDict = prepareEntrez()
	ensemblEntrezDict = ensemblToEntrez(entrezDict)
	return ensemblEntrezDict

def getPatientData(fileName,pid,ensemblEntrezDict): 
	"""
	Returns the patient's expression data specified by pid 
	INPUTS: 
		fileName - file specifying expression data matrix 
		pid - patient id specifying which patient from file 
		ensemblEntrezDict - created from createGeneIdMapping(), ensembl id to entrez id dictionary 
	OUTPUTS: 
		patientName - patient (sample) name 
		data - patient's exprsesion val dictionary {Entrez Gene Id: Expression Value}
	"""
	data = {}

	with open(fileName,'r+') as expData: 
		header = expData.readline().split('\t')
		patientName = header[pid]
		numPatients = len(header) - 1 
		if pid > numPatients: 
			print 'Invalid patient index!'
			return 

		for row in expData: 
			row = row.split('\t')
			gene = row[0]
			exp = row[pid]

			# Convert from Ensembl -> Entrez Gene Id 
			if gene in ensemblEntrezDict.keys(): 
				geneids = ensemblEntrezDict[gene]
				for g in geneids: #Several Gene Id Transcript Versions 
					data[g] = exp 

	return patientName,data


	#List of patients 
def patientRxnMappings(RNASeqFileName,patientIds,recon1RxnMappings,ensemblEntrezDict): 
	patientToRxns = {}

	for pid in patientIds: 
		patientName, pdata = getPatientData(RNASeqFileName,pid,ensemblEntrezDict)

		patientToRxns[patientName] = {}

		for rxnId in recon1RxnMappings.keys(): 
			rxn = recon1RxnMappings[rxnId]
			rxnVal = aggregate(rxn,pdata,pdata.keys())
			patientToRxns[patientName][rxnId] = rxnVal

	return patientToRxns

if __name__ == '__main__':
	ensemblEntrezDict = createGeneIdMapping()
	recon1RxnMappings = parse_gpr_mapping('./RECON1.json')
	patientToRxns = patientRxnMappings('../../../../../Desktop/RNASeq_Files/GSE81089_FPKM_cufflinks_nslc.tsv',[1],recon1RxnMappings,ensemblEntrezDict)
	
	#print patientToRxns

	# # print mappings[mappings.keys()[9]]
	# rxn = mappings[mappings.keys()[5]]
	# pName, pdata = getPatientData('../../../../../Desktop/RNASeq_Files/GSE81089_FPKM_cufflinks_nslc.tsv',3,ensemblEntrezDict)
	# #print pdata.keys()
	# aggValue = aggregate(rxn,pdata,pdata.keys())
	# print aggValue
	#aggregate()
	#aggregate(mapping,expression,genes)


	# 	"""
	# aggregate returns the biochemical expression value of a reaction given set of genes and corresponding express. value
	# INPUTS:
	# 	mapping - individual reaction 
	# 	expression - dictionary of gene:expressionval 
	# 	genes - list of genes
	# OUTPUT: 
	# 	biochemical exp. value 
	# """

	#patientDataFile = open('patient_expression.json','w')
	# patientDataFile.write(pdata)
	#json.dump(pdata,patientDataFile)


# for patient in range(1,numPatients+1): 
# 			patientData[header[patient]] = {}
# 			expData.seek(1)
# 			for row in expData: 
# 				row = row.split('\t')
# 				gene = row[0]
# 				exp = row[patient]
# 				patientData[header[patient]][gene] = exp 