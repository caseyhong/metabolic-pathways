import json 
from gpr_mapping import * #aggregate, parse_gpr_mapping 
from geneid_mapping import * #prepareEntrez,ensemblToEntrez
from diff_paths import * #difference(x,y), createNetwork
from make_hypergraph import * #make_hypergraph 
from igraph import * 
from colon_preprocessor import * #colon_translator

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
	if 'colon' in fileName: 
		colon_dict = colon_translator()

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
			if 'colon' in fileName and gene in colon_dict.keys(): 
				gene = colon_dict[gene.replace('\"','')]
			exp = row[pid]

			# Convert from Ensembl -> Entrez Gene Id 
			if gene in ensemblEntrezDict.keys(): 
				geneids = ensemblEntrezDict[gene]
				for g in geneids: #Several Gene Id Transcript Versions 
					data[g] = exp 

	return patientName,data


	#List of patients , rxnid, value 
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

def generateDigraph(patientToRxns): 

	#Initialize Hypergraph for all patients for tumor type 
	hg = make_hypergraph()

	#Add edge weights for each patient 
	for patient in patientToRxns.keys(): 
		rxn_expVal = patientToRxns[patient]
		hg = add_to_HG(hg,rxn_expVal)

	#Create final digraph
	dg = convert_to_DG(hg)
	return dg 

def mcf(fileName,cancer_ids,healthy_ids): 
	"""
	mcf returns the final differential expressed graph for a particular cancer file 
	INPUTS: 
		fileName - filename to parse from (needs preprocessing to have certain structure)
		cancer_ids - list of cancer sample ids 
		healthy_ids - list of healthy sample ids 
	"""
	ensemblEntrezDict = createGeneIdMapping()
	recon1RxnMappings = parse_gpr_mapping('./RECON1.json')
	rxnMetaboliteMapping = get_metabolite_associations('./RECON1.json')
	metabolite_pairs = get_metabolite_pairs(rxnMetaboliteMapping)

	# For Cancerous Patients 
	c_patient_rxns = patientRxnMappings(fileName,cancer_ids,recon1RxnMappings,ensemblEntrezDict)
	cancer_dg = generateDigraph(c_patient_rxns)

	# For Healthy Patients 
	h_patient_rxns = patientRxnMappings(fileName,healthy_ids,recon1RxnMappings,ensemblEntrezDict)
	healthy_dg = generateDigraph(h_patient_rxns)

	# Perform graph subtraction 
	diff_g = difference(cancer_dg,healthy_dg)
	return diff_g

 



if __name__ == '__main__':
	lung_file = '../../../../../Desktop/RNASeq_Files/GSE81089_FPKM_cufflinks_nslc.tsv'
	lung_diff_g = mcf(lung_file,[1,2],[39,44])

	colon_file = '../../../../../Desktop/RNASeq_Files/GSE41258_series_matrix_colon.txt'
	x = mcf(colon_file,[1],[1])
	print x 






	# ensemblEntrezDict = createGeneIdMapping()
	# recon1RxnMappings = parse_gpr_mapping('./RECON1.json')
	# rxnMetaboliteMapping = get_metabolite_associations('./RECON1.json')
	# metabolite_pairs = get_metabolite_pairs(rxnMetaboliteMapping)
	# patientToRxns = patientRxnMappings('../../../../../Desktop/RNASeq_Files/GSE81089_FPKM_cufflinks_nslc.tsv',[1,2],recon1RxnMappings,ensemblEntrezDict)
	# dg = generateDigraph(patientToRxns)
	# print dg















