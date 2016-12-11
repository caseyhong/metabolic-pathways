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

	#Invert the weights 
	hg,inverted_weights = invert_weights(hg)

	
	#Create final digraph
	#dg = convert_to_DG(hg)
	dg = convert_to_DG(hg,inverted_weights)
	return dg 

def mcf(lung_file,lung_ids,colon_file,colon_ids):
	"""
	mcf returns the final differential expressed graph for a particular cancer file 
	INPUTS: 
		lung_file - filename to parse for lung data 
		lung_ids - index values of which samples to take from file 
		colon_file - filename to parse for colon data 
		colon_ids - index values of which samples to take from file 
	"""
	ensemblEntrezDict = createGeneIdMapping()
	recon1RxnMappings = parse_gpr_mapping('./RECON1.json')
	rxnMetaboliteMapping = get_metabolite_associations('./RECON1.json')
	metabolite_pairs = get_metabolite_pairs(rxnMetaboliteMapping)

	# For Lung Patients 
	print 'Processing Lung DG'
	lung_patient_rxns = patientRxnMappings(lung_file,lung_ids,recon1RxnMappings,ensemblEntrezDict)
	lung_dg = generateDigraph(lung_patient_rxns)
	#print 'lungdg_',lung_dg
	# for e in lung_dg.es: 
	# 	print e['weight']

	# For Colon Patients 
	print 'Processing Colon DG'
	colon_patient_rxns = patientRxnMappings(colon_file,colon_ids,recon1RxnMappings,ensemblEntrezDict)
	colon_dg = generateDigraph(colon_patient_rxns)

	# Perform graph subtraction 
	print 'Performing Graph Subtraction'
	diff_g = difference(lung_dg,colon_dg)
	return diff_g

 
# source - a list containing the source vertex IDs which should be included in the result. If None, all vertices will be considered.
# target - a list containing the target vertex IDs which should be included in the result. If None, all vertices will be considered.
# weights - a list containing the edge weights. It can also be an attribute name (edge weights are retrieved from the given attribute) or None (all edges have equal weight).
# mode - the type of shortest paths to be used for the calculation in directed graphs. OUT means only outgoing, IN means only incoming paths. ALL means to consider the directed graph as an undirected one.
# Returns:
# the shortest path lengths for given vertices in a matrix

def run_dijkstra(graph):
	print 'Running Bellman Ford on Differential Graph'
	x = graph.shortest_paths()
	return x 
	#graph.shortest_paths_dijkstra(source=None, target=None, weights=None, mode=OUT)
	#print x 



if __name__ == '__main__':
	lung_file = '../../../../../Desktop/RNASeq_Files/GSE81089_FPKM_cufflinks_nslc.tsv'
	colon_file = '../../../../../Desktop/RNASeq_Files/GSE41258_series_matrix_colon.txt'

	diff_g = mcf(lung_file,[1,2],colon_file,[1,2])
	x = run_dijkstra(diff_g)
	print x 
	# print diff_g
	# for e in diff_g.es: 
	# 	print e['weight']

	#run_dijkstra(diff_g) 



















